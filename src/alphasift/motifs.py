"""
The motif engine: do unrelated partners put similar residues on the same bait residues, at the
same spacing along their own sequence, more often than shuffled letters would?

Extracted verbatim from viz_int.motif_engine (2026-10-02: tie-safe p-values, every interface,
one-motif-or-several grouping) so that its validation applies to this code unchanged. Changes:
`partner_3di` and `run_motif_engine` take each prediction's model file (`cif_of`) instead of
globbing a screen folder, since AlphaSIFT also reads models kept only as `.cif.zst`; the Foldseek
binary is found via `foldseek_bin`, $ALPHASIFT_FOLDSEEK or PATH.

Secondary structure is a choice (2026-10-06, user decision), `ss`:
  "dssp" (default)  as in viz_int: letters are shuffled only within a DSSP class (helix / strand /
                    other), so that "every partner binds with a helix" is not read as a motif. Needs
                    mkdssp. This is the validated engine; outputs are identical to viz_int's.
  "none"            letters are shuffled among ALL residues placed on the interface; no external
                    program, several times faster. Not validated: the viz_int numbers (known motifs
                    12/13 etc.) were measured with "dssp".
Validation: viz_int VALIDATION_LOG (motif stage, 2026-09-29 ..
2026-10-02), IMPROVEMENTS #34 and its addendum.

Per confidence interface (every one, significant or not; each result carries the flag):
  positions  the interface's bait residues. A partner residue counts at every bait residue
             it touches: confident pair (s = 1/(1+(PAE/d0)^2) >= tau) with Cb-Cb <= touch_cb
             (virtual Cb for Gly); one vote per partner per position.
  groups     partners >= id_group sequence identity are one group (paralogs count once).
  letters    engine "c": amino acid, similarity BLOSUM62.
             engine "3di_aa": (amino acid, Foldseek 3Di state of the partner chain on its
             own), similarity BLOSUM62 + w3di x 3Di matrix (the 3Di matrix's spread is ~3x
             BLOSUM62's; w3di = 1/3 gives the two equal weight).
  null       letters shuffled among the residues placed on the interface -- within DSSP class
             (ss="dssp") or all together (ss="none"), see above -- n_perm times; BH over the tests.
  Panel 1    per bait residue touched by >= min_groups groups: between-group letter agreement
             vs the null.
  Panel 2    rules (2-3 positions filled at a constant spacing <= max_span residues along the
             partner, >= min_groups groups, agreement vs the null) -> parts -> ELEMENTS
             (assembled when carried by the same groups at a constant short linker, same
             direction and spacing, coherent) -> near-duplicates dropped -> element FILTER:
             a partner stays only if it carries the element (>= 2 and >= half of its positions
             similar (>= 1) to the other groups' most frequent letter; fills >= half of the
             element's positions); elements left with < min_groups groups are dropped.
             Element pairs: bipartite (co-carried beyond chance, long linker) / alternative.
  controls   n_ctrl shuffled-letter Panel-1 runs (false columns) and one whole Panel-2 run on
             shuffled letters, filter included (false elements).
  motifs     (2026-10-02) one motif or several: elements are ONE motif when the same partner groups
             carry both beyond chance, or when they share bait residues with similar preferred amino acids
             there -- never when their carriers overlap less than chance. Motif pairs: same_pocket /
             separate_sites; status "different" only with that positive evidence, else "unresolved".
  scope      every confidence interface (significant or not, flagged): the bait-side Panel 1 exists
             whenever any partner touches the interface, whether or not a motif element is found.
Needs: mkdssp (DSSP 4.x) on PATH for ss="dssp"; engine "3di_aa" also a Foldseek binary (createdb) --
`foldseek_bin`, $ALPHASIFT_FOLDSEEK or `foldseek` on PATH.
"""

from __future__ import annotations

import collections
import itertools
import os
import shutil
import subprocess
import tempfile
import zlib
from pathlib import Path

import numpy as np
from Bio.Align import PairwiseAligner, substitution_matrices
from scipy.stats import false_discovery_control, hypergeom

from alphasift.structures import aa1

AA = "ARNDCQEGHILKMFPSTWYV"
AI = {a: i for i, a in enumerate(AA)}
_B = substitution_matrices.load("BLOSUM62")
BM_AA = np.array([[_B[a, b] for b in AA] for a in AA], float)

ENGINES = ("c", "3di_aa")
# A shuffled score counts as ">= observed" when it is within TIE_EPS of it. Equal letter sets give equal scores in
# exact arithmetic, but with a non-integer 3Di weight the float sums can differ in the last bits and break ties
# arbitrarily (port check 2026-10-01: weight 1/3 vs 0.333333 changed Panel 1 on 30/31 baits). Engine "c" sums
# exact values, so this never changes it (checked on all 31 baits).
TIE_EPS = 1e-9
ENGINE_DEFAULTS = dict(TAU=0.5, D0=10.0, TOUCH_CB=7.0, ID_GROUP=0.40, N_PERM=1000, MIN_GROUPS=3,
                       MAX_SPAN=10, Q=0.05, ASM_OVERLAP=0.7, ASM_LINK=4, BIP_LINK=5, ORIENT=0.8,
                       ASM_TOL=1, ASM_CONSIST=0.8, EL_CONSIST=0.9, REDUNDANT=0.7, N_CTRL=3,
                       CHUNK=250, W3DI=1 / 3, EL_FILTER=True, SS="dssp")


def _load_3di_matrix():
    rows = [l.split() for l in (Path(__file__).with_name("mat3di.txt")).read_text().splitlines()
            if l.strip() and not l.startswith("#")]
    L3 = [x for x in rows[0] if x != "X"]
    I3 = {c: i for i, c in enumerate(L3)}
    M3 = np.array([[float(r[1 + rows[0].index(b)]) for b in L3] for r in rows[1:] if r[0] in I3], float)
    return L3, I3, M3


# ------------------------------------------------------------------ inputs from model files
def dssp_classes(cif, bait_chain):
    """{partner residue: 'H' (DSSP H/G/I) | 'E' (E/B) | 'C'} for the non-bait chain of one model
    (mkdssp 4.x, classic output). None if mkdssp fails."""
    with tempfile.NamedTemporaryFile(suffix=".dssp") as tf:
        r = subprocess.run(["mkdssp", "--output-format", "dssp", str(cif), tf.name],
                           capture_output=True, timeout=300)
        if r.returncode != 0:
            return None
        lines = Path(tf.name).read_text().splitlines()
    k = next(i for i, l in enumerate(lines) if l.startswith("  #  RESIDUE"))
    out = {}
    for l in lines[k + 1:]:
        if l[13] == "!" or l[11] == bait_chain:
            continue
        ss = l[16]
        out[int(l[5:10])] = "H" if ss in "HGI" else "E" if ss in "EB" else "C"
    return out


def partner_3di(cif_of, bait, preds, foldseek_bin=None, threads=8):
    """{prediction: {partner residue: 3Di state}} -- Foldseek createdb 3Di of the partner chain
    (computed per chain, so the bait does not shape it), read from the same model CIFs the
    pipeline parses. Strings are mapped to the chain's residues with a CA, in order, and kept only
    when length and amino-acid string match (a mismatch is skipped, never guessed). Checked
    identical to the validation letters (80_, cropped dimers) on PARVA: 18,729 / 18,729."""
    from Bio.PDB import MMCIFParser
    fs = foldseek_bin or os.environ.get("ALPHASIFT_FOLDSEEK") or shutil.which("foldseek")
    if not fs:
        raise RuntimeError("engine '3di_aa' needs a Foldseek binary: pass foldseek_bin, set "
                           "ALPHASIFT_FOLDSEEK, or put foldseek on PATH")
    out = {}
    with tempfile.TemporaryDirectory() as td:
        td = Path(td); (td / "in").mkdir()
        cifs = {}
        for pr in preds:
            cif = Path(cif_of[pr])
            (td / "in" / f"{pr}.cif").symlink_to(cif.resolve()); cifs[pr] = cif
        run = lambda *a: subprocess.run([fs, *a, "-v", "0"], capture_output=True, check=True)
        run("createdb", str(td / "in"), str(td / "db"), "--threads", str(threads))
        run("lndb", str(td / "db_h"), str(td / "db_ss_h"))
        run("convert2fasta", str(td / "db_ss"), str(td / "ss.fa"))
        run("convert2fasta", str(td / "db"), str(td / "aa.fa"))

        def fasta(f):
            seqs, name = {}, None
            for l in f.read_text().splitlines():
                if l.startswith(">"):
                    name = l[1:].split()[0]; seqs[name] = ""
                else:
                    seqs[name] += l.strip()
            return seqs
        s3, sa = fasta(td / "ss.fa"), fasta(td / "aa.fa")
    parser = MMCIFParser(QUIET=True)
    for pr, cif in cifs.items():
        model = parser.get_structure(pr, str(cif))[0]
        ic = [c.id for c in model if c.id != bait["chain_ids"][pr]][0]
        res = [(r.id[1], aa1(r.resname)) for r in model[ic] if r.id[0] == " " and "CA" in r]
        a3, aa = s3.get(f"{pr}_{ic}"), sa.get(f"{pr}_{ic}")
        if a3 and aa and len(a3) == len(res) and aa == "".join(a for _, a in res):
            out[pr] = {r: c for (r, _), c in zip(res, a3)}
    return out


# ------------------------------------------------------------------ engine
class _Engine:
    def __init__(self, engine, seed_key, ifaces, pairs, seq, ss, bait_cb, part_cb, letters3, params):
        assert engine in ENGINES, engine
        self.engine, self.seed_key, self.P = engine, seed_key, params
        self.ifaces, self.pairs, self.seq, self.ss = ifaces, pairs, seq, ss
        self.stratify = params.get("SS", "dssp") == "dssp"
        self.bait_cb, self.part_cb, self.t3 = bait_cb, part_cb, letters3 or {}
        self.BM = BM_AA
        if engine == "3di_aa":
            self.L3, self.I3, self.M3 = _load_3di_matrix()
            self.BM = (BM_AA[:, None, :, None] + params["W3DI"] * self.M3[None, :, None, :]).reshape(400, 400)

    # -- partner groups, letters
    def seeded(self, *key):
        return np.random.default_rng(zlib.crc32("|".join(map(str, key)).encode()))

    def partner_groups(self, prs):
        al = PairwiseAligner(mode="global", substitution_matrix=_B, open_gap_score=-10, extend_gap_score=-0.5)
        seqs = {p: "".join(self.seq[p][k] for k in sorted(self.seq[p])) for p in prs}
        km = {p: {s[i:i + 3] for i in range(len(s) - 2)} for p, s in seqs.items()}
        par = {p: p for p in prs}

        def find(x):
            while par[x] != x:
                par[x] = par[par[x]]; x = par[x]
            return x
        for a, b in itertools.combinations(prs, 2):
            if not km[a] or not km[b] or len(km[a] & km[b]) / len(km[a] | km[b]) < 0.05:
                continue
            aln = al.align(seqs[a], seqs[b])[0]
            ident = sum(1 for (s1, e1), (s2, e2) in zip(*aln.aligned)
                        for x, y in zip(seqs[a][s1:e1], seqs[b][s2:e2]) if x == y)
            if ident / min(len(seqs[a]), len(seqs[b])) >= self.P["ID_GROUP"]:
                par[find(a)] = find(b)
        return {p: find(p) for p in prs}

    def letter_array(self, p):
        s = self.seq[p]; arr = np.full(max(s) + 1, -1)
        if self.engine == "c":
            for k, a in s.items():
                arr[k] = AI.get(a, -1)
        else:
            t3 = self.t3.get(p, {})
            for k, a in s.items():
                c3 = self.I3.get(t3.get(k, ""), -1)
                if c3 >= 0 and AI.get(a, -1) >= 0:
                    arr[k] = AI[a] * 20 + c3
        return arr

    def LET(self, code):
        return AA[code] if self.engine == "c" else AA[code // 20]

    def ss_of(self, p, j):
        return (self.ss.get(p) or {}).get(j, "C")

    # -- positions (C: bait residues)
    def obs_c(self, iface):
        out = []
        for p in iface["partners"]:
            b, j, pae = self.pairs[p]
            s = 1 / (1 + (pae / self.P["D0"]) ** 2)
            m = np.isin(b, iface["residues"]) & (s >= self.P["TAU"])
            b, j = b[m], j[m]
            d = np.array([np.linalg.norm(self.bait_cb[p][x] - self.part_cb[p][y])
                          if x in self.bait_cb[p] and y in self.part_cb[p] else np.inf for x, y in zip(b, j)])
            touch = d <= self.P["TOUCH_CB"]
            pk_of = collections.defaultdict(set)
            for x, y in zip(b[touch].tolist(), j[touch].tolist()):
                pk_of[y].add(x)
            per = collections.defaultdict(list)
            for y, pks in pk_of.items():
                for x in pks:
                    per[x].append([y, 1.0])
            for x, lst in per.items():
                tot = sum(w for _, w in lst)
                out += [(x, p, y, w / max(1.0, tot)) for y, w in lst]   # one vote per partner per pocket
        return out

    # -- null, agreement
    def make_null(self, obs, real, rng):
        """Letters shuffled among the residues placed on the interface: one pool per DSSP class
        (ss="dssp", as viz_int) or one pool for all (ss="none")."""
        used = sorted({(p, j) for _, p, j, _ in obs if j < len(real[p]) and real[p][j] >= 0})
        cls = [self.ss_of(p, j) if self.stratify else "all" for p, j in used]
        vals = np.array([real[p][j] for p, j in used])
        idx_by = collections.defaultdict(list)
        for i, c in enumerate(cls):
            idx_by[c].append(i)
        perm = {p: np.tile(real[p].astype(np.int16), (self.P["N_PERM"], 1))
                for p in {u[0] for u in used} | {o[1] for o in obs if o[1] in real}}

        def one():
            v = vals.copy()
            for ids in idx_by.values():
                v[ids] = rng.permutation(vals[ids])
            return v
        for t in range(self.P["N_PERM"]):
            for (p, j), x in zip(used, one()):
                perm[p][t, j] = x
        ctrl = {p: real[p].copy() for p in perm}
        for (p, j), x in zip(used, one()):
            ctrl[p][j] = x
        return perm, ctrl

    def agreement(self, obs, arr, grp):
        CH = self.P["CHUNK"]
        if arr[obs[0][0]].ndim == 2 and arr[obs[0][0]].shape[0] > CH:      # chunk permutations (memory)
            n = arr[obs[0][0]].shape[0]; ps = {o[0] for o in obs}
            return np.concatenate([self.agreement(obs, {p: arr[p][s:s + CH] for p in ps}, grp)
                                   for s in range(0, n, CH)])
        cols = []
        for p, j, _ in obs:
            a = arr[p]
            cols.append(np.array([a[j] if j < len(a) else -1]) if a.ndim == 1 else
                        (a[:, j] if j < a.shape[1] else np.full(a.shape[0], -1)))
        Lm = np.stack(cols, 1)
        w = np.array([o[2] for o in obs], float)
        gid = np.array([hash(grp[o[0]]) for o in obs])
        W = np.triu(w[:, None] * w[None, :] * (gid[:, None] != gid[None, :]), 1)
        V = Lm >= 0; Li = np.maximum(Lm, 0)
        M = W[None] * (V[:, :, None] & V[:, None, :])
        den = M.sum((1, 2))
        return (self.BM[Li[:, :, None], Li[:, None, :]] * M).sum((1, 2)) / np.maximum(den, 1e-12)

    # -- Panel 1
    def panel1(self, obs, real, nulls, grp):
        cols = collections.defaultdict(list)
        for pos, p, j, w in obs:
            cols[pos].append((p, j, w))
        order = sorted(cols, key=lambda c: (c, str(c)))
        out, pv = [], {k: [] for k in nulls}
        for c in order:
            o = cols[c]; ng = len({grp[p] for p, _, _ in o})
            row = dict(pos=c, n_groups=ng, obs=o, agree=None, p={}, q={})
            if ng >= self.P["MIN_GROUPS"]:
                a = float(self.agreement(o, real, grp)[0]); row["agree"] = a
                for k, (perm, _) in nulls.items():
                    an = self.agreement(o, perm, grp)
                    row["p"][k] = (1 + int((an >= a - TIE_EPS).sum())) / (1 + self.P["N_PERM"]); pv[k].append(row["p"][k])
            out.append(row)
        for k in nulls:
            if pv[k]:
                it = iter(false_discovery_control(pv[k], method="bh"))
                for r in out:
                    if k in r["p"]:
                        r["q"][k] = float(next(it))
        return out

    def marks_and_repeats(self, cols, real, grp):
        marks, rep = [], collections.Counter()
        for c1, c2 in zip(cols, cols[1:]):
            by1 = collections.defaultdict(set); by2 = collections.defaultdict(set)
            for p, j, _ in c1["obs"]:
                by1[p].add(j)
            for p, j, _ in c2["obs"]:
                by2[p].add(j)
            both = set(by1) & set(by2)
            allg = {grp[p] for p in set(by1) | set(by2)}
            shared = {grp[p] for p in both if by1[p] & by2[p]}
            consec = {grp[p] for p in both if not (by1[p] & by2[p]) and any(abs(x - y) == 1 for x in by1[p] for y in by2[p])}
            marks.append(dict(shared=len(shared) / max(1, len(allg)), consecutive=len(consec) / max(1, len(allg))))

            def top(c):
                w = collections.Counter()
                for p, j, ww in c["obs"]:
                    if j < len(real[p]) and real[p][j] >= 0:
                        w[real[p][j]] += ww
                return w.most_common(1)[0][0] if w else None
            t = top(c1)
            if t is not None and t == top(c2):
                rep["repeated_pairs"] += 1
                for p in set(by1) | set(by2):
                    a1 = {j for j in by1.get(p, ()) if j < len(real[p]) and real[p][j] == t}
                    a2 = {j for j in by2.get(p, ()) if j < len(real[p]) and real[p][j] == t}
                    if a1 and a2:
                        rep["same_residue" if a1 & a2 else "consecutive" if any(abs(x - y) == 1 for x in a1 for y in a2)
                            else "other_same_partner"] += 1
                    elif a1 or a2:
                        rep["different_partners"] += 1
        return marks, dict(rep)

    # -- Panel 2: rules -> parts -> elements
    def rules(self, obs, lett, perm, grp):
        P = self.P
        by_p = collections.defaultdict(list)
        for pos, p, j, w in obs:
            by_p[p].append((pos, j))
        order = {pos: (pos, str(pos)) for pos in {o[0] for o in obs}}
        cand = collections.defaultdict(list)

        def add(p, combo):
            poss = [c[0] for c in combo]; js = [c[1] for c in combo]
            d = [x - js[0] for x in js[1:]]
            if not (all(x > 0 for x in d) or all(x < 0 for x in d)):
                return
            cand[(tuple(poss), tuple(abs(x) for x in d))].append((p, tuple(js), 1 if d[0] > 0 else -1))
        near = {}
        for p, lst in by_p.items():
            lst = sorted(set(lst), key=lambda x: order[x[0]])
            near[p] = [(a, b) for a, b in itertools.combinations(lst, 2)
                       if a[0] != b[0] and a[1] != b[1] and abs(a[1] - b[1]) <= P["MAX_SPAN"]]
            for a, b in near[p]:
                add(p, (a, b))
        freq2 = {k for k, v in cand.items() if len({grp[x[0]] for x in v}) >= P["MIN_GROUPS"]}
        for p, prs in near.items():
            good = [(a, b) for a, b in prs if ((a[0], b[0]), (abs(b[1] - a[1]),)) in freq2]
            gset = set(good) | {(b, a) for a, b in good}
            nodes = sorted({x for ab in good for x in ab}, key=lambda x: order[x[0]])
            for a, b, c in itertools.combinations(nodes, 3):
                if len({a[0], b[0], c[0]}) == 3 and len({a[1], b[1], c[1]}) == 3 and \
                   max(a[1], b[1], c[1]) - min(a[1], b[1], c[1]) <= P["MAX_SPAN"] and \
                   (a, b) in gset and (a, c) in gset and (b, c) in gset:
                    add(p, (a, b, c))
        out, pv = [], []
        for key, sup in cand.items():
            best = {}
            for p, js, o in sup:
                best.setdefault(p, (js, o))
            if len({grp[p] for p in best}) < P["MIN_GROUPS"]:
                continue
            T = float(np.mean([self.agreement([(p, js[i], 1.0) for p, (js, _) in best.items()], lett, grp)[0]
                               for i in range(len(key[0]))]))
            Tn = np.mean([self.agreement([(p, js[i], 1.0) for p, (js, _) in best.items()], perm, grp)
                          for i in range(len(key[0]))], 0)
            pval = (1 + int((Tn >= T - TIE_EPS).sum())) / (1 + P["N_PERM"])
            out.append(dict(pos=list(key[0]), spacing=list(key[1]), agree=T, p=pval,
                            sup={p: (list(js), o) for p, (js, o) in best.items()})); pv.append(pval)
        if pv:
            for r, q in zip(out, false_discovery_control(pv, method="bh")):
                r["q"] = float(q)
        return len(out), [r for r in out if r.get("q", 1) < P["Q"]]

    @staticmethod
    def parts_of(sig):
        kept = []
        for r in sorted(sig, key=lambda r: (r["q"], -len(r["sup"]))):
            s = {(p, j) for p, (js, _) in r["sup"].items() for j in js}
            if all(len(s & k) < 0.5 * min(len(s), len(k)) for _, k in kept):
                kept.append((r, s))
        return [r for r, _ in kept]

    @staticmethod
    def linker(A, B, p):
        a, b = A["sup"][p][0], B["sup"][p][0]
        if max(a) < min(b):
            return "A>B", min(b) - max(a) - 1
        if max(b) < min(a):
            return "B>A", min(a) - max(b) - 1
        return "inter", 0

    @staticmethod
    def _el_consistency(el):
        sup = collections.defaultdict(dict)
        for part in el:
            for p, (js, _) in part["sup"].items():
                for pos, j in zip(part["pos"], js):
                    sup[p][pos] = j
        poss = sorted({pos for part in el for pos in part["pos"]}, key=str)
        anchor = max(poss, key=lambda pos: sum(1 for m in sup.values() if pos in m))
        dirs, offs = {}, collections.defaultdict(list)
        for p, m in sup.items():
            ks = [pos for pos in el[0]["pos"] if pos in m] or [pos for pos in poss if pos in m]
            js = [m[pos] for pos in ks]
            dirs[p] = -1 if len(js) >= 2 and js[-1] < js[0] else 1
            if anchor in m:
                for pos in m:
                    offs[pos].append(dirs[p] * (m[pos] - m[anchor]))
        med = {pos: int(np.median(v)) for pos, v in offs.items()}
        hit = tot = 0
        for p, m in sup.items():
            if anchor in m:
                for pos, j in m.items():
                    if pos in med:
                        tot += 1; hit += (m[anchor] + dirs[p] * med[pos]) == j
        return hit / tot if tot else 1.0

    def _compatible(self, A, B, grp):
        P = self.P
        gA = {grp[p] for p in A["sup"]}; gB = {grp[p] for p in B["sup"]}
        if len(gA & gB) < P["ASM_OVERLAP"] * min(len(gA), len(gB)):
            return False
        common = sorted(set(A["sup"]) & set(B["sup"]))
        lk = [self.linker(A, B, p) for p in common]
        if not lk or sum(1 for o, x in lk if o == "inter" or x <= P["ASM_LINK"]) < 0.8 * len(lk):
            return False
        same_dir = [p for p in common if A["sup"][p][1] == B["sup"][p][1]]
        if len(same_dir) < P["ASM_CONSIST"] * len(common):
            return False
        deltas = [(B["sup"][p][0][0] - A["sup"][p][0][0]) * A["sup"][p][1] for p in same_dir]
        mode = collections.Counter(deltas).most_common(1)[0][0]
        return sum(1 for d_ in deltas if abs(d_ - mode) <= P["ASM_TOL"]) >= P["ASM_CONSIST"] * len(common)

    def elements(self, parts, grp):
        els = []
        for part in sorted(parts, key=lambda r: (r["q"], -len(r["sup"]))):
            placed = False
            for el in els:
                if any(self._compatible(part, q, grp) for q in el) and \
                   self._el_consistency(el + [part]) >= self.P["EL_CONSIST"]:
                    el.append(part); placed = True; break
            if not placed:
                els.append([part])
        return els

    def nonredundant(self, els, grp):
        def key(el):
            return (min(pt["q"] for pt in el), -len({grp[p] for pt in el for p in pt["sup"]}),
                    -len({pos for pt in el for pos in pt["pos"]}))
        kept = []
        for el in sorted(els, key=key):
            pos = {pos for pt in el for pos in pt["pos"]}; gs = {grp[p] for pt in el for p in pt["sup"]}
            if any(len(pos & kp) >= self.P["REDUNDANT"] * len(pos) and len(gs & kg) >= self.P["REDUNDANT"] * len(gs)
                   for _, kp, kg in kept):
                continue
            kept.append((el, pos, gs))
        return [el for el, _, _ in kept]

    def trim(self, el, grp, lett):
        """Element filter: keep only the supporters that carry the element."""
        m = collections.defaultdict(dict)
        for part in el:
            for p, (js, _) in part["sup"].items():
                for pos, j in zip(part["pos"], js):
                    m[p][pos] = j
        npos = len({pos for part in el for pos in part["pos"]})
        three = self.engine == "3di_aa"

        def letter(p, j):
            a = lett[p]
            if j >= len(a) or a[j] < 0:
                return -1
            return int(a[j]) // 20 if three else int(a[j])

        def state(p, j):
            a = lett[p]
            return int(a[j]) % 20 if three and j < len(a) and a[j] >= 0 else -1
        keep = set()
        for p in sorted(m):
            mp = m[p]
            if len(mp) < 2 or len(mp) < 0.5 * npos:
                continue
            ok = 0
            for pos, j in mp.items():
                a = letter(p, j)
                if a < 0:
                    continue
                cnt = collections.Counter(letter(q, mq[pos]) for q, mq in m.items() if pos in mq and grp[q] != grp[p])
                cnt.pop(-1, None)
                if cnt:
                    top = max(cnt.items(), key=lambda kv: (kv[1], -kv[0]))[0]
                    if three:                   # amino acid + w3di x 3Di, each vs the other groups' top
                        c3 = collections.Counter(state(q, mq[pos]) for q, mq in m.items()
                                                 if pos in mq and grp[q] != grp[p])
                        c3.pop(-1, None)
                        t3 = max(c3.items(), key=lambda kv: (kv[1], -kv[0]))[0]
                        ok += BM_AA[a, top] + self.P["W3DI"] * self.M3[state(p, j), t3] >= 1
                    else:
                        ok += BM_AA[a, top] >= 1
            if ok >= 2 and ok >= 0.5 * len(mp):
                keep.add(p)
        if len({grp[p] for p in keep}) < self.P["MIN_GROUPS"]:
            return None
        out = [dict(part, sup={p: v for p, v in part["sup"].items() if p in keep}) for part in el]
        return [part for part in out if part["sup"]]

    def describe(self, el, real, grp):
        sup = collections.defaultdict(dict)
        ori = collections.Counter()
        for part in el:
            for p, (js, o) in part["sup"].items():
                for pos, j in zip(part["pos"], js):
                    sup[p][pos] = j
                ori[(grp[p], o)] += 1
        poss = sorted({pos for part in el for pos in part["pos"]}, key=lambda c: (c, str(c)))
        letters = []
        for pos in poss:
            cnt = collections.Counter(self.LET(real[p][m[pos]]) for p, m in sup.items()
                                      if pos in m and m[pos] < len(real[p]) and real[p][m[pos]] >= 0)
            letters.append({a: n for a, n in cnt.most_common(4)})
        gs = {grp[p] for p in sup}
        par = sum(1 for g_ in gs if ori[(g_, 1)] >= ori[(g_, -1)] and ori[(g_, 1)])
        anti = sum(1 for g_ in gs if ori[(g_, -1)] > ori[(g_, 1)])
        call = ("parallel" if par >= self.P["ORIENT"] * len(gs) else
                "antiparallel" if anti >= self.P["ORIENT"] * len(gs) else "mixed")
        anchor = max(poss, key=lambda pos: sum(1 for m in sup.values() if pos in m))
        dirs, offs = {}, collections.defaultdict(list)
        for p, m in sup.items():
            ks = [pos for pos in poss if pos in m]; js = [m[pos] for pos in ks]
            dirs[p] = -1 if len(js) >= 2 and js[-1] < js[0] else 1
            if anchor in m:
                for pos in ks:
                    offs[pos].append(dirs[p] * (m[pos] - m[anchor]))
        med = {pos: int(np.median(v)) for pos, v in offs.items()}
        hit = tot = 0
        for p, m in sup.items():
            if anchor not in m:
                continue
            for pos, j in m.items():
                if pos in med:
                    tot += 1; hit += (m[anchor] + dirs[p] * med[pos]) == j
        return dict(n_parts=len(el), positions=[str(x) for x in poss], touched=None, n_groups=len(gs),
                    letters=letters, orientation=call, parallel=par, antiparallel=anti,
                    q=min(pt["q"] for pt in el), residues={p: sorted(set(m.values())) for p, m in sup.items()},
                    map={p: {str(k): v for k, v in m.items()} for p, m in sup.items()},
                    consistency=hit / tot if tot else None)

    def element_pairs(self, els, descs, grp, groups_present):
        N = max(1, len(groups_present)); out = []
        for i, j in itertools.combinations(range(len(els)), 2):
            if set(descs[i]["positions"]) & set(descs[j]["positions"]):
                continue
            rA = {(p, r) for p, rs in descs[i]["residues"].items() for r in rs}
            rB = {(p, r) for p, rs in descs[j]["residues"].items() for r in rs}
            if rA & rB:
                continue
            gA = {grp[p] for p in descs[i]["residues"]}; gB = {grp[p] for p in descs[j]["residues"]}
            k = len(gA & gB)
            pe = float(hypergeom.sf(k - 1, N, len(gA), len(gB))); pd_ = float(hypergeom.cdf(k, N, len(gA), len(gB)))
            lk = []
            for p in sorted(set(descs[i]["residues"]) & set(descs[j]["residues"])):
                a, b = descs[i]["residues"][p], descs[j]["residues"][p]
                lk.append(("A>B", min(b) - max(a) - 1) if max(a) < min(b) else
                          ("B>A", min(a) - max(b) - 1) if max(b) < min(a) else ("inter", 0))
            out.append(dict(a=i, b=j, k=k, nA=len(gA), nB=len(gB), p_enr=pe, p_dep=pd_, linkers=lk))
        if out:
            for key, qk in (("p_enr", "q_enr"), ("p_dep", "q_dep")):
                for r, q in zip(out, false_discovery_control([r[key] for r in out], method="bh")):
                    r[qk] = float(q)
        bip = [r for r in out if r["q_enr"] < self.P["Q"] and r["linkers"] and
               sum(1 for o, x in r["linkers"] if o != "inter" and x >= self.P["BIP_LINK"]) >= 0.8 * len(r["linkers"])]
        for r in bip:
            L = [x for o, x in r["linkers"] if o != "inter"]
            r["fixed_linker"] = bool(L) and max(L) - min(L) <= 2
        alt = [r for r in out if r["q_dep"] < self.P["Q"]]
        return len(out), bip, alt

    def motif_groups(self, descs, grp, groups_present, bip):
        """One motif or several (2026-10-02, user; revised the same day after 92b_ showed the first version split one
        known motif into up to 4 "motifs": a 3-group element can never reach significance, and "not shown to be the
        same" was being read as "different").
        Two elements are LINKED (same motif) when either
          carriers  the same partner groups carry both more often than chance (hypergeometric over the groups placed
                    on the interface, BH q < Q), or
          letters   they share bait residues and, at >= half of the shared residues (>= 1), their most frequent
                    partner amino acids are similar (BLOSUM62 >= 1 -- the element filter's rule),
        unless their carriers overlap LESS than chance (depletion q < Q): positive evidence that different partners
        use them, which is never merged. Linked elements form a MOTIF (connected components). Each pair of motifs gets
          relation  same_pocket (shared bait residues) | separate_sites
          status    "different" only with positive evidence (carrier depletion q < Q, i.e. `exclusive`), otherwise
                    "unresolved" -- too few carriers to tell. No new threshold: Q and the filter's BLOSUM62 >= 1."""
        n, N, Q = len(descs), max(1, len(groups_present)), self.P["Q"]
        G = [{grp[p] for p in d["residues"]} for d in descs]

        def test(rows, B):
            for i, j in itertools.combinations(range(len(B)), 2):
                k = len(B[i] & B[j])
                rows.append(dict(a=i, b=j, k=k, p_enr=float(hypergeom.sf(k - 1, N, len(B[i]), len(B[j]))),
                                 p_dep=float(hypergeom.cdf(k, N, len(B[i]), len(B[j])))))
            for key, qk in (("p_enr", "q_enr"), ("p_dep", "q_dep")):
                if rows:
                    for r, q in zip(rows, false_discovery_control([r[key] for r in rows], method="bh")):
                        r[qk] = float(q)
            return rows

        def top_letter(d, pos):
            cnt = d["letters"][d["positions"].index(pos)]
            return min(cnt.items(), key=lambda kv: (-kv[1], kv[0]))[0] if cnt else None

        def letters_similar(a, b):
            shared = [x for x in descs[a]["positions"] if x in descs[b]["positions"]]
            if not shared:
                return False, shared
            ok = 0
            for x in shared:
                ta, tb = top_letter(descs[a], x), top_letter(descs[b], x)
                ok += ta is not None and tb is not None and _B[ta, tb] >= 1
            return ok >= 1 and ok >= 0.5 * len(shared), shared
        epairs = test([], G)
        par = list(range(n))

        def find(x):
            while par[x] != x:
                par[x] = par[par[x]]; x = par[x]
            return x
        for r in epairs:
            sim, shared = letters_similar(r["a"], r["b"])
            r["shared_positions"], r["letters_similar"] = shared, sim
            r["link"] = (None if r["q_dep"] < Q else "carriers" if r["q_enr"] < Q else "letters" if sim else None)
            if r["link"]:
                par[find(r["a"])] = find(r["b"])
        comps = collections.defaultdict(list)
        for i in range(n):
            comps[find(i)].append(i)
        bset = {(r["a"], r["b"]) for r in bip}
        motifs = []
        for idx in sorted(comps.values(), key=min):
            carriers = sorted({p for i in idx for p in descs[i]["residues"]})
            links = sorted({r["link"] for r in epairs if r["link"] and r["a"] in idx and r["b"] in idx})
            motifs.append(dict(id=f"M{len(motifs) + 1}", elements=idx,
                               positions=sorted({x for i in idx for x in descs[i]["positions"]}, key=lambda c: (int(c), c)),
                               n_groups=len({grp[p] for p in carriers}), carriers=carriers, linked_by=links,
                               bipartite=any((a, b) in bset for a in idx for b in idx)))
        rel = test([], [{grp[p] for p in m["carriers"]} for m in motifs])
        for r in rel:
            shared = sorted(set(motifs[r["a"]]["positions"]) & set(motifs[r["b"]]["positions"]), key=lambda c: (int(c), c))
            r.update(a=motifs[r["a"]]["id"], b=motifs[r["b"]]["id"], shared_positions=shared,
                     relation="same_pocket" if shared else "separate_sites", exclusive=r["q_dep"] < Q,
                     status="different" if r["q_dep"] < Q else "unresolved")
        return motifs, rel, epairs

    # -- one interface
    def run_iface(self, iface, real, grp):
        P = self.P
        obs = self.obs_c(iface)
        if not obs:
            return None
        rng = self.seeded(self.seed_key, iface["id"], "C")
        nulls = {"shuffle": self.make_null(obs, real, rng)}
        cols = self.panel1(obs, real, nulls, grp)
        marks, rep = self.marks_and_repeats(cols, real, grp)
        ctrl_p1 = []
        for k_ in range(P["N_CTRL"]):
            _, ctrl_letters = self.make_null(obs, real, self.seeded(self.seed_key, iface["id"], "C", "ctrl", k_))
            cfull = {p: ctrl_letters.get(p, real[p]) for p in real}
            ccols = self.panel1(obs, cfull, {"shuffle": nulls["shuffle"]}, grp)
            ctrl_p1.append(sum(1 for c in ccols if c["q"].get("shuffle") is not None and c["q"]["shuffle"] < 0.05))
        ncand, sig = self.rules(obs, real, nulls["shuffle"][0], grp)
        prt = self.parts_of(sig); els = self.nonredundant(self.elements(prt, grp), grp)
        if P["EL_FILTER"]:
            els = [x for x in (self.trim(e, grp, real) for e in els) if x]
        descs = [self.describe(e, real, grp) for e in els]
        gp = {grp[p] for _, p, _, _ in obs}
        npairs, bip, alt = self.element_pairs(els, descs, grp, gp)
        motifs, mrel, epairs = self.motif_groups(descs, grp, gp, bip)
        clett = nulls["shuffle"][1] | {p: real[p] for p in real if p not in nulls["shuffle"][1]}
        _, csig = self.rules(obs, clett, nulls["shuffle"][0], grp)
        c_el = self.nonredundant(self.elements(self.parts_of(csig), grp), grp)
        if P["EL_FILTER"]:
            c_el = [x for x in (self.trim(e, grp, clett) for e in c_el) if x]
        panel1_out = []
        for c in cols:
            w = collections.Counter()
            for p, j, ww in c["obs"]:
                if j < len(real[p]) and real[p][j] >= 0:
                    w[self.LET(real[p][j])] += ww
            tot = sum(w.values())
            panel1_out.append(dict(pos=str(c["pos"]), touched=None, n_groups=c["n_groups"], agree=c["agree"],
                                   q=c["q"], letters={a: round(v / tot, 3) for a, v in w.items()} if tot else {}))
        return dict(panel1=panel1_out, marks=marks, repeats=rep,
                    panel2={"shuffle": dict(candidates=ncand, rules=len(sig), parts=len(prt), elements=descs,
                                       pairs=npairs, bipartite=bip, alternative=alt,
                                       motifs=motifs, motif_relations=mrel, element_links=epairs)},
                    significant=iface.get("significant", True),
                    control_panel1=ctrl_p1, control=dict(rules=len(csig), elements=len(c_el)),
                    partners_used=len({p for _, p, _, _ in obs}), partners_total=len(iface["partners"]))

    def run(self):
        all_p = sorted({p for i in self.ifaces for p in i["partners"]})
        grp = self.partner_groups(all_p)
        real = {p: self.letter_array(p) for p in all_p}
        return grp, {i["id"]: self.run_iface(i, real, grp) for i in self.ifaces}


def _plain(o):
    """numpy scalars/arrays -> plain Python (JSON-safe), recursively."""
    if isinstance(o, dict):
        return {(k if isinstance(k, str) else str(k) if not isinstance(k, (int, np.integer)) else int(k)): _plain(v)
                for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [_plain(x) for x in o]
    if isinstance(o, np.integer):
        return int(o)
    if isinstance(o, np.floating):
        return float(o)
    if isinstance(o, np.bool_):
        return bool(o)
    return o


def run_motif_engine(engine, seed_key, ifaces, pairs, int_seq, cif_of, bait, bait_cb, part_cb,
                     params=None, foldseek_bin=None):
    """Run one motif engine on a bait's confidence interfaces.

    seed_key  fixes the permutations (alphasift.screen passes the screen folder's name, as viz_int)
    ifaces    [{'id', 'residues', 'significant', 'partners': [prediction, ...]}] (partners = the
              interface's confidently binding partners, alphasift.partners.build_partner_annotations)
    pairs     {prediction: (bait_res, partner_res, pae)} (alphasift.matrix.contact_pairs)
    int_seq   {prediction: {residue: three-letter name}} (alphasift.structures.load_binder_structures)
    cif_of    {prediction: plain-text model CIF path} (for DSSP and Foldseek)
    params    ENGINE_DEFAULTS overrides; SS = "dssp" (default) or "none" (see the module docstring)
    bait_cb / part_cb  {prediction: {residue: Cb xyz}} (load_binder_structures(..., with_cb=True))
    Returns {'engine', 'params', 'groups': {prediction: group id}, 'ifaces': {id: result | None}}.
    """
    P = {**ENGINE_DEFAULTS, **(params or {})}
    if P["SS"] not in ("dssp", "none"):
        raise ValueError(f"SS must be 'dssp' or 'none', not {P['SS']!r}")
    if P["SS"] == "dssp" and shutil.which("mkdssp") is None:
        raise RuntimeError("the motif engine's secondary-structure null needs mkdssp (DSSP 4.x) on PATH, "
                           "e.g. `conda install -c conda-forge dssp`; or run without it (motif_ss='none', "
                           "CLI --motif-ss none), which is not the validated null")
    used =sorted({p for i in ifaces for p in i["partners"] if pairs.get(p) is not None})
    ifaces = [dict(i, partners=[p for p in i["partners"] if pairs.get(p) is not None]) for i in ifaces]
    seq = {p: {r: aa1(n) for r, n in int_seq[p].items()} for p in used}
    ss = {p: dssp_classes(cif_of[p], bait["chain_ids"][p]) for p in used} if P["SS"] == "dssp" else {}
    t3 = partner_3di(cif_of, bait, used, foldseek_bin) if engine == "3di_aa" else None
    pr = {p: tuple(np.asarray(x) for x in pairs[p]) for p in used}
    eng = _Engine(engine, seed_key, ifaces, pr, seq, ss, bait_cb, part_cb, t3, P)
    grp, res = eng.run()
    return _plain(dict(engine=engine, params=P, groups=grp, ifaces=res))
