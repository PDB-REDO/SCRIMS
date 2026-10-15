# SCRIMS — statistical inference of functional interaction sites and shared motifs from structure-prediction screens

SCRIMS turns a set of AlphaFold predictions that share one **bait** protein into
three things:

1. a **confidence matrix** — one row per partner prediction, one column per bait
   residue, valued by how confident AlphaFold is about the contact there;
2. the **interaction sites** on that bait that the screen statistically supports,
   merged into 3-D interfaces; and
3. per interface, the **partners** that bind it confidently, with the partner
   residues involved, and the **sequence motifs** unrelated partners share there.

It answers "which parts of this protein do many partners confidently engage, who
binds there, and do they bind the same way?" — from predicted aligned error, not
from distance alone.

## What it computes

**The matrix.** Per bait residue *i* and prediction *p*, geometry only *selects*
candidate pairs — Cβ–Cβ ≤ `d0_contact` (8 Å) to any partner residue. The value is
then read from PAE alone: the direction-symmetrised **minimum PAE** over those
selected partner residues. Symmetrising is not cosmetic; using one PAE direction
collapses the signal on real data. That raw min-PAE (Å, undefined where the
prediction has no contact) is what the figure displays. For the statistics it
passes through the pTM-style pair term

    signal = 1 / (1 + (PAE / d0)²),    d0 = 10 Å, fixed

which flips direction (lower PAE → higher signal) and bounds to (0, 1], so `0`
means *no contact* and `0.5` means *exactly PAE = 10 Å*. `d0` is **one constant
for every screen and every partner**, deliberately: a per-pair, length-dependent
`d0` makes larger partners score higher for worse PAE. The value 10 Å follows
pDockQ2 (Zhu et al., *Bioinformatics* 2023).

**Recurrence — a test.** Sites where many partners converge. The null shifts each
partner's *whole row* by a random circular offset, randomising only each
partner's position relative to the others, which is precisely what recurrence
claims. Clusters are formed above the `pct`-th percentile of that null, each
observed cluster is compared against null clusters of **similar extent** (w/1.5
to 1.5w; when fewer than 50 fall in that band it widens ×1.25 per step until 50
do), and the p-values are BH-corrected. A site is `recurrent` at q < 0.05.

**Strength — a criterion, not a test.** A partner qualifies when *it alone* holds
at least `min_single` (4) residues at signal ≥ `tau_single` (0.80 ≈ PAE 5 Å)
anywhere on the bait; every sequence piece of it, allowing gaps of up to
`single_gap` (2), is then a strength region, whatever its length
(`--strength-scope group` restores the earlier rule: four residues per piece).
**It carries no p-value, by design.** A
value claim ("the model is very confident here") has no null in a matrix that
only supports positional nulls; five candidate tests were built and all failed
calibration. Report it as a declared criterion.

**The 3-D merge.** Two sites far apart in sequence can be one surface in space.
For each pair of sites, `support` is the fraction of partner models in which
their bait heavy atoms come within `d_merge` (8 Å). Two sites **link** when they
touch (support ≥ `merge_support`, 0.5), at least `merge_min_shared` (2) partners
bind both (signal ≥ `partner_tau_in`, 0.5, on a site residue), and the shared
partners are at least `merge_share_small` (0.7) of the smaller site's partners
and `merge_share_large` (0.2) of the larger's; interfaces are the connected
components. A site left alone that is short (< `attach_min_len`, 4 residues) or
weakly bound (< 2 partners at 0.5) first groups with other such sites it touches
when they share one of their own partners (signal ≥ `attach_tau_strong`, 0.8, on
it if any, else ≥ `attach_tau_weak`, 0.3); each group then joins one interface it
touches whose site one of its partners binds, or stands alone. A site no partner
binds at 0.5 groups and joins on closeness alone (support ≥ `merge_prox_support`,
0.9). No partner relatedness is used. `--merge-method average` restores the
earlier rule (average linkage on 1 − support). Every site reports the interface
it joined.

**Partners at each interface — a description.** Every partner with at least one
contact pair into the interface at signal ≥ `partner_tau_in` (0.5, PAE 10 Å) is
listed, with each of its residues in contact tiered by its best pair: **H** at
≥ `partner_tau_high` (0.8, PAE 5 Å), **M** at ≥ 0.5, **L** below, in exact runs of
one tier. The pairs are the matrix's own (Cβ–Cβ ≤ 8 Å, symmetrised PAE), seen from
the partner chain. Worked example at d0 = 10 Å: PAE 4 Å → s = 1/(1 + 0.16) = 0.86 → H;
PAE 7 Å → s = 0.67 → M; PAE 12 Å → s = 0.41 → L. The tier says how far to trust a
residue (gold-set precision H 72 %, M 37 %, L 14 %); it is not a test. These always
come from PAE, whichever `--source` built the matrix.

**Shared motifs — a test (the motif engine).** At every interface, significant or
not: do *unrelated* partners put *similar* residues on the *same* bait residues more
often than chance? A partner residue counts where it touches the bait (pair at
s ≥ 0.5 and Cβ–Cβ ≤ `motif_touch_cb`, 7 Å); partners ≥ 40 % identical form one group,
so paralogs count once. Residues are compared by BLOSUM62 (`--motif-engine c`, the
default) or by amino acid plus the partner's Foldseek 3Di state at weight 1/3
(`3di_aa`). The null shuffles the letters placed on the interface, 1000 times; a
shuffle within 1e-9 of the observed score counts as a tie, and p = (k + 1)/(n + 1),
BH-corrected. By default (`--motif-ss dssp`) letters are shuffled only **within
DSSP class** (helix / strand / other), so "every partner binds with a helix" is not
a motif; this is the default null and needs `mkdssp`. `--motif-ss none` shuffles
across the whole interface instead: no external program and several times faster.
Both were benchmarked in viz_int on 31 baits (VALIDATION_LOG 2026-10-06, "Motif null
WITHOUT secondary structure"): engine `c` with DSSP grouped known carriers better than
random for 12 of 13 known motif entries (8 of 9 on held-out baits), 1 false element in
shuffled controls, gold contact precision 62 %, known vs background 2.9×; without DSSP
11 of 13 (7 of 9), 0 false elements, 62 → 60 %, 3.0×. The entry lost (a MAPK8 D-site
with 3 carriers) is within what one partner flips, so the benchmark does not separate
the two nulls; DSSP stays the default because the published numbers rest on it. The **bait side** tests, per bait residue touched
by ≥ 3 groups, how much the groups agree. The **partner side** tests bait positions
filled at a constant spacing (≤ 10 residues) along the partner and assembles the
significant ones into *elements*; a declared filter then keeps in each element only
the partners that carry it. Elements carried by the same groups beyond chance, or
preferring similar residues on shared bait residues, are **one motif**; two motifs
are "different" only when their carriers overlap *less* than chance, otherwise
"unresolved". Validated in viz_int with the DSSP null: the default engine grouped
known carriers better than random for 12 of 13 known motif entries (8 of 9 on baits
added after the design), with 1 false element in shuffled controls; `3di_aa` 8 of 13,
with the higher contact precision (66 % vs 62 %). These numbers apply to
`--motif-ss dssp` only. `3di_aa` also needs Foldseek.

## Install

    pip install scrims                    # numpy, scipy, pandas, matplotlib, biopython
    pip install "scrims[zst]"             # + zstandard, for .zst AFDB downloads

Two programs are not on PyPI. **`mkdssp`** (DSSP 4.x) is needed by the default,
validated motif null (`conda install -c conda-forge dssp`); without it `scrims run`
stops and says so — pass `--motif-ss none` to run without secondary structure, or
`--motif-engine none` to skip motifs. **Foldseek** is needed only for
`--motif-engine 3di_aa` (`conda install -c bioconda foldseek`, or point
`--foldseek-bin` / `$SCRIMS_FOLDSEEK` at a binary).

### Or use the provided environment

If you would rather not touch your own environment, `environment.yml` builds a
complete one, test runner included:

    conda env create -f environment.yml
    conda activate scrims
    pip install -e .

### Running the tests

    pytest

The suite needs no prediction data: it builds a synthetic six-partner screen in a
temporary directory and checks the parts that could fail silently — that 0.5 of
the signal is exactly PAE = 10 Å, that 0 means *no contact* rather than weak
contact, that a planted convergent site is called recurrent, that a screen where
every partner binds somewhere different calls **nothing** recurrent, that the
strength criterion rejects three partners covering two residues each but accepts
one partner covering four, and counts the partner rather than each piece, that a
wide cluster on a screen with no convergence is not called, that the merge links
touching sites sharing partners but not when the sharing is one-sided, that short
sites attach, group or stand alone, that the bait is found without AlphaBridge,
that a single pair refuses to report a q-value, that the AlphaBridge source runs
end to end, and that the figure renders. For the partner and motif stages it checks
that the partner pairs reduce exactly to the matrix values, the H/M/L tiers, that
the motif engine runs on every interface, in both null modes, that DSSP is the
default and a missing `mkdssp` is a clear error, that `--motif-ss none` needs no
external program, the motif-logo alignment rules, the classic logo's exact small-sample
correction and gap-free letter stacks (and that the report's script computes the same), that partner names follow the
sequence rather than the file header's chain link, and that every output file is
written. Tests of the DSSP mode are skipped where `mkdssp` is not installed.

`pytest --run-slow` adds an integration test on three real AlphaFold3 screens (ezrin,
moesin, radixin; 61, 41 and 24 predictions) in `data/`, which is not in git: each runs
with `--source pae` on plain AlphaFold3 output, with `--source pae` and AlphaBridge files
present (results must be identical), and with `--source alphabridge`. It takes about
6 minutes and is skipped when `data/` is absent.

Note that `scipy >= 1.11` is required: the recurrence FDR uses
`scipy.stats.false_discovery_control`, which older versions do not have.

## Use

    scrims run path/to/screen -o out/

writes:

| file | what |
|---|---|
| `scrims.json` | everything below in one versioned document (schema `scrims/1`, described in `scrims/export.py`) |
| `report.html` | one self-contained page mirroring the companion website's protein page — matrix, interface cards, motif views (bait side and per-motif partner logos, each as a classic logo in bits or an evidence view; a per-carrier table; a motif lane filters the bait side to that motif), partners, all binders — with linked selection; opens offline (no 3-D: use the ChimeraX script) |
| `overview.png` | the matrix with RMSF and recurrence tracks and one interface lane |
| `interfaces/<I>.png` | per interface: bait-side letter columns (evidence view), motif lanes, and per motif its partner-sequence logo (classic, bits) over the per-carrier table |
| `interfaces/<I>_partners.png` | per interface: every partner with its H/M/L residue ranges, motif residues in bold |
| `scrims.cxc` (+ `models/`) | ChimeraX script: interfaces coloured, motif residues in cyan, every interface and motif a named selection |
| `matrix.png` | the matrix with the summed-signal track and the called sites |
| `sites.tsv`, `interfaces.tsv` | called sites; 3-D interfaces |
| `partners.tsv`, `bait_side.tsv`, `motifs.tsv` | partner ranges; bait-side columns; one row per motif |
| `matrix.npz`, `summary.json` | the raw matrices; bait, parameters, counts, skipped predictions |

The figures and the JSON follow the views of the companion website (same colours,
row order, letter colours and motif-logo alignment). Every parameter above is a flag
(`--pct`, `--d0-signal`, `--motif-engine`, `--motif-ss`, …). A 216-prediction screen
(RHOQ) takes about 6 minutes on one core with the DSSP null, of which mkdssp is about
3.7 (75 models × ~3 s).

```python
from scrims import analyze
res = analyze("path/to/screen")
res.M_pae        # partners × residues, raw min-PAE (Å); NaN = no contact
res.M_signal     # the same, transformed to (0, 1]; 0 = no contact
res.sites        # significance, channel, cores, interface
res.interfaces   # 3-D merged interfaces
res.partners     # interface -> partners with tiered residue ranges
res.motifs       # the motif engine's output, every interface
```

## Input formats

A **screen** is a directory of predictions of one bait against many partners. The
bait is identified by exact sequence identity — the chain whose full sequence
appears in *every* prediction. Chain letter and length are never used.

### Two confidence sources

`--source pae` (the default) works for **any** predictor and runs straight off
the prediction files. `--source alphabridge` is available for **AlphaFold3
only**, and takes the highest cut-off at which each bait residue survives
AlphaBridge's interface sweep — a signal that exists because AF3 carries a
distogram, which AlphaFold2 output does not.

It is not free. Using it means **installing AlphaBridge and running it over
every prediction in the screen first**, under AlphaBridge's own licence. A
*partial* run is worse than none: recurrence asks whether partners converge on
the same residues, so predictions missing from the matrix change the test rather
than merely shrinking it. SCRIMS lists every excluded prediction in
`summary.json` and warns when any are missing.

The two sources share a 0–1 scale but not its meaning: on the PAE channel 0.8 is
PAE 5 Å, while on the AlphaBridge channel it means the residue survives to
cut-off 0.8. Threshold parameters — including the merge's partner thresholds —
keep their numeric defaults across sources but not their interpretation.

| Predictor | Layout | Status |
|---|---|---|
| AlphaFold3 (server or local) | `<pred>/<pred>_model.cif` + `<pred>_confidences.json` | supported |
| AlphaFold2 / ColabFold, AFDB bulk packaging | `<id>-model_v1.cif[.zst]` + `<id>-predicted_aligned_error_v1.json[.zst]` | supported (the validation data: the NVIDIA heterodimer set from EMBL-EBI, AlphaFold-Multimer via ColabFold 1.6.0 — not AlphaFold3) |
| AlphaBridge run (AF3) | `<pred>/AlphaBridge/alphabridge_data.json` | supported as `--source alphabridge`, and as a PAE shortcut where a `confidence_heatmap` is stored and the raw PAE is not |
| AlphaFold2 official (`run_alphafold.py`) | `ranked_0.pdb` + `result_model_*.pkl` + `ranking_debug.json` | not yet |
| ColabFold local (`colabfold_batch`) | `*_unrelaxed_rank_001_*.pdb` + `*_scores_rank_001_*.json` | not yet |

If a screen has two shared sequences, or you are analysing a single pair (where
every chain is trivially "common"), name the bait yourself:

    scrims run path/to/pair -o out/ --bait-chain A

### Running AlphaBridge first (only for `--source alphabridge`)

AlphaBridge is a **separate tool**, not a dependency of SCRIMS and not
installed with it: <https://github.com/PDB-REDO/AlphaBridge> (web server:
<https://alpha-bridge.eu/>). It carries its own licence terms and its own
citation requirement — read them before use.

    git clone https://github.com/PDB-REDO/AlphaBridge
    conda env create -f AlphaBridge/environment.yml --name alphabridge
    conda activate alphabridge

It processes one AlphaFold3 prediction folder at a time, writing
`AlphaBridge/alphabridge_data.json` inside it:

    python3 define_interfaces.py -i /path/to/screen/<prediction>

SCRIMS needs that file for **every** prediction in the screen, so run it across
the whole directory before analysing:

    for d in /path/to/screen/*/; do python3 define_interfaces.py -i "$d"; done
    scrims run /path/to/screen -o out/ --source alphabridge

Predictions AlphaBridge fails on are listed in `summary.json` and excluded from
the matrix, and SCRIMS warns when any are — see the note above on why a partial
run changes the test rather than merely shrinking it.

**AlphaBridge is not required.** It is one optional source, not a dependency:
with no AlphaBridge output present, sequences and PAE are read from the
prediction files directly.

## Deliberately not here

From the original pipeline: the legacy motif panels (a consensus logo per pocket and
one anchored linear motif, replaced by the motif engine), the AlphaBridge-link-based
contacts, peptides and registers, the geometry-only channel, and the website's own
data export. Partner names come from the model files' headers (UniProt accession,
entry name, gene, organism) — no network lookups.

## Scope limit worth stating plainly

Recurrence needs a screen: it asks whether partners converge *relative to each
other*. A single pair cannot answer that, so with one prediction SCRIMS reports
the matrix and strength regions and explicitly declines to report a q-value.

## Provenance

The matrix construction, the statistical engine, the 3-D merge, the partner
annotation, the motif engine and the RMSF track are extracted **verbatim** from the
`viz_int` research pipeline (re-synced 2026-10-06), so its validation applies to this
code unchanged. Checked: on four screens (PARVA, RHOQ, VAPA, XPO1 — two layouts,
88 to 441 predictions) matrix, sites, interfaces, partners and the full motif-engine
output are identical to viz_int's, floats included. The motif-logo alignment is a
port of the companion website's own code and matches it on all 51 motifs of its 31
proteins. What is new here is the input dispatch layer, bait detection without
AlphaBridge, partner names read by sequence from the model headers, the figures,
the JSON and table outputs, the ChimeraX script, and the CLI.

## Licence and citation

**Licence: not yet fixed — see the note below before distributing this.** The
intended terms are free availability for non-commercial users, and continued
availability for at least two years after publication, as the journal requires.

**Relationship to AlphaBridge.** SCRIMS contains no AlphaBridge code, does not
import it, and does not require it to run; AlphaBridge output is one optional
input format among several. AlphaBridge is distributed by its authors under its
own, separate agreement, which restricts redistribution and derivative works and
requires that results obtained through its use cite its publication. If you
produced your screens with AlphaBridge, those terms apply to that use — cite it
accordingly.

Citation for SCRIMS: the Application Note and Zenodo DOI, once published.
