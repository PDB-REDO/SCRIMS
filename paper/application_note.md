# CONFINE: confidence matrices and statistically supported interaction sites from AlphaFold screens

<!-- TARGET: Bioinformatics, Application Note (≤4 pages ≈ 2,600 words, or 2,000 + one figure).
     Body below is ~2,050 words excluding abstract, caption, references.
     TODO before submission: author list, affiliations, corresponding author,
     funding, the AlphaBridge citation (ref. 12 — fill in the preprint details),
     the PDBe-KB citation (ref. 13 — verify), repository URL and Zenodo DOI. -->

**Daniel Álvarez-Salmoral**<sup>1,\*</sup>, *[co-authors TBC]*

<sup>1</sup> *[affiliation TBC]*

<sup>\*</sup> To whom correspondence should be addressed.

---

## Abstract

**Summary:** AlphaFold has made it practical to predict one protein against
hundreds of candidate partners, but the scores that come out of such a screen —
ipTM, pDockQ, ipSAE — rank whole complexes. They do not say *which part* of the
bait protein is engaged, nor which of the regions that appear engaged are
supported by more than one prediction. CONFINE fills that gap. It reduces a
screen to a confidence matrix, one row per partner and one column per bait
residue, valued by AlphaFold's predicted aligned error rather than by
interatomic distance, and it calls the interaction sites on that matrix under an
explicitly calibrated null. Two distinct claims are kept separate: *recurrence*,
that many partners converge on a region more than chance allows, which is tested
and FDR-corrected; and *strength*, that a single partner engages a region with
high confidence, which is reported as a declared criterion because no valid null
for it exists in this data. Sites are then merged in three dimensions, so
regions distant in sequence but adjacent in space are reported as one interface.
CONFINE reads AlphaFold2 and AlphaFold3 output directly and requires no upstream
interface-analysis step.

**Availability and implementation:** Source code, documentation and test suite
are freely available for non-commercial use at *[repository URL]* and archived
at *[Zenodo DOI]*. CONFINE is pure Python (≥3.10), depends only on NumPy, SciPy,
pandas and Matplotlib, installs with `pip`, and requires no network access,
compilation or registration.

**Contact:** *[email]*

**Supplementary information:** Supplementary data are available at *Bioinformatics*
online.

---

## 1 Introduction

Structure prediction has changed how protein–protein interactions are surveyed.
With AlphaFold2 and AlphaFold-Multimer (Jumper *et al.*, 2021; Evans *et al.*,
2022), ColabFold (Mirdita *et al.*, 2022) and AlphaFold3 (Abramson *et al.*,
2024), predicting one **bait** protein against hundreds of candidate partners is
now routine. Interpreting the result is not.

The established confidence measures score a *complex*. ipTM summarises a whole
interface in one number; pDockQ and pDockQ2 (Zhu *et al.*, 2023) estimate
docking quality per model; ipSAE (Dunbrack, 2025) corrects ipTM's well-described
length bias. All of them answer "is this pair likely to interact?" None answers
the question a screen actually poses: **which surface of the bait is being
engaged, and which of the regions that appear engaged would survive a null
model?**

That second question cannot be answered one prediction at a time. A region that
one partner touches may be an artefact; the same region touched by forty
partners, more often than random placement would produce, is a claim about the
protein rather than about a model. Making that claim requires (i) a per-residue
view of the whole screen and (ii) a null against which convergence can be
measured.

A natural first attempt is to build the per-residue view from distance: mark the
bait residues within some cutoff of a partner. We took that route initially and
abandoned it. Distance answers "how close did the model place these chains",
which is not the same as "how sure is the model that this placement is correct".
The two diverge exactly where it matters — in flexible regions and low-confidence
tails, where models place chains in contact with no confidence at all.

CONFINE therefore reads the per-residue signal from predicted aligned error
(PAE), uses geometry only to decide which residue pairs are *candidates*, and
supplies a calibrated statistic over the resulting matrix. It runs on AlphaFold2
and AlphaFold3 output directly.

## 2 Implementation

### 2.1 Input

A **screen** is a directory of predictions of one bait against many partners.
The bait is identified by exact sequence identity: the chain whose full sequence
occurs in every prediction. Chain letter and chain length are never used as
criteria, because neither is stable across predictions. Where that rule cannot
resolve a bait — a single pair, or a screen with two shared sequences — the
chain can be named explicitly.

CONFINE reads AlphaFold3 output (mmCIF plus `confidences.json`) and AlphaFold2
output in the AFDB/ColabFold bulk packaging (`-model_v1.cif[.zst]` plus
`-predicted_aligned_error_v1.json[.zst]`). Output from AlphaBridge
(Álvarez-Salmoral *et al.*), which carries a precomputed per-residue reduction,
can be used where the raw PAE matrices have not been retained. **AlphaBridge is
an optional input format, not a dependency**: with none present, sequences and
PAE are read from the prediction files directly. Adding a predictor means adding
a reader; nothing downstream changes.

### 2.2 The confidence matrix

For bait residue *i* in prediction *p*, geometry only **selects** candidate
pairs: partner residues *j* whose Cβ lies within 8 Å of residue *i*'s Cβ. The
value is then read from PAE alone, as the direction-symmetrised minimum,
min over *j* of min(PAE(*i*,*j*), PAE(*j*,*i*)). Symmetrising is not cosmetic;
using a single PAE direction collapses the signal on real data.

That raw minimum PAE, in Ångström and undefined where the prediction has no
contact at residue *i*, is what the figure displays. For the statistics it
passes through the pTM-style pair term

*s* = 1 / (1 + (PAE / *d*<sub>0</sub>)<sup>2</sup>),  *d*<sub>0</sub> = 10 Å,

which reverses direction — lower error, higher signal — and bounds the result to
(0, 1], so that **0 means no contact** rather than a weak one. The distinction
matters: the statistics below assume absence is encoded as zero.

*d*<sub>0</sub> is **one constant for every screen and every partner**. This is a
deliberate departure from the length-dependent *d*<sub>0</sub>(*L*) of TM-score
(Zhang and Skolnick, 2004), which ipTM inherits and which Dunbrack (2025)
identifies as the source of its bias: adding non-interacting sequence enlarges
*d*<sub>0</sub>, scaling PAE down and inflating the score. We reproduced that
artefact on our own data. Across 24 bait proteins, larger partners had
*worse* raw PAE (median Spearman ρ = +0.28, positive in 23 of 24 screens), yet a
length-dependent *d*<sub>0</sub> gave them *higher* signal (median ρ = +0.32,
positive in 24 of 24). A fixed *d*<sub>0</sub> removes the inversion. The value
10 Å follows pDockQ2, which fits the same pair term against DockQ, and is
interpretable: signal 0.5 is exactly PAE 10 Å, and 0.8 is PAE 5 Å.

### 2.3 Recurrence: a test

Recurrence asks whether partners converge on the same residues more than chance
allows. The null shifts **each partner's entire row** by an independent random
circular offset. This randomises each partner's position *relative to the
others* — precisely the claim being tested — while preserving each partner's own
footprint shape, its number of segments and their spacing. An earlier null that
relocated each contiguous segment independently was anti-conservative: real
partners bind compactly, scattered fragments overlap less often, and the
observed data cleared the bar too easily.

Clusters are formed where the summed signal exceeds the 95th percentile of the
null, and each observed cluster is compared against null clusters of **similar
extent** (±1.5×, falling back to the full pool below 50 members) rather than
against clusters of every width, because intact rows produce wide, heavy null
clusters that would otherwise mask narrow real sites. P-values are corrected by
Benjamini–Hochberg (1995); sites are called recurrent at *q* < 0.05.

### 2.4 Strength: a criterion, not a test

A region qualifies as a strength site when **one** partner holds at least four
residues at signal ≥ 0.80 (PAE ≤ 5 Å), allowing gaps of up to two positions.
Regions are formed and judged on the same evidence, so every strength region
qualifies by construction.

It carries **no p-value, and this is a result rather than an omission**. A value
claim — "the model is very confident here" — has no null in a matrix that
supports only positional nulls. Five candidate tests were built and all failed
calibration: a pooled-window null makes almost no calls, because 54 % of all
four-residue windows in the corpus carry some confident engagement; a
per-partner null fired in 63 % of synthetic negative screens; and once corrected
for every partner and position it makes zero calls anywhere, real or synthetic.
Testing a value claim would require control predictions of pairs known not to
interact. Until such controls exist, CONFINE reports the criterion and labels it
as one.

### 2.5 Merging sites into interfaces

Two sites far apart in sequence may form one surface. For each pair of sites,
*support* is the fraction of partner models in which their bait heavy atoms come
within 8 Å; sites are clustered on 1 − support and cut at 0.5. Each site reports
the interface it joined.

### 2.6 Output and performance

`confine run SCREEN -o OUT` writes the matrix figure, a table of sites (start,
end, interface, channel, significance, *q*, the significant clusters and
confident footprints inside each site, and the number of confident partners), a
table of interfaces, the matrix itself, and a run summary recording every
prediction that was skipped **with its reason**. A Python API returns the same
objects. On one core of an Intel i9-13900K, a screen of 216 predictions against
a 205-residue bait completes in 21.6 s including 5,000 permutations; the cost is
dominated by parsing prediction files, not by the statistics.

## 3 Validation

**Calibration.** Synthetic negative screens were built with the same bait length
and partner count as a real screen, but with rows taken from *other* proteins and
placed at random offsets, so recurrence's null hypothesis holds by construction
and any call is a false positive. The whole-row null with extent matching calls
15 clusters across 15 of 120 such screens (12 %), against 45 clusters in 17 of 24
under the earlier segment-scattering null. This is honestly still above the 5 %
that *q* < 0.05 implies under a complete null, and we report it as such.

**Agreement with known interfaces.** Called residues were compared against
interface residues annotated in experimentally determined structures (PDBe-KB),
restricted to heterotypic interfaces, over 19 baits and 2,006 annotated
residues. Recurrence calls 1,513 residues, of which 47 % are annotated interface
residues — 2.57× the density expected from the protein's own interface coverage
— at 42 % recall, while calling sites in only 3 of 24 synthetic negatives.
Including strength sites raises recall to 63 % at 45 % precision, with the two
channels recovering largely disjoint residue sets.<!-- NOTE: this combined
figure was measured at the previous cluster-forming percentile (90), not the
shipped default (95). Re-measure at 95 before submission. -->

**A solved-complex gold set.** The sharpest test uses 68 bait–partner pairs
whose exact complex is in the PDB and whose partner is in the screen, giving
true interface residues on both sides. **71.6 % of true bait interface residues
fall inside a called site, 70.4 % inside a significant one.** The corresponding
precision figure, 12.6 %, must be read as a **floor**: a bait residue called for
a different partner than the gold one is counted as a miss although it may be
perfectly correct. Six of the gold cases score zero — there the predicted
binding mode does not overlap the crystal interface at all, an error of
prediction that no downstream statistic can repair.

## 4 Conclusion

CONFINE turns an AlphaFold screen into two artefacts a structural biologist can
act on: a figure showing which parts of a bait its candidate partners engage and
how confidently, and a table of the regions that survive an explicitly
calibrated null. Its design commitments are that confidence and distance are
different quantities, that a fixed *d*<sub>0</sub> is preferable to one that
rewards large partners, and that a claim without a valid null should be labelled
a criterion rather than dressed as a test. It reads AlphaFold2 and AlphaFold3
output directly, adds no required upstream stage, and installs in one command.

---

**Figure 1.** *CONFINE output for a 216-partner screen against RHOQ (205
residues).* Top: the confidence matrix, one row per partner prediction, ordered
by profile similarity. Colour is the confidence signal, 1/(1 + (PAE/10)²); grey
marks residues with no contact in that prediction, which is a categorical state
rather than a low value. Bars beneath the matrix mark called sites — solid for
recurrent (*q* < 0.05), dashed for strength. Bottom: summed signal per bait
residue, the quantity the recurrence null acts on, with called sites shaded. The
two recurrent sites, 27–56 (*q* = 2.7 × 10⁻⁴, 26 confident partners) and 58–80,
merge into a single three-dimensional interface together with the sequence-distant
site 166–184.

---

## References

1. Jumper, J. *et al.* (2021) Highly accurate protein structure prediction with
   AlphaFold. *Nature*, **596**, 583–589.
2. Evans, R. *et al.* (2022) Protein complex prediction with AlphaFold-Multimer.
   *bioRxiv*, 2021.10.04.463034.
3. Mirdita, M. *et al.* (2022) ColabFold: making protein folding accessible to
   all. *Nat. Methods*, **19**, 679–682.
4. Abramson, J. *et al.* (2024) Accurate structure prediction of biomolecular
   interactions with AlphaFold 3. *Nature*, **630**, 493–500.
5. Zhu, W. *et al.* (2023) Evaluation of AlphaFold-Multimer prediction on
   multi-chain protein complexes. *Bioinformatics*, **39**, btad424.
6. Dunbrack, R.L. Jr. (2025) Rēs ipSAE loquuntur: what's wrong with AlphaFold's
   ipTM score and how to fix it. *bioRxiv*, 2025.02.10.637595.
7. Zhang, Y. and Skolnick, J. (2004) Scoring function for automated assessment
   of protein structure template quality. *Proteins*, **57**, 702–710.
8. Benjamini, Y. and Hochberg, Y. (1995) Controlling the false discovery rate.
   *J. R. Stat. Soc. B*, **57**, 289–300.
9. Álvarez-Salmoral, D. *et al.* AlphaBridge. *[citation to be completed]*
10. PDBe-KB consortium. *[citation to be verified]*
