# CONFINE: confidence matrices and statistically supported interaction sites from AlphaFold screens

<!-- TARGET: Bioinformatics, Application Note (≤4 pages ≈ 2,600 words, or 2,000 + one figure).
     Body is 2,330 words + 194-word abstract + 625 words of Supplementary
     Methods (S1-S4), excluding caption and references. Still ~300 over the
     2,000-plus-one-figure guidance; remaining levers are Fig. 1 to
     supplementary, or further cuts to 2.1/2.2.
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
length bias. All answer "is this pair likely to interact?" None answers what a
screen actually poses: **which surface of the bait is engaged, and which of the
regions that appear engaged would survive a null model?** That cannot be settled
one prediction at a time. A region one partner touches may be an artefact; the
same region touched by forty partners, more often than random placement would
produce, is a claim about the protein rather than about a model.

Building the per-residue view from distance is the obvious first attempt — mark
the bait residues within some cutoff of a partner — and we abandoned it.
Distance answers how close the model placed two chains, not how sure it is that
the placement is right, and the two diverge exactly where it matters, in
flexible regions and low-confidence tails where chains are placed in contact
with no confidence at all. CONFINE therefore reads the per-residue signal from
predicted aligned error (PAE), uses geometry only to decide which residue pairs
are *candidates*, and supplies a calibrated statistic over the resulting matrix.

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

*d*<sub>0</sub> is **one constant for every screen and every partner**, a
deliberate departure from the length-dependent *d*<sub>0</sub>(*L*) of TM-score
(Zhang and Skolnick, 2004) that ipTM inherits and that Dunbrack (2025)
identifies as the source of its bias: adding non-interacting sequence enlarges
*d*<sub>0</sub>, scaling PAE down and inflating the score. We reproduced the
artefact across 24 bait proteins — larger partners had *worse* raw PAE yet a
length-dependent *d*<sub>0</sub> gave them *higher* signal, in 24 of 24 screens
(Supplementary Methods S1). A fixed *d*<sub>0</sub> removes the inversion. The
value 10 Å follows pDockQ2, which fits the same pair term against DockQ, and is
interpretable: signal 0.5 is exactly PAE 10 Å, and 0.8 is PAE 5 Å.

### 2.3 Recurrence: a test

Recurrence asks whether partners converge on the same bait residues more than
chance allows. Let *M* be the signal matrix and *s*(*j*) = Σ<sub>*i*</sub>
*M*(*i*,*j*) the summed signal at residue *j*. The procedure is:

1. **Null.** In each of *P* = 5,000 realisations, every partner row is shifted by
   its own independent, uniformly drawn **circular** offset and the shifted rows
   are summed, giving a null profile of the bait's length.
2. **Cluster-forming threshold.** *t* is the 95th percentile of the pooled null
   values. Clusters are maximal runs above *t*, tolerating single-residue gaps,
   of at least four residues — applied identically to observed and null profiles.
3. **Cluster mass.** A cluster's statistic is Σ *s*(*j*) over its residues, so
   extent and height both contribute, rather than peak height alone.
4. **Extent-matched p-value.** Every cluster of every null realisation enters a
   pool as a (mass, width) pair. An observed cluster of width *w* is compared
   only against null clusters of width *w*/1.5 to 1.5*w*, falling back to the
   whole pool when fewer than 50 qualify, giving *p* = (*n*<sub>≥</sub> + 1)/(*n* + 1).
5. **Correction.** Benjamini–Hochberg (1995) across observed clusters; recurrent
   at *q* < 0.05.

Two choices carry the weight, and both were settled by measurement.

**Whole rows, not segments.** A null must randomise exactly what the alternative
hypothesis asserts and nothing else. Recurrence asserts that partners converge
*relative to one another*; it says nothing about the internal shape of any one
partner's footprint. A circular shift of the whole row leaves that shape
untouched — segment count, widths, spacing and the multiset of values are all
preserved — and randomises only position. Relocating each contiguous segment
independently, as an earlier version did, destroys the wrong thing: real
partners engaging one region contribute several nearby segments, while scattered
fragments overlap between partners less often, so the null under-produces
clusters and observed data clears the bar too easily. That null produced 45
spurious clusters in 17 of 24 synthetic negative screens. The shift wraps around
the sequence terminus; the non-wrapping alternative proved over-conservative,
discarding calls that annotated interfaces support (Supplementary Methods S1).

**Comparing like with like.** With rows intact, null clusters are wider and
heavier, so scoring a narrow observed cluster against null clusters of every
width penalises precisely the sites of interest. Extent matching is what makes
narrow sites detectable: on a 16-site regression set, 9 of 16 survived without
it and 16 of 16 with it. The cluster-forming percentile, in contrast, is not a
sensitive knob — swept from 80 to 95 it changes the number of calls by under 5 %
(Supplementary Methods S1), which is why the conservative end is the default.

### 2.4 Strength: a criterion, not a test

Recurrence is blind to a real and common case: a single partner engaging a
region with high confidence. That is a claim about one prediction rather than
about convergence, and CONFINE reports it separately.

**Formation.** For each partner independently, take its residues at signal
≥ 0.80 (PAE ≤ 5 Å); group residues separated by at most two non-qualifying
positions; keep groups of at least four residues; merge overlapping groups
across partners. Forming and judging use the same evidence, so every strength
region qualifies by construction.

Each element was forced by the failure of a simpler one. Judging on a
per-residue maximum cannot separate a region held by 27 partners at PAE 1.2 Å
from one held by a single partner at 6.2 Å — under that rule 179 of 298 regions
qualified. Forming regions by a loose union across partners and then judging
them strictly on one partner guarantees candidates that can never pass; 63 such
regions existed in an earlier release. Requiring four *consecutive* residues
loses most real sites, because confident contacts are rarely contiguous: RHOQ
164–180 has nine partners with at least four confident residues and none with
four in a row. Parameter sweeps for the threshold, the minimum size and the gap
are in Supplementary Methods S2.

CONFINE also reports, without gating on it, the total confident engagement per
region — the summed Ångström margin below the contact bar over all partners —
which rises with both the number of engaged partners and their confidence, and
so separates the two cases a maximum conflates.

**Why there is no p-value, and why that is a result.** A value claim — "the
model is very confident here" — has no null in a matrix that supports only
positional nulls. Five designs were built and all failed calibration: a pooled
window null makes almost no calls, since 54 % of all four-residue windows in the
corpus carry some confident engagement; a per-partner null looked convincing on
real data but fired in 63 % of synthetic negative screens; corrected honestly
for every partner and position it makes zero calls anywhere, real or synthetic
(Supplementary Methods S2). Testing this claim needs control predictions of
pairs known not to interact — an experiment, not a reweighting. Until those
exist, CONFINE states the criterion, reports its parameters, and prints no
q-value beside it.

### 2.5 Merging sites into interfaces

Sites are called in sequence space, but a binding surface is three-dimensional:
two regions 100 residues apart can be adjacent on the folded protein, and
reporting them separately would be an artefact of the coordinate system rather
than a statement about the protein.

For every pair of called sites and **every partner model** in the screen,
CONFINE computes the minimum heavy-atom distance between the two sites' bait
residues in that model. *Support* is the fraction of models in which that
distance is below 8 Å. Sites are clustered by average linkage on 1 − support and
cut at 0.5, so two sites join an interface when they are adjacent in at least
half the models in which both are resolved.

Using the **ensemble** rather than one reference structure is deliberate: the
bait's conformation differs between predictions, sometimes substantially in
flexible regions, so a single model would make the merge an accident of which
prediction was chosen. Requiring agreement across half the screen makes the
result reference-free and, in practice, stable.

The rule was checked against experimental data. Among site pairs where both
members are covered by at least three residues of some known partner interface —
the pairs on which the question is answerable — 54 of 57 have the **same**
partner covering both, so the merge groups regions that genuinely serve one
binding function; the shipped criterion retains 46 of those 57 at 98 %
precision, and an alternative rule based on co-engagement by the same partners
added nothing (Supplementary Methods S3). The caveat is that the coordinates
being merged are predicted, so geometric agreement with a crystal conformation
is a weaker test than the partner-identity agreement quoted here.

An interface in this output is a **spatial grouping of bait-side regions**, not
an assertion that a single partner occupies all of it at once. Which partners
engage which part of an interface is a separate question CONFINE does not answer.

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

Three instruments were used, in increasing order of directness (Supplementary
Methods S4).

| instrument | what it tests | result |
|---|---|---|
| 120 synthetic negative screens | false-positive rate, null true by construction | 15 clusters in 15 screens (12 %); the earlier segment-scattering null gave 45 in 17 of 24 |
| annotated interfaces, 19 baits, 2,006 residues | agreement with experiment | recurrence calls 1,513 residues, 47 % annotated, **2.57×** enrichment, 42 % recall, and fires in 3 of 24 negatives; adding strength gives 63 % recall at 45 % precision<!-- combined figure measured at pct 90; re-measure at 95 --> |
| 68 solved complexes, truth on both sides | recovery of true interface residues | **71.6 %** of true bait interface residues fall in a called site, **70.4 %** in a significant one |

Three caveats belong with those numbers. The 12 % false-positive rate is a
real improvement on 45 clusters in 17 of 24, but it remains above the 5 % that
*q* < 0.05 implies under a complete null. The gold set's precision, 12.6 %, is a
**floor** only: a bait residue called for a partner other than the gold one
counts as a miss although it may be correct. And six gold cases score zero,
where the predicted binding mode does not overlap the crystal interface at all —
a prediction error that no downstream statistic can repair.

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

## Supplementary Methods

**S1 — Fixed *d*<sub>0</sub>, null variants and the cluster-forming
percentile.** With a length-dependent *d*<sub>0</sub>, partner length correlated
positively with raw PAE (median Spearman ρ = +0.28, positive in 23 of 24
screens) *and* positively with the resulting signal (median ρ = +0.32, positive
in 24 of 24) — larger partners scored better despite worse error. At a fixed
*d*<sub>0</sub> the same PAE gives the same signal in every screen. The
circular shift was compared against a non-wrapping variant, which places each
row at a random offset without wrap-around. Of the clusters on which the two
disagree, 42 are called only by the circular null; 45 % of their residues are
annotated interface residues, at 2.10× enrichment, against 44 % and 2.25× for
the 84 clusters both nulls call, and the non-wrapping null recovered nothing the
circular one missed. It is therefore over-conservative, discarding calls that
external data supports. Cluster-forming percentile sweep (24 baits, real screens
and one synthetic negative each):

| percentile | clusters | residues called | precision | enrichment | recall | negatives with a call |
|---|---|---|---|---|---|---|
| 80 | 154 | 1,723 | 44 % | 2.05× | 40 % | 2/24 |
| 85 | 147 | 1,607 | 44 % | 2.03× | 37 % | 3/24 |
| 90 | 150 | 1,588 | 45 % | 2.08× | 37 % | 4/24 |
| **95** | **151** | **1,513** | **47 %** | **2.57×** | **42 %** | **3/24** |

**S2 — Strength: rejected tests and parameter sweeps.** The five designs tested
were a pooled window null, a per-partner null, a per-partner null with an
ipTM-matched decoy set, a two-facet formulation separating amount from
single-partner confidence, and a Monte-Carlo cluster null over both facets.
Gap tolerance sweep: gap 0 gives 163 regions and 12 of 16 known sites; gap 2,
206 regions and 15 of 16; gap 4, 177 and 15 of 16; gap 8, 160 and 15 of 16.
Threshold sweep against annotated interfaces: the shipped 0.80 / 4 residues /
gap 2 gives 49 % precision, 2.21× enrichment, 56 % recall; 0.92 / 4 / gap 2
gives 56 %, 2.40×, 40 %; 0.92 / 4 / gap 0 gives 68 %, 2.74×, 19 %; 0.92 / 3 /
gap 0 gives 65 %, 2.84×, 23 %; loosening to 0.50 / 4 / gap 2 gives 40 %, 1.87×,
78 %. The looser shipped setting is a declared trade-off: the strength channel
exists to surface candidates recurrence structurally cannot see.

**S3 — Merge alternatives.** A first comparison of ensemble support against
partner co-engagement appeared to favour co-engagement, but scored unannotated
pairs as negatives. Restricted to evaluable pairs — both members covered by ≥ 3
residues of some known partner interface — the shipped support criterion keeps
46 of 57 pairs at 98 % precision and co-engagement adds nothing.

**S4 — Validation instruments.** *Synthetic negative screens* carry the bait
length and partner count of a real screen, but their rows are real partner rows
taken from the other proteins in the corpus, each cropped or embedded at a
random offset; recurrence's null hypothesis is therefore true by construction
and every call is a false positive. Cropping truncates donors, so these rows
carry roughly two-thirds of real signal and the false-positive counts are
conservative. The generator is seeded deterministically — an earlier harness
using a per-process hash gave 0 of 24 and 4 of 24 on successive runs of the same
test. *Annotated interfaces* are the interface residues aggregated over all PDB
structures for each bait, in UniProt numbering, restricted to heterotypic
interfaces; enrichment is the share of called residues that are annotated,
divided by that protein's own interface density, so 1.0 is chance. Absence of
annotation is not evidence of non-interaction, which is why the *gold set* — 68
bait–partner pairs whose exact complex is solved and whose partner is in the
screen, giving true residues on both sides — is preferred wherever it can answer
the question.

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
