# AlphaSIFT — statistical inference of functional interaction sites from structure-prediction screens

AlphaSIFT turns a set of AlphaFold predictions that share one **bait** protein into
two things:

1. a **confidence matrix** — one row per partner prediction, one column per bait
   residue, valued by how confident AlphaFold is about the contact there; and
2. the **interaction sites** on that bait that the screen statistically supports,
   merged into 3-D interfaces.

It answers "which parts of this protein do many partners confidently engage?" and
it does so from predicted aligned error, not from distance alone.

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

## Install

    pip install alphasift                    # numpy, scipy, pandas, matplotlib
    pip install "alphasift[zst]"             # + zstandard, for .zst AFDB downloads

### Or use the provided environment

If you would rather not touch your own environment, `environment.yml` builds a
complete one, test runner included:

    conda env create -f environment.yml
    conda activate alphasift
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
end to end, and that the figure renders.

Note that `scipy >= 1.11` is required: the recurrence FDR uses
`scipy.stats.false_discovery_control`, which older versions do not have.

## Use

    alphasift run path/to/screen -o out/

writes `matrix.png`, `sites.tsv`, `interfaces.tsv`, `matrix.npz` and
`summary.json`. Every parameter above is a flag (`--pct`, `--d0-signal`, …).

```python
from alphasift import analyze
res = analyze("path/to/screen")
res.M_pae        # partners × residues, raw min-PAE (Å); NaN = no contact
res.M_signal     # the same, transformed to (0, 1]; 0 = no contact
res.sites        # significance, channel, cores, interface
res.interfaces   # 3-D merged interfaces
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
than merely shrinking it. AlphaSIFT lists every excluded prediction in
`summary.json` and warns when any are missing.

The two sources share a 0–1 scale but not its meaning: on the PAE channel 0.8 is
PAE 5 Å, while on the AlphaBridge channel it means the residue survives to
cut-off 0.8. Threshold parameters — including the merge's partner thresholds —
keep their numeric defaults across sources but not their interpretation.

| Predictor | Layout | Status |
|---|---|---|
| AlphaFold3 (server or local) | `<pred>/<pred>_model.cif` + `<pred>_confidences.json` | supported |
| AlphaFold2 / ColabFold, AFDB bulk packaging | `<id>-model_v1.cif[.zst]` + `<id>-predicted_aligned_error_v1.json[.zst]` | supported |
| AlphaBridge run (AF3) | `<pred>/AlphaBridge/alphabridge_data.json` | supported as `--source alphabridge`, and as a PAE shortcut where a `confidence_heatmap` is stored and the raw PAE is not |
| AlphaFold2 official (`run_alphafold.py`) | `ranked_0.pdb` + `result_model_*.pkl` + `ranking_debug.json` | not yet |
| ColabFold local (`colabfold_batch`) | `*_unrelaxed_rank_001_*.pdb` + `*_scores_rank_001_*.json` | not yet |

If a screen has two shared sequences, or you are analysing a single pair (where
every chain is trivially "common"), name the bait yourself:

    alphasift run path/to/pair -o out/ --bait-chain A

### Running AlphaBridge first (only for `--source alphabridge`)

AlphaBridge is a **separate tool**, not a dependency of AlphaSIFT and not
installed with it: <https://github.com/PDB-REDO/AlphaBridge> (web server:
<https://alpha-bridge.eu/>). It carries its own licence terms and its own
citation requirement — read them before use.

    git clone https://github.com/PDB-REDO/AlphaBridge
    conda env create -f AlphaBridge/environment.yml --name alphabridge
    conda activate alphabridge

It processes one AlphaFold3 prediction folder at a time, writing
`AlphaBridge/alphabridge_data.json` inside it:

    python3 define_interfaces.py -i /path/to/screen/<prediction>

AlphaSIFT needs that file for **every** prediction in the screen, so run it across
the whole directory before analysing:

    for d in /path/to/screen/*/; do python3 define_interfaces.py -i "$d"; done
    alphasift run /path/to/screen -o out/ --source alphabridge

Predictions AlphaBridge fails on are listed in `summary.json` and excluded from
the matrix, and AlphaSIFT warns when any are — see the note above on why a partial
run changes the test rather than merely shrinking it.

**AlphaBridge is not required.** It is one optional source, not a dependency:
with no AlphaBridge output present, sequences and PAE are read from the
prediction files directly.

## Deliberately out of scope

The partner-side steps that followed site calling in the original pipeline —
per-partner contact sets, binding peptides, structure-anchored registers and
sequence motifs. Those depend on AlphaBridge link confidence and on motif
machinery, and are a different paper.

## Scope limit worth stating plainly

Recurrence needs a screen: it asks whether partners converge *relative to each
other*. A single pair cannot answer that, so with one prediction AlphaSIFT reports
the matrix and strength regions and explicitly declines to report a q-value.

## Provenance

The matrix construction, the statistical engine and the 3-D merge are extracted
**verbatim** from the `viz_int` research pipeline (last re-synced 2026-10-01), so
published validation of those numbers applies to this code unchanged. What is new here is the input
dispatch layer, bait detection without AlphaBridge, the figure, and the CLI.

## Licence and citation

**Licence: not yet fixed — see the note below before distributing this.** The
intended terms are free availability for non-commercial users, and continued
availability for at least two years after publication, as the journal requires.

**Relationship to AlphaBridge.** AlphaSIFT contains no AlphaBridge code, does not
import it, and does not require it to run; AlphaBridge output is one optional
input format among several. AlphaBridge is distributed by its authors under its
own, separate agreement, which restricts redistribution and derivative works and
requires that results obtained through its use cite its publication. If you
produced your screens with AlphaBridge, those terms apply to that use — cite it
accordingly.

Citation for AlphaSIFT: the Application Note and Zenodo DOI, once published.
