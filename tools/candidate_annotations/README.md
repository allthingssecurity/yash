# Candidate annotations

For a readable view, [explore the lead Design #3 on the main site](https://rnafoundry.com/yash/#vaccine) or read the [annotation summary](../../designs/design-3/review/annotations.md).
The interactive page belongs to the RNA Foundry website repository. To export its verified data, use:

```sh
python tools/candidate_annotations/site_data.py --output /path/to/rnafoundry-www/site/yash/antigens.json
# Add --verify to check the exported asset and summary without rewriting them.
```

This exports recorded annotations; it performs no new inference or site build.

Add mutant/reference binding comparisons and MHCflurry processing/presentation
scores without changing a recorded selection or construct.

The [results](../../designs/design-3/review/candidate-annotations/README.md) annotate
the agreed 30-region plan and all 167 currently eligible splice windows.
These annotations do not change Design #3 selection or assembly.

## Run or verify

From the repository root, in the existing Python 3.11/3.12 environment:

```sh
# Read-only plan; no network or compute starts.
python tools/candidate_annotations/run.py

# New inference on a Modal A10 using the existing private reference volume.
# Choose a fresh directory for a new run; prior measurements are not overwritten.
python tools/candidate_annotations/run.py --execute --output work/candidate-annotations-new

# Offline replay of every annotation and ranking from collected predictions.
python tools/candidate_annotations/run.py --verify
```

Remote execution requires the Modal SDK, account access to the documented image
and reference volume, and the provider's model bundle already installed there.
The model is external, not redistributed. Offline verification needs only the
standard library and the included input/search/prediction evidence.

## Inputs and target peptides

`analysis.py` builds a request from 18 mutation regions (including the explicitly
approved H2BC12), five CTA/self regions, two fusion regions, and 167 splice
windows. The 167 splice windows have already passed the recorded reconstruction,
reading-frame, normal-junction ND, strand, broad self, affinity, expression and
hydropathy rules, including full-window reference-join amino-acid identity.
No new aligner or read-pattern gate is introduced.

The current plan marks five of those splice windows, giving 30 regions in total.
These inputs do not assemble a protein, optimize linkers or encode new RNA.
Splice annotations use the original window ID as `annotation_id`; other sources
use `region_id`. Two TSPY2 window records share a sequence/event region ID but
remain separate evidence records. Their peptide counts must not be combined.
There are 192 annotated region/window records and 191 distinct region IDs.

Only the recorded eligible target peptides are annotated. Mutation peptides must
overlap the altered interval. The in-frame fusion uses its recorded junction
codon and self-excluded peptide list, including its explicit 508.95-nM admission
exception. Ordinary flanks are not promoted into fusion-specific targets.
CTA/self peptides remain intentionally reference-derived.

## Mutant/reference comparison

For a single-residue substitution, restore the focal annotated reference residue
at the recorded mutation position and predict the counterpart peptide on the
same HLA allele. Report **mutant nM / reference nM**, without an exclusion cutoff.
Lower values mean a predicted binding gain for the mutant. TGFBR2's frameshift
has no single-residue reversion in this analysis and is explicitly not scored.
This implements the mutant/reference comparison proposed in
[PR #4](https://github.com/popvax/yash/pull/4), with fresh predictions and no
agretopicity cutoff.

These are focal-reference counterfactuals, not reconstructed fully phased normal
proteins. Other RNA-supported background residues remain unchanged, including
PPIC's additional background difference. The independent self screens remain in
place. Neither a large binding gain nor a ratio near one establishes T-cell
recognition, tolerance or safety.

## Processing and presentation

`predict.py` uses MHCflurry 2.3.14 and the exact 44-file model manifest already
pinned in this repository. Fresh affinities are compared with each region's
recorded best eligible affinity; relative disagreement above 0.0001 stops replay.

Each query is keyed by **peptide, N flank, C flank and context mode**. Identical
peptides from different contexts are predicted separately. This extends the
presentation analysis proposed in [PR #5](https://github.com/popvax/yash/pull/5)
without reusing a score across different contexts. The worker scores every
peptide without flanks and, where both sides are available, with up to five
recorded residues per side.
Missing residues are never invented.

Mutation contexts come from the recorded Isovar proteins; Guardant contexts are
explicitly reference-derived. The in-frame fusion has its recorded longer source
protein. Splice, CTA/self and the other fusion use the available 25-aa windows.
The output marks context basis, length, mode and complete five-residue flanks.

These are source-context annotations. The
[thirty-region assembly](../../designs/design-3/README.md) uses the
same selection; processing/presentation scores using its final vaccine linkers
and SEC/MITD boundaries are not claimed here. A model score is not a patient-specific probability of display or
a measurement of immunogenicity.

## Splice comparisons

All comparisons use the same 167 admitted windows and one best window per gene.
The original weighted-score/event/region ordering supplies deterministic ties.
Selection is unchanged.

- **Best presentation:** best no-flank score among a window's eligible
  peptide/HLA pairs with fresh affinity ≤500 nM. This reproduces the type of
  summary proposed in [PR #5](https://github.com/popvax/yash/pull/5), but drops
  expression weighting and summed binding breadth. Differences therefore do not
  isolate the effect of processing alone.
- **Expression-weighted presentation sum:** sum the no-flank presentation scores
  over target peptide occurrences and HLA pairs with affinity ≤500 nM, then
  multiply by √supplied gene value. This retains the previous splice expression
  multiplier and occurrence-count semantics. It is an explicitly exploratory
  comparison we defined, not a native Vaxrank score or a validated selection rule.
- **Available-context sensitivity:** repeat the comparisons with recorded flanks.
  Scores near window edges can depend on how much context is available, so these
  are annotations rather than a reason to select a one-residue shift as optimal.

Neither expression nor RNA counts enter MHCflurry. The expression multiplier is
applied afterward, where named explicitly. Summing model scores does not estimate
the number or abundance of presented peptides, or independent immune responses.

## Class II with CapHLA

`class_ii.py` adds class-II annotations to the same current plan and eligible
splice pool. It preserves the class-I inference files and every recorded design.

```sh
python tools/candidate_annotations/class_ii.py           # Plan only.
python tools/candidate_annotations/class_ii.py --execute # Start fresh Modal A10 inference.
python tools/candidate_annotations/class_ii.py --collect # Collect after completion.
python tools/candidate_annotations/class_ii.py --verify  # Offline replay.
```

For a new run use `--output work/class-II-new` on both execution and collection;
existing requests are not overwritten. Model inputs and per-target predictions
are saved on the persistent Modal volume. Repeating `--collect` while the job
is running only reports its state. Model definitions and weights are obtained from the provider's
MIT-licensed [CapHLA repository](https://github.com/changyunjian/CapHLA), pinned to
revision `33ebdd6ce6dadbbb1c66b026ce4b5d81dbf3a831`. The external adapter and all
14 model/source/license files in `caphla-manifest.json` must match their recorded
identities. No provider code, weights or training data are copied into this repo.
Provider example predictions must reproduce before patient scoring.

The input is the recorded 25-aa region, with every overlapping 15–20-aa peptide
that touches its altered interval or join. CTA/self regions contribute all
peptides. No sequence is extended or padded. For 17 single-residue mutation
regions, a focal-reference counterpart is also predicted on the same allele;
other background residues are retained. The TGFBR2 frameshift has no such
single-residue comparison.

[Class-II typing evidence](../../data/hla/class-II.json) records our earlier
tumour-DNA/RNA T1K calls and separately attributed, externally contributed
matched-normal WGS corroboration from [PR #7](https://github.com/popvax/yash/pull/7). The contributed caller was not
rerun here. Primary summaries use DRB1*03:01 and DRB1*08:02; additional summaries
keep the eight observed DQ/DP alpha-by-beta combinations exploratory. DRB3 is
excluded from active annotation while its typing evidence is retained.
DQ/DP pairs are not labelled measured cis/trans, and their functional expression
is not established. Exact model identities are required, with no nearest-allele
fallback. The observed DRA*01:01 and DRA*01:02 peptide-binding domains are
identical; their four possible DRA/DRB1 combinations map to two model targets.
Their recorded difference is in the cytoplasmic tail. This does not establish
equal cellular expression or trafficking. The reference alignment identity and
pairings are in the typing file. CapHLA's DR alpha-chain assumption and X tokens
in its DQ pseudosequences remain visible in the model-coverage table.

The five-fold BA and EL ensembles yield separate binding and presentation scores.
Predicted nM is `50000 ** (1 - mean BA score)`; EL is the mean positive-class
score. No percentile, clinical probability, expression multiplier, hard score
cutoff, combined class-I/class-II rank or selection change is introduced.

Every long peptide and focal-reference counterpart receives a fresh exact search
in the same pinned Ensembl 103/115 and UniProt canonical/additional catalogs,
including alternate translations. The existing personal-normal sequence library
is searched too, but **it was constructed for ≤11-aa targets**: its local
haplotype clustering and flanks do not establish complete 15–20-aa coverage.
Positive matches are witnesses; negative searches remain explicitly scoped and
are not called complete class-II personal-normal exclusion. No normal tissue
expression, translation or presentation is measured.

Standard CapHLA output does not assign an HLA-II binding core. Whole-peptide
overlap with a mutation/join and absence of an exact full-peptide reference match
do not prove that the actual binding core is altered or that CD4 cells recognize
it. All predictions remain available, alongside exact-search hit flags and
separate summaries of peptides without a hit in the searched catalogs. Intentional
reference CTA/self antigens are not silently discarded.

The [class-II results](../../designs/design-3/review/class-II-annotations/README.md)
join these annotations to the existing class-I table for the planned 30 regions.
The separate final build does not change this annotation run: no vaccine-junction
or full-cassette class-II processing result is claimed here.

## Reference-versus-candidate display

`reference_comparisons.py` reconstructs the five selected splice comparisons
from their original GENCODE 37 CDS anchors and transcript models. It compares
annotated protein sequence with the observed candidate, not measured normal RNA.
F5 and KMT2D retain their explicit deleted segments; HSF1 retains both distinct
reference-isoform contexts. TGFBR2 uses its Isovar upstream context to locate the
reference frameshift counterpart. Single substitutions restore only the focal
reference residue when exporting website data. These display annotations do not
change selection, peptides, affinities or the construct.

```sh
# Check the bundled comparisons and their selected-window identities.
python tools/candidate_annotations/reference_comparisons.py --verify
# Rebuild, or verify against the prepared upstream GENCODE reference models.
python tools/candidate_annotations/reference_comparisons.py --models /path/to/transcript-models.json.gz
python tools/candidate_annotations/reference_comparisons.py --verify --models /path/to/transcript-models.json.gz
```
