# Binding and presentation annotations

Fresh annotations for the **agreed 30-region plan**: 18 mutation, five
junction-derived splice, five CTA/self and two fusion regions. This plan includes
H2BC12 and the current splice shortlist. The
[Design #3 construct](../../README.md) assembles these thirty regions.
Scores here use source contexts, not the final linker contexts.

## Explore the results

| File | Contents |
|---|---|
| [Current plan](current-plan.tsv) | Thirty region rows with mutant/reference binding ratios and best presentation scores. |
| [Splice gene comparison](splice-gene-comparison.tsv) | Ten genes compared by original rank, best presentation and expression-weighted presentation. |
| [All splice windows](splice-window-comparison.tsv) | All 167 eligible 25-aa windows, their scores, sequences and ranking positions. |
| [All region annotations](region-annotations.tsv) | The 30-region plan plus the other eligible splice windows: 192 window/region records, 191 distinct region IDs. |
| [Peptide annotations](peptide-annotations.tsv.gz) | Every scored target occurrence/HLA pair, reference counterpart where applicable, and both context modes. |
| [Summary](summary.json) | Counts, scoring definitions, affinity parity and selection status. |
| [Measured predictions](predictions.tsv.gz) | Fresh, context-keyed MHCflurry outputs: 16,500 rows. |
| [Inference record](inference.json) · [verification](verification.json) | Model/runtime identities, input/output checksums and offline replay checks. |

## What changes under presentation ranking?

If we retain the existing splice expression multiplier and replace the binding
contribution with an exploratory presentation sum, **the same five genes and
exact 25-aa windows remain selected**. Only their ordering changes:

| Gene | Original rank | Expression-weighted presentation rank | Exploratory score |
|---|---:|---:|---:|
| F5 | 1 | 1 | 118.79 |
| HLA-B | 2 | 2 | 35.90 |
| NBPF14 | 4 | 3 | 23.24 |
| KMT2D | 5 | 4 | 21.80 |
| HSF1 | 3 | 5 | 17.93 |

The exploratory score is the sum of no-flank presentation scores over target
peptide occurrences/HLA pairs passing ≤500 nM, multiplied by √supplied gene
value. It is not a native Vaxrank score, a calibrated probability or measured
peptide abundance. MUC5AC is next at 9.67. Using available source flanks retains
the same five genes, but changes the preferred F5 and HLA-B window boundaries.

Selecting solely by the highest single presentation score would instead give
**TSPY2, MUC5AC, F5, HLA-B and PLK5**. That comparison also removes the previous
expression weighting and summed evidence across target peptides. It is not
enough, by itself, to justify replacing three genes. No selection rule or
construct was changed in this analysis.

## Scope

- Seventeen missense regions have paired mutant/reference affinity annotations.
  The TGFBR2 frameshift is marked outside that single-residue comparison.
- The 13 paired-DNA/RNA mutation regions and five reference-derived Guardant
  additions retain their different sequence-evidence bases. Predicted
  presentation does not supply missing Guardant RNA evidence.
- Source-context and no-flank predictions are separate. Splice context is limited
  to recorded 25-aa windows; incomplete context is marked. The proposed vaccine
  has not yet supplied assembled-context predictions.
- A mutant/reference ratio describes binding change, not T-cell discrimination.
  Presentation scores predict ligand-like properties, not a measured response.

Splice records are keyed by window ID (`annotation_id`), not only region ID.
Two TSPY2 windows share a region ID and are annotated separately, without
combining or double-counting their peptide contributions.

All 192 window/region best target affinities reproduce their recorded values within a maximum
relative difference of 0.00000166. Inference used MHCflurry 2.3.14 and the pinned
44-file bundle on a Modal NVIDIA A10; measured model time was 10.92 seconds.
The sandbox was terminated after collection. No external models are bundled.

[Code and complete scoring rules](../../../../tools/candidate_annotations/README.md)
explain the context handling, target restrictions and replay commands.

```sh
python tools/candidate_annotations/run.py --verify
```
