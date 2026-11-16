# Patient-RNA splice evidence

The [current five-gene RNA shortlist](junction-only/README.md) is **F5, HLA-B,
HSF1, NBPF14 and KMT2D**. It adds full-window observed/reference-join amino-acid
identity to the earlier reconstruction/ND/strand and broad self-peptide rules.
The selected windows are assembled in the [current Design #3](../../../designs/design-3/README.md),
which has 30 regions and an 885-aa protein.
The [Design #3 change summary](../../../designs/design-3/README.md#rna-candidate-pipeline-changes)
and [current selection methods](../../../designs/design-3/methods.md#current-rna-selection-rules)
explain the added checks and their effect on the selected genes and windows.

This shared evidence reconstructs splice-source candidates from Yash's actual
RNA bases before HLA prediction or selection. Reference annotation supplies a
coding-frame anchor, rather than replacement bases for his sequence.

**Splice-source is a discovery label, not a confirmed explanation of every
sequence difference.** Earlier selected windows contain RNA-supported amino-acid
changes relative to their nominated reference: F5 S→F, MUC5B S→I, MUC3A I→F and
MUC6 T→M. These are earlier sequence versions, not the current five windows.
Their translations retain useful sequence evidence even when repeats
make the exact splice assignment uncertain. They are not direct protein
measurements or confirmed RNA-editing events.

## Find the evidence

For the chosen approach, start with the [current five-gene selection](junction-only/README.md)
and its [manual review](junction-only/manual-review.md). The default
[normal-peptide gate](../../../tools/rna_evidence/README.md#mandatory-normal-peptide-exclusion)
requires full-window reference-join identity and the broad self check.

The [recorded-alignment sequence screen](sequence-screen/README.md) and runs
below remain intact for reproducing the designs. Its 505 primary windows precede
the added identity rule, which leaves 167 eligible windows. The original MUC6
window fails the broader reference check; its substitution-dependent proposals
also fail the current identity rule.

The [expanded protein search](independent-assessment/protein-screen/protein-screen.json)
provides exact-reference and bounded personal-normal matches for the active
novelty gate. No remapping or alternative-selection workflow is required.

| Files | Contents |
|---|---|
| [Current selection](junction-only/README.md), [five selected windows](junction-only/selected.tsv) | Active reference-join identity rule: 167 eligible windows, 18 events, 10 genes; select one window from each of the top five genes |
| [policy.json](policy.json) | Read, frame, reference, binding and expression rules |
| [Input events](reconstruction/input-events.tsv), [outcomes](reconstruction/event-audit.tsv) | Every normal-filtered event, including rejected and unresolved events |
| [Patient windows](reconstruction/patient-windows.jsonl.gz) | Coding sequences, amino acids, alignment paths, fragment names and stop-free upstream frame witnesses |
| [Reference transcripts](reconstruction/reference-transcripts.tsv) | Complete reference CDS-to-protein translation checks |
| [Native review](review/receipt.json) | LENS review of actual patient coding sequences |
| [Exact-self check](reference-check/receipt.json) | Full annotated-reference inventory and two complete exact-match scans |
| [Affinities](scoring/affinities.tsv.gz), [all windows](scoring/all-windows.tsv.gz) | Fresh six-allele class-I predictions and scored candidates |
| [Event ranking](scoring/event-ranking.tsv), [gene ranking](scoring/gene-ranking.tsv) | Best qualifying window per event and gene |
| [Top five](scoring/selected.tsv), [sequences](scoring/selected.fasta) | Splice regions used by Design #1 |
| [Independent BAM check](read-verification/verification.json) | Every window and named fragment checked through an independent alignment-array implementation |
| [Inference receipt](scoring/receipt.json) | Inputs, model hashes, package versions and executed worker |
| [Sequence screening](sequence-screen/README.md), [admitted windows](sequence-screen/admitted.json) | Additional preselection checks across 752 previously eligible windows; 505 primary-tier windows across 15 genes, with review/excluded decisions retained |

Start with all 4,929 junction UIDs passing both native MSNAF normal-control
screens, before reference translation or HLA filtering. Require actual 75-nt
haplotypes in at least two qualifying fragment names. Overlapping mates must
agree; different fragments are never assembled into a synthetic consensus.
A verified GENCODE 37 upstream CDS anchors the frame, and observed codons from
that anchor through the join must be high-quality and stop-free.

Native LENS event and peptide review precedes admission. Exact reference-matching
join peptides do not contribute to scoring. The complete native Ensembl 103
inventory includes annotated protein isoforms and their source identities.
Fresh MHCflurry affinity predictions are checked against the OpenVax adapter.
An admitted region needs at least one nonself native-PASS peptide at <=500 nM.
Native Vaxrank affinity sums are multiplied by the square root of the supplied
gene value; no splice-fraction or raw-read abundance multiplier is used.

Rank windows, then the best per event and gene. Select five positive-expression
genes, one region per gene, with no pairwise sequence-similarity deduplication. A window is a
25-aa region, an event is an RNA join, and a gene may contribute several events.
These are the original-stage results; current selection additionally applies
the reference-join identity rule described above.

## Verify offline

After installing the repository, run from its root:

```sh
python workflow/stages/04_region_integration/patient_splice/replay.py \
  --results data/shared-rna/splicing --output work/splice-verification.json
```

This checks the input partition and hashes, translations, fragment/path/frame
witnesses, independent BAM-check receipt, exact-self decisions, complete allele
coverage, affinity contributions, expression multiplication and gene ranking.
It starts no cloud job, rereads no BAM and makes no new model predictions.
The normal initial-design reproduction includes the same verification.

The independent BAM check confirms bases and gaps in the recorded alignments.
It does not realign the reads against competing repeat placements. Junction
catalog non-detection, normal-DNA allele evidence and reference-peptide matches
answer separate questions and must not be substituted for one another.

For raw-read processing and fresh inference, follow the
[upstream instructions](../../../workflow/stages/04_region_integration/patient_splice/README.md).
Use matching staged inputs and authenticated Modal access. Explicit launchers
print plans unless --execute is provided.

A supported local sequence and stop-free frame do not establish a complete
translated isoform, tumour specificity, peptide presentation or a T-cell response.
Fragment names are not UMI molecule counts. Linear CIGAR gaps cannot establish
circular/back-splice joins; unresolved and nonlinear coordinates remain in the audit.
