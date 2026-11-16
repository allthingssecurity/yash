# Design #3 inputs

These files supply the current thirty-region selection. The selected output is
[selection/selected.tsv](../selection/selected.tsv); there is one Design #3 build.

| Input | Contents |
|---|---|
| [Source region catalog](source-regions.json) | Recorded candidate rows used to retain the mutation, CTA/self and fusion regions. The selector ignores splice rows in this upstream catalog and uses the current shared RNA selection instead. Source AG labels are local; compare stable `region_id` values. |
| [H2BC12 source catalog](h2bc12/source-regions.json) · [reference check](h2bc12/result.json) | The explicitly included H2BC12 E114K region, with missing expression. Only its row is used from this source catalog. Reference request, worker and receipt are kept alongside the result. |
| [Guardant nominations](guardant/reported-variants.tsv) · [all finding decisions](guardant/finding-decisions.tsv) | Report labels, plasma fractions and the reasons for the five missense additions. |
| [Matched-normal and tumour measurements](guardant/summary.tsv) | Normal genotypes/read counts and January tumour DNA/RNA counts. The Guardant sequence hypotheses are not full-window RNA reconstructions. |
| [Reference contexts](contexts/contexts.json) · [FASTA](contexts/contexts.fasta) | Ensembl reference flanks with each reported substitution applied. |
| [Predictions](predictions/peptide-HLA-affinity.tsv) · [selected additions](predictions/selected-additions.tsv) | Fresh class-I predictions, window scores and evidence exceptions. |
| [Reference matches](predictions/reference-matches.json) · [personal-normal matches](predictions/patient-normal-matches.json) | Scoped exact short-peptide search results for the Guardant additions. |
| [Sequence/normal checks](sequence-review/README.md) | Measured normal-DNA and broad reference evidence used by the RNA workflow. |
| [Current RNA selection](../../../data/shared-rna/splicing/junction-only/selected.json) | Five reference-join-identity regions from the original supplied RNA alignment. |
| [Policy](../policy.json) | Pinned input identities, source counts and explicit exceptions. |

Source catalogs preserve upstream measurements; their unselected rows are not
additional constructs. [Data and reports](../../../data/README.md) ·
[External resources](../../../docs/external-resources.md).
