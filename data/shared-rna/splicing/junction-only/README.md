# Junction-derived RNA selection

The selected windows contain no additional amino-acid substitutions: actual RNA translation equals the reference bases translated across the same proposed join on every recorded path. Synonymous nucleotide differences are allowed. The proposed join is retained in that comparison; this is not a claim that the RNA is a normal transcript.

The earlier pool has 752 windows; 505 pass the preceding sequence/self checks. The added identity rule leaves **167 windows, 18 events and 10 genes**. We select one window from each of the top 5 genes.

| Rank | Gene | Selected 25-aa region | Best predicted target | HLA | Affinity, nM | Supplied gene value | Weighted score |
|---:|---|---|---|---|---:|---:|---:|
| 1 | F5 | `HTTLSLDFSQTNLSPELSHTTLSLD` | `LSPELSHTTL` | HLA-C*01:02 | 23.58 | 470.598 | 157.61 |
| 2 | HLA-B | `EAAREAEQRSLPGGHVRGVAPQIPG` | `EQRSLPGGHVR` | HLA-A*33:03 | 66.16 | 185.958 | 44.19 |
| 3 | HSF1 | `FSVDTSALLDVSGAPPPRLPAPPPR` | `VSGAPPPRL` | HLA-C*01:02 | 44.45 | 21.840 | 39.63 |
| 4 | NBPF14 | `QEQLAEGCRLAQHLVQKLSPGKVAI` | `KLSPGKVAI` | HLA-C*01:02 | 44.53 | 133.821 | 33.61 |
| 5 | KMT2D | `HLSPQPEEPHLSPQAEGPHLSPQPE` | `LSPQPEEPHL` | HLA-C*01:02 | 85.57 | 40.502 | 25.33 |

Selection and manual review use the original supplied RNA alignment and recorded rules. Scores and native event/region tie-breaks are unchanged. Equal-score windows can have different flanks; the selected IDs, not gene names alone, identify this selection.

[Manual review](manual-review.md) · [Selected JSON](selected.json) · [Selected TSV](selected.tsv) · [All ten genes](gene-ranking.tsv) · [All window decisions](all-windows.tsv.gz) · [Policy](policy.json) · [Summary](summary.json) · [Verification](verification.json)

## Reproduce

From the repository root:

```sh
python tools/rna_evidence/normal_peptides.py --top-genes 5 --results data/shared-rna/splicing/junction-only
python tools/rna_evidence/normal_peptides.py --verify --results data/shared-rna/splicing/junction-only
```

Verification recomputes admission and ranking from checksum-verified measured searches and recorded scores, then compares every generated selection artifact byte for byte. It does not reread BAMs, repeat protein searches or run HLA inference. Unknown search results stop the process. The selected windows remain computational hypotheses; translated tumour protein, HLA presentation and T-cell recognition have not been measured. No construct is rebuilt by this selection stage.
