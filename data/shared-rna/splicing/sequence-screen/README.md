# RNA candidate screening before selection

The stage screened **752 reconstructed windows, 42 events and 22 genes**. These already passed the recorded RNA/frame/native-peptide checks, coordinate-normal ND and limited strand checks. It examined all recorded alignment paths, every window nucleotide difference and matched-normal evidence before reranking. The construct was not changed.

After the new checks: **15 genes across 33 events** have a primary-tier window. Window decisions: PRIMARY: 505, EXCLUDE: 212, REVIEW_NORMAL_EVIDENCE: 35.

## Primary-tier candidates

| Rank | Gene | Selected 25-aa region | Best predicted target, nM | HLA | Full-window fragments | Supplied gene value | Weighted score |
|---:|---|---|---:|---|---:|---:|---:|
| 1 | F5 | `ALGQMPIFPDLSHTTLSLDFSQTNL` | 24.24 | HLA-C*01:02 | 2 | 470.598 | 272.94 |
| 2 | MUC5B | `LTTTATTPAATISTVTPSSALGTTH` | 44.87 | HLA-C*01:02 | 2 | 284.038 | 86.20 |
| 3 | MUC6 | `LPTHLPPFSTSLVTPSTHTVIITTH` | 11.63 | HLA-C*01:02 | 2 | 65.815 | 57.32 |
| 4 | MUC3A | `TSQPTTITPSSVGFSGSLPMMTDLT` | 260.85 | HLA-C*01:02 | 3 | 656.582 | 45.34 |
| 5 | HLA-B | `ITQRKWEAAREAEQRSLPGGHVRGV` | 66.16 | HLA-A*33:03 | 2 | 185.958 | 44.19 |
| 6 | HSF1 | `FSVDTSALLDVSGAPPPRLPAPPPR` | 44.45 | HLA-C*01:02 | 2 | 21.840 | 39.63 |
| 7 | NBPF14 | `RLAQHLVQKLSPGKVAIGPDDPKSQ` | 44.53 | HLA-C*01:02 | 2 | 133.821 | 33.61 |
| 8 | KMT2D | `HLSPQPEEPHLSPQAEGPHLSPQPE` | 85.57 | HLA-C*01:02 | 3 | 40.502 | 25.33 |
| 9 | POLE | `NCQRKMAWQWRGEFSECRPPFSTRW` | 48.09 | HLA-A*33:03 | 11 | 41.727 | 18.85 |
| 10 | MUC5AC | `TRTTSAPTTSTTLAPTTSTTSAPTT` | 28.85 | HLA-C*01:02 | 4 | 27.522 | 10.45 |
| 11 | MUC17 | `ASSPPPTAEGTSIPTSTSSEGNTPL` | 215.72 | HLA-C*01:02 | 2 | 133.620 | 9.01 |
| 12 | TSPY2 | `GRLEARRMRPEGSLTYRVPERLRQG` | 42.70 | HLA-C*07:02 | 7 | 4.458 | 7.77 |
| 13 | MDC1 | `DQSVTSEPTYQATRGRKNRSSVKTP` | 268.43 | HLA-B*15:08 | 2 | 98.145 | 7.28 |
| 14 | JAKMIP3 | `TPCTPDDDLEEVTRVSCRVLGLEEV` | 158.52 | HLA-A*33:03 | 21 | 0.324 | 2.98 |
| 15 | PLK5 | `KYHPRPHLPNTCPRLPPAHTSHLLT` | 50.29 | HLA-C*01:02 | 2 | 1.767 | 2.97 |

The table chooses the best eligible window per event, then per gene; no quota or automatic additions are applied. Fragment names support sequence and are not UMI-confirmed independent molecules. The score uses cached class-I affinity contributions × √supplied gene value, not raw-read abundance.

## What the checks mean

- Exact target-peptide matches in all Ensembl translations, UniProt canonical/additional isoforms or the bounded patient-normal coding/HLA models are excluded from scoring.
- Every nucleotide difference is compared with normal DNA, including synonymous/background differences. If novelty relies on extra RNA-inferred amino-acid changes, the primary tier requires resolved normal evidence and at least one reference-supported change. A normal-carried background does not by itself invalidate a separate RNA change. Missing/conflicting or solely normal-carried changes remain in the review tier.
- A same-alignment reference-bases counterfactual distinguishes junction-context novelty from novelty requiring extra RNA changes. It is an attribution check, not a proposed antigen. When that junction context is already nonself, no somatic-DNA requirement is invented for the RNA junction.
- All recorded paths must agree on primary admission. These are the alignments already present in the BAM; the stage does not establish their unique correctness or realign repeats.
- Admission retains a native-PASS peptide with predicted affinity ≤500 nM, positive supplied expression and the existing hydropathy rule. A model prediction is not measured translation, presentation, immunogenicity or tumour specificity.

## Evidence and reruns

The run examined **102 unique nominated alleles** and **2362 reference/counterfactual peptide queries**. It reused 12,204 checksum-pinned peptide/HLA predictions; it made no new HLA-model inference. The original native score is verified for every input window before adding exclusions.

- [Primary gene ranking](gene-ranking.tsv) · [event ranking](event-ranking.tsv) · [all windows](windows.tsv.gz).
- [Review-inclusive gene ranking](review-gene-ranking.tsv) · [per-peptide tiers and novelty attribution](peptides.tsv.gz).
- [Normal DNA and fresh focal tumour/RNA counts](normal-sites.json.gz) · [all-path sequence differences](window-evidence.json.gz).
- [Reference witnesses](reference-matches.json.gz) · [scope and summary](summary.json) · [completion receipt](complete.json) · [offline replay](verification.json).
- [Frozen inputs and source hashes](request.json.gz) · [pipeline code and commands](../../../../tools/rna_evidence/README.md).

## Use the checked candidates

[admitted.json](admitted.json) contains only primary-tier windows, with source IDs, sequences, target masks and nonself target peptides. [exports.json](exports.json) binds those exports to the measured screening inputs. Select additional regions from this pool and run the selection gate before assembly; review-tier windows require a separate documented decision.
