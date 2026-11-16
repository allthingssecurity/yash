# Design #3 — current antigen annotations

**30 selected regions:** 18 mutation, 5 splice, 5 CTA/self and 2 fusion regions. The [assembled construct](../README.md) has 885 amino acids; assembly, junction checks and RNA encoding passed offline replay. Design #3 is the lead candidate.

[Explore on the main site](https://rnafoundry.com/yash/#vaccine) · [Combined source table](class-II-annotations/current-plan.tsv) · [Class-I methods](../../../tools/candidate_annotations/README.md) · [Class-II methods](class-II-annotations/README.md)

Binding is predicted affinity in nM (lower is stronger). Presentation is a separate model score (higher is more ligand-like; not a probability of a response). Each column is the best within the region and can refer to a different peptide or allele. No new cutoff, combined score or selection change is applied.

DR summaries use DRB1*03:01/*08:02. Eight DQ/DP pairing hypotheses are kept separately in the interactive view; DRB3 is excluded. Both observed DRA peptide-binding domains are identical, so four possible DRA/DRB1 combinations map to two model targets.

## Mutation

| Antigen | Supplied TPM | Class-I binding, nM | Class-I presentation | DR binding, nM | DR presentation | Notes |
|---|---:|---:|---:|---:|---:|---|
| AQP3 p.C11Y | 19.64 | 29.42 | 0.973 | 84.83 | 0.174 | Paired DNA + tumour RNA |
| AR p.P555S | 3.59 | 44.21 | 0.838 | 468.93 | 0.204 | Reference-derived window; Normal coverage exception: 0/19 ALT |
| ARHGAP40 p.R438W | 6.27 | 138.22 | 0.551 | 1,593.92 | 0.029 | Paired DNA + tumour RNA |
| DHX36 p.R160K | 30.39 | 100.76 | 0.775 | 76.41 | 0.775 | Reference-derived window |
| ERBB4 p.D113V | 0.51 | 41.65 | 0.957 | 207.93 | 0.195 | Reference-derived window |
| H2BC12 p.E114K | — | 136.46 | 0.396 | 685.25 | 0.994 | Explicit inclusion; supplied expression missing |
| HSPG2 p.R3952L | 167.61 | 81.11 | 0.658 | 77.71 | 0.997 | Paired DNA + tumour RNA |
| KCTD19 p.G920R | 9.91 | 36.94 | 0.901 | 149.74 | 0.996 | Paired DNA + tumour RNA |
| KLF5 p.D418N | 1,297.24 | 37.04 | 0.925 | 368.08 | 0.417 | Paired DNA + tumour RNA |
| KMT2D p.R2789Q | 40.50 | 41.45 | 0.780 | 725.72 | 0.616 | Reference-derived window |
| PLEC p.R2171W | 92.29 | 413.26 | 0.237 | 508.04 | 0.455 | Paired DNA + tumour RNA |
| PPIC p.P187S | 72.92 | 36.51 | 0.834 | 230.70 | 0.352 | Paired DNA + tumour RNA |
| RAF1 p.G535S | 19.83 | 118.47 | 0.829 | 1,051.32 | 0.256 | Reference-derived window |
| RXRA p.A241S | 23.01 | 58.08 | 0.896 | 889.68 | 0.740 | Paired DNA + tumour RNA |
| SF1 p.R218Q | 157.41 | 41.47 | 0.854 | 2,953.30 | 0.975 | Paired DNA + tumour RNA |
| SYT12 p.P391L | 7.92 | 243.74 | 0.574 | 50.67 | 0.309 | Paired DNA + tumour RNA |
| TGFBR2 p.N515fs | 55.74 | 38.03 | 0.785 | 960.90 | 0.490 | Paired DNA + tumour RNA |
| VPS11 p.D151A | 15.71 | 115.60 | 0.714 | 74.12 | 0.547 | Paired DNA + tumour RNA |

## Splice

| Antigen | Supplied TPM | Class-I binding, nM | Class-I presentation | DR binding, nM | DR presentation | Notes |
|---|---:|---:|---:|---:|---:|---|
| F5  | 470.60 | 23.58 | 0.933 | 230.71 | 0.126 | Best DR binding peptide has an exact searched self match |
| HLA-B  | 185.96 | 66.16 | 0.918 | 1,375.61 | 0.243 | RNA · reference-join translation |
| HSF1  | 21.84 | 44.45 | 0.724 | 837.58 | 0.142 | RNA · reference-join translation |
| KMT2D  | 40.50 | 85.57 | 0.770 | 3,557.52 | 0.181 | READ_END_BIAS |
| NBPF14  | 133.82 | 44.53 | 0.820 | 736.84 | 0.349 | RNA · reference-join translation |

## CTA/self

| Antigen | Supplied TPM | Class-I binding, nM | Class-I presentation | DR binding, nM | DR presentation | Notes |
|---|---:|---:|---:|---:|---:|---|
| LRRC36  | 9.15 | 21.01 | 0.975 | 256.65 | 1.000 | Expressed gene · reference protein |
| PRAME  | 17.01 | 15.30 | 0.956 | 175.70 | 0.595 | Expressed gene · reference protein |
| SPATA31A6  | 13.61 | 33.56 | 0.928 | 4,497.27 | 0.097 | Expressed gene · reference protein |
| TEKT5  | 4.84 | 23.74 | 0.946 | 145.33 | 0.831 | Expressed gene · reference protein |
| YBX2  | 15.48 | 29.75 | 0.965 | 148.50 | 0.720 | Expressed gene · reference protein |

## Fusion

| Antigen | Supplied TPM | Class-I binding, nM | Class-I presentation | DR binding, nM | DR presentation | Notes |
|---|---:|---:|---:|---:|---:|---|
| CDH1--NRG1  | — | 70.20 | 0.915 | 437.54 | 0.515 | RNA fusion-derived region |
| CDH1::NRG1  | — | 508.95 | 0.120 | 40.44 | 0.738 | Recorded class-I affinity exception: 508.95 nM |

## How to read these annotations

- Supplied TPM is the vendor-transformed gene value, shown as an annotation; it is not an input to either prediction model.
- Class I uses MHCflurry 2.3.14. CapHLA separately scores 15–20-aa peptides for class II.
- Reference comparisons restore only the nominated missense residue. Other background residues stay unchanged; these are not fully phased normal proteins.
- Long-peptide exact matches are reported explicitly. CTA/self matches are intentional; for altered sources, a match prevents treating that particular peptide as unique.
- Whole-peptide overlap with an alteration does not establish that the class-II binding core contains it. CapHLA does not assign that core here.
- Personal-normal sequence models were built for ≤11-aa targets; their negative results are incomplete for class II. Reference searches do not measure presentation in healthy cells.

For allele-specific winners, counterpart predictions and the remaining splice pool, use the [interactive view on the main site](https://rnafoundry.com/yash/#vaccine).
