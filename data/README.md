# Data for Yash

Start with the datasets below. Reports, expression values, class-I typing and the linked candidate tables are included in the repository.

[Samples](samples.tsv) · [Clinical timeline](clinical-timeline.md) · [File catalog](files.tsv) · [Lead design](../designs/design-3/README.md)

## Sequencing data

Sizes are decimal GB/MB. Download the index alongside a BAM or indexed VCF when your software needs random access.

<!-- datasets:start -->
| Dataset | When collected | Collected by | Files | Size | What it contains |
|---|---|---|---|---|---|
| Tumour DNA | 15 Jan 2025 | Caris | [BAM](https://yash-data.samarth-3f0.workers.dev/raw/TN25-112103/DNA_TN25-112103.bam) · [index](https://yash-data.samarth-3f0.workers.dev/raw/TN25-112103/DNA_TN25-112103.bam.bai) | 13.26 GB + 4.40 MB | Caris-aligned tumour DNA reads. |
| Caris variant calls | 15 Jan 2025 | Caris | [VCF](https://yash-data.samarth-3f0.workers.dev/raw/TN25-112103/DNA_TN25-112103.vcf) | 487.47 kB | Original vendor calls; separate from our paired tumour–normal calls below. |
| Tumour RNA reads | 15 Jan 2025 | Caris | [FASTQ R1](https://yash-data.samarth-3f0.workers.dev/raw/TN25-112103/original/RNA_TN25-112103_S45.R1.fastq.gz) · [FASTQ R2](https://yash-data.samarth-3f0.workers.dev/raw/TN25-112103/original/RNA_TN25-112103_S45.R2.fastq.gz) | 678.54 MB + 703.01 MB | Original paired RNA reads. |
| Tumour RNA alignment | 15 Jan 2025 | Caris | [BAM](https://yash-data.samarth-3f0.workers.dev/raw/TN25-112103/RNA_TN25-112103.bam) · [index](https://yash-data.samarth-3f0.workers.dev/raw/TN25-112103/RNA_TN25-112103.bam.bai) | 1.32 GB + 2.99 MB | Supplied tumour RNA alignment. |
| Matched-normal DNA | 30 Oct 2026 |  | [BAM](https://yash-data.samarth-3f0.workers.dev/raw/TN25-112103/normal/normal.recal.bam) · [index](https://yash-data.samarth-3f0.workers.dev/raw/TN25-112103/normal/normal.recal.bam.bai) | **143.29 GB + 9.06 MB** | White-blood-cell genomic DNA; sorted, duplicates marked and base qualities recalibrated. |
<!-- datasets:end -->

The normal BAM is the full final alignment. Original normal FASTQs and tumour DNA FASTQs are not included in the public file catalog. Dates refer to specimen collection, including derived files. The tumour date is recorded on page 1 of the [Caris report](reports/caris-tissue-2025-02-22.pdf). See the [sample table](samples.tsv) for sources.

Alignments use GRCh38 coordinates. Caris's exact custom-masked reference is unavailable; the [reference guide](../docs/external-resources.md) and [normal-alignment methods](../designs/design-2/methods.md#1-normal-alignment-and-sample-identity) identify the compatible analysis reference and its checksum.

## Expression and HLA typing

| Dataset | Open in the repository | Description |
|---|---|---|
| Gene expression | [vendor-transformed-expression.tsv](expression/vendor-transformed-expression.tsv) | Supplied gene-level expression, preserved byte for byte. [Original download](https://yash-data.samarth-3f0.workers.dev/raw/TN25-112103/vendor-transformed-expression.tsv). |
| HLA class I | [class-I.tsv](hla/class-I.tsv) | The six recorded class-I alleles used in our binding analyses, with their source. This is a summary of existing typing. |

The [methods](../designs/design-3/methods.md) explain how expression and HLA typing enter the design. Other typing or expression analyses can be contributed with their sources.

## Clinical reports

| Report | Date | Contents |
|---|---|---|
| [Stanford pathology and IHC](reports/stanford-pathology-ihc-2025-01.pdf) | January 2025 | Liver biopsy pathology, diagnostic IHC and the Claudin-18.2 addendum. |
| [Stanford STAMP](reports/stanford-stamp-2025-01-31.pdf) | 31 January 2025 | Tumour genetic panel. |
| [Caris tumour profiling](reports/caris-tissue-2025-02-22.pdf) | 22 February 2025 | Tumour DNA/RNA profiling, fusion, HLA and IHC results. |
| [Guardant360](reports/guardant360-2025-10-21.pdf) | 21 October 2025 | Liquid-biopsy cfDNA results. |

## Processed data and candidate evidence

These are analysis outputs. Use the specimen-level data above to develop an independent workflow, or start from the processed results here.

| Dataset | Files | Description |
|---|---|---|
| Paired tumour–normal calls | [All filtered alleles](https://yash-data.samarth-3f0.workers.dev/analysis/TN25-112103/design-2/paired/filter/all-unique.vcf.gz) · [index](https://yash-data.samarth-3f0.workers.dev/analysis/TN25-112103/design-2/paired/filter/all-unique.vcf.gz.tbi) · [PASS alleles](https://yash-data.samarth-3f0.workers.dev/analysis/TN25-112103/design-2/paired/filter/pass.vcf.gz) · [index](https://yash-data.samarth-3f0.workers.dev/analysis/TN25-112103/design-2/paired/filter/pass.vcf.gz.tbi) | 9.95 MB all-allele VCF; 357 kB PASS VCF. Native filter labels are retained; vaccine eligibility uses further criteria. |
| Normal coding genotypes | [GVCF](https://yash-data.samarth-3f0.workers.dev/analysis/TN25-112103/design-2/normal-coding/gather/normal.coding.g.vcf.gz) · [index](https://yash-data.samarth-3f0.workers.dev/analysis/TN25-112103/design-2/normal-coding/gather/normal.coding.g.vcf.gz.tbi) | 30.30 MB GVCF. Covers annotated coding/start/stop intervals plus 100-bp flanks; not a whole-genome genotype callset. |
| Paired mutation evidence | [DNA and normal evidence](../designs/design-2/inputs/mutations/variant-evidence.tsv) · [RNA evidence](../designs/design-2/inputs/rna/rna-evidence.tsv) · [RNA-supported proteins](../designs/design-2/inputs/rna/primary-proteins.tsv) | Paired-DNA assessment and RNA reconstruction for Design #2, including unsuccessful outcomes. |
| Paired-source candidate pool | [All eligible source regions](../designs/design-2/inputs/integrated/candidate-pool.tsv) · [selected regions](../designs/design-2/selection/selected.tsv) | Starting points for an alternative selection. |
| Paired-source binding predictions | [Mutation peptide/HLA summary](../designs/design-2/inputs/predictions/mutations/best-peptide-per-HLA.tsv) · [fusion affinities](../designs/design-2/inputs/predictions/fusions/affinities.tsv) | Recorded HLA-I results; detailed measurements remain with their sources. |
| Shared splice evidence | [Evidence guide](shared-rna/splicing/README.md) · [normal-screened gene ranking](shared-rna/splicing/normal-screen/gene-ranking.tsv) | RNA sequence reconstruction, read support, normal junction screening and binding predictions reused by the recorded designs. |
| RNA sequence screening | [Candidate table and evidence](shared-rna/splicing/sequence-screen/README.md) · [mandatory self-peptide gate](../tools/rna_evidence/README.md#mandatory-normal-peptide-exclusion) | Recorded-alignment evidence and broad self exclusion before the added junction-only identity rule. The broad self check is part of the current RNA gate. |
| Junction-derived RNA shortlist | [Selected five genes](shared-rna/splicing/junction-only/README.md) · [manual review](shared-rna/splicing/junction-only/manual-review.md) · [selected regions](shared-rna/splicing/junction-only/selected.tsv) | Current selection: F5, HLA-B, HSF1, NBPF14 and KMT2D. Observed/reference-join amino-acid identity throughout 25 aa; broad self exclusion retained. 167 windows/18 events/10 genes qualify. Used in the final Design #3 build. |
| Design #3 | [Build and methods](../designs/design-3/README.md) · [selection](../designs/design-3/selection/selected.tsv) · [protein](../designs/design-3/construct/protein/cassette.fasta) | Eighteen mutation, five splice, five CTA/self and two fusion regions; new joint assembly and PopRNA encoding, with independent junction and offline replay checks. The canonical Design #3 directory contains the current build. |
| Candidate annotations | [Current 30-region plan](../designs/design-3/review/candidate-annotations/current-plan.tsv) · [binding/presentation analysis](../designs/design-3/review/candidate-annotations/README.md) · [splice gene comparison](../designs/design-3/review/candidate-annotations/splice-gene-comparison.tsv) | Fresh context-keyed processing/presentation predictions and focal-reference affinity comparisons; all 167 eligible splice windows annotated. Selection unchanged. |
| Class-II annotations | [Combined class-I/class-II plan](../designs/design-3/review/class-II-annotations/current-plan.tsv) · [CapHLA results](../designs/design-3/review/class-II-annotations/README.md) · [typing evidence](hla/class-II.json) | Fresh 15–20-aa binding/presentation predictions, explicit DRB1/exploratory grouping and scoped long-peptide reference checks. No selection change. |
| Shared CTA/fusion evidence | [Scored source windows](shared-rna/checkpoints/lens-scored-windows.tsv.gz) · [in-frame fusion evidence](shared-rna/fusions/README.md) | Recorded RNA-source contexts and predictions. |
| Earlier mutation analysis | [Design #1 evidence](../designs/initial/inputs/mutations/README.md) | Tumour-only analysis used for the earlier design, without matched-normal DNA. |

Selected antigens and construct sequences live with the recorded lead
[Design #3](../designs/design-3/README.md), or the earlier
[Design #2](../designs/design-2/README.md) and
[Design #1](../designs/initial/README.md). The [file catalog](files.tsv) includes
public downloads and repository data files; its dataset, design and role columns
distinguish results from execution records.

## Download with checksum verification

For a verified download, run this from the repository root. This example downloads only the small expression file:

```sh
python tools/r2/restore.py --destination ../yash-data \
  --dataset expression --execute
```
To inspect or download a named dataset after installing the CLI:

```sh
yash-vaccine data list --dataset tumour-rna
yash-vaccine data restore --dataset tumour-rna --destination ../yash-data
```

The second command prints a plan. Add `--execute` to download. File IDs in
[files.tsv](files.tsv) can be selected individually with `--file`. The catalog
also generates this page's sequencing table and the root README's download table;
`python tools/catalog.py --write` refreshes those local views.
