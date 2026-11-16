# Design #3 methods

Design #3 contains 30 regions: thirteen paired-DNA/RNA mutations, five Guardant reference-derived mutation hypotheses, five splice regions, five CTA/self regions and two fusion hypotheses. The resulting SEC–antigen–MITD protein is 885 aa. Source measurements, selected sequences and known nucleotide segments are preserved by checksum.

## Paired-DNA/RNA mutation regions

The paired tumour–normal calling, DNA filters and Isovar reconstruction are described in [Design #2's methods](../design-2/methods.md). Twelve selected regions are reused with their measured evidence. H2BC12 E114K is explicitly included from the same results: normal ALT/total 0/54, RNA ALT/total 11/57 and eight longer-context fragment names. Its supplied expression is missing and not imputed. The [additional Ensembl 115 check](inputs/h2bc12/result.json) found no hit among its 38 scored targets. [Issue #1](https://github.com/popvax/yash/issues/1) highlighted its earlier missing-expression exclusion. FES remains outside the fixed 25-aa recipe; no reference padding is invented for its 15-aa RNA context.

## CTA/self and fusion regions

PRAME, SPATA31A6, LRRC36, YBX2 and TEKT5 are reference-protein hypotheses with positive supplied gene expression. They are intentionally self-derived and do not require full-window patient-RNA reconstruction. Their discovery, normal-tissue and mTEC annotations remain in the [shared source evidence](../../data/shared-rna/checkpoints/lens-scored-windows.tsv.gz).

The two CDH1::NRG1 regions represent distinct sequence hypotheses, not two independent driver discoveries. The in-frame product is an explicit ~508.95-nM exception; the alternative is a predicted frameshift-tail product. [Fusion evidence](../../data/shared-rna/fusions/README.md). Neither binding predictions nor RNA joins establish the actual tumour-cell ligand.

## 1. Guardant nomination and coordinate mapping

The [Guardant360 report](../../data/reports/guardant360-2025-10-21.pdf) describes
blood collected on **14 October 2025**, reported on 21 October. It labels the
five calls below as variants of uncertain significance. The earlier tumour
material was collected on **15 January 2025**; the matched-normal white-blood-cell
sample on **30 October 2026**. Differences between these time points are not
automatically sequencing errors or proof of acquired tumour evolution.

As [Yannick Buccella noted](https://x.com/YannickBuccella/status/2120612621921702280),
intervening treatment makes the age of the tumour profile relevant. This analysis
supplies no newer 2026 tumour/plasma profile; the normal-cell sample cannot
substitute for one.

| Guardant protein label | Plasma allele fraction | Inferred GRCh38 allele |
|---|---:|---|
| AR P555S | 0.4% | chrX:67643302 C>T |
| DHX36 R160K | 0.2% | chr3:154315170 C>T |
| ERBB4 D113V | 0.3% | chr2:211947513 T>A |
| KMT2D R2789Q | 3.2% | chr12:49039222 C>T |
| RAF1 G535S | 0.2% | chr3:12585187 C>T |

### Why the expansion has five additions

We accounted for the report’s current variant list before defining the additions.
Its **11 current VUS** comprise **five new missense calls, two mutations already
included, two previously assessed stop-gains and two synonymous changes**.
The separate fusion, promoter, potential clonal-haematopoiesis and historical
response-map entries do not add another missense region to this design.

| Finding | Treatment in Design #3 | Rationale |
|---|---|---|
| AR P555S; DHX36 R160K; ERBB4 D113V; KMT2D R2789Q; RAF1 G535S | Five reference-derived additions | Additional current missense calls without matching entries in the parent protein-changing review. Each produced a nonself predicted binder ≤500 nM, a hydropathy-compatible window and a positive supplied gene value. Normal ALT was not detected; AR retains its coverage exception. |
| RXRA A214S; TGFBR2 N515Kfs*26 | Already included | Parent regions represent these events, using p.A241S and p.N515fs annotations respectively. No duplicate regions were added; transcript-dependent numbering is preserved. |
| JUN R276* | Previous review retained | The stop allele has RNA support, but Isovar’s `protein_sequence_contains_mutation` check was false and no eligible altered-protein region was selected. |
| MERTK E192* | Previous review retained | The earlier RNA check found 0 alternate fragments among 25 covering fragments. No new stop-terminal target analysis was performed. |
| DHX36 L147L; ROS1 F1699F | Not added by this sequence-change screen | Synonymous labels do not directly change the encoded amino acid. Possible RNA/splicing consequences were not evaluated from these labels. |
| CDH1–NRG1 fusion | Already targeted | Both parent fusion regions are retained. The report supplies no additional distinct fusion-junction sequence to reconstruct. |
| TERT promoter SNV | Outside the coding-mutation branch | The promoter change does not directly encode an altered amino-acid region. This does not evaluate TERT as an expression/self-antigen target. |
| DNMT3A C351Afs*56 | Not nominated as tumour-derived | Guardant flags potential clonal haematopoiesis. We did not independently confirm blood-cell origin or perform a targeted normal check for this call in the five-variant review. |
| BCOR L1632P | Historical entry, not nominated | The response map shows 0.2% in March and ND in October. ND is not proof of biological absence; no separate target assessment was performed. |

The RNA-reconstruction requirement was deliberately waived for the five missense
additions. **JUN and MERTK were not excluded by pretending the new five had
stronger RNA evidence.** Both introduce premature stops, which do not encode a
replacement amino acid at the stop position. Their ordinary upstream reference
sequence does not supply the changed-residue window used for this expansion;
new terminal or altered-processing targets would need a separate analysis.
The prior JUN/MERTK RNA outcomes above are retained evidence, not a uniformly
reapplied RNA admission rule for this new branch.

All five additions also passed the reported affinity, exact-self and hydropathy
checks. **There was no fixed five-antigen quota and no claim that these are the
only possible vaccine targets in the report.** This is a missense-focused
expansion of the existing design, with explicit nomination and evidence
exceptions. VUS status was not a rejection criterion.

The [decision log](inputs/guardant/finding-decisions.tsv) records every genomic
finding in the current summary/VUS/CH lists and the historical BCOR response-map
entry. TMB, MSI/HRD, tumour fraction, pharmacogenetics and HLA typing are contextual
measurements, not additional altered peptide sequences nominated by this screen.
This accounting documents the scope and recorded decisions; it is not a new
executable filter or a retrospective claim that all findings underwent the same
new prediction pipeline.

### Coordinate mapping

The report supplies protein labels without genomic coordinates or transcript
accessions. We searched compatible Ensembl 103 transcripts, checked the stated
reference amino acid, mapped coding nucleotides through exon coordinates and
verified reference codons against GRCh38. Candidate single-nucleotide changes
had to produce the stated amino-acid substitution and round-trip through Varcode.
Each label had one compatible allele in this search. These coordinates are
**inferred interpretations**, not coordinates supplied by Guardant.

One ERBB4 transcript lacked a complete cDNA coding string and began with an
unresolved amino acid. A separate check translated its genomic CDS from the
first complete codon and matched the 112 complete amino acids through position
113. It gave the same D113V allele as the two cDNA-mapped transcripts; missing
first-codon nucleotides were not invented.

## 2. Direct matched-normal and tumour review

The existing duplicate-marked, recalibrated normal BAM and native coding-region
GVCF were read on Modal. The GVCF/index hashes were checked against the completed
normal run. The BAM’s recorded hash and measured file size identify it; its
143-GB full-file hash was not recomputed. The original sample-identity receipt
reports a tumour–normal match; identity processing was not rerun here.

Read evidence required MQ≥20 and focal BQ≥20, excluded unmapped, secondary,
supplementary, duplicate and QC-fail alignments, and required NH=1 when present.
Paired alignments were grouped by read name; these are not UMI-defined independent
molecules. Conflicting mates were reported separately. An independent direct
alignment-position counter agreed with the normal REF/ALT counts at all five sites.

| Variant | Normal ALT/total fragments | Normal genotype GQ | January tumour DNA ALT/total | January RNA ALT/total |
|---|---:|---:|---:|---:|
| AR P555S | 0/19 | 81 | 0/465 | 1/19 |
| DHX36 R160K | 0/34 | 99 | 0/53 | 0/7 |
| ERBB4 D113V | 0/34 | 99 | 1/871 | 0/22 |
| KMT2D R2789Q | 0/59 | 99 | 3/2,563 | 0/139 |
| RAF1 G535S | 0/36 | 99 | 0/1,442 | 1/218 |

All five have reference-genotype evidence and no qualifying normal ALT fragments.
Four meet the parent’s additional ≥20-fragment normal-reference rule. **AR has
19 fragments and is explicitly included despite that coverage exception.**
Normal non-detection does not establish tumour origin, rule out every low-level
blood-cell clone or prove that the altered peptide is absent from all healthy cells.

RXRA A214S, JUN R276* and MERTK E192* served as controls. Their inferred coordinates
match the prior records and recover the expected tumour DNA support. None of the
five new inferred alleles has an exact record in the existing paired all-allele
callset; direct read measurements, rather than missing VCF records, support the
table above.

Evidence: [complete measurements](inputs/guardant/results.json),
[summary](inputs/guardant/summary.tsv), [execution receipt](inputs/guardant/complete.json)
and [incomplete-transcript resolution](inputs/guardant/incomplete-transcript-resolution.json).

## 3. Reference-derived altered protein contexts

For every compatible mapped transcript, we took up to 24 reference amino acids
on each side of the nominated position and substituted the reported alternate
amino acid. Identical mutant contexts at the same local target position were
deduplicated while retaining their transcript IDs. Every context had to contain
at least 25 standard amino acids.

These are reference-based protein models. They do not incorporate unmeasured
local tumour haplotypes and are not claimed to be full-window patient-RNA
reconstructions. The January read counts remain measurements alongside them,
not confirmation of these complete sequence models.

## 4. HLA-I prediction and 25-aa selection

We enumerated 8–11-aa peptides overlapping the substituted amino acid and ran
fresh MHCflurry 2.3.14 affinity predictions using the same pinned 2.3.0 affinity
model bundle and six HLA-I alleles as the parent. The small antigen screen used
the CPU analysis runtime and native mhctools affinity adapter. Larger junction
requests used an A10 GPU with batch size 16,384 and complete adapter-parity checks.
Model file identities were checked in both stages.
No class-II, presentation or immunogenicity model entered selection.

Mutation-target scoring excluded exact matches in the native annotated reference
proteome and the parent’s explicitly scoped normal coding/HLA sequence models.
No hit in those references is not proof of absence from all healthy tissue.

Every possible 25-aa window containing the changed residue was scored using
native Vaxrank epitope logic. The sum was independently checked against the
recorded affinity contribution:

```text
contribution(a) = (1 + exp(-350/150)) / (1 + exp((a - 350)/150))
```

Here `a` is predicted affinity in nM. Contributions at affinity ≥5,000 nM, or
below 0.00001, are zero. Scores sum eligible peptide/HLA contributions, not just
the strongest peptide. Within each mutation, windows with any 19-aa mean
Kyte–Doolittle score >1.6 were excluded. The highest-scoring remaining window
was selected, breaking ties by mutation centering, context identifier and start.

**All five nominated variants are included.** The parent’s ≤500-nM gate is
reported for each addition; it is not silently described as an admission rule
for this explicitly nominated tier. Any weak-binding exception is recorded in
`admission_exception` alongside the absence of RNA reconstruction and AR’s
normal-coverage exception.
In the completed antigen screen, **all five meet the ≤500-nM affinity gate**;
no weak-binding exception was needed.

Supplied gene values are annotations. Guardant cfDNA allele fractions are not
RNA fractions and are not multiplied by expression. Neither expression magnitude
nor a top-N quota removes one of these five additions.

## Current RNA selection rules

The change addresses two separate questions: **is a scored peptide already in
the searched normal/reference sequences, and does its sequence require the
nominated join without additional amino-acid substitutions?** Read support,
normal-junction detection and protein-reference membership answer different
parts of that assessment. This revision concerns the splice-derived pool;
fusion and CTA/self selection is unchanged.
The specificity-attribution finding and MUC6 alternate-reference match are
described in [Issue #1](https://github.com/popvax/yash/issues/1), with a subsequent
[independent reproduction](https://github.com/popvax/yash/issues/1#issuecomment-6804436820).
[PR #6](https://github.com/popvax/yash/pull/6) proposes a related alternative;
the admission rules below describe our current implementation, rather than
adopting that PR’s recipe or construct.

1. **Retain the original reconstruction and frame requirements.** Start with
   LENS-nominated events and windows reconstructed from actual bases in the
   supplied tumour RNA alignment. Every base of the complete 75-nt/25-aa window
   needs the recorded support from at least two qualifying fragment names;
   overlapping mates must agree. A verified GENCODE 37 upstream CDS anchors a
   stop-free observed translation. Fragment names are not UMI-confirmed molecules.
2. **Retain coordinate-aware normal ND and strand admission.** Look up the exact
   genomic junction coordinates across gene assignments and both strands in the
   processed GTEx and TCGA-adjacent-normal catalogs. ND requires no coordinate
   entry in either catalog; entries with zero counts or opposite-strand labels
   remain review findings. Missing identifiers are not imputed as measured zeros.
   Full-window witnesses must retain a strict expected-orientation majority under
   the recorded library-orientation rule. Motif disagreement remains report-only.
   These checks leave 752 of the 3,804 previously eligible windows.
3. **Apply the sequence/normal-evidence and broad peptide checks.** The preceding
   [sequence screen](../../data/shared-rna/splicing/sequence-screen/README.md)
   examines recorded paths, nucleotide differences and their normal-DNA evidence.
   Exact 8–11-aa target matches in the broad reference/personal catalogs are
   removed before recalculating affinity contributions and admission. This stage
   has 505 primary windows across 33 events and 15 genes. These primary labels
   alone do not imply passing the next rule.
4. **Require full-window reference-join identity.** For each recorded supporting
   path, retain its proposed join and frame, replace observed nucleotide
   differences with their nominated genomic reference bases, and translate.
   Every resulting 25-aa translation must equal the observed RNA window.
   Synonymous differences are allowed. A window requiring an extra amino-acid
   substitution fails even if that substitution has RNA support. The reference
   comparison still contains the proposed join; it is not a known normal
   transcript. This leaves 167 windows across 18 events and 10 genes.
5. **Retain the existing class-I admission and ranking.** A window needs a remaining
   native-PASS nonself target with affinity ≤500 nM, positive supplied expression
   and the existing hydropathy admission. Rank by the native Vaxrank affinity
   contribution sum × √supplied transformed gene value, using the original
   event/region tie-breaks; choose the best passing window from each of the top
   five genes. Raw RNA abundance and splice fractions are not expression weights.

This selects **F5, HLA-B, HSF1, NBPF14 and KMT2D**. The broad self gate applies to scored altered targets, while ordinary self flanks may remain. CTA/self regions are intentionally outside that novelty gate. Negative reference searches are scoped to the catalogs and bounded genotyped models; they do not measure healthy-tissue presentation.

The exploratory class-I presentation-sum × √gene-value comparison retains the same five windows. Presentation and class-II results remain annotations, not selection criteria.

```sh
python tools/rna_evidence/normal_peptides.py --verify --results data/shared-rna/splicing/junction-only
```

[All decisions](../../data/shared-rna/splicing/junction-only/all-windows.tsv.gz) · [Selected regions](../../data/shared-rna/splicing/junction-only/selected.tsv).

## Additional annotations

Fresh focal-reference affinity comparisons follow the question raised in [PR #4](https://github.com/popvax/yash/pull/4); no agretopicity cutoff is applied. Processing/presentation annotations extend [PR #5](https://github.com/popvax/yash/pull/5), keying predictions by peptide and flanking context. [Class-I implementation](../../tools/candidate_annotations/README.md).

CapHLA separately predicts 15–20-aa class-II binding and presentation. The two DRB1 targets are primary; eight unphased DQ/DP combinations remain exploratory and DRB3 is excluded. Matched-normal typing contributed in [PR #7](https://github.com/popvax/yash/pull/7) corroborates the earlier T1K calls; that contributed caller was not rerun. Class-II core novelty and complete personal-normal exclusion are not established. [Class-II scope and results](review/class-II-annotations/README.md).

## Joint order and linker optimization

Use the unchanged patient-supported SEC and MITD segments and the eleven recorded
linker choices. Consider every ordered pair of distinct antigens, every SEC-entry
pair and every MITD-exit pair. With thirty antigens this is 930 directed pairs ×
eleven linkers = 10,230 options. Every 8–11-aa peptide touching a linker or direct
join is scored against the six recorded HLA-I alleles by fresh MHCflurry 2.3.14
inference. Identical peptide sequences are deduplicated for prediction only;
edge costs preserve occurrences and allele-specific contributions.

Reject internal options failing the existing 19-aa Kyte–Doolittle threshold.
Use the unchanged four-phase OR-Tools CP-SAT objective: minimize boundary
peptide/HLA occurrences at percentile ≤0.5, then ≤2; maximize the worst boundary
percentile; minimize inserted residues. Each phase must be proven optimal before
export. SEC/MITD edges are included in prediction and optimization; fixed routing
hydropathy overlaps are reported separately. Full-chain enumeration must agree
with edge-wise enumeration and contain every antigen exactly once.

Class-I presentation, focal-reference affinity ratios and CapHLA class-II scores
remain [source-region annotations](review/annotations.md). They are not added
to this optimization objective. Their predictions do not describe the final
linker context, and no class-II whole-cassette screen is claimed.

## Independent junction and hydropathy checks

Run fresh inference again on only the selected junction peptides. Check exact
enumeration and the same six alleles, identical ≤50/≤500-nM and ≤0.5/≤2-percentile
classifications, relative affinity disagreement below 0.0001 and absolute
percentile disagreement below 0.0001. Independently rescan the placed full chain,
verify all component sequences and apply the existing hydropathy policy.
These are computational checks, not measured processing, topology or presentation.

## RNA encoding and annotated exports

Use the pinned PopRNA implementation with its recorded human codon table,
10% codon cutoff, motif list, random seed 0 and folding/alternative-frame-stop
heuristics. Protein order and linker sequences are fixed before encoding. Keep
the existing 5′ UTR, three stop codons and 3′ UTR. Only the coding sequence changes
synonymously. PolyA and any downstream extension remain unspecified; the output
FASTA and GenBank represent known segments, not a complete molecular sequence.

Verify the exact protein translation, each antigen/linker/routing component,
codon/motif constraints, fixed elements, GenBank feature coordinates and source
identities. The GenBank splice notes reflect the actual reference-join admission
rule. Measure the final folding heuristic after all synonymous steps and report
it separately: an intermediate folding target is not assumed to survive later
codon/repeat edits. No expression, stability or global-optimality claim is made.

## Reproduction

[run.py](run.py) explicitly separates read-only planning, new Modal execution,
collection and offline replay. Each stage freezes its source/input/output
identities on the persistent volume. [replay.py](replay.py) verifies collected
receipts, repeats selection and recorded-cost assembly, checks chosen junctions
and repeats nucleotide encoding/annotation with pinned dependencies. Large
uncollected native all-option prediction files are reported as hash-only;
their absent bytes are not claimed to have been verified.
