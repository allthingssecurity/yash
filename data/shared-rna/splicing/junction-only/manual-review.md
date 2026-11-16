# Manual review of the five junction-derived RNA regions

Reviewed 15 November 2026 using the original supplied RNA alignment and recorded
reconstruction, support, reading-frame, native review and strand rules. The
current selection is [selected.json](selected.json). This review adds no
remapping, distinct-pattern or contextual-error requirements.

The added sequence rule remains: the entire observed **25-aa sequence** must
match reference bases translated across the **same proposed join on every
originally recorded supporting path**. Synonymous nucleotide changes are allowed.
The broad reference/personal-self peptide exclusion also remains mandatory.

There are **167 passing windows across 18 events and 10 genes**. Existing weighted
scores and event/region tie-breaks select one window from each of the top five
genes. This is an RNA shortlist; no new vaccine construct has been assembled.

## Original-alignment evidence and normal-peptide check

| Gene | Selected 25-aa region | Best predicted target | HLA | Affinity, nM | Original full-window fragment names | Expected/opposing names | Scored target peptides | Exact self matches among scored targets |
|---|---|---|---|---:|---:|---:|---:|---:|
| F5 | `HTTLSLDFSQTNLSPELSHTTLSLD` | `LSPELSHTTL` | C*01:02 | 23.58 | 2 | 2 / 0 | 6 | 0 |
| HLA-B | `EAAREAEQRSLPGGHVRGVAPQIPG` | `EQRSLPGGHVR` | A*33:03 | 66.16 | 3 | 3 / 0 | 8 | 0 |
| HSF1 | `FSVDTSALLDVSGAPPPRLPAPPPR` | `VSGAPPPRL` | C*01:02 | 44.45 | 2 | 2 / 0 | 19 | 0 |
| NBPF14 | `QEQLAEGCRLAQHLVQKLSPGKVAI` | `KLSPGKVAI` | C*01:02 | 44.53 | 6 | 4 / 2 | 2 | 0 |
| KMT2D | `HLSPQPEEPHLSPQAEGPHLSPQPE` | `LSPQPEEPHL` | C*01:02 | 85.57 | 3 | 3 / 0 | 4 | 0 |

**All five pass the original rules plus the new sequence-identity rule and broad
self exclusion.** Across the five regions, all **39 distinct scored 8–11-aa
target peptides** have no exact match in the searched Ensembl 103/115, UniProt
canonical/additional isoforms or bounded patient-normal coding/HLA models.
Alternative contigs and pseudogene translations are included, as in the check
that excluded the original MUC6 region. Matching peptides elsewhere in the
windows are excluded from scoring; ordinary self flanks can remain in a 25-aa
region.

The nominated junction coordinates are also ND in both processed normal
junction catalogs. That is a separate coordinate-level check. These searches
do not measure presentation in every healthy tissue.

## Individual review

**F5.** Window `RNA_SPLICE_de3577be74252545fecc`, region
`REG_6a492ad25bc739969e40`, coordinates `chr1:169541117-169541199(-)`.
The target crosses the new adjacency as `LSPE | LSHTTL`. Both original fragment
names support the complete framed window in the expected orientation. Observed
and reference-join translations are identical, with no additional nucleotide
changes on the recorded path. Native review is PASS/100. Score: 157.61;
supplied gene value: 470.59834. It passes with two full-window supporting names.

**HLA-B.** Window `RNA_SPLICE_de2328085e66ea5701c0`, region
`REG_43c17f5844fb05456835`, coordinates `chr6:31271142-31356237(-)`.
Three original names support the complete framed window in the expected
orientation. Its nonself target is the joined tail, including `EQRSLPGGHVR`;
ordinary self-like sequence is not counted as novel. The source gene is HLA-B,
while the best predicted presenting allele is A*33:03. Native review is PASS/100.
Score: 44.19; supplied gene value: 185.95777.

**HSF1.** Window `RNA_SPLICE_dfffafab89cf697bf26f`, region
`REG_463147fb26e0938a4cc4`, coordinates `chr8:144313640-144313717(+)`.
The two original supporting names contain the complete framed window and agree
with the expected orientation. Reference-join and observed translations are
identical; the scored junction target remains nonself. Native review is
PASS/100. Score: 39.63; supplied gene value: 21.84034. It passes with two
full-window supporting names.

**NBPF14.** Window `RNA_SPLICE_63ad26a57e2e2b74ce43`, region
`REG_0a914a772ba5bd2a7738`, coordinates `chr1:148586249-148592547(-)`.
Six original names support the complete sequence: four expected and two
opposing. This passes the original strict-majority strand rule; the opposing
minority is retained in the evidence. Reference-join and observed translations
are identical and scored targets pass the self search. Native review is
PASS/100. Score: 33.61; supplied gene value: 133.82139.

**KMT2D.** Window `RNA_SPLICE_335695e305c0bc6b7cfb`, region
`REG_43525196a04b74490f29`, coordinates `chr12:49051335-49051363(-)`.
Three original names support the complete framed sequence in the expected
orientation. The original reference-join translation equals the observed
sequence. Native review remains PASS/85 with its existing `READ_END_BIAS` flag;
that flag is not hidden or turned into a new exclusion. Scored targets pass the
self search. Score: 25.33; supplied gene value: 40.50225.

All five meet the native-PASS <=500-nM target requirement, have positive supplied
expression and pass the original internal-region hydropathy rule. The original
motif-direction annotations remain report-only, as before. Fragment names are
not UMI-confirmed independent molecule counts. These passing labels establish
computational eligibility under the chosen rules, not measured tumour protein,
HLA presentation or T-cell recognition.

## Sources and reproduction

- [Original RNA window/path witnesses](original-rna-windows.jsonl.gz) and [review evidence](manual-evidence.json).
- [Reference-identity filter, scores and commands](README.md), [selected TSV](selected.tsv), [all ten genes](gene-ranking.tsv), [selection verification](verification.json).
- [Original complete-window BAM verification](../read-verification/verification.json).
- [Original coordinate-normal and strand checks](../normal-screen/README.md).
- [Reference-path translations](../sequence-screen/window-evidence.json.gz) and [broad reference/personal-self decisions](../sequence-screen/reference-matches.json.gz).
