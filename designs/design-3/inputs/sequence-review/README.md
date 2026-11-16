# RNA-supported amino-acid changes and normal-sequence screening

The RNA reads contain nucleotide sequences predicted to encode amino-acid changes relative to the nominated reference. Splice analysis discovered these candidates, but repeats leave their exact genomic placement and splice origin uncertain. The bases are useful sequence evidence; they are not direct protein measurements or confirmed RNA-editing events.

## Matched-normal DNA

Fresh focal-base pileups and native normal genotypes were checked at every nucleotide difference in the four selected RNA windows. Reference and alternate counts use fragment names, not UMI-confirmed independent molecules. Both the Isovar collector and a direct alignment-position implementation agree on alternate counts.

| Candidate | Reference → RNA-inferred amino acid | Reference → observed RNA codon | Normal ALT / REF fragments | Native normal GQ |
|---|---|---|---:|---:|
| F5 | S→F | `UCU→UUU` | 0 / 50 | 99 |
| MUC5B | S→I | `AGC→AUC` | 0 / 52 | 99 |
| MUC3A | I→F | `AUC→UUC` | 0 / 74 | 99 |
| MUC6 | T→M | `ACG→AUG` | 0 / 150 | 99 |

No alternate fragments were detected at these four positions. This is reassuring at the nominated loci; it does not establish the origin of a repetitive read, exclude every inherited repeat copy or prove tumour-specific translation. Normal DNA also cannot establish whether an unusual RNA join occurs in healthy tissue.

The full-window check additionally records MUC6 R→P at the first residue: 24 alternate normal fragments, but a reference native genotype, giving unresolved evidence rather than a confident inherited call. A third MUC6 nucleotide difference is synonymous and present in normal DNA. Neither changes the separate T→M result above.

## Broader reference and patient-self search

All 8–11-aa peptides in the 24 altered-source regions of the pre-filter pool were searched against all **115,262 Ensembl 103 peptide translations**, the **85,968 primary protein-coding subset**, and the recorded bounded patient-normal coding/HLA models. Target counts below use the existing target intervals. Intentional CTA/self regions are outside this novelty requirement.

| RNA candidate | Primary-nonself target peptides | Remain after broad-reference + personal-self exclusion | Decision |
|---|---:|---:|---|
| F5 | 26 | 26 | Retain as an RNA-sequence hypothesis |
| MUC5B | 22 | 22 | Retain as an RNA-sequence hypothesis |
| MUC3A | 18 | 18 | Retain as an RNA-sequence hypothesis |
| MUC6 | 23 | 0 | Exclude from altered-target selection |

MUC6’s entire selected 25-aa region matches **ENSP00000487059.1 / ENST00000627256.2**, on alternative contig **CHR_HSCHR11_2_CTG1**, annotated **polymorphic_pseudogene**. All 35 target-overlapping peptides match this broader reference; all 23 that were nonself against the primary coding subset therefore lose novelty under the expanded rule. Design #3 excludes this region.

This is a predicted alternative human translation, not a measurement that a healthy tissue expresses or displays it. The conservative policy excludes that exact sequence match even though the focal normal-DNA test supports the reference allele.

No additional target-peptide exclusions appeared in the other 23 altered-source regions relative to the primary-reference plus personal-self search. No negative search establishes absence from every normal protein: the reference release, predicted translations and bounded patient-normal models define the scope. No healthy-tissue immunopeptidome search or new HLA-binding inference was performed by this audit.

## Second protein catalog: UniProt

A supplementary exact search covered **286,017 sequences** from UniProt's human reference proteome **UP000005640**, including canonical and additional/isoform FASTA files last modified 4 November 2026. All 1,584 queried peptides were searched. It adds no target-peptide exclusions beyond the Ensembl/patient-normal result; the retained F5, MUC5B and MUC3A target sets are unchanged.

MUC6's match is in the extended Ensembl inventory; neither its full window nor `SMAKTSTSL` matched this UniProt download. Its exclusion is a conservative response to a predicted alternative human translation, not a claim of measured healthy-tissue protein or a reviewed UniProt match.

[Query-level decisions, provider URLs and file hashes](uniprot-reference-check.json) retain the complete search. Reference FASTA files come from [UniProt](https://www.uniprot.org/help/proteome), not from this repository.

## Files and reruns

- [Summary table](summary.tsv) · [all altered-region decisions](summary.json).
- [Normal-DNA genotypes and fresh normal/tumour/RNA base counts](normal-dna.json).
- [RNA windows, codons and alignment paths](windows.json) · [reference witnesses](reference-matches.json).
- [Reference/sample scope](reference-scope.json) · [audit request](request.json) · [candidate pool before exclusion](candidate-pool.json) · [execution receipt](complete.json).
- [Runnable audit tools](../../../../tools/rna_evidence/README.md).
