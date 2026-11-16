# Class II annotations with CapHLA

Class-II binding and presentation annotations for the **agreed 30-region plan**
and the rest of the currently eligible splice pool. The 30 regions comprise
18 mutation, five splice, five CTA/self and two fusion regions. No region,
selection rule or recorded construct is changed.

The completed run scored **2,567 distinct target peptides**, plus 857 additional
focal-reference counterpart queries, against ten class-II targets:
**34,240 peptide/target pairs**, each with separate five-fold BA and EL results.
All 192 window/region records and all thirty planned regions have complete
length and target accounting. Fresh inference used an NVIDIA A10 for 14.05 seconds; reference
searches followed. Provider controls, model identities and numerical checks
passed. Outputs remain on the persistent volume, and compute was terminated.

## Explore the results

| File | Contents |
|---|---|
| [Current plan](current-plan.tsv) | Thirty regions, joining the existing class-I annotations to the new class-II summaries. |
| [All region annotations](region-annotations.tsv) | The planned regions and all 167 eligible splice windows: 192 window/region records, 191 distinct region IDs. |
| [Region by class-II target](region-target-annotations.tsv) | Separate best binding and presentation peptides for every region and each of the ten targets. |
| [Peptide annotations](peptide-annotations.tsv.gz) | Every 15–20-aa target occurrence, its model scores, exact-search flags and focal-reference counterpart where defined. |
| [HLA/model coverage](hla-coverage.tsv) | Requested and actual model identities, pseudosequence coverage, typing confidence and pairing assumptions. |
| [Long-peptide reference matches](reference-matches.json.gz) · [search scope](normal-search.json) | Fresh exact searches of all four reference-protein catalogs and the separately scoped personal-normal models. |
| [Summary](summary.json) · [length accounting](peptide-length-coverage.tsv) | Complete prediction and occurrence counts; no threshold or selection change. |
| [Measured predictions](predictions.tsv.gz) · [inference](inference.json) · [verification](verification.json) | Fresh BA/EL predictions, pinned model/runtime identities and offline replay. |

## Read the columns

The `classII_DRB1_` columns summarize **DRB1*03:01 and DRB1*08:02**. The separate
`classII_exploratory_` columns include the eight observed DQ/DP alpha-by-beta
combinations. DRB3 is excluded from active annotation; its typing evidence is
retained. Pairings are hypotheses, not measured cis/trans assignments or ten
confirmed functional molecules.

Each splice row is keyed by its original window ID (`annotation_id`). Two TSPY2
records share a region ID but are scored and summarized separately; their
peptide occurrences are not pooled or counted twice within a window.

`best_affinity_nm` is the strongest predicted binding; `best_presentation_score`
is the highest CapHLA EL score. These may identify different peptides or
restrictions. Higher EL and lower nM are favorable in their respective models,
but the EL number is not a patient-specific probability or a class-I-equivalent
score. No percentile, hard cutoff, expression weighting or combined class-I/II
rank is supplied.

All long peptides are searched against the pinned Ensembl 103/115 translations
and UniProt canonical/additional proteins, including alternate sequences. Raw
best predictions are retained even when an exact reference hit exists. Separate
`best_no_exact_search_hit_` columns summarize peptides without a hit in the
searched reference/personal catalogs. Intentional reference CTA/self regions
remain visible rather than being removed.

**The existing personal-normal library was constructed for ≤11-aa peptides.**
It is searched for long-peptide hits, but its local haplotype clusters and
flanks do not guarantee complete 15–20-aa coverage. Positive hits remain useful;
negative searches are explicitly incomplete for class II. None of these searches
measures healthy-cell translation or presentation.

**CapHLA does not assign a binding core here.** A long peptide can touch a
mutation or join while its actual HLA-II core remains unchanged. Full-peptide
novelty and predicted binding therefore do not establish mutation-specific CD4
recognition. Focal-reference comparisons are supplied for missense candidates;
they restore only that residue and retain other background differences.

## Source and replay

[Typing evidence](../../../../data/hla/class-II.json) records native tumour-DNA/RNA
T1K calls and separately attributed matched-normal corroboration from
[PR #7](https://github.com/popvax/yash/pull/7).
The contributed caller was not rerun. Both observed DRA proteins pair with DR
beta chains. Their peptide-binding domains are identical: the complete reference
representatives differ only at precursor residue 242 (mature residue 217), in
the cytoplasmic tail. The four possible DRA/DRB1 combinations therefore map to
two distinct CapHLA predictions, one per DRB1 allele. This model equivalence does
not establish equal expression, stability or trafficking. The alignment source
and four combinations are recorded in the typing file; model assumptions and
the X token in upstream DQ pseudosequences are explicit in `hla-coverage.tsv`.

The MIT-licensed [CapHLA code and weights](https://github.com/changyunjian/CapHLA/tree/33ebdd6ce6dadbbb1c66b026ce4b5d81dbf3a831)
remain external. Five BA and five EL folds are used through the existing
OpenVax/mhctools adapter. Every pinned provider file and the adapter are checked;
provider example predictions must reproduce before patient scoring.
Scientific reference: [Chang and Wu, CapHLA](https://doi.org/10.1093/bib/bbae595).

The [Design #3 build](../../README.md) assembles these thirty regions.
These annotations use the individual regions, with no sequence extension or invented flanks; no full
vaccine-cassette class-II processing prediction is claimed.

```sh
python tools/candidate_annotations/class_ii.py --verify
```

The [code and execution guide](../../../../tools/candidate_annotations/README.md#class-ii-with-caphla)
explains fresh inference, persistent outputs and collection. Verification replays
the tables from included predictions without downloading models or starting jobs.
