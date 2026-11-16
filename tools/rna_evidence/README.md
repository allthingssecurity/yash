# Screen RNA candidates before selection

The active selection is **F5, HLA-B, HSF1, NBPF14 and KMT2D**, assembled in
[Design #3](../../designs/design-3/README.md). The [methods](../../designs/design-3/methods.md#current-rna-selection-rules)
explain the reconstruction, frame, native review, coordinate-normal ND, strand,
expression-weighted ranking and broad exact-self rules.

Each observed 25-aa window must equal reference bases translated across the
same proposed join on every recorded supporting path. Synonymous nucleotide
differences are allowed; additional amino-acid substitutions fail this rule.
The resulting pool has 167 eligible windows across 18 events and 10 genes.

```sh
python tools/rna_evidence/normal_peptides.py --verify --results data/shared-rna/splicing/junction-only
```

## Mandatory normal-peptide exclusion

[normal_peptides.py](normal_peptides.py) checks every target-overlapping 8–11-aa
peptide against all Ensembl 103 and 115 translations, UniProt canonical and
additional isoforms, and the bounded patient-normal coding/HLA models. This
includes alternative contigs and pseudogene translations, which caught MUC6.
Catalog membership is checked by exact amino-acid sequence, across all genes.

Remove matching peptides from the affinity sum, recompute the expression-weighted
score and binding admission, then rank the eligible windows. A region with no
remaining qualifying nonself binder is excluded. The assembly gate checks exact
window identity, supplied target lists and scores, preventing a stale list or
score from restoring a removed peptide. Missing search decisions stop the gate.
The original MUC6 region fails the self search. Both substitution-dependent
MUC6 proposals also fail the current reference-join identity rule.

```sh
python tools/rna_evidence/normal_peptides.py \
  --selection data/shared-rna/splicing/junction-only/selected.json

# Export score-ready windows for another selection; no quota is imposed here.
python tools/rna_evidence/normal_peptides.py \
  --export work/rna-nonself-candidates.json
```

The default gate and `pipeline.py --check-selection` apply the new identity rule.
Recorded Design #3's controller/direct replay explicitly retain their historical
substitution-permitting admission while requiring broad self exclusion; the CLI
equivalent is `--include-substitution-dependent`. It is not the default for new
selections. The gate verifies and reuses completed protein-search
measurements, including the expanded exact-protein search;
it does not apply that assessment's mapper, pattern or contextual-error gates.
Use `--screen` and `--proteins` for compatible new measured outputs. Unsearched
peptides need a new measured catalog search before they can enter a selection.

This excludes **self-matching targets from scoring**, not every ordinary peptide
in the surrounding 25-aa region. Intentional CTA/self targets remain outside
this novelty rule. A negative catalog search is not evidence that a peptide is
absent from every healthy tissue or is presented on tumour cells.

## What it does

1. **Keep the existing sequence evidence.** Start with actual 75-nt RNA windows
   already passing frame/read reconstruction, native peptide review,
   coordinate-aware normal-junction ND and the recorded strand checks. Check
   every alignment path recorded for each window.
2. **Identify every additional base change.** Translate the reference bases along
   each recorded path and compare them with the observed RNA sequence. Record all
   differences, including synonymous changes, rather than assuming every amino-acid
   difference came from the splice join.
3. **Read matched-normal evidence.** Check nominated differences against normal
   genotypes and qualifying fragments. An independent direct counter checks ALT
   counts. Record tumour-DNA and RNA counts separately. Missing or conflicting
   evidence stays unresolved.
4. **Search for self peptides.** Search target-overlapping 8–11-aa peptides in all
   Ensembl 103 translations, UniProt canonical/additional isoforms and the existing
   bounded patient-normal coding/HLA models. Exclude exact matches from scoring.
   This includes alternative-contig and pseudogene translations, not just the
   primary coding-reference subset.
5. **Separate the reason for novelty.** Ask whether the same alignment path,
   translated with reference bases, already gives a nonself peptide. If not,
   novelty depends on extra RNA amino-acid changes and needs the explicit normal
   evidence rule below. All recorded paths must support primary admission.
6. **Rescore and gate selection.** Verify each original native Vaxrank score,
   remove excluded/review-only contributions from the primary score, then rank
   the best window per event and gene. Admission retains a native-PASS target
   peptide ≤500 nM, positive supplied expression and the existing hydropathy rule.
   Scores use affinity sum × √supplied gene value; no read-abundance or splice-fraction
   multiplier. Require an exact primary checked window and the mandatory broader
   normal-peptide gate before using the score or assembling a selection.

The [policy](pool-policy.json) pins the rules. For change-dependent novelty,
primary admission needs resolved normal evidence for the peptide-overlapping
amino-acid changes and at least one reference-supported change. A known normal
background does not invalidate a separate RNA change. Reference support requires
at least 20 usable normal fragments, native GQ≥20 and fewer than three ALT
fragments; MQ/BQ≥20 apply to the normal read evidence. These are explicit research
rules, not proof of RNA editing or normal-tissue absence.

## Replay the included results offline

From the repository root, using Python 3.11/3.12:

```sh
python tools/rna_evidence/pipeline.py --verify \
  --output data/shared-rna/splicing/sequence-screen \
  --check-selection designs/design-3/selection/selected.json
```

No Modal access, BAM download, model weights or additional packages are needed for
this replay. It checks output hashes, independently recomputes normal-evidence
decisions, native affinity contributions, peptide tiers and window/gene rankings,
and regenerates the readable table and checked export. It does not reread BAMs or
repeat the protein-catalog search.

`admitted.json` contains only primary-tier windows. The selection gate checks the
exact window ID, region ID and amino-acid sequence, plus coding sequence/target
bounds when supplied. Unassessed, changed, review-only or excluded RNA regions
stop the gate. Other source types are explicitly outside its scope. Design #3’s
selector and replay run the reference-join and broad self gate; Design #2’s recorded outputs
remain unchanged.

## Run a new measured screen on Modal

The launcher prints a plan by default. Use a new run ID and output directory:

```sh
python tools/rna_evidence/pipeline.py \
  --run-id rna-screen-my-run --output work/rna-screen-my-run
python tools/rna_evidence/pipeline.py --execute \
  --run-id rna-screen-my-run --output work/rna-screen-my-run
python tools/rna_evidence/pipeline.py --collect \
  --run-id rna-screen-my-run --output work/rna-screen-my-run
```

`--execute` starts one CPU worker and returns its handle. `--collect` polls that
handle, then downloads and verifies completed outputs; it never starts another
job. Repeat collection if still running. Completed native measurements are not
overwritten. Add `--verify --output work/rna-screen-my-run` for offline replay.

Use `--candidates /path/to/windows.tsv` to screen a subset of the eligible pool;
the TSV needs a `window_id` column. A selection JSON with `source: "SPLICE"` and
`window_id` is also accepted. Newly discovered windows first need the documented
[upstream reconstruction and QC stages](../../data/shared-rna/splicing/README.md).
This launcher is not an independent RNA variant-discovery caller.

New measured execution uses authenticated Modal access, app
`research-yash-matched-normal`, volume `research-lens-workspace`, CPU image
`im-hX1h6HknDMvasm9Z6q530H` and the existing completed workspace
`/data/results/TN25-112103/matched-normal-20261107`. It needs the tumour-DNA/RNA
BAMs, indexed normal BAM, native coding GVCF, sample-identity receipt, compatible
genome/PyEnsembl 103 reference and scoped personal-normal models staged there.
These are the recorded case runtime contracts, not a portable fresh-cloud setup.
The normal BAM’s recorded hash and file size identify it; this stage does not
rehash all 143 GB. Reference catalogs are obtained from providers and hash checked;
no reference sequences or model weights are bundled.

The run reuses **12,204** checksum-pinned peptide/HLA affinity predictions; it
performs no fresh HLA inference. Its worker calls native Vaxrank 3.46.0 to recover
contributions and checks all original scores before exclusions. Package versions,
source hashes and native outputs are in the completion receipt and summary.

Code: [driver](pipeline.py), [measured worker](pool_worker.py),
[offline replay and selection gate](pool_report.py). The worker also uses the
unchanged [focal-evidence helpers](worker.py) and Design #2’s normal/self helpers.

---

# Recorded selected-window audit

This read-only audit separates an observed RNA nucleotide sequence, its inferred
amino-acid changes, the recorded splice alignment and its normal-sequence matches.
It does not alter raw data, call an RNA-editing mechanism or measure proteins.

The [results and decisions](../../designs/design-3/inputs/sequence-review/README.md)
cover focal normal-DNA checks for F5, MUC5B, MUC3A and MUC6, plus reference searches
for all 24 altered-source regions in the audit pool. The pool includes MUC6,
which is excluded from Design #3 after this check. The five intentional CTA/self
regions are outside the altered-peptide novelty requirement.

## Verify the collected evidence

From the repository root, with the project's Python environment:

```sh
python tools/rna_evidence/report.py
```

This verifies collected output hashes, source-window identities and the recorded
normal/reference decisions, then regenerates the readable summary. It does not
reread a BAM or make a new inference. To repeat the reference sequence search,
obtain the Ensembl 103 protein FASTA from the
[provider guide](../../docs/external-resources.md), then run:

```sh
python tools/rna_evidence/report.py \
  --reference-proteins /path/to/Homo_sapiens.GRCh38.pep.all.fa.gz
```

All reference queries are checked against the primary protein-coding subset and
all provider translations. Exact matches on alternative contigs and pseudogene
translations are conservatively excluded; their presence is not proof of normal
tissue expression. Patient-normal coding/HLA-model decisions remain recorded
measurements in this local check, not an exhaustive normal-proteome reconstruction.

## Repeat the normal-DNA and personal-self checks

This needs authenticated access to the existing Modal reference/data volume and
the documented analysis image; `modal` is required locally. The source BAMs,
reference genome and native normal genotypes are already staged there. No new
alignment or mutation discovery is launched.

```sh
python tools/rna_evidence/run.py
python tools/rna_evidence/run.py --execute \
  --run-id rna-normal-self-rerun --output work/rna-normal-self-rerun
python tools/rna_evidence/run.py --collect \
  --run-id rna-normal-self-rerun --output work/rna-normal-self-rerun
```

The first command prints a plan. `--execute` starts a worker and records its
handle; `--collect` polls that same handle and verifies completed output hashes.
Use a new run ID and output directory for a fresh rerun. Completed measurements
are never overwritten. The frozen candidate pool, RNA windows and provider
FASTA checksum define the inputs. Every full-window nucleotide difference is
examined, including synonymous changes; direct pileups and the existing Isovar
collector independently check normal alternate counts. The original candidate
pool remains in the results so MUC6's exclusion can be reproduced.

The worker source is [worker.py](worker.py). Normal genotype interpretation uses
[the existing helper](../../designs/design-2/pipeline/normal_evidence.py), and
personal-self search uses [the existing scoped models](../../designs/design-2/pipeline/self_search.py).

## Check the second protein catalog

The supplementary [UniProt result](../../designs/design-3/inputs/sequence-review/uniprot-reference-check.json)
uses 286,017 human reference-proteome sequences, including canonical and additional
isoforms. The provider files were last modified 4 November 2026. Download the
[canonical FASTA](https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Eukaryota/UP000005640/UP000005640_9606.fasta.gz)
and [additional FASTA](https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Eukaryota/UP000005640/UP000005640_9606_additional.fasta.gz),
then concatenate the two gzip files in that order. Python's gzip reader supports
the resulting concatenated members.

```sh
cat UP000005640_9606.fasta.gz UP000005640_9606_additional.fasta.gz > human-with-isoforms.fasta.gz
python tools/rna_evidence/check_uniprot.py \
  --reference-proteins human-with-isoforms.fasta.gz \
  --download-receipt designs/design-3/inputs/sequence-review/uniprot-reference-check.json
```

The recorded component and combined hashes must match; acquire the corresponding
archived provider release if `current_release` has changed. No reference sequences
or weights are redistributed. [validate.py](validate.py) checks the collected
decisions quickly without downloading references, re-reading BAMs or running
models. The design controller applies this check before replay or new execution;
newly matching target peptides stop it for renewed selection review.
