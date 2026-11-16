# Start here

Read README.md, then data/README.md for datasets or docs/usage.md for code. This is a research handoff: shared evidence,
our runnable implementation and vaccine designs. **Design #3 is the lead** in
designs/design-3/ and uses matched-normal DNA. Design #1 in designs/initial/ was
generated without matched-normal DNA and is retained for comparison. Independent methods,
antigen selections, scoring models and construct designs are welcome.

## Make a new design

- `yash-vaccine init NAME` creates designs/NAME/README.md and a starter run.py.
- Use any code, tools, inputs and output formats appropriate to your question.
  No design.json, fixed quotas, region length or prescribed scoring method is
  required for an independent design. Reuse or change src/ helpers when useful.
- Put the question, dependencies, data sources, commands and findings in the
  design's README. Keep useful outputs alongside it, and large scratch files in work/.
- A design can be an antigen list, an analysis or a construct. A cassette is not
  required if the work has not reached assembly. Explain what was actually done.

## Lead design and earlier implementation

- Start with designs/design-3/README.md and methods.md for the lead design.
  Verify collected results with `yash-vaccine verify`;
  replay selection, recorded-cost assembly, junction checks and RNA with
  `yash-vaccine reproduce --output work/design-3-replay`. These default to the lead.
  These are offline checks. Raw-read/prediction reruns require the documented
  Modal inputs and explicit --execute. Keep recorded evidence and checksums intact. Historical project source/docs are archived
  under provenance/snapshots/. Three full upstream RNA modules are not bundled;
  their expected hashes and minimal provider-source patches are explicit.
  Verification reports hash-only external records separately; never claim their
  absent bytes were checked. Maintained-code changes are checked by replay.
- The dated matched-normal-20261107 identifiers in Modal/work locations and
  historical execution receipts identify the actual run; the current repository
  directory is designs/design-2/. Do not rename recorded remote data or rewrite
  native receipts just to change the design's display name.

## Design #3

- `designs/design-3/` is the only Design #3 directory: 30 regions and an 885-aa
  protein. The canonical selection, construct, policy, methods and run/replay
  entry points are directly in that directory. Design #3 is the default lead.
- Thirteen mutation regions have paired-DNA/RNA evidence, including H2BC12 with
  expression left missing. Five Guardant additions are reference-derived VUS
  hypotheses, not full-window RNA reconstructions. Retain AR's normal-coverage
  exception and the in-frame fusion's binding exception.
- RNA selection uses `data/shared-rna/splicing/junction-only/selected.json`:
  F5, HLA-B, HSF1, NBPF14 and KMT2D. Preserve the original alignment, support,
  frame/native/normal-ND/strand rules and broad self exclusion. Every 25-aa window
  must match reference bases translated across the same join on every recorded
  supporting path; synonymous differences are allowed. No additional mapper,
  distinct-pattern or contextual-error gate is part of this recipe.
- Verify with `python designs/design-3/selection.py --verify`; replay with
  `python designs/design-3/run.py --replay --output work/design-3-replay` using a
  fresh output. Keep source measurements, selected sequence bytes and native
  receipts intact. Native paths are resolved to current inputs only on read.
- `tools/rna_evidence/normal_peptides.py` gates RNA selection. Intentional CTA/self
  regions are outside the altered-peptide novelty gate. Unknown targets or stale
  scores stop selection. All 167 eligible windows remain available for analysis.
- `tools/candidate_annotations/` records focal-reference affinity ratios and
  context-keyed class-I presentation; `class_ii.py` adds CapHLA 15–20-aa annotations.
  Both have offline `--verify` modes and explicit fresh-run modes. Never reuse
  a processing score across different peptide/flank contexts. Annotation scores
  do not select regions or change the assembly objective.
- Class II uses two primary DRB1 targets and eight exploratory unphased DQ/DP
  combinations. DRB3 is excluded as requested; typing evidence is retained.
  Equivalent DRA binding domains map four chain pairings to two model targets.
  Do not claim measured phase, binding-core novelty, a complete long-peptide
  personal-normal exclusion or a T-cell response. `data/hla/class-II.json` records
  T1K calls and separately attributed PR #7 normal-WGS corroboration.
- The interactive view lives on the main website at https://rnafoundry.com/yash/#vaccine,
  not in this research repo. `tools/candidate_annotations/site_data.py --output PATH`
  exports verified annotations for that site; add `--verify` to check freshness.
  Keep website HTML/CSS/JS in the RNA Foundry website repository.
- Keep private reviews, superseded builds and unused alternatives under ignored
  `work/`. Do not add them to the catalog, documentation links or future commits.
  Preserve Designs #1/#2 and the 145 published data objects unchanged.

## Design #1 and optional helpers

- `yash-vaccine reproduce designs/initial` checks the complete 29-region, 852-aa protein and
  known RNA/DNA segments before unspecified PolyA using included evidence;
  install `.[mrna]`. Selection, joint SEC/MITD assembly, junction checks and
  RNA encoding are stages of designs/initial/.
  Keep its outputs and the original evidence as references.
- Splice selection requires every base of the full 75-nt window in actual patient
  RNA, with at least two qualifying fragment names, consistent overlapping mates
  and a stop-free observed path from a verified GENCODE 37 upstream CDS anchor.
  It starts before HLA filtering, ranks all qualifying windows, uses one region
  per splice gene and applies no identity filter. data/shared-rna/splicing/ records
  the complete pool, native QC, fresh affinities and independent BAM verification.
- Our recipe applies `kyte_doolittle_19_gt_1_6` before source quotas and rejects
  hydrophobic internal joins before optimization. Every unmasked full-chain score
  is exported; fixed SEC/MITD overlaps are separately reported, not topology-validated.
  See designs/initial/review/hydrophobicity/README.md. Preserve original DNA/RNA evidence.
- Methods specific to a design belong in that design's directory. Design #1's
  scientific decisions are in designs/initial/methods.md; Design #2's are in
  designs/design-2/methods.md. Keep docs/ for shared
  usage and workflow guides.
- Python 3.11/3.12 and `pip install -e '.[mrna]'` are supported for our implementation.
  Independent approaches can document their own environments.
- To vary Design #1's recipe, use `yash-vaccine init NAME --from designs/initial/selection`.
  This helper copies the selection/internal-cassette stage, before SEC/MITD and RNA.
  It retains our case context, fixed 25-aa inputs and class-I objective;
  see docs/recipe-runner.md. These are limitations of that helper, not rules for
  other designs. Do not silently bypass its evidence/cache checks.
  This helper uses Design #1's tumour-only evidence, not the matched-normal lead.
  For changes to Design #3, create a plain workspace with `init NAME`, start from
  its selection/selected.tsv and methods.md, plus the shared and Design #2 source evidence,
  and implement/document the changed analysis and any new predictions yourself.
- `compare` reads the three complete recorded designs or compatible recipe-runner outputs. Independent analyses can use
  their own comparisons. When reusing our regions, compare region_id values;
  AG01-style labels are local to one selection.

## Evidence and execution

- data/README.md explains shared files and public read-only downloads;
  data/catalog/access.json records the endpoint. Keep the bucket private and
  preserve the exact patient/results allowlist. Deployment or access changes
  require an explicit user instruction. Models, software
  and general references come from upstream providers via docs/external-resources.md. List inputs with
  `yash-vaccine data list`; restore defaults to a plan. Keep the bucket private.
  Do not reintroduce withdrawn third-party mirrors or full source in snapshots.
- Preserve original measurements, flags and source identities in shared evidence.
  New models, HLA types or expression policies should use appropriate inputs and
  be documented as new analysis, rather than relabelling cached predictions.
- Keep assumptions and the methods used explicit. Separate local replay,
  new inference, raw-read processing and experimental results in reporting.
- workflow/stages contains source-specific upstream workers. Read docs/upstream.md
  before use; some workers execute at import. Do not bulk-import them.
- Remote prediction and upstream helpers require --execute. Missing cache entries
  stop our recipe helper; never silently drop queries or claim new inference.
- Run existing checks appropriate to changes. Do not add tests or CI unless the
  user changes that instruction. tools/check_reproduction.py checks the reference;
  tools/verify_fresh_clone.py exercises independent and optional-helper paths.
- Keep credentials outside Git/logs. Maintain the recorded candidate and frozen
  evidence while developing alternatives on a branch.

## Data catalog maintenance

- `data/files.tsv` is the maintained index of repository data and public downloads.
  Run `python tools/catalog.py` to validate it, or add `--write` to regenerate
  marked README tables and local download JSON. These commands never publish.
- Use `--dataset` or `--file` with `yash-vaccine data list/restore`; restore is
  still a plan until `--execute`. Reports and HLA summaries are bundled files.
- Original paths inside immutable records are resolved by
  `src/yash_vaccine/paths.py` and `provenance/path-map.json`. Keep measurement
  bytes and native receipts intact. Do not recreate duplicate old directories.
