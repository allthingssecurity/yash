# Use the data and build another approach

Design #3 is the lead: 30 regions in an 885-aa protein, using matched-normal DNA
and additional plasma/RNA candidates. Designs #1 and #2 remain available for comparison. Independent methods and formats are welcome.

## Verify the lead

With Python 3.11 or 3.12, install from the repository root:

```sh
python -m pip install -e '.[mrna]'
yash-vaccine verify
```

This read-only command checks recorded evidence, native completion receipts,
selected antigens, protein and GenBank annotations. Large native outputs absent
locally are explicitly reported as hash-only. It starts
no jobs and downloads nothing. `make verify` uses the same lead-design check.

Verification and replay refuse `python -O` / `PYTHONOPTIMIZE`, which disable
assertion-based checks. Recorded subprocess workers run without an inherited
`PYTHONOPTIMIZE` flag, addressing the verification gap raised in
[PR #2](https://github.com/popvax/yash/pull/2). Splice replay compares regenerated
gzip payloads by their decompressed contents, as proposed in
[PR #3](https://github.com/popvax/yash/pull/3); recorded files still require their
exact original hashes.

## Replay the lead locally

```sh
yash-vaccine reproduce --output work/design-3-replay
```

This reselects regions from the included candidate pool, re-solves the recorded
joint order/linker cost graph, independently recounts chosen junction predictions
and reruns synonymous encoding. Selected TSV/FASTA, protein/components/order and
nucleotide outputs must match byte for byte. It uses recorded model measurements;
it does not repeat variant calling, raw-read processing or HLA model inference.

Omit `--output` to use a temporary directory, removed after verification. Existing
output directories are never overwritten. `make reproduce` also targets the lead.

```sh
yash-vaccine reproduce designs/design-3 --output work/explicit-lead
python designs/design-3/run.py --replay --output work/scripted-lead
```

## Earlier design and comparisons

```sh
yash-vaccine reproduce designs/initial --output work/design-1-replay
yash-vaccine compare designs/initial designs/design-2 --output work/design-comparison
```

The first replays Design #1's full 29-region/852-aa protein and known nucleotide
segments. Its selection uses the maintained implementation; protein assembly
uses verified recorded source in a temporary input tree so historical code hashes
remain accurate. The second compares the two complete proteins by `region_id`, reporting
22 retained, two added and seven removed regions. `AG01` labels are local to each
selection. DNA/self evidence and prediction queries differ between these designs;
the comparison is not a controlled measurement of vaccine performance.

All three designs leave PolyA and downstream extensions unspecified. No complete
molecular RNA/template length is assigned.

## Design #3

[Design #3](../designs/design-3/README.md) is the 30-region, 885-aa Guardant-augmented
candidate, including H2BC12 and the five reference-join splice regions.
The default verify/reproduce commands above target this design.
[Explore its interactive annotations on the main site](https://rnafoundry.com/yash/#vaccine).

```sh
python designs/design-3/selection.py --verify
python designs/design-3/run.py                         # Plan only.
python designs/design-3/run.py --replay --output work/design-3-replay
```

The [methods](../designs/design-3/methods.md) distinguish paired-DNA/RNA mutations
from reference-derived Guardant hypotheses. Replay checks selection, order/linker
optimization, fresh selected-junction measurements and RNA encoding. Large native
files not bundled locally are explicitly hash-only. The `compare` helper accepts all three complete recorded designs and compatible
recipe outputs. Use stable `region_id` values when comparing selections.

## Start independent work

The chosen RNA approach keeps the earlier reconstruction/frame/native/normal-ND/
strand and broad self rules, and requires reference-join amino-acid identity
throughout each observed 25-aa window. The
[selected five genes](../data/shared-rna/splicing/junction-only/README.md) are
F5, HLA-B, HSF1, NBPF14 and KMT2D; their manual review is linked there. Follow the
[RNA screening guide](../tools/rna_evidence/README.md); verify the included
recorded-alignment screen and gate your selection:

```sh
python tools/rna_evidence/pipeline.py --verify \
  --output data/shared-rna/splicing/sequence-screen \
  --check-selection designs/my-candidate/selection.json

# Obtain score-ready windows with exact self targets removed.
python tools/rna_evidence/normal_peptides.py \
  --export work/rna-nonself-candidates.json

# Select or replay the current five-gene junction-only shortlist.
python tools/rna_evidence/normal_peptides.py --top-genes 5 --results data/shared-rna/splicing/junction-only
python tools/rna_evidence/normal_peptides.py --verify --results data/shared-rna/splicing/junction-only
```

This checks exact sequences, target bounds, supplied target lists and scores.
Self-matching targets lose their contribution before ranking; unknown queries
stop selection. The original MUC6 window stays excluded. Ordinary flanks and
intentional CTA/self targets are outside this altered-target novelty rule.
The chosen RNA shortlist is assembled in [Design #3](../designs/design-3/README.md).
H2BC12 is included with missing expression; FES remains excluded under the 25-aa
rule. No alternative mapper or pattern gate is required by this recipe.

The [candidate annotations](../designs/design-3/review/candidate-annotations/README.md)
add focal-reference affinity comparisons and fresh context-keyed processing and
presentation predictions to the agreed 30-region plan. All 167 eligible splice
windows are compared, without changing the selection or assembly objective.
See the [scoring and execution guide](../tools/candidate_annotations/README.md).

```sh
python tools/candidate_annotations/run.py          # Read-only inference plan.
python tools/candidate_annotations/run.py --verify # Offline annotation replay.
```

The [class-II annotations](../designs/design-3/review/class-II-annotations/README.md)
use CapHLA on 15–20-aa peptides, with typing and pairing uncertainties explicit.
The combined thirty-region table keeps class-I and class-II scores separate.

```sh
python tools/candidate_annotations/class_ii.py --verify
# New inference: choose a fresh output directory, then collect that same run.
python tools/candidate_annotations/class_ii.py --execute --output work/class-II-new
python tools/candidate_annotations/class_ii.py --collect --output work/class-II-new
```

```sh
yash-vaccine init my-candidate
python designs/my-candidate/run.py
```

`init` creates a README and starter script, not an analysis or vaccine. Replace
the script with your work, or document another language/workflow. You can create
the directory manually instead. No design.json, quota, HLA class, window length
or construct format is required for an independent contribution.

Describe your question, methods, inputs, dependencies, commands, results and
unfinished work in that directory. Keep useful plain-text, TSV, FASTA or other
outputs alongside it. A design can be an antigen list or analysis without a cassette.

To change the lead Design #3, start your script from its
[input guide](../designs/design-3/inputs/README.md),
[selected regions](../designs/design-3/selection/selected.tsv) and
[methods](../designs/design-3/methods.md). The recorded lead replay checks a fixed
design; changing its inputs is not a supported way to create another design.
Document any new predictions or assembly checks your changes require.

The optional [Design #1 recipe helper](recipe-runner.md) varies the earlier
tumour-only selection and internal cassette. `init --from` does not copy the
matched-normal lead or add its SEC/MITD and RNA stages.

## Download source data

```sh
yash-vaccine data list
yash-vaccine data list --dataset tumour-rna
python tools/r2/restore.py --destination ../yash-data --dataset tumour-rna
```

The last command prints an offline plan. Add `--execute` to download, or use
`--file` for catalog file IDs (`--key` still accepts original keys). Public downloads need no credentials. The endpoint and availability are
recorded in data/catalog/access.json. The raw category includes the large
matched-normal BAM; inspect the byte total. Original normal FASTQs are excluded.

[Data](../data/README.md), [normal-data scope](../data/README.md#sequencing-data)
and [download instructions](../tools/r2/README.md) explain files, hashes and access.
New discovery/model runs have separate [runtime contracts](upstream.md).
Their third-party dependencies are obtained directly from
[upstream providers](external-resources.md), not this project's data endpoint.
Verification reports three upstream source files as hash references when they
have not been obtained and byte-checked locally; offline replay does not run them.

## Contribute

Fork the repository, create a branch and add your work under `designs/`. Include
antigens/rationale, sequences where produced, supporting results, code, dependencies
and reproduction commands. Label unspecified sequence components explicitly.
When comparing with our approach, use the lead Design #3 unless studying an earlier design.
See [CONTRIBUTING](../CONTRIBUTING.md).
