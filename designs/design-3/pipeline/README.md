# Joint assembly and RNA encoding

The entry point is [run.py](../run.py). It uploads a checksum-pinned input bundle
to a separate persistent Modal run root. Scientific inputs are fixed by
[policy.json](../policy.json); [selection.py](../selection.py) verifies all thirty
regions against the current annotated plan and replays the splice admission gate.

1. **Prepare.** Enumerate 30 × 29 internal directed pairs plus 30 SEC-entry and
   30 MITD-exit pairs, each with eleven linker choices: **10,230 options**.
   Enumerate every overlapping 8–11-aa peptide touching each linker or direct
   boundary; deduplicate sequences only for model inference, retaining occurrence
   counts when calculating edge costs.
2. **Predict.** Run fresh MHCflurry 2.3.14 affinity/percentile predictions against
   all six recorded HLA-I alleles on a Modal A10. Check every provider model file
   and the recorded runtime. Missing queries stop assembly.
3. **Solve.** Reject hydrophobic internal options. Jointly choose a path containing
   every antigen once and one linker per edge, including SEC and MITD boundaries.
   The four existing lexicographic objectives minimize counts at ≤0.5 percentile,
   then ≤2 percentile; maximize the worst boundary percentile; then minimize added
   linker residues. Each phase must be proven optimal. Scan the full chain again
   for hydropathy and independently enumerate the placed junction peptides.
4. **Predict selected.** Make a separate fresh prediction run for the selected
   junction peptides. This does not reuse optimization measurements.
5. **Verify.** Require identical boundary enumeration, threshold classifications
   and numerical agreement with the optimization predictions. Verify all thirty
   sequences appear exactly once, the routing elements are unchanged and the
   full-chain hydropathy policy passes. Fixed SEC/MITD overlaps are reported
   separately; this does not establish cellular topology.
6. **Encode.** Synonymously encode the fixed protein using the same pinned PopRNA
   implementation/settings, human codon table and recorded seed. Preserve the
   existing UTRs and three stops. Check exact translation, motif/codon constraints,
   component identities and GenBank annotations. Report the final folding
   heuristic separately from hard sequence checks. PolyA and downstream extension
   have no assigned sequence or length.

[cassette.py](cassette.py) adapts our earlier assembly worker to this selection;
the objective and inference engine are unchanged. The controller freezes our
maintained source files in its uploaded bundle; provider code and models remain
external. Computational checks do not establish tumour-cell presentation or
immunogenicity.

The recorded authenticated Modal workspace is required for fresh inference.
Offline replay uses collected edge costs and independent selected-junction
predictions, verifies present files exactly, and reports any large uncollected
native prediction table as hash-only. It does not claim to verify absent bytes.
