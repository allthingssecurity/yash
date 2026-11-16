# Execution records

This build uses the separate Modal root
`/data/results/TN25-112103/design-3-final-20261116/attempt-01` on the persistent
`research-lens-workspace` volume. Original sequencing inputs, other runs and
public download settings are unchanged.

`inputs.manifest.json` records uploaded file sizes/hashes and their repository
paths. Each completed stage writes a `complete.json` with actual output
identities. `modal.json` records stage handles after the build finishes. Local
run state/logs remain under ignored `work/design-3-final-20261116/attempt-01/`.

`offline-replay.json`, when present, records successful local selection,
assembly, junction, encoding and annotation replay from the included evidence.
Native receipts remain unchanged. Failed stages require inspection before a
new attempt; observation timeouts alone are not a reason to restart computation.

The layout has one canonical Design #3 directory. Native receipts retain the paths
used at execution; replay resolves their former `final/` prefix to this directory.
`recorded-source/` holds the original selection policy and selector named in the
selection receipt. The maintained selector replays the same selection against
its current input paths. These source records are not alternative designs.
