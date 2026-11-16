"""Exclude exact self targets before RNA ranking and gate the scored selection.

This reuses checksum-verified protein-search measurements only. It does not
apply the independent assessment's mapper, pattern or contextual-error gates.
"""
from collections import Counter, defaultdict
from pathlib import Path
import argparse
import csv
import gzip
import hashlib
import io
import json
import math

from pool_report import contribution, load, table_rows

ROOT = Path(__file__).resolve().parents[2]
DEFAULT_SCREEN = ROOT / 'data/shared-rna/splicing/sequence-screen'
DEFAULT_PROTEINS = ROOT / 'data/shared-rna/splicing/independent-assessment/protein-screen'
CATALOGS = {'ensembl103', 'ensembl115', 'uniprot_canonical', 'uniprot_additional', 'personal_normal'}


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def verify_native(directory, request):
    receipt = load(directory / 'complete.json')
    if receipt['status'] != 'COMPLETED' or receipt['request_sha256'] != sha(request):
        raise ValueError(f'Unverified measured search: {directory}')
    for name, record in receipt['outputs'].items():
        path = directory / name
        if path.stat().st_size != record['bytes'] or sha(path) != record['sha256']:
            raise ValueError(f'Changed measured output: {path}')
    return receipt


def checked_pool(screen=DEFAULT_SCREEN, proteins=DEFAULT_PROTEINS, *, reference_join_only=True):
    """Return checked windows; by default require reference-join AA identity."""
    screen, proteins = Path(screen), Path(proteins)
    verify_native(screen, screen / 'request.json.gz')
    verify_native(proteins, proteins.parent / 'request.json.gz')
    request = load(screen / 'request.json.gz')
    reference = load(screen / 'reference-matches.json.gz')
    broader = load(proteins / 'protein-matches.json.gz')
    metadata = load(proteins / 'protein-screen.json')
    if {s['name'] for s in metadata['sources']} != CATALOGS:
        raise ValueError('All configured protein catalogs must be searched')
    if not set(request['affinities']) <= set(reference) & set(broader):
        raise ValueError('Missing protein-search decisions; unknown is not nonself')
    sources = {s['name']: s for s in metadata['sources']}
    if sources['ensembl103']['sha256'] != request['reference_ensembl_sha256']:
        raise ValueError('Incompatible Ensembl search measurements')
    components = request['reference_uniprot']['component_files']
    if {s['sha256'] for s in components} != {sources[k]['sha256'] for k in ['uniprot_canonical', 'uniprot_additional']}:
        raise ValueError('Incompatible UniProt search measurements')

    witnesses = {}
    for peptide in request['affinities']:
        old, new = reference[peptide], broader[peptide]
        old_hit = bool(old['ensembl_all'] or old['uniprot'] or old['personal']['occurs_in_searched_models'])
        new_hit = any(bool(new.get(k)) for k in CATALOGS)
        if old_hit != old['self_match'] or new_hit != new['self_match']:
            raise ValueError(f'Inconsistent self-search witness: {peptide}')
        witnesses[peptide] = dict(self_match=old_hit or new_hit,
            catalogs=sorted({k for k in CATALOGS if new.get(k)} |
                ({'personal_normal'} if old['personal']['occurs_in_searched_models'] else set())))

    by_window = defaultdict(list)
    for row in table_rows(screen / 'peptides.tsv.gz'):
        by_window[row['window_id']].append(row)
    window_metrics = {r['window_id']: r for r in table_rows(screen / 'windows.tsv.gz')}
    path_evidence = {r['window_id']: r for r in load(screen / 'window-evidence.json.gz')}
    decisions = {}
    changed = 0
    for window in request['windows']:
        wid = window['window_id']
        peptide_rows = by_window[wid]
        lo, hi = int(window['target_start_in_window']), int(window['target_end_in_window'])
        expected = [(i, i + k, window['sequence'][i:i + k])
            for k in request['policy']['peptide_lengths']
            for i in range(len(window['sequence']) - k + 1)
            if (i < hi and i + k > lo if hi > lo else i < lo < i + k)]
        if Counter(expected) != Counter((int(p['start']), int(p['end']), p['peptide']) for p in peptide_rows):
            raise ValueError(f'Incomplete target enumeration: {wid}')
        retained = []
        for row in peptide_rows:
            peptide = row['peptide']
            values = request['affinities'][peptide]
            if set(values) != set(request['alleles']):
                raise ValueError(f'Missing peptide–HLA predictions: {peptide}')
            score = sum(contribution(nm) for nm in values.values())
            if not math.isclose(score, float(row['native_contribution']), rel_tol=1e-8, abs_tol=1e-8):
                raise ValueError(f'Changed affinity score: {peptide}')
            # Retain the earlier normal-evidence tier; add no remapping gates.
            if row['tier'] == 'PRIMARY' and not witnesses[peptide]['self_match']:
                retained.append(row)
            elif row['tier'] == 'PRIMARY' and witnesses[peptide]['self_match']:
                changed += 1
        score = sum(float(p['native_contribution']) for p in retained)
        seeds = [p for p in retained if p['native_PASS_seed'] == 'True'
                 and float(p['best_affinity_nm']) <= request['policy']['admission_affinity_nm']]
        expression = float(window['supplied_transformed_gene_tpm'])
        from_rows = window_metrics[wid]
        earlier_eligible = bool(seeds) and score > 0 and expression > 0 and float(from_rows['max_KD19']) <= 1.6
        paths = path_evidence[wid]['path_evidence']
        reference_translations = [p['reference_translation_on_recorded_path'] for p in paths]
        reference_matches = bool(paths) and all(
            len(ref) == len(window['sequence']) == 25 and ref == window['sequence']
            for ref in reference_translations)
        eligible = earlier_eligible and (reference_matches or not reference_join_only)
        best = min(retained, key=lambda p: float(p['best_affinity_nm']), default=None)
        decisions[wid] = dict(window_id=wid, context_id=wid, region_id=window['region_id'],
            event_id=window['event_id'], gene_id=window['gene_id'], gene=window['gene'],
            source='SPLICE', sequence=window['sequence'], coding_sequence=window['coding_sequence'],
            window_start=0, window_end=len(window['sequence']), fragment_count=window['fragment_count'],
            supplied_transformed_gene_tpm=expression, expression_proxy_left=expression,
            target_start_in_window=lo, target_end_in_window=hi,
            target_peptides=sorted({p['peptide'] for p in retained if float(p['native_contribution']) > 0}),
            native_affinity_sum=score, weighted_score_left=score * math.sqrt(expression),
            best_target_peptide=best['peptide'] if best else None,
            best_target_allele=best['best_allele'] if best else None,
            best_target_affinity_nm=float(best['best_affinity_nm']) if best else None,
            eligible_before_reference_join_filter=earlier_eligible,
            reference_join_translation_matches=reference_matches,
            reference_join_path_count=len(paths),
            reference_join_translations=reference_translations,
            primary_window_eligible=eligible,
            reason='PASS' if eligible else 'REFERENCE_JOIN_TRANSLATION_MISMATCH'
                if earlier_eligible and reference_join_only and not reference_matches
                else 'NO_QUALIFYING_TARGET_AFTER_SELF_AND_RECORDED_EVIDENCE_CHECKS')
    result = dict(status='BROAD_SELF_PEPTIDE_SCREEN_PASSED', windows=len(decisions),
        eligible_windows=sum(w['primary_window_eligible'] for w in decisions.values()),
        reference_join_translation_required=reference_join_only,
        eligible_windows_before_reference_join_filter=sum(w['eligible_before_reference_join_filter'] for w in decisions.values()),
        reference_join_matching_windows=sum(w['reference_join_translation_matches'] for w in decisions.values()),
        eligible_events=len({w['event_id'] for w in decisions.values() if w['primary_window_eligible']}),
        eligible_genes=len({w['gene_id'] for w in decisions.values() if w['primary_window_eligible']}),
        queried_target_peptides=len(witnesses), exact_self_target_peptides=sum(w['self_match'] for w in witnesses.values()),
        additional_primary_occurrences_removed=changed, catalogs=[{k: v for k, v in s.items() if k != 'path'} for s in metadata['sources']],
        measured_search_receipts={str(p.relative_to(ROOT)) if p.is_relative_to(ROOT) else str(p): sha(p)
            for p in [screen / 'complete.json', proteins / 'complete.json']},
        independent_mapper_or_pattern_gates_applied=False, new_reference_search=False,
        new_model_inference=False, normal_tissue_expression_or_presentation_measured=False,
        scope='Exact altered-target peptide matches, including alternate contigs, pseudogenes, isoforms and bounded personal-normal models. CTA/self targets are outside this gate.')
    return decisions, witnesses, result


def require_selection(selection, screen=DEFAULT_SCREEN, proteins=DEFAULT_PROTEINS, *, reference_join_only=True):
    """Reject self targets, stale scores or unassessed RNA windows before assembly."""
    decisions, _, result = checked_pool(screen, proteins, reference_join_only=reference_join_only)
    selected = load(selection) if str(selection).endswith('.json') else table_rows(selection)
    candidates = [r for r in selected if r.get('source') == 'SPLICE']
    if not candidates:
        raise ValueError('No SPLICE-source regions; other branches are outside this gate')
    failures = []
    for row in candidates:
        checked = decisions.get(row.get('window_id'))
        reasons = []
        if not checked:
            reasons.append('UNSEARCHED_WINDOW')
        elif not checked['primary_window_eligible']:
            reasons.append(checked['reason'])
        if checked:
            for key in ['region_id', 'sequence', 'coding_sequence', 'target_start_in_window', 'target_end_in_window']:
                if key in row and str(row[key]) != str(checked[key]):
                    reasons.append('CHANGED_' + key.upper())
            targets = row.get('target_peptides')
            if targets is not None:
                targets = set(filter(None, targets.split(';'))) if isinstance(targets, str) else set(targets)
                if targets != set(checked['target_peptides']):
                    reasons.append('STALE_OR_SELF_TARGET_LIST')
            for key in ['native_affinity_sum', 'weighted_score_left', 'best_target_affinity_nm']:
                if key in row and (checked[key] is None or not math.isclose(float(row[key]), checked[key], rel_tol=1e-8, abs_tol=1e-8)):
                    reasons.append('STALE_' + key.upper())
        if reasons:
            failures.append(dict(gene=row.get('gene'), window_id=row.get('window_id'), reasons=reasons))
    if failures:
        raise ValueError(json.dumps(dict(status='RNA_SELF_PEPTIDE_GATE_FAILED', failures=failures)))
    result.update(status='RNA_SELF_PEPTIDE_GATE_PASSED', selected_RNA_regions=len(candidates),
        genes=[r['gene'] for r in candidates], other_source_regions_not_assessed=len(selected) - len(candidates))
    return result


def selection_outputs(decisions, result, count):
    """Apply the recorded weighted-score/event/region tie-break and gene quota."""
    if count < 1:
        raise ValueError('The gene quota must be positive')
    ranked = sorted((r for r in decisions.values() if r['primary_window_eligible']),
        key=lambda r: (-r['weighted_score_left'], r['event_id'], r['region_id']))
    genes = {}
    for row in ranked:
        genes.setdefault(row['gene_id'], row)
    ranking = [dict(row, gene_rank=i + 1) for i, row in enumerate(genes.values())]
    if len(ranking) < count:
        raise ValueError(f'Only {len(ranking)} passing genes; cannot select {count}')
    selected = ranking[:count]
    policy = dict(schema_version=1, region_length_aa=25, top_genes=count,
        reference_join_identity='Observed amino-acid sequence must equal the reference-base translation across the same proposed join on every recorded supporting path.',
        synonymous_changes_allowed=True, reference_translation_is_not_a_normal_transcript=True,
        read_evidence_basis='Original supplied RNA alignment and recorded reconstruction/frame/native/strand rules only.',
        unchanged='Recorded RNA support/frame/native review, coordinate-normal ND, strand, exact self exclusions, <=500 nM native-PASS admission, positive supplied expression and hydropathy rules.',
        score='Native Vaxrank affinity contribution sum * sqrt(supplied transformed gene value)',
        ranking='Descending weighted score, then ascending event_id and region_id; best passing window per gene, one region per gene.',
        mapper_pattern_contextual_error_gates_required=False,
        no_automatic_construct_rebuild=True)
    summary = dict(result, status='REFERENCE_JOIN_ONLY_RNA_SELECTION_COMPLETE',
        selected_genes=[r['gene'] for r in selected], selected_window_ids=[r['window_id'] for r in selected],
        selected_regions=len(selected), policy=policy,
        source_sha256=sha(Path(__file__)), construct_changed=False)

    def encoded(value):
        return (json.dumps(value, indent=2) + '\n').encode()

    def tsv(records):
        stream = io.StringIO(newline='')
        fields = list(dict.fromkeys(k for row in records for k in row))
        writer = csv.DictWriter(stream, fieldnames=fields, delimiter='\t', lineterminator='\n')
        writer.writeheader()
        for row in records:
            writer.writerow({k: json.dumps(v, separators=(',', ':')) if isinstance(v, (dict, list)) else v for k, v in row.items()})
        return stream.getvalue().encode()

    def compressed(body):
        stream = io.BytesIO()
        # GzipFile fixes the portable header; gzip.compress(mtime=0) has
        # different OS-header bytes across supported Python versions.
        with gzip.GzipFile(filename='', mode='wb', fileobj=stream, mtime=0) as handle:
            handle.write(body)
        return stream.getvalue()

    lines = ['# Junction-derived RNA selection', '',
        'The selected windows contain no additional amino-acid substitutions: actual RNA translation equals the reference bases translated across the same proposed join on every recorded path. Synonymous nucleotide differences are allowed. The proposed join is retained in that comparison; this is not a claim that the RNA is a normal transcript.', '',
        f"The earlier pool has {result['windows']} windows; {result['eligible_windows_before_reference_join_filter']} pass the preceding sequence/self checks. The added identity rule leaves **{result['eligible_windows']} windows, {result['eligible_events']} events and {result['eligible_genes']} genes**. We select one window from each of the top {count} genes.", '',
        '| Rank | Gene | Selected 25-aa region | Best predicted target | HLA | Affinity, nM | Supplied gene value | Weighted score |',
        '|---:|---|---|---|---|---:|---:|---:|']
    for row in selected:
        lines.append(f"| {row['gene_rank']} | {row['gene']} | `{row['sequence']}` | `{row['best_target_peptide']}` | {row['best_target_allele']} | {row['best_target_affinity_nm']:.2f} | {row['supplied_transformed_gene_tpm']:.3f} | {row['weighted_score_left']:.2f} |")
    lines += ['', 'Selection and manual review use the original supplied RNA alignment and recorded rules. Scores and native event/region tie-breaks are unchanged. Equal-score windows can have different flanks; the selected IDs, not gene names alone, identify this selection.', '',
        '[Manual review](manual-review.md) · [Selected JSON](selected.json) · [Selected TSV](selected.tsv) · [All ten genes](gene-ranking.tsv) · [All window decisions](all-windows.tsv.gz) · [Policy](policy.json) · [Summary](summary.json) · [Verification](verification.json)', '',
        '## Reproduce', '', 'From the repository root:', '', '```sh',
        f'python tools/rna_evidence/normal_peptides.py --top-genes {count} --results data/shared-rna/splicing/junction-only',
        'python tools/rna_evidence/normal_peptides.py --verify --results data/shared-rna/splicing/junction-only', '```', '',
        'Verification recomputes admission and ranking from checksum-verified measured searches and recorded scores, then compares every generated selection artifact byte for byte. It does not reread BAMs, repeat protein searches or run HLA inference. Unknown search results stop the process. The selected windows remain computational hypotheses; translated tumour protein, HLA presentation and T-cell recognition have not been measured. No construct is rebuilt by this selection stage.', '']
    outputs = {'policy.json': encoded(policy), 'summary.json': encoded(summary),
        'selected.json': encoded(selected), 'selected.tsv': tsv(selected),
        'gene-ranking.tsv': tsv(ranking),
        'selected.fasta': ''.join(f">{r['gene']} {r['window_id']} {r['region_id']}\n{r['sequence']}\n" for r in selected).encode(),
        'all-windows.tsv.gz': compressed(tsv([decisions[k] for k in sorted(decisions)])),
        'README.md': '\n'.join(lines).encode()}
    receipt = dict(status='REFERENCE_JOIN_SELECTION_ARTIFACTS_VERIFIED',
        files={name: dict(bytes=len(body), sha256=hashlib.sha256(body).hexdigest()) for name, body in outputs.items()},
        source_sha256=sha(Path(__file__)), selected_regions=len(selected),
        exact_translation_and_all_paths_checked=True, quota_and_original_tie_break_applied=True,
        source_measurements_verified=True, new_inference=False, BAMs_reread=False, construct_changed=False)
    outputs['verification.json'] = encoded(receipt)
    return outputs, summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--selection', type=Path, help='JSON/TSV selection to gate before assembly')
    parser.add_argument('--screen', type=Path, default=DEFAULT_SCREEN)
    parser.add_argument('--proteins', type=Path, default=DEFAULT_PROTEINS)
    parser.add_argument('--output', type=Path, help='Write a verification summary')
    parser.add_argument('--export', type=Path, help='Write score-ready, eligible windows after exact-self exclusion')
    parser.add_argument('--include-substitution-dependent', action='store_true',
        help='Historical replay only: omit the new full-window reference-join identity rule')
    parser.add_argument('--top-genes', type=int, help='Select this many genes using the existing score and tie-break')
    parser.add_argument('--results', type=Path, help='Directory for the reproducible junction-only selection artifacts')
    parser.add_argument('--verify', action='store_true', help='Recompute and byte-compare an existing --results selection')
    args = parser.parse_args()
    if (args.top_genes is not None or args.verify) and not args.results:
        parser.error('--top-genes/--verify requires --results')
    if args.results and (args.include_substitution_dependent or args.selection or args.export or args.output):
        parser.error('--results is the junction-only selection mode; do not combine it with the other output/historical modes')
    reference_join_only = not args.include_substitution_dependent
    decisions, _, result = checked_pool(args.screen, args.proteins, reference_join_only=reference_join_only)
    if args.results:
        count = args.top_genes if args.top_genes is not None else (load(args.results / 'policy.json')['top_genes'] if args.verify else 5)
        outputs, result = selection_outputs(decisions, result, count)
        for name, body in outputs.items():
            path = args.results / name
            if path.exists():
                if path.read_bytes() != body:
                    raise ValueError(f'Selection artifact differs: {path}; preserve recorded results and use a new directory')
            elif args.verify:
                raise ValueError(f'Missing selection artifact: {path}')
            else:
                args.results.mkdir(parents=True, exist_ok=True)
                path.write_bytes(body)
        require_selection(args.results / 'selected.json', args.screen, args.proteins)
        print(json.dumps(result, indent=2))
        return
    if args.selection:
        result = require_selection(args.selection, args.screen, args.proteins, reference_join_only=reference_join_only)
    if args.export:
        args.export.parent.mkdir(parents=True, exist_ok=True)
        args.export.write_text(json.dumps([w for w in decisions.values() if w['primary_window_eligible']], indent=2) + '\n')
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
