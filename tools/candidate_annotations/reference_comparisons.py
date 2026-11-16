"""Explain selected splice/frameshift regions against anchored reference transcripts.

No selection, prediction or sequence optimization is performed. Reference models
are obtained from the documented GENCODE/GRCh38 preparation, not redistributed.
"""
from pathlib import Path
from collections import defaultdict
import argparse
import csv
import gzip
import hashlib
import json

ROOT = Path(__file__).resolve().parents[2]
DESIGN = ROOT / 'designs/design-3'
WINDOWS = ROOT / 'data/shared-rna/splicing/junction-only/original-rna-windows.jsonl.gz'
RNA = ROOT / 'designs/design-2/inputs/rna/primary-proteins.tsv'
OUTPUT = DESIGN / 'review/reference-comparisons.json'


def require(value, message):
    if not value:
        raise ValueError(message)


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def load(path):
    with (gzip.open(path, 'rt') if str(path).endswith('.gz') else Path(path).open()) as stream:
        return json.load(stream)


def selected_windows():
    with gzip.open(WINDOWS, 'rt') as stream:
        return [json.loads(line) for line in stream]


def coding_positions(model):
    # These GTF model intervals are one-based, inclusive; read witnesses are zero-based.
    parts = sorted(model['parts'], reverse=model['strand'] == '-')
    positions = []
    for start, end, phase in parts:
        positions.extend(range(end-1, start-2, -1) if model['strand'] == '-' else range(start-1, end))
    require(parts[0][2] == 0, 'Partial leading CDS requires a separately assessed comparison')
    return positions


def generate(model_path):
    models = load(model_path)['models']
    selections = {r['region_id']: r for r in load(DESIGN/'selection/selected.json')}
    by_window = {r['window_id']:r for r in selections.values() if r['source']=='SPLICE'}
    with (ROOT/'data/shared-rna/splicing/reconstruction/reference-transcripts.tsv').open() as stream:
        exact = {r['transcript_id'] for r in csv.DictReader(stream, delimiter='\t') if r['reference_translation_exact'] == 'True'}
    entries = {}
    for window in selected_windows():
        candidate = window['sequence']
        region_id = by_window[window['window_id']]['region_id']
        require(selections[region_id]['sequence'] == candidate, 'Stale selected RNA window')
        variants = {}
        for anchor in window['upstream_frame_anchors']:
            tid = anchor['transcript_id']; model = models[tid]
            require(tid in exact and model['gene_id'] == window['gene_id'], 'Unverified reference transcript')
            positions = coding_positions(model); index = {p: i for i, p in enumerate(positions)}
            start_nt = anchor['upstream_CDS_index_0based'] + window['start_nt_relative_to_join'] - anchor['upstream_query_offset']
            require(start_nt >= 0 and start_nt % 3 == 0, 'Reference start is not an anchored complete codon')
            require(index.get(window['genomic_positions_0based'][0]) == start_nt, 'Reference start coordinate mismatch')
            first_anchor = anchor['upstream_coding_positions_0based'][0]
            require(index.get(first_anchor) == anchor['upstream_CDS_index_0based'] + anchor['first_complete_codon_query_offset'] - anchor['upstream_query_offset'], 'Reference anchor mismatch')
            join = window['junction_nt_offset']
            donor = index.get(window['genomic_positions_0based'][join-1])
            acceptor = index.get(window['genomic_positions_0based'][join])
            deleted = acceptor-donor-1 if donor is not None and acceptor is not None and acceptor > donor else None
            effect = 'in_frame_deletion' if deleted is not None and deleted % 3 == 0 else 'different_downstream_sequence'
            length = 25 + deleted//3 if effect == 'in_frame_deletion' else 25
            reference = model['reference_protein'][start_nt//3:start_nt//3+length]
            require(len(reference) == length, 'Truncated reference interval')
            variant = dict(reference=reference, effect=effect, deleted_coding_nt=deleted,
                           transcripts=[], transcript_intervals=[])
            if effect == 'in_frame_deletion':
                lo = window['target_start_in_window']; split = int(join % 3 != 0)
                right = 25-lo-split
                prefix = reference[:lo]; ref_codon = reference[lo:lo+split]; alt_codon = candidate[lo:lo+split]
                removed = reference[lo+split:len(reference)-right]
                suffix = reference[-right:] if right else ''
                require(prefix + alt_codon + suffix == candidate, 'Deletion does not reconstruct candidate')
                require(len(removed) == deleted//3, 'Wrong removed reference interval')
                variant.update(prefix=prefix, reference_codon=ref_codon, candidate_codon=alt_codon,
                               removed=removed, suffix=suffix, removed_aa=len(removed))
            key = (reference, effect, deleted)
            if key not in variants:
                variants[key] = variant
            variants[key]['transcripts'].append(tid)
            variants[key]['transcript_intervals'].append(dict(transcript_id=tid, start_aa=start_nt//3+1, end_aa=start_nt//3+length))
        entries[region_id] = dict(source='SPLICE', gene=window['gene'], candidate=candidate,
            window_id=window['window_id'], join_nt_offset=window['junction_nt_offset'],
            target_start=window['target_start_in_window'], target_end=window['target_end_in_window'],
            variants=list(variants.values()))

    # A frameshift needs its transcript sequence, not a single-residue reversion.
    with RNA.open() as stream:
        contexts = {r['variant_id']:r for r in csv.DictReader(stream, delimiter='\t')}
    for selected in selections.values():
        if selected['source'] != 'MUTATION' or selected.get('protein_change') != 'p.N515fs':
            continue
        context = contexts[selected['variant_id']]
        prefix = context['sequence'][:int(context['mutation_start_0based'])]
        ids = set(json.loads(context['transcript_ids'])); variants = {}
        for tid, model in models.items():
            if tid.split('.')[0] not in ids:
                continue
            protein = model['reference_protein']; offset = protein.find(prefix)
            require(offset >= 0 and protein.find(prefix, offset+1) < 0, 'Frameshift upstream anchor is ambiguous')
            start = offset + int(selected['window_start_0based'])
            reference = protein[start:start+25]
            require(len(reference) == 25 and reference[:int(selected['target_start_in_window'])] == selected['sequence'][:int(selected['target_start_in_window'])], 'Frameshift reference prefix mismatch')
            variant = variants.setdefault(reference, dict(reference=reference, effect='frameshift', transcripts=[], transcript_intervals=[]))
            variant['transcripts'].append(tid)
            variant['transcript_intervals'].append(dict(transcript_id=tid, start_aa=start+1, end_aa=start+25))
        require(variants, 'Missing frameshift reference transcripts')
        entries[selected['region_id']] = dict(source='MUTATION', gene=selected['gene'], candidate=selected['sequence'], variants=list(variants.values()))
    return dict(schema_version=1, reference='GENCODE 37 / GRCh38 annotated protein transcripts',
        reference_models_sha256=sha(model_path), selected_sha256=sha(DESIGN/'selection/selected.json'),
        original_windows_sha256=sha(WINDOWS), RNA_contexts_sha256=sha(RNA),
        normal_RNA_measured=False, selection_changed=False, entries=entries)


def verify(path=OUTPUT):
    result = load(path)
    require(result['selected_sha256'] == sha(DESIGN/'selection/selected.json') and
            result['original_windows_sha256'] == sha(WINDOWS) and result['RNA_contexts_sha256'] == sha(RNA), 'Comparison inputs changed')
    selected = {r['region_id']:r for r in load(DESIGN/'selection/selected.json')}
    windows = {r['window_id']:r for r in selected_windows()}
    require(len(result['entries']) == 6 and not result['normal_RNA_measured'] and not result['selection_changed'], 'Wrong comparison scope')
    for rid, entry in result['entries'].items():
        require(entry['candidate'] == selected[rid]['sequence'] and entry['variants'], 'Wrong comparison candidate')
        if entry['source'] == 'SPLICE':
            window = windows[entry['window_id']]
            require(entry['candidate'] == window['sequence'], 'Wrong observed splice window')
        for variant in entry['variants']:
            require(variant['transcripts'] and len(variant['reference']) >= 25, 'Incomplete reference witness')
            if variant['effect'] == 'in_frame_deletion':
                require(variant['prefix'] + variant['reference_codon'] + variant['removed'] + variant['suffix'] == variant['reference'], 'Wrong reference deletion blocks')
                require(variant['prefix'] + variant['candidate_codon'] + variant['suffix'] == entry['candidate'], 'Wrong spliced deletion blocks')
                require(len(variant['removed']) == variant['removed_aa'] == variant['deleted_coding_nt']//3, 'Wrong deletion size')
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--models', type=Path, help='Prepared GENCODE 37 transcript-models.json.gz; not bundled')
    parser.add_argument('--output', type=Path, default=OUTPUT)
    parser.add_argument('--verify', action='store_true')
    args = parser.parse_args()
    if args.verify:
        result = verify(args.output)
        if args.models:
            require(generate(args.models) == result, 'Reference model replay differs')
    else:
        if not args.models:
            parser.error('--models is required to regenerate reference transcript comparisons')
        result = generate(args.models)
        args.output.write_text(json.dumps(result, indent=2, sort_keys=True)+'\n')
        verify(args.output)
    print('Reference comparisons verified: five splice regions and the TGFBR2 frameshift; no selection or prediction changes.')


if __name__ == '__main__':
    main()
