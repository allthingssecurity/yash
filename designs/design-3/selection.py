"""Freeze or verify the agreed thirty-region selection; no new inference."""
from pathlib import Path
from collections import Counter
import argparse
import csv
import hashlib
import io
import json
import sys

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
sys.path.insert(0, str(REPO / 'src'))
sys.path.insert(0, str(REPO / 'tools/rna_evidence'))
from yash_vaccine.hydrophobicity import summarize
from yash_vaccine.io import require_assertions
from normal_peptides import require_selection


def require(condition, message):
    if not condition:
        raise ValueError(message)


def load(path):
    return json.loads(Path(path).read_text())


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def encoded(value):
    return (json.dumps(value, indent=2, allow_nan=False) + '\n').encode()


def selection():
    policy = load(HERE / 'policy.json')
    for path, digest in policy['scientific_inputs'].items():
        require(sha(REPO / path) == digest, 'Scientific input changed: ' + path)
    parent = load(HERE / 'inputs/source-regions.json')
    retained = [r for r in parent if r['source'] != 'SPLICE']
    require(len(retained) == 24, 'Expected 24 unchanged non-splice parent regions')
    h2 = [r for r in load(HERE / 'inputs/h2bc12/source-regions.json')
          if r['source'] == 'MUTATION' and r['gene'] == 'H2BC12']
    require(len(h2) == 1, 'Expected the explicitly approved H2BC12 region')
    require(h2[0]['supplied_gene_value'] is None and h2[0]['manual_inclusion'],
            'H2BC12 expression must remain missing and explicitly included')
    newer = load(HERE / 'inputs/h2bc12/result.json')
    require(newer['status'] == 'PASS_NO_NEWER_REFERENCE_HITS' and not newer['matches'],
            'H2BC12 newer-reference check failed')
    require(newer['region_id'] == h2[0]['region_id'] and
            set(newer['peptides']) == set(h2[0]['scored_target_peptides']),
            'H2BC12 reference search does not cover the selected targets')
    splice_path = REPO / 'data/shared-rna/splicing/junction-only/selected.json'
    gate = require_selection(splice_path)
    splices = load(splice_path)
    require([r['gene'] for r in splices] == policy['RNA_genes'], 'RNA selection changed')
    selected = ([r for r in retained if r['source'] == 'MUTATION'] + h2 + splices +
                [r for r in retained if r['source'] in {'CTA/SELF', 'FUSION'}])
    annotated = [r for r in load(REPO / 'designs/design-3/review/candidate-annotations/request.json')['regions']
                 if r['current_plan']]
    require({(r['region_id'], r['sequence']) for r in selected} ==
            {(r['region_id'], r['sequence']) for r in annotated},
            'Assembly selection differs from the annotated thirty-region plan')
    for index, row in enumerate(selected, 1):
        row['antigen_id'] = f'AG{index:02d}'
        row['hydrophobicity'] = summarize(row['sequence'])
        require(len(row['sequence']) == 25 and row['hydrophobicity']['status'] == 'PASS',
                'Region fails the fixed length/hydropathy rule: ' + row['region_id'])
    require(len(selected) == len({r['region_id'] for r in selected}) ==
            len({r['sequence'] for r in selected}) == 30, 'Expected thirty distinct regions')
    require(dict(Counter(r['source'] for r in selected)) == policy['source_counts'],
            'Source counts disagree with the approved plan')
    require(not any(r['gene'] in {'MUC6', 'FES'} for r in selected), 'Excluded region reintroduced')
    summary = dict(name=policy['name'], regions=30, source_counts=policy['source_counts'],
        retained_parent_nonsplice_regions=24, added_H2BC12=True, H2BC12_expression_missing=True,
        RNA_genes=policy['RNA_genes'], RNA_admission=gate, class_I_II_annotations_do_not_select=True,
        lead_changed=False, new_variant_calling=False, new_raw_read_processing=False)
    return selected, summary


def artifacts(selected, summary):
    fields = list(dict.fromkeys(k for row in selected for k in row))
    stream = io.StringIO()
    writer = csv.DictWriter(stream, fieldnames=fields, delimiter='\t', lineterminator='\n')
    writer.writeheader()
    for row in selected:
        writer.writerow({k: json.dumps(v) if isinstance(v, (list, dict)) else v for k, v in row.items()})
    return {'selected.json': encoded(selected), 'selected.tsv': stream.getvalue().encode(),
        'selected.fasta': ''.join('>' + r['antigen_id'] + ' ' + r['gene'] + ' ' + r['region_id'] +
                                  '\n' + r['sequence'] + '\n' for r in selected).encode(),
        'summary.json': encoded(summary)}


def main():
    require_assertions('Final thirty-region selection')
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--verify', action='store_true')
    parser.add_argument('--output', type=Path, default=HERE / 'selection')
    args = parser.parse_args()
    selected, summary = selection()
    outputs = artifacts(selected, summary)
    if args.verify:
        for name, body in outputs.items():
            require((args.output / name).read_bytes() == body, 'Selection changed: ' + name)
        receipt = load(args.output / 'complete.json')
        require(receipt['policy_sha256'] == sha(HERE / 'execution/recorded-source/policy.json') and
                receipt['selector_sha256'] == sha(HERE / 'execution/recorded-source/selection.py'), 'Selection source identity changed')
        for name, metadata in receipt['outputs'].items():
            require(sha(args.output / name) == metadata['sha256'], 'Selection hash changed: ' + name)
    else:
        require(not args.output.exists(), 'Choose a fresh output; selections are immutable')
        args.output.mkdir(parents=True)
        for name, body in outputs.items():
            (args.output / name).write_bytes(body)
        receipt = dict(status='COMPLETED', policy_sha256=sha(HERE / 'policy.json'),
            selector_sha256=sha(__file__), outputs={name: dict(bytes=len(body),
            sha256=hashlib.sha256(body).hexdigest()) for name, body in outputs.items()})
        (args.output / 'complete.json').write_bytes(encoded(receipt))
    print(json.dumps(dict(status='SELECTION_VERIFIED' if args.verify else 'SELECTION_CREATED',
                         regions=30, sources=summary['source_counts'])))


if __name__ == '__main__':
    main()
