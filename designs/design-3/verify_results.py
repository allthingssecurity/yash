"""Verify recorded lead inputs and outputs without optimization or inference."""
from pathlib import Path
import csv
import hashlib
import json
import subprocess
import sys

from Bio import SeqIO
from Bio.Seq import Seq

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
sys.path.insert(0, str(REPO / 'src'))
from yash_vaccine.io import require_assertions, child_env
from replay import recorded_input


def require(value, message):
    if not value:
        raise ValueError(message)


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def load(path):
    return json.loads(Path(path).read_text())


def verify():
    require_assertions('Design #3 verification')
    selection = subprocess.run(
        [sys.executable, str(HERE / 'selection.py'), '--verify'],
        cwd=REPO, capture_output=True, text=True, env=child_env(),
    )
    require(selection.returncode == 0, selection.stdout + selection.stderr)
    inputs = load(HERE / 'execution/inputs.manifest.json')
    for record in inputs.values():
        path = recorded_input(record['repository_path'])
        require(path.stat().st_size == record['bytes'] and sha(path) == record['sha256'],
                'Changed input ' + str(path))

    worker = sha(HERE / 'pipeline/cassette.py')
    checked, absent = 0, []
    for folder in ['prepare', 'all-inference', 'protein', 'selected-inference', 'verify', 'mrna']:
        base = HERE / 'construct' / folder
        receipt = load(base / 'complete.json')
        require(receipt['status'] == 'COMPLETED' and receipt['worker_sha256'] == worker,
                'Unverified stage ' + folder)
        for name, record in receipt['outputs'].items():
            path = base / name
            if path.is_file():
                require(path.stat().st_size == record['bytes'] and sha(path) == record['sha256'],
                        'Changed output ' + str(path))
                checked += 1
            else:
                require(folder == 'all-inference' and name == 'predictions.tsv' and record['bytes'] > 40_000_000,
                        'Missing output ' + str(path))
                absent.append(dict(file=str(path.relative_to(REPO)), **record,
                                   verification='HASH_ONLY_NATIVE_OUTPUT_NOT_DOWNLOADED'))

    selected = load(HERE / 'selection/selected.json')
    protein = str(SeqIO.read(HERE / 'construct/protein/cassette.fasta', 'fasta').seq)
    with (HERE / 'construct/protein/components.tsv').open() as stream:
        parts = list(csv.DictReader(stream, delimiter='\t'))
    require(len(selected) == 30 and len(protein) == 885, 'Wrong lead candidate')
    require(''.join(p['aa'] for p in parts) == protein, 'Component/protein disagreement')
    require(sorted((p['region_id'], p['aa']) for p in parts if p['kind'] == 'antigen') ==
            sorted((r['region_id'], r['sequence']) for r in selected), 'Antigen identity disagreement')

    cds = str(SeqIO.read(HERE / 'construct/mrna/cds.fasta', 'fasta').seq)
    require(cds[-3:] == 'TAA' and str(Seq(cds[:-3]).translate()) == protein, 'Translation mismatch')
    fixed = load(HERE / 'construct/mrna/backbone.json')['fixed_elements']
    rna = str(SeqIO.read(HERE / 'construct/mrna/mrna-segments.fasta', 'fasta').seq)
    expected = fixed['utr5_dna'] + cds[:-3] + fixed['stop_block_dna'] + fixed['utr3_body_dna']
    require(rna == expected.replace('T', 'U'), 'RNA fixed elements mismatch')
    gb = SeqIO.read(HERE / 'construct/mrna/template-segments.gb', 'genbank')
    require(len(gb.features) == 76, 'GenBank feature count')
    encoding = load(HERE / 'construct/mrna/verification.json')
    require(encoding['all_constraints_pass'] and encoding['polyA_status'] == 'unspecified',
            'Encoding constraints failed')
    require(load(HERE / 'construct/final-check.json')['protein_sha256'] ==
            sha(HERE / 'construct/protein/cassette.fasta'), 'Junction/protein mismatch')
    return dict(
        status='VERIFIED', design='design-3', selected_regions=30, protein_aa=885,
        known_rna_nt=len(rna), included_inputs_byte_verified=len(inputs),
        included_outputs_byte_verified=checked, hash_only_native_outputs=absent,
        PolyA='unspecified', new_inference=False, BAMs_reread=False,
    )


if __name__ == '__main__':
    print(json.dumps(verify()))
