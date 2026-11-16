"""Check the added H2BC12 target peptides against the pinned newer protein catalog."""
from pathlib import Path
import argparse
import gzip
import hashlib
import json


def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda:f.read(8*1024*1024),b''):h.update(block)
    return h.hexdigest()


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--request',type=Path,required=True)
    p.add_argument('--reference',type=Path);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    q=json.loads(a.request.read_text());reference=a.reference or Path(q['reference'])
    assert sha(reference)==q['reference_sha256']
    assert not a.output.exists(),'Choose a new reference-review output directory'
    peptides=sorted(set(q['peptides']));assert len(peptides)==38
    hits={};records=0;name=None;parts=[]
    def inspect():
        nonlocal records
        if name is None:return
        records+=1;sequence=''.join(parts)
        for peptide in peptides:
            if peptide not in hits and peptide in sequence:
                hits[peptide]=dict(record=name.split()[0],start_aa=sequence.index(peptide),protein_sequence_sha256=hashlib.sha256(sequence.encode()).hexdigest())
    with gzip.open(reference,'rt') as f:
        for line in f:
            if line.startswith('>'):inspect();name=line[1:].strip();parts=[]
            else:parts.append(line.strip())
    inspect()
    result=dict(status='PASS_NO_NEWER_REFERENCE_HITS' if not hits else 'REVIEW_NEWER_REFERENCE_HITS',gene=q['gene'],region_id=q['region_id'],
        peptides=peptides,protein_records=records,matches=hits,reference=str(reference),reference_sha256=sha(reference),
        request_sha256=sha(a.request),worker_sha256=sha(__file__),
        scope='Additional exact Ensembl115 check for H2BC12; complements the recorded Ensembl103/UniProt/personal-normal review. No hit is not measured normal-tissue absence.')
    a.output.mkdir(parents=True);target=a.output/'result.json';target.write_text(json.dumps(result,indent=2)+'\n')
    (a.output/'complete.json').write_text(json.dumps(dict(status='COMPLETED',outputs={'result.json':dict(bytes=target.stat().st_size,sha256=sha(target))}),indent=2)+'\n')
    print(json.dumps(dict(status=result['status'],peptides=len(peptides),records=records,matches=len(hits))))


if __name__=='__main__':main()
