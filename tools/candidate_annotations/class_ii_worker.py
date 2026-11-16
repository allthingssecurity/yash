"""Fresh CapHLA inference and scoped exact long-peptide protein searches."""
from pathlib import Path
import argparse
import csv
import gzip
import hashlib
import importlib.metadata
import inspect
import io
import json
import math
import os
import subprocess
import time


def require(condition,message):
    if not condition:raise ValueError(message)


def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for chunk in iter(lambda:f.read(8*1024*1024),b''):h.update(chunk)
    return h.hexdigest()


def write_json(path,value,compressed=False):
    body=(json.dumps(value,indent=2,sort_keys=True,allow_nan=False)+'\n').encode()
    if compressed:
        with Path(path).open('wb') as f:
            with gzip.GzipFile(filename='',fileobj=f,mode='wb',mtime=0) as z:z.write(body)
    else:Path(path).write_bytes(body)


def write_table(path,rows):
    stream=io.StringIO();w=csv.DictWriter(stream,fieldnames=list(rows[0]),delimiter='\t',lineterminator='\n')
    w.writeheader();w.writerows(rows);body=stream.getvalue().encode()
    if str(path).endswith('.gz'):
        with Path(path).open('wb') as f:
            with gzip.GzipFile(filename='',fileobj=f,mode='wb',mtime=0) as z:z.write(body)
    else:Path(path).write_bytes(body)


def fasta(path):
    with (gzip.open(path,'rt') if str(path).endswith('.gz') else Path(path).open()) as f:
        header=None;parts=[]
        for line in f:
            if line.startswith('>'):
                if header is not None:yield header,''.join(parts)
                header=line[1:].strip();parts=[]
            else:parts.append(line.strip())
        if header is not None:yield header,''.join(parts)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--request',type=Path,required=True);p.add_argument('--models',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    os.environ['OMP_NUM_THREADS']='4';os.environ['MKL_NUM_THREADS']='4'
    import numpy as np
    import pandas as pd
    import torch
    import ahocorasick
    import mhctools.caphla as cp
    from mhctools import CapHLA
    from mhctools.pred import Kind
    torch.set_num_threads(4);np.random.seed(0);torch.manual_seed(0)
    torch.set_float32_matmul_precision('highest')
    torch.backends.cuda.matmul.allow_tf32=False;torch.backends.cudnn.allow_tf32=False
    device='cuda' if torch.cuda.is_available() else 'cpu'
    if device=='cuda':torch.cuda.manual_seed_all(0)
    req=json.loads(a.request.read_text());manifest=req['model_manifest'];C=a.models
    require(cp.CAPHLA_REVISION==manifest['revision'],'Wrapper/upstream revision differs')
    require(sha(Path(inspect.getfile(cp)))==manifest['mhctools_caphla_wrapper_sha256'],'Changed external CapHLA adapter')
    require(subprocess.check_output(['git','-C',str(C),'rev-parse','HEAD'],text=True).strip()==manifest['revision'],
        'Wrong upstream checkout')
    for name,meta in manifest['model_files'].items():
        path=C/name
        require(path.stat().st_size==meta['bytes'] and sha(path)==meta['sha256'],'Changed model/source file: '+name)
    library=pd.read_csv(C/'HLA_library.csv').set_index('Allele Name')['MHC pseudo-seq'].to_dict()
    targets=[]
    for t in req['targets']:
        requested='HLA-'+t['id'];canonical,key,kind=cp.normalize_caphla_allele(requested)
        require(kind=='II' and key in library,'No exact model identity: '+requested)
        pseudoseq=library[key]
        require(len(pseudoseq)==34,'Unexpected HLA pseudosequence')
        targets.append(dict(id=t['id'],requested_allele=requested,canonical_allele=canonical,model_allele_key=key,
            model_pseudosequence=pseudoseq,model_pseudosequence_X_tokens=pseudoseq.count('X'),
            evidence_group='DRB1_primary' if t['id'] in req['policy']['primary_targets'] else 'exploratory',
            typing_status=t['status'],pairing_status=t['phase_status'],
            DR_alpha_policy='DRA*01:01 implicit; observed DRA*01:02 not independently modeled' if t['locus']=='DR' else 'Explicit alpha-beta hypothesis'))
    require(len({t['canonical_allele'] for t in targets})==len(req['targets']),'Duplicate normalized targets')
    predictor=CapHLA(alleles=[t['requested_allele'] for t in targets],mode='both',caphla_path=str(C),device=device,batch_size=1024)

    def unpack(result):
        components={x.kind:x for x in result.preds}
        require(len(components)==2,'Missing BA or EL result')
        ba,el=components[Kind.pMHC_affinity],components[Kind.pMHC_presentation]
        return dict(peptide=ba.peptide,canonical_allele=ba.allele,affinity_score=float(ba.score),
            affinity_nm=float(ba.value),presentation_score=float(el.score))

    # Use the provider's own examples as implementation controls; no patient
    # predictions are reused from the previous mutation-pool analysis.
    reference=pd.read_csv(C/'test_out.csv').groupby('Allele Name',sort=False).head(2).reset_index(drop=True)
    controls=[unpack(x) for x in predictor.predict_pairs(list(zip(reference.peptide,reference['Allele Name'])))]
    errors=[]
    for observed,expected in zip(controls,reference.to_dict('records')):
        for field in ['affinity_score','presentation_score']:
            err=abs(observed[field]-float(expected[field]));errors.append(err)
            require(math.isclose(observed[field],float(expected[field]),rel_tol=1e-4,abs_tol=1e-5),'Provider example mismatch')
    print('REFERENCE_CONTROLS_PASSED',len(controls),flush=True)
    a.output.mkdir(parents=True,exist_ok=False)
    write_table(a.output/'hla-coverage.tsv',targets)
    started=time.perf_counter();rows=[]
    for target in targets:
        results=[unpack(x) for x in predictor.predict_pairs([(p,target['requested_allele']) for p in req['queries']])]
        require(len(results)==len(req['queries']),'Missing model rows')
        require([x['peptide'] for x in results]==req['queries'],'Changed peptide order')
        for row in results:
            require(row['canonical_allele']==target['canonical_allele'],'Wrong target identity')
            require(all(math.isfinite(row[k]) for k in ['affinity_score','affinity_nm','presentation_score']), 'Non-finite model value')
            require(row['affinity_nm']>0 and 0<=row['presentation_score']<=1,'Invalid BA/EL value')
            row.update(target_id=target['id'],model_allele_key=target['model_allele_key'],evidence_group=target['evidence_group'])
            rows.append(row)
        print('CAPHLA_PREDICTED',target['id'],len(results),flush=True)
        write_table(a.output/('checkpoint-'+str(targets.index(target))+'.tsv.gz'),results)
    inference_seconds=time.perf_counter()-started
    write_table(a.output/'predictions.tsv.gz',rows)
    automaton=ahocorasick.Automaton()
    for peptide in req['queries']:automaton.add_word(peptide,peptide)
    automaton.make_automaton();hits={p:{} for p in req['queries']};sources=[]
    for source in req['normal_reference_sources']:
        path=Path(source['path'])
        require(sha(path)==source['sha256'],'Changed reference catalog: '+source['name'])
        count=0
        for header,seq in fasta(path):
            count+=1
            for end,peptide in automaton.iter(seq):
                hits[peptide].setdefault(source['name'],dict(reference_header=header,protein_offset_0based=end-len(peptide)+1))
        require(count==source['expected_records'],'Reference record-count mismatch')
        sources.append(dict(**source,records_searched=count,query_peptides_with_hits=sum(source['name'] in x for x in hits.values())))
        print('LONG_PEPTIDE_REFERENCE_SEARCH',source['name'],count,flush=True)
    write_json(a.output/'reference-matches.json.gz',hits,True)
    write_json(a.output/'normal-search.json',dict(algorithm='Exact substring search, first witness per catalog',
        query_lengths=[15,16,17,18,19,20],queries=len(req['queries']),sources=sources,
        personal_normal_catalog_designed_max_peptide_length=11,
        personal_normal_negative_is_complete_for_classII=False,
        personal_normal_caveat='Existing local haplotype clusters and flanks were constructed for <=11-aa targets. Positive long-peptide matches are witnesses; negative long-peptide searches do not establish complete genotype-compatible class-II coverage.',
        healthy_translation_or_presentation_measured=False))
    outputs=['hla-coverage.tsv','predictions.tsv.gz','reference-matches.json.gz','normal-search.json']
    meta=dict(status='COMPLETED',new_inference=True,request_sha256=sha(a.request),worker_sha256=sha(__file__),
        caphla_revision=manifest['revision'],all_upstream_model_files_verified=True,upstream_files=len(manifest['model_files']),
        reference_smoke_passed=True,reference_smoke_pairs=len(controls),max_reference_smoke_absolute_score_difference=max(errors),
        model_modes=['BA','EL'],folds_per_mode=5,targets=len(targets),queries=len(req['queries']),prediction_rows=len(rows),
        inference_elapsed_seconds=inference_seconds,total_elapsed_seconds=time.perf_counter()-started,GPU_used=device=='cuda',
        device=device,GPU_name=torch.cuda.get_device_name() if device=='cuda' else None,
        upstream_code_or_weights_redistributed=False,
        versions={k:importlib.metadata.version(k) for k in ['mhctools','mhcgnomes','torch','numpy','pandas','pyahocorasick']},
        output_sha256={name:sha(a.output/name) for name in outputs})
    write_json(a.output/'inference.json',meta);print(json.dumps(meta),flush=True)


if __name__=='__main__':main()
