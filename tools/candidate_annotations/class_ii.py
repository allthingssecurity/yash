"""CapHLA class-II annotation: explicit plan, fresh Modal inference or offline replay."""
from pathlib import Path
from collections import defaultdict
import argparse
import datetime
import json
import math
import re

from analysis import ROOT, DEFAULT as CLASS_I, require, sha, load, table, json_bytes, tsv_bytes, zipped

DEFAULT = ROOT/'designs/design-3/review/class-II-annotations'
MANIFEST = Path(__file__).with_name('caphla-manifest.json')
TYPING = ROOT/'data/hla/class-II.json'
MODEL_DIR = '/data/results/TN25-112103/class2-caphla-predictions-20261030/CapHLA'


def request():
    from analysis import finish
    finish(CLASS_I,verify=True)
    parent = load(CLASS_I/'request.json')
    require(len({r['annotation_id'] for r in parent['regions']})==len(parent['regions']), 'Duplicate window annotation identity')
    typing = load(TYPING)
    chains = typing['observed_class_II_DNA_chains']
    expected_targets = set(chains['DRB1']) | {
        a+'/'+b for alpha,beta in [('DQA1','DQB1'),('DPA1','DPB1')]
        for a in chains[alpha] for b in chains[beta]}
    require(len(typing['targets']) == len(expected_targets) and
            {t['id'] for t in typing['targets']} == expected_targets,
            'Active targets must be DRB1 and the observed DQ/DP hypotheses; DRB3 is excluded')
    occurrences,queries,coverage = [],set(),[]
    for region in parent['regions']:
        seq,lo,hi = region['sequence'],region['target_start'],region['target_end']
        change = re.fullmatch(r'p\.([A-Z])(\d+)([A-Z])',region['protein_change'] or '')
        reference = None
        if region['source'] == 'MUTATION' and change:
            require(hi == lo+1 and seq[lo] == change[3], 'Wrong mutation interval')
            reference = seq[:lo]+change[1]+seq[lo+1:]
        for length in range(15,21):
            count = 0
            for start in range(len(seq)-length+1):
                end = start+length
                overlaps = start < hi and end > lo if hi > lo else start < lo < end
                if region['source'] != 'CTA/SELF' and not overlaps:
                    continue
                peptide = seq[start:end]
                row = dict(annotation_id=region['annotation_id'],region_id=region['region_id'],gene=region['gene'],source=region['source'],
                    current_plan=region['current_plan'],peptide=peptide,length=length,start=start,end=end,
                    altered_start_in_peptide=max(0,lo-start),altered_end_in_peptide=min(length,hi-start),
                    reference_peptide=reference[start:end] if reference is not None else '',
                    binding_core_alteration_status='NOT_ASSIGNED_BY_MODEL_WHOLE_PEPTIDE_OVERLAP_ONLY')
                if region['source'] == 'CTA/SELF':
                    row['altered_start_in_peptide'] = row['altered_end_in_peptide'] = ''
                    row['binding_core_alteration_status'] = 'INTENTIONALLY_REFERENCE_CTA_SELF'
                occurrences.append(row)
                queries.add(peptide)
                if row['reference_peptide']:
                    require(row['reference_peptide'] != peptide, 'Unchanged focal counterfactual')
                    queries.add(row['reference_peptide'])
                count += 1
            coverage.append(dict(annotation_id=region['annotation_id'],region_id=region['region_id'],length=length,occurrences=count))
    source_meta = load(ROOT/'data/shared-rna/splicing/independent-assessment/protein-screen/protein-screen.json')
    sources = []
    for record in source_meta['sources']:
        name = record['name']
        path = record.get('path')
        if name == 'ensembl103':
            path = '/data/raft/references/homo_sapiens/pyensembl/GRCh38/ensembl103/Homo_sapiens.GRCh38.pep.all.fa.gz'
        elif name == 'personal_normal':
            path = '/data/results/TN25-112103/matched-normal-20261107/normal-self/normal-sequence-contexts.fasta'
        require(path, 'Missing normal/reference catalog path')
        sources.append(dict(name=name,path=path,sha256=record['sha256'],expected_records=record['records'],
            scope=record.get('scope',{})))
    paths = [CLASS_I/'request.json',CLASS_I/'current-plan.tsv',CLASS_I/'region-annotations.tsv',
        CLASS_I/'verification.json',TYPING,MANIFEST,
        ROOT/'data/shared-rna/splicing/independent-assessment/protein-screen/protein-screen.json']
    return dict(schema_version=2,regions=parent['regions'],occurrences=occurrences,
        queries=sorted(queries),length_coverage=coverage,targets=typing['targets'],
        normal_reference_sources=sources,model_manifest=load(MANIFEST),
        source_sha256={str(p.relative_to(ROOT)):sha(p) for p in paths},
        policy=dict(peptide_lengths=list(range(15,21)),ensemble_folds=5,model_modes=['BA','EL'],
            primary_targets=['DRB1*03:01','DRB1*08:02'],
            exploratory_targets='Unphased observed DQ/DP alpha-by-beta combinations',
            excluded_targets=[t['id'] for t in typing['excluded_targets']],
            DR_model_equivalence=typing['DR_model_equivalence'],
            matched_normal_corroboration='Contributed PR #7 calls agree; that caller was not rerun here',
            counterfactual='Focal annotated reference residue restored; other background differences unchanged',
            binding_core_assigned=False,RNA_count_or_expression_weighting=False,hard_score_threshold=None,
            existing_personal_normal_catalog_max_peptide_length=11,
            personal_normal_negative_is_complete_for_classII=False,
            no_hit_is_not_proof_of_normal_tissue_absence=True,
            selection_changed=False,construct_rebuilt=False,vaccine_context_prediction=False))


def products(req,predictions,coverage,hits):
    require(len({r['annotation_id'] for r in req['regions']})==len(req['regions']), 'Duplicate annotation identities')
    expected={(p,t['id']) for p in req['queries'] for t in req['targets']}
    values={(p['peptide'],p['target_id']):p for p in predictions}
    require(len(values) == len(predictions) and set(values) == expected, 'Missing/duplicate class-II predictions')
    targets={x['id']:x for x in coverage}
    require(set(targets) == {t['id'] for t in req['targets']}, 'Incomplete model allele coverage')
    require(set(hits) == set(req['queries']), 'Incomplete long-peptide self search')
    for x in predictions:
        require(x['canonical_allele'] == targets[x['target_id']]['canonical_allele'], 'Wrong allele attribution')
        require(x['model_allele_key'] == targets[x['target_id']]['model_allele_key'], 'Wrong model key')
        nm,ba,el = (float(x[k]) for k in ['affinity_nm','affinity_score','presentation_score'])
        require(all(math.isfinite(v) for v in [nm,ba,el]) and nm>0 and 0<=el<=1, 'Invalid class-II prediction')
        require(math.isclose(nm,50000**(1-ba),rel_tol=1e-8), 'Changed BA transform')
    peptide_rows=[]
    grouped=defaultdict(list)
    for occurrence in req['occurrences']:
        peptide=occurrence['peptide']
        source_hits=hits[peptide]
        reference_hit=any(source_hits.get(n) for n in ['ensembl103','ensembl115','uniprot_canonical','uniprot_additional'])
        personal_hit=bool(source_hits.get('personal_normal'))
        for target in req['targets']:
            pred=values[peptide,target['id']]
            row=dict(**occurrence,target_id=target['id'],evidence_group=pred['evidence_group'],
                canonical_allele=pred['canonical_allele'],model_allele_key=pred['model_allele_key'],
                model_pseudosequence_X_tokens=targets[target['id']]['model_pseudosequence_X_tokens'],
                typing_status=target['status'],pairing_status=target['phase_status'],
                affinity_nm=float(pred['affinity_nm']),affinity_score=float(pred['affinity_score']),
                presentation_score=float(pred['presentation_score']),exact_reference_match=reference_hit,
                scoped_personal_normal_match=personal_hit,
                exact_self_search_hit=reference_hit or personal_hit,
                personal_normal_negative_is_complete_for_classII=False,
                reference_affinity_nm='',mutant_over_reference_affinity='',reference_presentation_score='')
            if occurrence['reference_peptide']:
                ref=values[occurrence['reference_peptide'],target['id']]
                row.update(reference_affinity_nm=float(ref['affinity_nm']),
                    mutant_over_reference_affinity=row['affinity_nm']/float(ref['affinity_nm']),
                    reference_presentation_score=float(ref['presentation_score']))
            peptide_rows.append(row)
            grouped[row['annotation_id'],row['target_id']].append(row)
    by_target=[]
    for region in req['regions']:
        for target in req['targets']:
            rows=grouped[region['annotation_id'],target['id']]
            require(rows, 'Region has no class-II peptides')
            ba=min(rows,key=lambda x:(x['affinity_nm'],x['peptide']))
            el=max(rows,key=lambda x:(x['presentation_score'],-x['affinity_nm'],x['peptide']))
            no_hits=[x for x in rows if not x['exact_self_search_hit']]
            nh=max(no_hits,key=lambda x:(x['presentation_score'],-x['affinity_nm']),default=None)
            by_target.append(dict(annotation_id=region['annotation_id'],region_id=region['region_id'],window_id=region['window_id'],gene=region['gene'],source=region['source'],
                current_plan=region['current_plan'],target_id=target['id'],evidence_group=ba['evidence_group'],
                peptides_scored=len(rows),exact_reference_or_scoped_personal_matching_peptides=sum(x['exact_self_search_hit'] for x in rows),
                best_affinity_nm=ba['affinity_nm'],best_affinity_peptide=ba['peptide'],
                best_presentation_score=el['presentation_score'],best_presentation_peptide=el['peptide'],
                best_presentation_peptide_reference_matching=el['exact_reference_match'],
                best_presentation_peptide_scoped_personal_matching=el['scoped_personal_normal_match'],
                best_no_exact_search_hit_presentation=nh['presentation_score'] if nh else '',
                best_no_exact_search_hit_peptide=nh['peptide'] if nh else '',
                binding_core_alteration_assessed=False))
    summaries=[]
    for region in req['regions']:
        row={k:region[k] for k in ['annotation_id','region_id','window_id','gene','source','protein_change','sequence','current_plan']}
        for group,prefix in [('DRB1_primary','classII_DRB1'),('exploratory','classII_exploratory')]:
            subset=[x for x in by_target if x['annotation_id']==region['annotation_id'] and x['evidence_group']==group]
            ba=min(subset,key=lambda x:(x['best_affinity_nm'],x['target_id']))
            el=max(subset,key=lambda x:(x['best_presentation_score'],-x['best_affinity_nm']))
            no_hits=[x for x in subset if x['best_no_exact_search_hit_presentation']!='']
            nh=max(no_hits,key=lambda x:x['best_no_exact_search_hit_presentation'],default=None)
            row.update({prefix+'_best_affinity_nm':ba['best_affinity_nm'],prefix+'_best_affinity_peptide':ba['best_affinity_peptide'],
                prefix+'_best_affinity_target':ba['target_id'],prefix+'_best_presentation_score':el['best_presentation_score'],
                prefix+'_best_presentation_peptide':el['best_presentation_peptide'],prefix+'_best_presentation_target':el['target_id'],
                prefix+'_best_presentation_exact_reference_match':el['best_presentation_peptide_reference_matching'],
                prefix+'_best_no_exact_search_hit_presentation':nh['best_no_exact_search_hit_presentation'] if nh else '',
                prefix+'_best_no_exact_search_hit_peptide':nh['best_no_exact_search_hit_peptide'] if nh else '',
                prefix+'_best_no_exact_search_hit_target':nh['target_id'] if nh else ''})
        summaries.append(row)
    classI={r['annotation_id']:r for r in table(CLASS_I/'region-annotations.tsv')}
    current=[{**classI[r['annotation_id']],**r} for r in summaries if r['current_plan']]
    require(len(current)==30,'Wrong current-plan region count')
    actual={o['peptide'] for o in req['occurrences']}
    altered={o['peptide'] for o in req['occurrences'] if o['source']!='CTA/SELF'}
    matching={p for p in actual if any(hits[p].values())}
    summary=dict(status='CLASSII_ANNOTATIONS_COMPLETE_SELECTION_UNCHANGED',regions=len(summaries),current_plan_regions=30,
        distinct_region_ids=len({r['region_id'] for r in summaries}),annotation_identity='window_id for SPLICE, otherwise region_id',
        eligible_splice_windows=sum(r['source']=='SPLICE' for r in summaries),peptide_occurrences=len(req['occurrences']),
        unique_target_peptides=len(actual),focal_reference_counterpart_queries=len(set(req['queries'])-actual),
        distinct_peptides_inferred=len(req['queries']),classII_targets=len(req['targets']),peptide_target_predictions=len(predictions),
        primary_DRB1_targets=len(req['policy']['primary_targets']),
        exploratory_targets=len(req['targets'])-len(req['policy']['primary_targets']),
        excluded_targets=req['policy']['excluded_targets'],peptide_occurrence_target_rows=len(peptide_rows),
        target_peptides_with_exact_reference_or_scoped_personal_hits=len(matching),
        altered_target_peptides_with_exact_search_hits=len(altered & matching),
        current_plan_altered_regions_without_any_no_exact_search_hit_peptide=[r['gene'] for r in current
            if r['source']!='CTA/SELF' and r['classII_DRB1_best_no_exact_search_hit_peptide']==''],
        binding_core_assigned=False,hard_score_cutoff=None,expression_or_RNA_counts_used=False,
        classI_and_classII_scores_jointly_calibrated=False,selection_changed=False,construct_rebuilt=False,
        personal_normal_negative_is_complete_for_classII=False)
    return {'peptide-annotations.tsv.gz':zipped(tsv_bytes(peptide_rows)),
        'region-target-annotations.tsv':tsv_bytes(by_target),'region-annotations.tsv':tsv_bytes(summaries),
        'current-plan.tsv':tsv_bytes(current),'peptide-length-coverage.tsv':tsv_bytes(req['length_coverage']),
        'summary.json':json_bytes(summary)},summary


def finish(out,verify=False):
    out=Path(out)
    req=load(out/'request.json')
    require(req==request(),'Class-II inputs changed; use a fresh output directory')
    meta=load(out/'inference.json')
    require(meta['status']=='COMPLETED' and meta['new_inference'],'Unverified class-II inference')
    require(meta['request_sha256']==sha(out/'request.json'),'Inference/request mismatch')
    require(meta['worker_sha256']==sha(Path(__file__).with_name('class_ii_worker.py')),'Changed measured worker')
    require(meta['all_upstream_model_files_verified'] and meta['reference_smoke_passed'],'Unverified model execution')
    for name,expected in meta['output_sha256'].items():
        require(sha(out/name)==expected,'Changed measured class-II output: '+name)
    require(load(out/'normal-search.json')['personal_normal_catalog_designed_max_peptide_length']==11,
        'Personal-normal search scope changed')
    outputs,summary=products(req,table(out/'predictions.tsv.gz'),table(out/'hla-coverage.tsv'),load(out/'reference-matches.json.gz'))
    if verify:
        receipt=load(out/'verification.json')
        require(receipt['source_sha256']==req['source_sha256'],'Changed source identities')
        for name,body in outputs.items():
            require((out/name).read_bytes()==body,'Class-II annotation replay differs: '+name)
        for name,digest in receipt['output_sha256'].items():
            require(sha(out/name)==digest,'Changed class-II output: '+name)
    else:
        for name,body in outputs.items():(out/name).write_bytes(body)
        files=[*outputs,'request.json','inference.json',*meta['output_sha256']]
        receipt=dict(status='PASS',source_sha256=req['source_sha256'],
            output_sha256={name:sha(out/name) for name in files},analysis_sha256=sha(__file__),
            every_occurrence_target_replayed=True,selection_changed=False)
        (out/'verification.json').write_bytes(json_bytes(receipt))
    return summary


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output',type=Path,default=DEFAULT)
    mode=p.add_mutually_exclusive_group()
    mode.add_argument('--execute',action='store_true')
    mode.add_argument('--collect',action='store_true')
    mode.add_argument('--verify',action='store_true')
    a=p.parse_args();out=a.output.resolve()
    if a.verify:
        print(json.dumps(finish(out,True),indent=2));return
    req=request()
    print(json.dumps(dict(regions=len(req['regions']),current_plan=30,
        peptide_queries=len(req['queries']),targets=len(req['targets']),
        peptide_target_pairs=len(req['queries'])*len(req['targets']),
        execute=a.execute,selection_changed=False),indent=2),flush=True)
    if not a.execute and not a.collect:return
    import modal
    worker=Path(__file__).with_name('class_ii_worker.py')
    volume=modal.Volume.from_name('research-lens-workspace')
    if a.execute:
        require(not (out/'request.json').exists(),'Choose a fresh class-II output directory')
        out.mkdir(parents=True,exist_ok=True);(out/'request.json').write_bytes(json_bytes(req))
        stamp=datetime.datetime.now(datetime.timezone.utc).strftime('%Y%m%dT%H%M%SZ')
        remote='/data/results/TN25-112103/class-II-annotations-'+sha(out/'request.json')[:12]+'-'+sha(worker)[:12]+'-'+stamp
        with volume.batch_upload() as upload:
            upload.put_file(str(out/'request.json'),remote.removeprefix('/data/')+'/request.json')
            upload.put_file(str(worker),remote.removeprefix('/data/')+'/worker.py')
        image=modal.Image.from_id('im-DBMghmAQfB3Oz7McdNZeic').pip_install('pyahocorasick==2.2.0')
        sandbox=modal.Sandbox.create('python','-u',remote+'/worker.py','--request',remote+'/request.json',
            '--models',MODEL_DIR,'--output',remote+'/output',
            app=modal.App.lookup('research-yash-vaccine',create_if_missing=True),image=image,
            gpu='A10',cpu=4,memory=16384,timeout=3600,volumes={'/data':volume})
        state=dict(status='RUNNING',started_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
            remote=remote,sandbox_id=sandbox.object_id,request_sha256=sha(out/'request.json'),worker_sha256=sha(worker),terminated=False)
        (out/'execution.json').write_bytes(json_bytes(state));print(json.dumps(state,indent=2),flush=True);return
    state=load(out/'execution.json')
    require(state['request_sha256']==sha(out/'request.json') and state['worker_sha256']==sha(worker),'Changed live run inputs')
    sandbox=modal.Sandbox.from_id(state['sandbox_id']);rc=sandbox.poll()
    print('POLL',state['sandbox_id'],rc,flush=True)
    if rc is None:return
    state['exit_code']=rc
    (out/'stdout.log').write_text(sandbox.stdout.read());(out/'stderr.log').write_text(sandbox.stderr.read())
    try:
        if rc:raise RuntimeError('Class-II inference failed: '+(out/'stderr.log').read_text()[-5000:])
        prefix=state['remote'].removeprefix('/data/')+'/output/'
        meta=json.loads(b''.join(volume.read_file(prefix+'inference.json')))
        require(meta['worker_sha256']==sha(worker),'Wrong remote worker')
        for name in ['inference.json',*meta['output_sha256']]:
            (out/name).write_bytes(b''.join(volume.read_file(prefix+name)))
        print(json.dumps(finish(out),indent=2),flush=True);state['status']='COMPLETED'
    finally:
        sandbox.terminate(wait=True);state['terminated']=True
        state['finished_utc']=datetime.datetime.now(datetime.timezone.utc).isoformat()
        (out/'execution.json').write_bytes(json_bytes(state))


if __name__=='__main__':main()
