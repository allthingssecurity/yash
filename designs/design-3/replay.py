"""Verify and replay the final 30-region construct locally; no cloud or new inference."""
from pathlib import Path
import argparse
import csv
import hashlib
import importlib.util
import json
import shutil
import subprocess
import sys

HERE=Path(__file__).resolve().parent
REPO=HERE.parents[1]
sys.path.insert(0,str(HERE/'pipeline'))
sys.path.insert(0,str(REPO/'src'))
from yash_vaccine.hydrophobicity import protein_screen
from yash_vaccine.io import require_assertions


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def load(path):return json.loads(Path(path).read_text())
def rows(path):
    with Path(path).open() as stream:return list(csv.DictReader(stream,delimiter='\t'))
def module(path,name):
    spec=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m


def recorded_input(name):
    # Native receipts retain the paths used at execution. Only local reads move.
    prefix = 'designs/design-3/final/'
    if name == prefix + 'policy.json':
        return HERE / 'execution/recorded-source/policy.json'
    return HERE / name[len(prefix):] if name.startswith(prefix) else REPO / name


def main():
    require_assertions('Final thirty-region replay')
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    output=a.output.resolve();assert not output.exists(),'Choose a new replay output directory'
    subprocess.run([sys.executable,str(HERE/'selection.py'),'--verify'],check=True)
    manifest=load(HERE/'execution/inputs.manifest.json')
    for _,r in manifest.items():
        path=recorded_input(r['repository_path']);assert path.stat().st_size==r['bytes'] and sha(path)==r['sha256'],str(path)
    first_manifest=HERE/'execution/source-v1/inputs.manifest.json'
    if first_manifest.exists():
        for key,r in load(first_manifest).items():
            path=HERE/'execution/source-v1/cassette.py' if key=='code/cassette.py' else recorded_input(r['repository_path'])
            assert path.stat().st_size==r['bytes'] and sha(path)==r['sha256'],str(path)
    worker_hashes={sha(HERE/'pipeline/cassette.py')}
    if (HERE/'execution/source-v1/cassette.py').exists():worker_hashes.add(sha(HERE/'execution/source-v1/cassette.py'))
    hash_only=[];checked=0
    for folder in ['prepare','all-inference','protein','selected-inference','verify','mrna']:
        path=HERE/'construct'/folder;receipt=load(path/'complete.json');assert receipt['status']=='COMPLETED'
        assert receipt['worker_sha256'] in worker_hashes
        for name,r in receipt['outputs'].items():
            f=path/name
            if f.is_file():assert f.stat().st_size==r['bytes'] and sha(f)==r['sha256'],str(f);checked+=1
            else:
                assert folder=='all-inference' and name=='predictions.tsv' and r['bytes']>40_000_000
                hash_only.append(dict(file=str(f.relative_to(REPO)),**r,verification='HASH_ONLY_NATIVE_OUTPUT_NOT_DOWNLOADED'))
    selected=load(HERE/'selection/selected.json');n=len(selected);assert n==30
    h2=next(r for r in selected if r['gene']=='H2BC12')
    h2_source=HERE/'inputs/h2bc12'
    newer=h2_source;check=load(newer/'result.json');done=load(newer/'complete.json')
    assert check['status']=='PASS_NO_NEWER_REFERENCE_HITS' and not check['matches']
    assert check['region_id']==h2['region_id'] and set(check['peptides'])==set(h2['scored_target_peptides'])
    assert check['request_sha256']==sha(h2_source/'request.json')
    assert check['worker_sha256']==sha(h2_source/'check_reference.py')
    assert sha(newer/'result.json')==done['outputs']['result.json']['sha256']
    cassette=module(HERE/'pipeline/cassette.py','reviewed_cassette')
    optimizer=module(REPO/'designs/initial/construct/optimize.py','reviewed_joint_optimizer')
    costs=rows(HERE/'construct/protein/edge-costs.tsv')
    for e in costs:
        for k in ['left_index','right_index','strong_rank','mild_rank','affinity_50','affinity_500','added_aa','peptide_occurrences']:e[k]=int(e[k])
        for k in ['min_rank','min_affinity_nm','max_pair_mean_kd']:e[k]=float(e[k])
        e['hydrophobicity_pass']=e['hydrophobicity_pass']=='True'
    assert len(costs)==n*(n+1)*11 and len({(e['left_index'],e['right_index'],e['linker']) for e in costs})==len(costs)
    path,solver=optimizer.solve(n,[e for e in costs if e['hydrophobicity_pass']],seconds=600)
    assert solver['all_phases_proven_optimal']
    architecture=load(HERE/'construct/architecture.json')
    sequence,parts,joins=cassette.assemble(path,selected,architecture['elements']['SEC']['aa'],architecture['elements']['MITD']['aa'])
    output.mkdir(parents=True);protein=output/'construct/protein';protein.mkdir(parents=True)
    (protein/'cassette.fasta').write_text('>Yash_design_3_final_SEC_antigens_MITD\n'+sequence+'\n')
    cassette.table(protein/'components.tsv',parts)
    for name in ['cassette.fasta','components.tsv']:assert (protein/name).read_bytes()==(HERE/'construct/protein'/name).read_bytes(),name
    _,hydro=protein_screen(sequence,parts);assert hydro['status']=='PASS'
    request=load(HERE/'construct/protein/verification-request.json');alleles=request['alleles']
    fresh={(r['peptide'],r['allele']):(float(r['affinity_nm']),float(r['affinity_percentile'])) for r in rows(HERE/'construct/selected-inference/predictions.tsv')}
    assert set(fresh)=={(p,a) for p in request['peptides'] for a in alleles}
    inference=load(HERE/'construct/selected-inference/inference.json')
    assert inference['input_sha256']==sha(HERE/'construct/protein/verification-request.json')
    assert inference['predictions_sha256']==sha(HERE/'construct/selected-inference/predictions.tsv')
    _,totals,_=cassette.assess(sequence,joins,path,alleles,lambda p,h:fresh[p,alleles[h]])
    recorded=load(HERE/'construct/final-check.json');assert recorded['protein_sha256']==sha(protein/'cassette.fasta')
    for k in ['peptide_occurrences','peptide_hla_rows','strong_rank','mild_rank','affinity_50','affinity_500']:assert totals[k]==recorded['totals'][k],k
    mrna=output/'construct/mrna';mrna.mkdir()
    for name in ['run.py','genbank.py','optimizer.py','backbone.json','annotation-sources.json','human-codon-usage.json','requirements.txt']:shutil.copyfile(HERE/'construct/mrna'/name,mrna/name)
    (output/'selection').mkdir();shutil.copyfile(HERE/'selection/selected.tsv',output/'selection/selected.tsv')
    for name in ['architecture.json','final-check.json']:shutil.copyfile(HERE/'construct'/name,output/'construct'/name)
    with (output/'encoding.stdout.log').open('w') as stdout,(output/'encoding.stderr.log').open('w') as stderr:
        subprocess.run([sys.executable,str(mrna/'run.py'),'--output',str(mrna)],stdout=stdout,stderr=stderr,check=True)
    identical=['cds.fasta','mrna-segments.fasta','template-segments.fasta','template-segments.gb','components.tsv']
    for name in identical:assert (mrna/name).read_bytes()==(HERE/'construct/mrna'/name).read_bytes(),name
    v=load(mrna/'verification.json');assert v['all_constraints_pass'] and v['polyA_status']=='unspecified' and not v['complete_construct_verified']
    result=dict(status='DESIGN3_FINAL_DESIGN_OFFLINE_REPLAY_PASSED',regions=n,protein_length_aa=len(sequence),
        selection_and_RNA_admission_verified=True,all_solver_phases_optimal=True,protein_and_components_identical=True,
        known_nucleotide_and_GenBank_outputs_identical=identical,included_output_files_byte_verified=checked,
        hash_only_native_outputs=hash_only,new_inference_performed=False,BAMs_reread=False,
        retained_parent_regions=24,added_H2BC12=True,H2BC12_newer_reference_check_verified=True,
        RNA_genes=load(HERE/'policy.json')['RNA_genes'],PolyA='unspecified',hydrophobicity=hydro)
    (output/'replay.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))


if __name__=='__main__':main()
