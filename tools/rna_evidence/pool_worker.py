"""Screen any provided reconstructed RNA windows before antigen selection.

The historic four-window worker remains unchanged. This stage adds caching,
all recorded alignment paths, two protein catalogs, native-score parity and
explicit normal-evidence tiers. It does not realign reads or infer new HLA scores.
"""
from pathlib import Path
from collections import Counter,defaultdict
import argparse,csv,datetime,gzip,hashlib,importlib.metadata,io,json,math,sys,urllib.request,uuid
from concurrent.futures import ThreadPoolExecutor
from Bio.Seq import Seq
from worker import sha,save,translation_records,library,exact_hit,direct_counts,primary_coding
AA=set('ACDEFGHIKLMNPQRSTVWY')


def read_json(path):
    with (gzip.open(path,'rt') if str(path).endswith('.gz') else open(path)) as f:return json.load(f)


def write_json_gz(path,value):
    with gzip.open(path,'wt') as f:json.dump(value,f,default=str)


def table(path,rows):
    fields=list(dict.fromkeys(k for r in rows for k in r))
    with (gzip.open(path,'wt') if str(path).endswith('.gz') else open(path,'w')) as f:
        w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',lineterminator='\n');w.writeheader()
        for r in rows:w.writerow({k:json.dumps(v) if isinstance(v,(list,dict)) else v for k,v in r.items()})


def fetch_reference(meta,cache):
    cache.mkdir(parents=True,exist_ok=True)
    def get(component):
        path=cache/Path(component['url']).name
        if not path.exists():
            temporary=path.with_name(path.name+'.'+uuid.uuid4().hex+'.partial')
            with urllib.request.urlopen(component['url'],timeout=60) as response,temporary.open('wb') as f:
                while True:
                    data=response.read(1024*1024)
                    if not data:break
                    f.write(data)
            assert temporary.stat().st_size==component['bytes'] and sha(temporary)==component['sha256']
            temporary.replace(path)
        assert path.stat().st_size==component['bytes'] and sha(path)==component['sha256']
        return path
    with ThreadPoolExecutor(max_workers=2) as pool:paths=list(pool.map(get,meta['component_files']))
    # Hash the concatenated gzip members without copying the reference data.
    h=hashlib.sha256()
    for path in paths:
        with path.open('rb') as f:
            for block in iter(lambda:f.read(1024*1024),b''):h.update(block)
    assert h.hexdigest()==meta['sha256']
    records=[]
    for path in paths:records.extend(translation_records(path))
    return records


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root',type=Path,required=True);p.add_argument('--request',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();root=a.root;out=a.output
    if (out/'complete.json').exists():raise FileExistsError('Preserve a completed run; choose a new run ID')
    out.mkdir(parents=True,exist_ok=True);request=read_json(a.request);policy=request['policy']
    assert {k:policy[k] for k in ['normal_min_fragments','normal_min_GQ','normal_ALT_below','minimum_mapping_quality','minimum_base_quality']} == dict(normal_min_fragments=20,normal_min_GQ=20,normal_ALT_below=3,minimum_mapping_quality=20,minimum_base_quality=20), 'Change the evidence helper explicitly for different normal/read thresholds'
    def progress(stage,**kw):save(out/'status.json',dict(stage=stage,updated_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),**kw));print(stage,kw,flush=True)
    import pysam,varcode
    from pyensembl import Genome
    from normal_evidence import NormalEvidence
    from self_search import patient_self_matches
    from mhctools.pred import Kind,Prediction
    from vaxrank.candidate_epitope import CandidateEpitope
    from vaxrank.epitope_config import EpitopeConfig
    from vaxrank.epitope_dsl import attach_per_allele_scores,epitopes_for_ranking
    sys.path.insert(0,str(root/'frozen/src'))
    import yash_vaccine.hydrophobicity as hydropathy
    assert importlib.metadata.version('vaxrank')=='3.46.0'
    identity=read_json(root/'identity/complete.json');assert identity['identity']=='MATCH'
    normal=NormalEvidence(root);plan=read_json(root/'normal-coding-plan.json');done=read_json(root/'normal-coding/gather/complete.json')
    assert Path(plan['normal_bam']).stat().st_size==143294724201
    assert plan['normal_bam_sha256']==done['normal_bam_sha256']
    for name in ['normal.coding.g.vcf.gz','normal.coding.g.vcf.gz.tbi']:assert sha(root/'normal-coding/gather'/name)==done['outputs'][name]['sha256']
    genome=Genome(**read_json(root/'analysis-genome-config.json'))
    tumour=pysam.AlignmentFile('/data/inputs/TN25-112103/DNA_TN25-112103.bam','rb')
    rna=pysam.AlignmentFile('/data/inputs/TN25-112103/RNA_TN25-112103.bam','rb');tumour.check_index();rna.check_index()
    progress('reference_loading',windows=len(request['windows']))
    fasta=Path('/data/raft/references/homo_sapiens/pyensembl/GRCh38/ensembl103/Homo_sapiens.GRCh38.pep.all.fa.gz')
    assert sha(fasta)==request['reference_ensembl_sha256']
    ensembl_records=list(translation_records(fasta));broad=library(ensembl_records);primary=library([r for r in ensembl_records if primary_coding(r[0])])
    uni_records=fetch_reference(request['reference_uniprot'],Path('/data/results/TN25-112103/provider-cache/uniprot')/request['reference_uniprot']['sha256'][:16]);uni=library(uni_records)
    windows=[];site_cache={};wanted=set();comp=str.maketrans('ACGTN','TGCAN')
    for n,w in enumerate(request['windows'],1):
        assert len(w['sequence'])==25 and len(w['coding_sequence'])==75 and str(Seq(w['coding_sequence']).translate())==w['sequence']
        chrom=w['coordinates'].split(':')[0];minus=w['coordinates'].endswith('(-)');path_evidence=[]
        paths=w['genomic_path_support'];assert len(paths)==w['genomic_path_count']
        for path_index,path in enumerate(paths):
            positions=path['genomic_positions_0based'];assert len(positions)==75
            bases=[]
            for pos in positions:
                base=normal.fasta.fetch(chrom,pos,pos+1).upper() if pos is not None else 'N'
                bases.append(base.translate(comp) if minus else base)
            reference=''.join(bases);reference_aa=str(Seq(reference).translate());differences=[]
            for i,(ref_rna,alt_rna) in enumerate(zip(reference,w['coding_sequence'])):
                if ref_rna==alt_rna:continue
                aa=i//3;position=positions[i]
                if position is None:
                    differences.append(dict(window_nt_position_1based=i+1,aa_position_1based=aa+1,nonsynonymous=True,normal_classification='UNRESOLVED_ALLELE_PROJECTION',site_key=None));continue
                ref=ref_rna.translate(comp) if minus else ref_rna;alt=alt_rna.translate(comp) if minus else alt_rna
                key=':'.join(map(str,[chrom,position+1,ref,alt]))
                if key not in site_cache:
                    variant=varcode.Variant(chrom,position+1,ref,alt,ensembl=genome);measured=normal.assess([chrom,position+1,ref,alt],variant)
                    counts={name:direct_counts(bam,chrom,position+1) for name,bam in [('normal',normal.bam),('tumour_DNA',tumour),('tumour_RNA',rna)]}
                    agreement=counts['normal']['fragment_bases'].get(alt,0)==measured['read_evidence']['fragment_counts'].get('ALT',0)
                    site_cache[key]=dict(allele=[chrom,position+1,ref,alt],normal=measured,independent_counts=counts,independent_ALT_count_agreement=agreement,audit_normal_classification=measured['classification'] if agreement else 'UNRESOLVED_NORMAL_EVIDENCE')
                differences.append(dict(window_nt_position_1based=i+1,aa_position_1based=aa+1,reference_codon=reference[3*aa:3*aa+3],rna_codon=w['coding_sequence'][3*aa:3*aa+3],reference_amino_acid=reference_aa[aa],RNA_inferred_amino_acid=w['sequence'][aa],nonsynonymous=reference_aa[aa]!=w['sequence'][aa],normal_classification=site_cache[key]['audit_normal_classification'],site_key=key))
            path_evidence.append(dict(path_index=path_index,reference_translation_on_recorded_path=reference_aa,differences=differences,fragment_names=path['fragment_names']))
        lo,hi=int(w['target_start_in_window']),int(w['target_end_in_window']);occurrences=[]
        for k in policy['peptide_lengths']:
            for i in range(len(w['sequence'])-k+1):
                if i<hi and i+k>lo:
                    peptide=w['sequence'][i:i+k];occurrences.append(dict(peptide=peptide,start=i,end=i+k));wanted.add(peptide)
                    for path in path_evidence:
                        counterfactual=path['reference_translation_on_recorded_path'][i:i+k]
                        if set(counterfactual)<=AA:wanted.add(counterfactual)
        windows.append(dict(**w,path_evidence=path_evidence,occurrences=occurrences))
        if n%50==0:progress('normal_DNA_comparison',windows_processed=n,unique_sites=len(site_cache))
    progress('exact_reference_search',queries=len(wanted),unique_normal_sites=len(site_cache))
    personal=patient_self_matches(root,wanted)
    matches={}
    for n,peptide in enumerate(sorted(wanted),1):
        hit_primary=exact_hit(primary,peptide);hit_broad=exact_hit(broad,peptide);hit_uni=exact_hit(uni,peptide)
        is_self=bool(hit_broad or hit_uni or personal[peptide]['occurs_in_searched_models'])
        matches[peptide]=dict(primary=hit_primary,ensembl_all=hit_broad,uniprot=hit_uni,personal=personal[peptide],self_match=is_self)
        if n%500==0:progress('exact_reference_search',queries_completed=n,total=len(wanted))
    # Use recorded predictions; missing pairs or changed input hashes stop here.
    affinities=request['affinities'];alleles=request['alleles'];actual_peptides={o['peptide'] for w in windows for o in w['occurrences']}
    assert set(affinities)==actual_peptides
    unscored=[]
    for peptide in sorted(actual_peptides):
        assert set(affinities[peptide])==set(alleles)
        preds=tuple(Prediction(kind=Kind.pMHC_affinity,peptide=peptide,allele=allele,value=affinities[peptide][allele],score=affinities[peptide][allele],predictor_name='mhcflurry',predictor_version='2.3.14+models2.3.0') for allele in alleles)
        unscored.append(CandidateEpitope(sequence=peptide,source_sequence=peptide,offset=0,predictions=preds,overlaps_targetable=True,patient_alleles=tuple(alleles)))
    ec=EpitopeConfig();native=epitopes_for_ranking(attach_per_allele_scores(unscored,ec),ec);contributions={e.sequence:float(e.epitope_score) for e in native}
    assert set(contributions)<=actual_peptides
    # Native ranking omits peptides with exactly zero affinity contribution.
    # Verify that case independently rather than treating it as missing inference.
    for peptide in actual_peptides-set(contributions):
        values=[0 if nm>=5000 else (1+math.exp(-350/150))/(1+math.exp((nm-350)/150)) for nm in affinities[peptide].values()]
        expected=sum(value if value>=1e-5 else 0 for value in values)
        assert expected==0, ('Native ranking omitted a positive-score peptide',peptide,expected)
        contributions[peptide]=0.0
    decisions=[];peptide_rows=[];evidence=[]
    for w in windows:
        baseline=w['baseline'];old=sum(contributions[o['peptide']] for o in w['occurrences'] if not request['original_native_reference_matches'][o['peptide']])
        assert math.isclose(old,float(baseline['native_affinity_sum']),rel_tol=1e-8,abs_tol=1e-8),(w['window_id'],old,baseline['native_affinity_sum'])
        seeds=set(baseline['seed_peptides'].split(';'))- {''};review_occ=[];primary_occ=[];normal_flags=Counter();basis_counts=Counter()
        for o in w['occurrences']:
            peptide=o['peptide'];path_tiers=[];path_bases=[]
            if matches[peptide]['self_match']:tier='EXCLUDE_SELF';basis='REFERENCE_OR_PATIENT_SELF_MATCH'
            else:
                for path in w['path_evidence']:
                    p_ref=path['reference_translation_on_recorded_path'][o['start']:o['end']]
                    differences=[d for d in path['differences'] if d['nonsynonymous'] and o['start']<=d['aa_position_1based']-1<o['end']]
                    if set(p_ref)<=AA and not matches[p_ref]['self_match']:
                        path_tiers.append('PRIMARY');path_bases.append('JUNCTION_CONTEXT_ALREADY_NONSELF')
                    elif differences and all(d['normal_classification'] in {'NORMAL_REFERENCE_SUPPORTED','PRESENT_IN_NORMAL'} for d in differences) and any(d['normal_classification']=='NORMAL_REFERENCE_SUPPORTED' for d in differences):
                        path_tiers.append('PRIMARY');path_bases.append('RNA_DIFFERENCE_DEPENDENT_NORMAL_REFERENCE_SUPPORTED')
                    else:
                        path_tiers.append('REVIEW');path_bases.append('RNA_DIFFERENCE_DEPENDENT_NORMAL_PRESENT_OR_UNRESOLVED')
                    for d in differences:normal_flags[d['normal_classification']]+=1
                tier='PRIMARY' if path_tiers and all(t=='PRIMARY' for t in path_tiers) else 'REVIEW'
                basis=';'.join(sorted(set(path_bases)))
                review_occ.append(o)
                if tier=='PRIMARY':primary_occ.append(o)
            basis_counts[basis]+=1
            best=min(affinities[peptide],key=affinities[peptide].get)
            peptide_rows.append(dict(window_id=w['window_id'],region_id=w['region_id'],gene=w['gene'],peptide=peptide,start=o['start'],end=o['end'],tier=tier,novelty_basis=basis,native_contribution=contributions[peptide],native_PASS_seed=peptide in seeds,best_allele=best,best_affinity_nm=affinities[peptide][best]))
        value=float(w['supplied_transformed_gene_tpm']);hydro=hydropathy.summarize(w['sequence'])
        def scored(occ):
            best=min((dict(peptide=o['peptide'],allele=a,affinity_nm=affinities[o['peptide']][a]) for o in occ for a in alleles),key=lambda r:r['affinity_nm'],default=None)
            score=sum(contributions[o['peptide']] for o in occ)
            admitted=any(o['peptide'] in seeds and min(affinities[o['peptide']].values())<=policy['admission_affinity_nm'] for o in occ)
            eligible=bool(admitted and score>0 and value>0 and hydro['status']=='PASS')
            return dict(score=score,weighted_score=score*math.sqrt(value) if value>0 else None,best=best,eligible=eligible,seed_pass=admitted)
        primary_score=scored(primary_occ);review_score=scored(review_occ)
        if primary_score['eligible']:tier='PRIMARY'
        elif review_score['eligible']:tier='REVIEW_NORMAL_EVIDENCE'
        else:tier='EXCLUDE'
        reasons=[]
        if not review_occ:reasons.append('NO_NONSELF_TARGET_PEPTIDE')
        if not primary_score['seed_pass']:reasons.append('NO_PRIMARY_NATIVE_PASS_BINDER_LE500NM')
        if normal_flags and not primary_score['eligible']:reasons.append('DEPENDENT_NORMAL_DNA_EVIDENCE_REQUIRES_REVIEW')
        if hydro['status']!='PASS':reasons.append('HYDROPHOBIC_REGION')
        best=primary_score['best'] if primary_score['eligible'] else review_score['best']
        decisions.append(dict(window_id=w['window_id'],region_id=w['region_id'],event_id=w['event_id'],gene=w['gene'],gene_id=w['gene_id'],sequence=w['sequence'],coordinates=w['coordinates'],fragment_count=w['fragment_count'],recorded_alignment_paths=w['genomic_path_count'],supplied_gene_value=value,
          tier=tier,reasons=';'.join(reasons),primary_native_affinity_sum=primary_score['score'],primary_weighted_score=primary_score['weighted_score'],review_native_affinity_sum=review_score['score'],review_weighted_score=review_score['weighted_score'],best_target_peptide=best['peptide'] if best else '',best_target_allele=best['allele'] if best else '',best_target_affinity_nm=best['affinity_nm'] if best else '',
          primary_nonself_target_occurrences=len(primary_occ),all_nonself_target_occurrences=len(review_occ),normal_evidence_flags=dict(normal_flags),novelty_basis=dict(basis_counts),max_KD19=hydro['max_mean_kd'],native_old_score_verified=old,baseline_weighted_score=float(baseline['weighted_score_left'])))
        evidence.append(dict(window_id=w['window_id'],region_id=w['region_id'],gene=w['gene'],path_evidence=w['path_evidence']))
    def rank_key(r,review=False):return (-float(r['review_weighted_score'] if review else r['primary_weighted_score']),-float(r['review_native_affinity_sum'] if review else r['primary_native_affinity_sum']),r['event_id'],r['window_id'])
    def rankings(rows,review=False):
        events={};genes={}
        for r in sorted(rows,key=lambda r:rank_key(r,review)):events.setdefault(r['event_id'],r)
        for r in sorted(events.values(),key=lambda r:rank_key(r,review)):genes.setdefault(r['gene_id'],r)
        return ([dict(event_rank=i,**r) for i,r in enumerate(sorted(events.values(),key=lambda r:rank_key(r,review)),1)],
                [dict(gene_rank=i,**r) for i,r in enumerate(sorted(genes.values(),key=lambda r:rank_key(r,review)),1)])
    primary_events,primary_genes=rankings([r for r in decisions if r['tier']=='PRIMARY'])
    review_events,review_genes=rankings([r for r in decisions if r['tier']!='EXCLUDE'],True)
    for name,rows in [('windows.tsv.gz',decisions),('peptides.tsv.gz',peptide_rows),('event-ranking.tsv',primary_events),('gene-ranking.tsv',primary_genes),('review-gene-ranking.tsv',review_genes)]:table(out/name,rows)
    write_json_gz(out/'normal-sites.json.gz',site_cache);write_json_gz(out/'reference-matches.json.gz',matches);write_json_gz(out/'window-evidence.json.gz',evidence)
    summary=dict(status='RNA_POOL_SCREEN_COMPLETED',windows=len(windows),events=len({w['event_id'] for w in windows}),genes=len({w['gene_id'] for w in windows}),window_tiers=dict(Counter(r['tier'] for r in decisions)),primary_genes=len(primary_genes),primary_events=len(primary_events),unique_normal_sites=len(site_cache),reference_queries=len(wanted),original_affinity_peptides=len(actual_peptides),cached_prediction_pairs=len(actual_peptides)*len(alleles),new_model_inference=False,
      reference_ensembl_translations=len(ensembl_records),reference_uniprot_sequences=len(uni_records),all_native_old_scores_verified=True,all_recorded_alignment_paths_examined=True,normal_sites_classifications=dict(Counter(v['audit_normal_classification'] for v in site_cache.values())),normal_count_implementation_disagreements=sum(not v['independent_ALT_count_agreement'] for v in site_cache.values()),policy=policy,
      source_hashes=request['source_hashes'],normal_bam_sha256=plan['normal_bam_sha256'],hydrophobicity_source_sha256=sha(hydropathy.__file__),versions={n:importlib.metadata.version(n) for n in ['vaxrank','isovar','mhctools','pysam']},construct_changed=False)
    save(out/'summary.json',summary)
    paths=[p for p in out.iterdir() if p.is_file() and p.name!='status.json']
    save(out/'complete.json',dict(status='COMPLETED',completed_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),worker_sha256=sha(__file__),request_sha256=sha(a.request),outputs={p.name:dict(bytes=p.stat().st_size,sha256=sha(p)) for p in paths}))
    print('COMPLETED',json.dumps({k:v for k,v in summary.items() if k not in ['source_hashes','policy']}),flush=True)


if __name__=='__main__':main()
