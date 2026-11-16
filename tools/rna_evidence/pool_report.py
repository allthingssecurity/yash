"""Replay the collected screening decisions and produce a concise candidate table."""
from pathlib import Path
from collections import Counter,defaultdict
import csv,gzip,hashlib,json,math,os
AA=set('ACDEFGHIKLMNPQRSTVWY')


def load(path):
    with (gzip.open(path,'rt') if str(path).endswith('.gz') else open(path)) as f:return json.load(f)


def table_rows(path):
    with (gzip.open(path,'rt') if str(path).endswith('.gz') else open(path)) as f:return list(csv.DictReader(f,delimiter='\t'))


def contribution(nm):
    if nm>=5000:return 0.0
    score=(1+math.exp(-350/150))/(1+math.exp((nm-350)/150))
    return score if score>=1e-5 else 0.0


def replay_and_report(out):
    out=Path(out);request=load(out/'request.json.gz');summary=load(out/'summary.json')
    decisions=table_rows(out/'windows.tsv.gz');peptides=table_rows(out/'peptides.tsv.gz')
    evidence={r['window_id']:r for r in load(out/'window-evidence.json.gz')};matches=load(out/'reference-matches.json.gz');sites=load(out/'normal-sites.json.gz')
    requested={w['window_id']:w for w in request['windows']};by_window=defaultdict(list)
    assert len(decisions)==len(requested)==summary['windows']
    for key,row in sites.items():
        n=row['normal'];reads=n['read_evidence'];status=n['genotype_evidence']['status'];alt=reads['fragment_counts'].get('ALT',0)
        adequate=reads['usable_fragment_names']>=request['policy']['normal_min_fragments']
        expected='PRESENT_IN_NORMAL' if status=='NORMAL_TARGET_ALLELE_PRESENT' and adequate and alt>0 else 'NORMAL_REFERENCE_SUPPORTED' if status=='NORMAL_TARGET_ALLELE_NOT_CALLED' and adequate and alt<request['policy']['normal_ALT_below'] else 'UNRESOLVED_NORMAL_EVIDENCE'
        assert expected==n['classification'],key
        agreement=row['independent_counts']['normal']['fragment_bases'].get(row['allele'][3],0)==alt
        assert agreement==row['independent_ALT_count_agreement']
        assert row['audit_normal_classification']==(expected if agreement else 'UNRESOLVED_NORMAL_EVIDENCE')
    for r in peptides:
        w=requested[r['window_id']];p=r['peptide'];start,end=int(r['start']),int(r['end'])
        assert w['sequence'][start:end]==p
        assert start<int(w['target_end_in_window']) and end>int(w['target_start_in_window'])
        assert math.isclose(float(r['native_contribution']),sum(contribution(x) for x in request['affinities'][p].values()),rel_tol=1e-8,abs_tol=1e-8)
        if matches[p]['self_match']:expected='EXCLUDE_SELF'
        else:
            path_results=[]
            for path in evidence[w['window_id']]['path_evidence']:
                ref=path['reference_translation_on_recorded_path'][start:end]
                differences=[d for d in path['differences'] if d['nonsynonymous'] and start<=d['aa_position_1based']-1<end]
                if set(ref)<=AA and not matches[ref]['self_match']:path_results.append(True)
                else:
                    kinds=[d['normal_classification'] for d in differences]
                    path_results.append(bool(kinds) and all(k in {'NORMAL_REFERENCE_SUPPORTED','PRESENT_IN_NORMAL'} for k in kinds) and 'NORMAL_REFERENCE_SUPPORTED' in kinds)
                for d in differences:
                    if d['site_key'] is not None:assert d['normal_classification']==sites[d['site_key']]['audit_normal_classification']
            expected='PRIMARY' if path_results and all(path_results) else 'REVIEW'
        assert expected==r['tier'],(r['window_id'],p,expected,r['tier'])
        by_window[r['window_id']].append(r)
    for r in decisions:
        w=requested[r['window_id']];prs=by_window[r['window_id']]
        assert len(prs)==sum(1 for k in request['policy']['peptide_lengths'] for i in range(len(w['sequence'])-k+1) if i<int(w['target_end_in_window']) and i+k>int(w['target_start_in_window']))
        primary=[p for p in prs if p['tier']=='PRIMARY'];review=[p for p in prs if p['tier']!='EXCLUDE_SELF']
        for name,subset in [('primary',primary),('review',review)]:
            expected=sum(float(p['native_contribution']) for p in subset)
            assert math.isclose(float(r[name+'_native_affinity_sum']),expected,rel_tol=1e-8,abs_tol=1e-8)
            assert math.isclose(float(r[name+'_weighted_score']),expected*math.sqrt(float(r['supplied_gene_value'])),rel_tol=1e-8,abs_tol=1e-8)
        def admit(prs):return any(p['native_PASS_seed']=='True' and float(p['best_affinity_nm'])<=request['policy']['admission_affinity_nm'] for p in prs) and sum(float(p['native_contribution']) for p in prs)>0 and float(r['supplied_gene_value'])>0 and float(r['max_KD19'])<=1.6
        tier='PRIMARY' if admit(primary) else 'REVIEW_NORMAL_EVIDENCE' if admit(review) else 'EXCLUDE'
        assert tier==r['tier']
    # Independently recover the deterministic event/gene ranking.
    def key(r,review=False):return (-float(r['review_weighted_score' if review else 'primary_weighted_score']),-float(r['review_native_affinity_sum' if review else 'primary_native_affinity_sum']),r['event_id'],r['window_id'])
    for filename,tier_filter,review in [('gene-ranking.tsv',lambda r:r['tier']=='PRIMARY',False),('review-gene-ranking.tsv',lambda r:r['tier']!='EXCLUDE',True)]:
        events={};genes={}
        for r in sorted(filter(tier_filter,decisions),key=lambda r:key(r,review)):events.setdefault(r['event_id'],r)
        for r in sorted(events.values(),key=lambda r:key(r,review)):genes.setdefault(r['gene_id'],r)
        expected=sorted(genes.values(),key=lambda r:key(r,review));recorded=table_rows(out/filename)
        assert [r['window_id'] for r in expected]==[r['window_id'] for r in recorded]
    result=dict(status='RNA_POOL_OFFLINE_REPLAY_PASSED',windows=len(decisions),peptide_occurrences=len(peptides),unique_normal_sites=len(sites),window_tiers=dict(Counter(r['tier'] for r in decisions)),native_score_curve_checked=True,normal_classifications_recomputed=True,all_peptide_tiers_recomputed=True,window_scores_and_gene_ranking_recomputed=True,new_model_inference=False,new_BAM_reading=False,construct_changed=False)
    (out/'verification.json').write_text(json.dumps(result,indent=2)+'\n')
    primary=table_rows(out/'gene-ranking.tsv');review=table_rows(out/'review-gene-ranking.tsv')
    lines=['# RNA candidate screening before selection','',
      f"The stage screened **{summary['windows']} reconstructed windows, {summary['events']} events and {summary['genes']} genes**. These already passed the recorded RNA/frame/native-peptide checks, coordinate-normal ND and limited strand checks. It examined all recorded alignment paths, every window nucleotide difference and matched-normal evidence before reranking. The construct was not changed.",'',
      f"After the new checks: **{summary['primary_genes']} genes across {summary['primary_events']} events** have a primary-tier window. Window decisions: "+', '.join(f"{k}: {v}" for k,v in summary['window_tiers'].items())+'.','',
      '## Primary-tier candidates','',
      '| Rank | Gene | Selected 25-aa region | Best predicted target, nM | HLA | Full-window fragments | Supplied gene value | Weighted score |','|---:|---|---|---:|---|---:|---:|---:|']
    for r in primary:
        lines.append(f"| {r['gene_rank']} | {r['gene']} | `{r['sequence']}` | {float(r['best_target_affinity_nm']):.2f} | {r['best_target_allele']} | {r['fragment_count']} | {float(r['supplied_gene_value']):.3f} | {float(r['primary_weighted_score']):.2f} |")
    lines+=['','The table chooses the best eligible window per event, then per gene; no quota or automatic additions are applied. Fragment names support sequence and are not UMI-confirmed independent molecules. The score uses cached class-I affinity contributions × √supplied gene value, not raw-read abundance.','',
      '## What the checks mean','',
      '- Exact target-peptide matches in all Ensembl translations, UniProt canonical/additional isoforms or the bounded patient-normal coding/HLA models are excluded from scoring.',
      '- Every nucleotide difference is compared with normal DNA, including synonymous/background differences. If novelty relies on extra RNA-inferred amino-acid changes, the primary tier requires resolved normal evidence and at least one reference-supported change. A normal-carried background does not by itself invalidate a separate RNA change. Missing/conflicting or solely normal-carried changes remain in the review tier.',
      '- A same-alignment reference-bases counterfactual distinguishes junction-context novelty from novelty requiring extra RNA changes. It is an attribution check, not a proposed antigen. When that junction context is already nonself, no somatic-DNA requirement is invented for the RNA junction.',
      '- All recorded paths must agree on primary admission. These are the alignments already present in the BAM; the stage does not establish their unique correctness or realign repeats.',
      '- Admission retains a native-PASS peptide with predicted affinity ≤500 nM, positive supplied expression and the existing hydropathy rule. A model prediction is not measured translation, presentation, immunogenicity or tumour specificity.','',
      '## Evidence and reruns','',
      f"The run examined **{summary['unique_normal_sites']} unique nominated alleles** and **{summary['reference_queries']} reference/counterfactual peptide queries**. It reused {summary['cached_prediction_pairs']:,} checksum-pinned peptide/HLA predictions; it made no new HLA-model inference. The original native score is verified for every input window before adding exclusions.",'',
      '- [Primary gene ranking](gene-ranking.tsv) · [event ranking](event-ranking.tsv) · [all windows](windows.tsv.gz).',
      '- [Review-inclusive gene ranking](review-gene-ranking.tsv) · [per-peptide tiers and novelty attribution](peptides.tsv.gz).',
      '- [Normal DNA and fresh focal tumour/RNA counts](normal-sites.json.gz) · [all-path sequence differences](window-evidence.json.gz).',
      '- [Reference witnesses](reference-matches.json.gz) · [scope and summary](summary.json) · [completion receipt](complete.json) · [offline replay](verification.json).',
      '- [Frozen inputs and source hashes](request.json.gz) · [pipeline code and commands]('+os.path.relpath(Path(__file__).resolve().parent/'README.md',out)+').','']
    admitted=[]
    for r in decisions:
        if r['tier']!='PRIMARY':continue
        w=requested[r['window_id']]
        selected_targets=sorted({p['peptide'] for p in by_window[r['window_id']] if p['tier']=='PRIMARY' and float(p['native_contribution'])>0})
        admitted.append(dict(**r,source='SPLICE',coding_sequence=w['coding_sequence'],target_start_in_window=w['target_start_in_window'],target_end_in_window=w['target_end_in_window'],target_peptides=selected_targets,evidence_basis='RNA_SUPPORTED_SEQUENCE_HYPOTHESIS',unique_splice_origin_proven=False))
    (out/'admitted.json').write_text(json.dumps(admitted,indent=2)+'\n')
    exports=dict(status='VERIFIED_RNA_CANDIDATE_EXPORT',native_receipt_sha256=hashlib.sha256((out/'complete.json').read_bytes()).hexdigest(),request_sha256=hashlib.sha256((out/'request.json.gz').read_bytes()).hexdigest(),admitted_windows=len(admitted),
        files={name:dict(bytes=(out/name).stat().st_size,sha256=hashlib.sha256((out/name).read_bytes()).hexdigest()) for name in ['admitted.json','gene-ranking.tsv','event-ranking.tsv','review-gene-ranking.tsv','verification.json']},replay_source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    (out/'exports.json').write_text(json.dumps(exports,indent=2)+'\n')
    lines += ['## Use the checked candidates','',
        '[admitted.json](admitted.json) contains only primary-tier windows, with source IDs, sequences, target masks and nonself target peptides. [exports.json](exports.json) binds those exports to the measured screening inputs. Select additional regions from this pool and run the selection gate before assembly; review-tier windows require a separate documented decision.','']
    (out/'README.md').write_text('\n'.join(lines))
    return result


def require_selection(screen,selection):
    """Fail before assembly if a proposed RNA region has not passed this screen."""
    screen=Path(screen);selection=Path(selection)
    native=load(screen/'complete.json');assert native['status']=='COMPLETED'
    assert native['request_sha256']==hashlib.sha256((screen/'request.json.gz').read_bytes()).hexdigest()
    for name,r in native['outputs'].items():
        p=screen/name;assert p.stat().st_size==r['bytes'] and hashlib.sha256(p.read_bytes()).hexdigest()==r['sha256'],name
    exported=load(screen/'exports.json')
    assert exported['status']=='VERIFIED_RNA_CANDIDATE_EXPORT'
    assert exported['replay_source_sha256']==hashlib.sha256(Path(__file__).read_bytes()).hexdigest(), 'Replay the screen with the maintained report code before selection'
    assert exported['native_receipt_sha256']==hashlib.sha256((screen/'complete.json').read_bytes()).hexdigest()
    assert exported['request_sha256']==hashlib.sha256((screen/'request.json.gz').read_bytes()).hexdigest()
    for name,r in exported['files'].items():
        p=screen/name;assert p.stat().st_size==r['bytes'] and hashlib.sha256(p.read_bytes()).hexdigest()==r['sha256'],name
    selected=load(selection) if selection.suffix=='.json' else table_rows(selection)
    candidates=[r for r in selected if r.get('source')=='SPLICE']
    assert candidates,'No SPLICE-source regions were supplied; this gate does not assess other source types'
    admitted={r['window_id']:r for r in load(screen/'admitted.json')}
    failures=[]
    for r in candidates:
        w=admitted.get(r.get('window_id'))
        mismatch=not w or w['region_id']!=r['region_id'] or w['sequence']!=r['sequence']
        if w:
            for key in ['coding_sequence','target_start_in_window','target_end_in_window']:
                if key in r and str(r[key])!=str(w[key]):mismatch=True
        if mismatch:failures.append(dict(gene=r.get('gene'),window_id=r.get('window_id'),reason='NOT_IN_PRIMARY_CHECKED_POOL_OR_SEQUENCE_MISMATCH'))
    if failures:raise ValueError(json.dumps(dict(status='RNA_SELECTION_REQUIRES_REVIEW',failures=failures)))
    return dict(status='RNA_SELECTION_GATE_PASSED',RNA_regions=len(candidates),other_source_regions_not_assessed=len(selected)-len(candidates))
