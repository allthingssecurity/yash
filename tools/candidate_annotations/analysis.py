"""Prepare and replay annotations without changing recorded antigen selections."""
from pathlib import Path
import csv
import gzip
import hashlib
import io
import json
import math
import re
import sys

ROOT = Path(__file__).resolve().parents[2]
DEFAULT = ROOT / 'designs/design-3/review/candidate-annotations'
LENGTHS = (8, 9, 10, 11)


def require(condition, message):
    if not condition:
        raise ValueError(message)


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def load(path):
    with (gzip.open(path, 'rt') if str(path).endswith('.gz') else Path(path).open()) as f:
        return json.load(f)


def table(path):
    with (gzip.open(path, 'rt') if str(path).endswith('.gz') else Path(path).open()) as f:
        return list(csv.DictReader(f, delimiter='\t'))


def json_bytes(value):
    return (json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + '\n').encode()


def tsv_bytes(rows):
    stream = io.StringIO()
    fields = list(dict.fromkeys(k for r in rows for k in r))
    writer = csv.DictWriter(stream, fieldnames=fields, delimiter='\t', lineterminator='\n')
    writer.writeheader()
    for row in rows:
        writer.writerow({k: json.dumps(v, separators=(',', ':')) if isinstance(v, (list, dict)) else v
                         for k, v in row.items()})
    return stream.getvalue().encode()


def zipped(body):
    stream = io.BytesIO()
    with gzip.GzipFile(filename='', fileobj=stream, mode='wb', mtime=0) as f:
        f.write(body)
    return stream.getvalue()


def occurrences(sequence, lo, hi):
    return [(s, s + k, sequence[s:s+k]) for k in LENGTHS
            for s in range(len(sequence)-k+1)
            if (s < hi and s+k > lo if hi > lo else s < lo < s+k)]


def request():
    """Use the agreed 30-region plan and all 167 eligible splice windows."""
    sys.path.insert(0, str(ROOT / 'tools/rna_evidence'))
    from normal_peptides import checked_pool, require_selection
    selected_path = ROOT / 'data/shared-rna/splicing/junction-only/selected.json'
    require_selection(selected_path)
    decisions, _, checks = checked_pool()
    splices = [r for r in decisions.values() if r['primary_window_eligible']]
    selected_splices = load(selected_path)
    parent = load(ROOT / 'designs/design-3/inputs/source-regions.json')
    h2 = [r for r in load(ROOT / 'designs/design-3/inputs/h2bc12/source-regions.json')
          if r['gene'] == 'H2BC12' and r['source'] == 'MUTATION']
    require(len(h2) == 1, 'Expected the explicitly approved H2BC12 region')
    nonsplice = [r for r in parent if r['source'] != 'SPLICE'] + h2
    planned = nonsplice + selected_splices
    planned_ids = {r['region_id'] for r in planned}
    require(len(planned_ids) == len(planned) == 30, 'The current plan must contain 30 distinct regions')
    proteins = {r['variant_id']: r for r in table(ROOT / 'designs/design-2/inputs/rna/primary-proteins.tsv')}
    guardant = {r['variant_id']: r for r in table(ROOT / 'designs/design-3/inputs/contexts/contexts.tsv')}
    guardant_ref = load(ROOT / 'designs/design-3/inputs/predictions/reference-matches.json')
    guardant_normal = load(ROOT / 'designs/design-3/inputs/predictions/patient-normal-matches.json')
    normal_rna = load(ROOT / 'designs/design-2/inputs/integrated/rna-source-self-matches.json')
    inframe = load(ROOT / 'data/shared-rna/fusions/inframe-request.json')
    inframe_targets = set(inframe['peptides']) - set(inframe['reference_matching_peptides'])
    mutation_receipt = load(ROOT / 'designs/design-2/execution/rank/complete.json')
    require(mutation_receipt['exact_reference_matching_peptides'] == 0 and
            mutation_receipt['patient_normal_matching_peptides'] == 0,
            'Mutation annotation requires the recorded no-hit receipt')

    query_map, regions, records = {}, [], []

    def query(peptide, n='', c='', mode='without_flanks'):
        require(mode in {'without_flanks', 'with_flanks'}, 'Unknown query context mode')
        key = (peptide, n, c, mode)
        qid = 'Q_' + hashlib.sha256(json.dumps(key).encode()).hexdigest()[:24]
        value = dict(query_id=qid, peptide=peptide, n_flank=n, c_flank=c, mode=mode)
        require(qid not in query_map or query_map[qid] == value, 'Query-ID collision')
        query_map[qid] = value
        return qid

    for row in sorted(nonsplice + splices, key=lambda x: x['region_id']):
        source, seq, rid = row['source'], row['sequence'], row['region_id']
        annotation_id = row['window_id'] if source == 'SPLICE' else rid
        require(len(seq) == 25, f'Unsupported window length: {rid}')
        context, offset, basis = seq, 0, 'recorded_25aa_window_only'
        counterpart, agretopicity_status = None, 'NOT_A_SINGLE_RESIDUE_SUBSTITUTION'
        lo = int(row.get('target_start_in_window') or 0)
        hi = int(row.get('target_end_in_window') or 25)
        if source == 'SPLICE':
            lo, hi = int(row['target_start_in_window']), int(row['target_end_in_window'])
            targets = set(row['target_peptides'])
            self_scope = 'broad_reference_and_personal_normal_checked'
        elif source == 'MUTATION':
            ctx = proteins.get(row['variant_id']) or guardant.get(row['variant_id'])
            require(ctx is not None, f'Missing source context: {rid}')
            context = ctx['sequence']
            start_key = 'window_start_0based' if row['variant_id'] in proteins else 'window_start'
            offset = int(row[start_key])
            require(context[offset:offset+25] == seq, f'Wrong source-context placement: {rid}')
            focal_start = int(ctx.get('mutation_start_0based', ctx.get('mutation_start')))
            focal_end = int(ctx.get('mutation_end_0based_exclusive', ctx.get('mutation_end')))
            require((lo, hi) == (max(0, focal_start-offset), min(25, focal_end-offset)),
                    f'Mutation interval disagrees with reconstructed context: {rid}')
            basis = 'patient_RNA_Isovar_context' if row['variant_id'] in proteins else 'reference_derived_Guardant_context'
            targets = {p for _, _, p in occurrences(seq, lo, hi)}
            self_scope = 'recorded_mutation_reference_and_personal_normal_no_hit'
            if row['variant_id'] in guardant:
                require(targets <= guardant_ref.keys() & guardant_normal.keys(), 'Missing Guardant self decisions')
                targets = {p for p in targets if not guardant_ref[p]['occurs'] and
                           not guardant_normal[p]['occurs_in_searched_models']}
            change = row.get('protein_change') or row.get('DNA_annotation_protein_change')
            match = re.fullmatch(r'p\.([A-Z])(\d+)([A-Z])', change or '')
            if match:
                require(hi == lo+1 and seq[lo] == match[3], f'Wrong mutation position: {rid}')
                counterpart = seq[:lo] + match[1] + seq[lo+1:]
                agretopicity_status = 'FOCAL_REFERENCE_RESIDUE_RESTORED_OTHER_BACKGROUND_UNCHANGED'
        elif source == 'FUSION' and rid == 'REG_1f7babf47bbafeebeedc':
            context, offset, basis = inframe['source_sequence'], 558, 'recorded_inframe_fusion_source_protein'
            require(context[offset:offset+25] == seq, 'Wrong in-frame fusion placement')
            lo = inframe['junction_codon_0based'] - offset
            hi = lo+1
            targets = inframe_targets
            require(targets <= normal_rna.keys(), 'Missing in-frame personal-normal decisions')
            targets = {p for p in targets if not normal_rna[p]['occurs_in_searched_models']}
            self_scope = 'recorded_fusion_reference_and_personal_normal_checked'
        else:
            targets = set(filter(None, row.get('target_peptides', '').split(';')))
            self_scope = 'intentionally_reference_CTA_self' if source == 'CTA/SELF' else 'recorded_fusion_reference_and_personal_normal_checked'
        require(targets, f'No scored targets: {rid}')
        possible = occurrences(seq, lo, hi)
        require(targets <= {p for _, _, p in possible}, f'Target outside altered span: {rid}')
        region = dict(annotation_id=annotation_id, region_id=rid, window_id=row.get('window_id', ''), event_id=row['event_id'],
            gene=row['gene'], source=source, protein_change=row.get('protein_change') or row.get('DNA_annotation_protein_change', ''),
            sequence=seq, current_plan=rid in planned_ids, target_start=lo, target_end=hi,
            supplied_gene_value=row.get('supplied_gene_value', row.get('supplied_transformed_gene_tpm', '')),
            original_weighted_score=row.get('weighted_score_left', row.get('weighted_score', '')),
            original_best_affinity_nm=row.get('best_target_affinity_nm'), context_basis=basis,
            context_length=len(context), agretopicity_status=agretopicity_status,
            self_screen_scope=self_scope, target_peptides=sorted(targets))
        regions.append(region)
        for start, end, peptide in possible:
            if peptide not in targets:
                continue
            pos = offset+start
            n, c = context[max(0, pos-5):pos], context[pos+len(peptide):pos+len(peptide)+5]
            has_both = bool(n and c)
            record = dict(annotation_id=annotation_id, region_id=rid, start=start, end=end, peptide=peptide,
                n_flank=n, c_flank=c, full_5_5_flanks=len(n) == len(c) == 5,
                no_flank_query=query(peptide),
                source_context_query=query(peptide, n, c, 'with_flanks') if has_both else query(peptide),
                source_context_mode='with_flanks' if has_both else 'without_flanks')
            if counterpart is not None:
                refpep = counterpart[start:end]
                require(refpep != peptide, f'Counterfactual does not touch mutation: {rid}')
                record.update(reference_peptide=refpep, reference_query=query(refpep))
            records.append(record)
    require(len(regions) == 192 and len(splices) == 167, 'Unexpected eligible pool size; explicitly review changed inputs')
    require(len({r['annotation_id'] for r in regions}) == len(regions), 'Duplicate window annotation identities')
    paths = ['config/case.json', 'designs/design-3/inputs/source-regions.json',
        'designs/design-3/inputs/h2bc12/source-regions.json',
        'data/shared-rna/splicing/junction-only/selected.json',
        'data/shared-rna/splicing/junction-only/verification.json',
        'data/shared-rna/splicing/sequence-screen/complete.json',
        'data/shared-rna/splicing/independent-assessment/protein-screen/complete.json',
        'designs/design-2/inputs/rna/primary-proteins.tsv', 'designs/design-2/execution/rank/complete.json',
        'designs/design-3/inputs/contexts/contexts.tsv', 'designs/design-3/inputs/predictions/reference-matches.json',
        'designs/design-3/inputs/predictions/patient-normal-matches.json',
        'designs/design-2/inputs/integrated/rna-source-self-matches.json',
        'data/shared-rna/fusions/inframe-request.json',
        'designs/initial/inputs/checkpoints/provenance/construct-model-files.json']
    return dict(schema_version=2, alleles=load(ROOT/'config/case.json')['hla_class_I'],
        regions=regions, occurrences=records, queries=sorted(query_map.values(), key=lambda q:q['query_id']),
        source_sha256={p:sha(ROOT/p) for p in paths}, RNA_eligibility_check=checks,
        policy=dict(selection_changed=False, construct_rebuilt=False,
            splice_comparison='Best no-flank presentation among target peptide/allele pairs with fresh affinity <=500 nM; best window per gene. Existing weighted score, event_id and region_id break ties.',
            presentation_uses_expression=False, expression_used_only_as_existing_tie_break=True,
            flanked_scores='Sensitivity annotations with available recorded source context, never synthetic padding.',
            query_identity=['peptide','n_flank','c_flank','mode'],
            agretopicity='Mutant nM / focal-reference nM on the same allele; no cutoff.',
            counterpart_is_not_fully_phased_patient_normal=True,
            current_plan_not_assembled=True, vaccine_context_prediction_performed=False))


def products(req, predictions):
    values = {(p['query_id'], p['allele']): p for p in predictions}
    require(len(values) == len(predictions) == len(req['queries'])*len(req['alleles']), 'Missing/duplicate predictions')
    require(set(values) == {(q['query_id'], a) for q in req['queries'] for a in req['alleles']}, 'Unexpected inference queries')
    queries = {q['query_id']:q for q in req['queries']}
    for value in predictions:
        query = queries[value['query_id']]
        require(all(query[k] == value[k] for k in ['peptide','n_flank','c_flank','mode']), 'Wrong context attached to prediction')
        require(float(value['affinity_nm']) > 0 and all(math.isfinite(float(value[k])) for k in
            ['affinity_nm','processing_score','presentation_score']), 'Invalid model value')
        require(all(0 <= float(value[k]) <= 1 for k in ['processing_score','presentation_score']), 'Invalid model score')
    by_region = {r['annotation_id']: [] for r in req['regions']}
    require(len(by_region) == len(req['regions']), 'Duplicate annotation identities')
    peptide_rows = []
    for occurrence in req['occurrences']:
        for allele in req['alleles']:
            bare = values[occurrence['no_flank_query'], allele]
            flanked = values[occurrence['source_context_query'], allele]
            row = dict(annotation_id=occurrence['annotation_id'], region_id=occurrence['region_id'], peptide=occurrence['peptide'], allele=allele,
                start=occurrence['start'], end=occurrence['end'], n_flank=occurrence['n_flank'], c_flank=occurrence['c_flank'],
                full_5_5_flanks=occurrence['full_5_5_flanks'], source_context_mode=occurrence['source_context_mode'],
                affinity_nm=float(bare['affinity_nm']),
                processing_no_flanks=float(bare['processing_score']), presentation_no_flanks=float(bare['presentation_score']),
                processing_source_context=float(flanked['processing_score']), presentation_source_context=float(flanked['presentation_score']),
                reference_peptide=occurrence.get('reference_peptide',''), reference_affinity_nm='', mutant_over_reference_affinity='')
            if occurrence.get('reference_query'):
                ref = values[occurrence['reference_query'], allele]
                row['reference_affinity_nm'] = float(ref['affinity_nm'])
                row['mutant_over_reference_affinity'] = row['affinity_nm']/row['reference_affinity_nm']
            peptide_rows.append(row)
            by_region[row['annotation_id']].append(row)
    region_rows = []
    max_affinity_difference = 0.0
    for region in req['regions']:
        rows = by_region[region['annotation_id']]
        best_affinity = min(rows, key=lambda r:(r['affinity_nm'], r['peptide'], r['allele']))
        admitted = [r for r in rows if r['affinity_nm'] <= 500]
        # The explicit in-frame fusion exception remains visible rather than silently disappearing.
        best_bare = max(admitted or rows, key=lambda r:(r['presentation_no_flanks'], -r['affinity_nm']))
        best_source = max(admitted or rows, key=lambda r:(r['presentation_source_context'], -r['affinity_nm']))
        original = float(region['original_best_affinity_nm'])
        diff = abs(best_affinity['affinity_nm']/original-1)
        max_affinity_difference = max(max_affinity_difference, diff)
        require(diff < 1e-4, f'Fresh affinity differs from recorded target set: {region["region_id"]}: {diff}')
        row = {k:v for k,v in region.items() if k != 'target_peptides'}
        row.update(target_peptide_count=len(set(r['peptide'] for r in rows)),
            best_affinity_nm=best_affinity['affinity_nm'], best_affinity_peptide=best_affinity['peptide'],
            best_affinity_allele=best_affinity['allele'], affinity_gate_exception=not bool(admitted),
            best_presentation_no_flanks=best_bare['presentation_no_flanks'],
            best_presentation_no_flanks_peptide=best_bare['peptide'], best_presentation_no_flanks_allele=best_bare['allele'],
            best_presentation_source_context=best_source['presentation_source_context'],
            best_presentation_source_context_peptide=best_source['peptide'], best_presentation_source_context_allele=best_source['allele'],
            best_presentation_source_full_flanks=best_source['full_5_5_flanks'],
            reference_affinity_of_best_binder_nm=best_affinity['reference_affinity_nm'],
            mutant_over_reference_of_best_binder=best_affinity['mutant_over_reference_affinity'],
            lowest_mutant_over_reference_affinity=min((r['mutant_over_reference_affinity'] for r in admitted
                if r['mutant_over_reference_affinity'] != ''), default=''))
        # Exploratory replacement for the binding contribution, retaining the
        # existing splice expression multiplier and occurrence-count semantics.
        # Sum only target peptide/allele pairs meeting the 500-nM affinity gate.
        # This is our explicit comparison rule, not a native Vaxrank score.
        row.update(presentation_sum_no_flanks=sum(r['presentation_no_flanks'] for r in admitted),
            presentation_sum_source_context=sum(r['presentation_source_context'] for r in admitted),
            expression_weighted_presentation_no_flanks='', expression_weighted_presentation_source_context='')
        if region['source'] == 'SPLICE':
            expression = float(region['supplied_gene_value'])
            row.update(expression_weighted_presentation_no_flanks=row['presentation_sum_no_flanks']*math.sqrt(expression),
                expression_weighted_presentation_source_context=row['presentation_sum_source_context']*math.sqrt(expression))
        region_rows.append(row)

    splices = [r for r in region_rows if r['source'] == 'SPLICE']
    def original_key(r):
        return (-float(r['original_weighted_score']), r['event_id'], r['region_id'])
    def rank(column):
        ordered = sorted(splices, key=lambda r:(-r[column], *original_key(r))) if column else sorted(splices, key=original_key)
        genes = {}
        for r in ordered:
            genes.setdefault(r['gene'], r)
        return {r['annotation_id']:i+1 for i,r in enumerate(ordered)}, {g:(i+1,r) for i,(g,r) in enumerate(genes.items())}
    affinity_windows, affinity_genes = rank(None)
    bare_windows, bare_genes = rank('best_presentation_no_flanks')
    source_windows, source_genes = rank('best_presentation_source_context')
    weighted_windows, weighted_genes = rank('expression_weighted_presentation_no_flanks')
    weighted_source_windows, weighted_source_genes = rank('expression_weighted_presentation_source_context')
    comparisons = []
    for gene,(oldrank,old) in affinity_genes.items():
        nr,new = bare_genes[gene]
        sr,sensitivity = source_genes[gene]
        wr,weighted = weighted_genes[gene]
        wsr,weighted_source = weighted_source_genes[gene]
        comparisons.append(dict(gene=gene, original_gene_rank=oldrank,
            original_region_id=old['region_id'], original_window_id=old['window_id'], original_sequence=old['sequence'],
            original_weighted_score=old['original_weighted_score'],
            original_best_presentation_no_flanks=old['best_presentation_no_flanks'],
            presentation_gene_rank=nr, presentation_region_id=new['region_id'], presentation_window_id=new['window_id'],
            presentation_sequence=new['sequence'], presentation_event_id=new['event_id'],
            best_presentation_no_flanks=new['best_presentation_no_flanks'],
            presentation_best_peptide=new['best_presentation_no_flanks_peptide'],
            presentation_best_allele=new['best_presentation_no_flanks_allele'],
            presentation_original_weighted_score=new['original_weighted_score'],
            source_context_gene_rank=sr, source_context_region_id=sensitivity['region_id'],
            best_presentation_source_context=sensitivity['best_presentation_source_context'],
            source_context_best_has_full_flanks=sensitivity['best_presentation_source_full_flanks'],
            original_top5=oldrank <= 5, presentation_top5=nr <= 5, source_context_top5=sr <= 5,
            same_window=old['annotation_id'] == new['annotation_id'],
            weighted_presentation_gene_rank=wr, weighted_presentation_region_id=weighted['region_id'], weighted_presentation_window_id=weighted['window_id'],
            weighted_presentation_sequence=weighted['sequence'],
            expression_weighted_presentation_no_flanks=weighted['expression_weighted_presentation_no_flanks'],
            weighted_presentation_same_window=old['annotation_id'] == weighted['annotation_id'],
            weighted_presentation_top5=wr <= 5,
            weighted_source_context_gene_rank=wsr, weighted_source_context_region_id=weighted_source['region_id'],
            expression_weighted_presentation_source_context=weighted_source['expression_weighted_presentation_source_context']))
    for r in splices:
        r.update(original_window_rank=affinity_windows[r['annotation_id']],
            presentation_window_rank=bare_windows[r['annotation_id']], source_context_window_rank=source_windows[r['annotation_id']],
            weighted_presentation_window_rank=weighted_windows[r['annotation_id']],
            weighted_source_context_window_rank=weighted_source_windows[r['annotation_id']])
    current = [r for r in region_rows if r['current_plan']]
    require(len(current) == 30, 'Incomplete current-plan annotations')
    original_top = [g for g,(i,_) in affinity_genes.items() if i <= 5]
    new_top = [g for g,(i,_) in bare_genes.items() if i <= 5]
    source_top = [g for g,(i,_) in source_genes.items() if i <= 5]
    weighted_top = [g for g,(i,_) in weighted_genes.items() if i <= 5]
    require(set(original_top) == {r['gene'] for r in current if r['source'] == 'SPLICE'}, 'Original gene ranking drift')
    summary = dict(status='ANNOTATIONS_COMPLETE_SELECTION_UNCHANGED', current_plan_regions=len(current),
        regions_annotated=len(region_rows), splice_windows_compared=len(splices), splice_genes=len(affinity_genes),
        distinct_region_ids=len({r['region_id'] for r in region_rows}),annotation_identity='window_id for SPLICE, otherwise region_id',
        unique_queries=len(req['queries']), peptide_allele_inference_rows=len(predictions),
        peptide_occurrence_annotation_rows=len(peptide_rows),
        mutation_regions=sum(r['source'] == 'MUTATION' for r in current),
        mutation_regions_with_reference_comparison=sum(r['source'] == 'MUTATION' and
            r['reference_affinity_of_best_binder_nm'] != '' for r in current),
        original_splice_top5=original_top, presentation_no_flanks_top5=new_top,
        presentation_source_context_top5=source_top, shared_top5_genes=sorted(set(original_top)&set(new_top)),
        expression_weighted_presentation_top5=weighted_top,
        expression_weighted_presentation_source_context_top5=[g for g,(i,_) in weighted_source_genes.items() if i <= 5],
        weighted_presentation_keeps_all_five_gene_windows=all(
            affinity_genes[g][1]['annotation_id'] == weighted_genes[g][1]['annotation_id'] for g in weighted_top)
            and set(weighted_top) == set(original_top),
        exploratory_weighted_score='Sum no-flank presentation scores over eligible target peptide occurrences and HLA pairs with affinity <=500 nM, times sqrt(supplied gene value). One best window per gene. Not a native Vaxrank score or a validated selection rule.',
        max_relative_fresh_vs_recorded_best_affinity_difference=max_affinity_difference,
        selection_changed=False, construct_rebuilt=False, vaccine_context_prediction_performed=False)
    return {'peptide-annotations.tsv.gz':zipped(tsv_bytes(peptide_rows)),
        'region-annotations.tsv':tsv_bytes(region_rows), 'current-plan.tsv':tsv_bytes(current),
        'splice-window-comparison.tsv':tsv_bytes(sorted(splices,key=lambda r:r['presentation_window_rank'])),
        'splice-gene-comparison.tsv':tsv_bytes(sorted(comparisons,key=lambda r:r['weighted_presentation_gene_rank'])),
        'summary.json':json_bytes(summary)}, summary


def current_source_paths(record):
    """Resolve moved source catalogs without changing native inference records."""
    mapping = {
        'designs/design-3/selection/selected.json': 'designs/design-3/inputs/source-regions.json',
        'designs/design-3-rna-reviewed/selection/selected.json': 'designs/design-3/inputs/h2bc12/source-regions.json',
    }
    value = dict(record)
    value['source_sha256'] = {mapping.get(name, name): digest for name, digest in record['source_sha256'].items()}
    for name, digest in value['source_sha256'].items():
        require(sha(ROOT / name) == digest, 'Changed annotation source: ' + name)
    return value


def finish(out, *, verify=False):
    out = Path(out)
    req = load(out/'request.json')
    require(current_source_paths(req) == request(), 'Current inputs changed; use a new annotation directory')
    inference = load(out/'inference.json')
    require(inference['request_sha256'] == sha(out/'request.json'), 'Inference/request mismatch')
    require(inference['model_files_verified'] and inference['new_model_inference'], 'Unverified inference')
    require(inference['worker_sha256'] == sha(Path(__file__).with_name('predict.py')), 'Changed measured inference worker')
    require(inference['model_manifest_sha256'] == sha(ROOT/'designs/initial/inputs/checkpoints/provenance/construct-model-files.json'), 'Changed model manifest')
    require(inference['prediction_sha256'] == sha(out/'predictions.tsv.gz'), 'Changed measured predictions')
    outputs, summary = products(req, table(out/'predictions.tsv.gz'))
    if verify:
        receipt = load(out/'verification.json')
        require(receipt['source_sha256'] == req['source_sha256'], 'Changed source identities')
        for name,body in outputs.items():
            require((out/name).read_bytes() == body, 'Annotation replay differs: '+name)
        for name,digest in receipt['output_sha256'].items():
            require(sha(out/name) == digest, 'Changed output: '+name)
    else:
        for name,body in outputs.items():
            (out/name).write_bytes(body)
        receipt = dict(status='PASS', source_sha256=req['source_sha256'],
            output_sha256={name:sha(out/name) for name in [*outputs,'request.json','predictions.tsv.gz','inference.json']},
            analysis_sha256=sha(Path(__file__)), every_region_and_occurrence_replayed=True,
            fresh_affinity_parity_passed=True, frozen_designs_modified=False)
        (out/'verification.json').write_bytes(json_bytes(receipt))
    return summary
