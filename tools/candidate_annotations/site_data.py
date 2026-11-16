"""Export verified scientific data for the main website; no inference or site build."""
from pathlib import Path
import argparse
import json
import re

from analysis import ROOT, load, table, sha, require
from class_ii import DEFAULT, finish
from reference_comparisons import verify as verify_reference_comparisons

REVIEW = ROOT / 'designs/design-3/review/annotations.md'
LABELS = {'MUTATION': 'Mutation', 'SPLICE': 'Splice', 'CTA/SELF': 'CTA/self', 'FUSION': 'Fusion'}


def number(value):
    return None if value in ('', None) else float(value)


def prediction(row):
    return dict(peptide=row['peptide'], allele=row['target_id'],
                affinity=number(row['affinity_nm']), presentation=number(row['presentation_score']),
                reference_match=row['exact_reference_match'] == 'True',
                personal_match=row['scoped_personal_normal_match'] == 'True',
                self_match=row['exact_self_search_hit'] == 'True',
                reference_affinity=number(row['reference_affinity_nm']),
                reference_presentation=number(row['reference_presentation_score']))


def final_construct(current):
    directory = ROOT / 'designs/design-3'
    replay = directory / 'execution/offline-replay.json'
    if not replay.exists():
        return None, []
    files = [replay, directory/'selection/selected.json',
        directory/'construct/protein/construct.json', directory/'construct/protein/antigen-order.tsv',
        directory/'construct/protein/cassette.fasta', directory/'construct/final-check.json',
        directory/'construct/mrna/verification.json']
    checked, selected, construct = load(replay), load(files[1]), load(files[2])
    require(checked['status'] == 'DESIGN3_FINAL_DESIGN_OFFLINE_REPLAY_PASSED', 'Final build has no successful replay')
    require({(r['region_id'], r['sequence']) for r in selected} ==
            {(r['region_id'], r['sequence']) for r in current}, 'Final construct differs from the annotated plan')
    require(load(files[5])['protein_sha256'] == sha(files[4]), 'Final protein identity changed')
    require(load(files[6])['all_constraints_pass'], 'Final encoding checks failed')
    by_id = {r['region_id']: r for r in selected}
    order = table(files[3])
    require(len(order) == 30, 'Incomplete final construct order')
    return dict(regions=30, protein_length=construct['protein_length_aa'],
        linker_counts=construct['linker_counts'], totals=construct['totals'],
        antigens=[dict(position=int(r['position']), gene=r['gene'], sequence=r['sequence'],
            change=by_id[r['region_id']].get('protein_change', ''),
            source=by_id[r['region_id']]['source']) for r in order]), files


def build(destination):
    finish(DEFAULT, verify=True)
    req = load(DEFAULT/'request.json')
    reference_comparisons = verify_reference_comparisons()['entries']
    primary = set(req['policy']['primary_targets'])
    class_i = {r['annotation_id']: r for r in table(DEFAULT.parent/'candidate-annotations/region-annotations.tsv')}
    by_allele = {}
    for row in table(DEFAULT.parent/'candidate-annotations/peptide-annotations.tsv.gz'):
        key = (row['annotation_id'], row['allele'])
        candidate = dict(allele=row['allele'], peptide=row['peptide'],
            affinity_nm=float(row['affinity_nm']), start=int(row['start']))
        previous = by_allele.get(key)
        if previous is None or (candidate['affinity_nm'], candidate['peptide'], candidate['start']) < (previous['affinity_nm'], previous['peptide'], previous['start']):
            by_allele[key] = candidate
    # Retain separate BA/EL winners. They can be different peptides and alleles.
    by_target = {}
    for r in table(DEFAULT/'peptide-annotations.tsv.gz'):
        key = r['annotation_id'], r['target_id']
        p = prediction(r)
        group = by_target.setdefault(key, {})
        for label, eligible in [('all', True), ('no_hit', not p['self_match'])]:
            if not eligible:
                continue
            ba, el = label+'_binding', label+'_presentation'
            if ba not in group or (p['affinity'], p['peptide']) < (group[ba]['affinity'], group[ba]['peptide']):
                group[ba] = p
            if el not in group or (p['presentation'], -p['affinity'], p['peptide']) > (
                    group[el]['presentation'], -group[el]['affinity'], group[el]['peptide']):
                group[el] = p
    expression = {}
    with (ROOT/'data/expression/vendor-transformed-expression.tsv').open() as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            gene, value = line.rstrip().split('\t')
            require(gene not in expression, 'Duplicate supplied gene value')
            expression[gene] = float(value)
    manual = {r['selection']['window_id']: r['original_alignment_evidence'] for r in
              load(ROOT/'data/shared-rna/splicing/junction-only/manual-evidence.json')['entries']}
    guardant = {r['gene']: r for r in table(ROOT/'designs/design-3/inputs/guardant/summary.tsv')}
    rna = {r['variant_id']:r for r in table(ROOT/'designs/design-2/inputs/rna/primary-proteins.tsv')}
    recorded_order = table(ROOT/'designs/design-3/construct/protein/antigen-order.tsv')
    recorded = {r['region_id']: int(r['position']) for r in recorded_order}
    recorded_selection = {r['region_id']:r for r in load(ROOT/'designs/design-3/selection/selected.json')}
    protein_parts = {r['region_id']:r for r in table(ROOT/'designs/design-3/construct/protein/components.tsv') if r['kind']=='antigen'}
    class_i_alleles = [r['allele'] for r in table(ROOT/'data/hla/class-I.tsv')]
    rows = []
    for r in req['regions']:
        i = class_i[r['annotation_id']]
        region = dict(id=r['annotation_id'], region_id=r['region_id'], window_id=r['window_id'],
            gene=r['gene'], source=r['source'], source_label=LABELS[r['source']],
            change=r['protein_change'], sequence=r['sequence'], target_start=r['target_start'],
            target_end=r['target_end'], planned=r['current_plan'], expression=expression.get(r['gene']),
            recorded_position=recorded.get(r['region_id']), evidence=[], notes=[], class_ii=[])
        if r['source'] == 'MUTATION':
            if r['gene'] in guardant and r['context_basis'].startswith('reference'):
                g = guardant[r['gene']]
                region['basis'] = 'Guardant · reference-derived window'
                region['evidence'] = [
                    'Explicitly included from the October Guardant plasma report; the full window is not reconstructed from tumour RNA.',
                    f"Matched-normal DNA: {g['normal_ALT']} alternate / {int(g['normal_REF'])+int(g['normal_ALT'])+int(g['normal_OTHER'])} qualifying fragments.",
                    f"January tumour RNA at the nominated DNA site: {g['tumour_RNA_ALT']} alternate / {int(g['tumour_RNA_REF'])+int(g['tumour_RNA_ALT'])} REF+ALT fragments. This is a different specimen/time point."]
                region['notes'].append('Reference-derived window')
                if r['gene'] == 'AR':
                    region['notes'].append('Normal coverage exception: 0/19 ALT')
            else:
                region['basis'] = 'Paired DNA + tumour RNA'
                region['evidence'] = ['Included from matched-normal-supported DNA calls and an altered protein context reconstructed from tumour RNA.']
                context = rna.get(r['event_id'])
                if context:
                    region['evidence'].append(f"Isovar altered-protein context support: {context['supporting_fragments']} RNA fragments. This is sequence support, not an expression estimate.")
                if r['gene'] == 'H2BC12':
                    region['notes'].append('Explicit inclusion; supplied expression missing')
        elif r['source'] == 'SPLICE':
            region['basis'] = 'RNA · reference-join translation'
            region['evidence'] = ['All 25 amino acids match reference bases translated across the nominated join on every recorded supporting path; no additional amino-acid substitution is required.',
                'Scored class-I target peptides pass the broad reference/personal-self search. Nominated coordinates are not detected in the searched normal-junction catalogs.']
            if r['window_id'] in manual:
                m = manual[r['window_id']]
                region['evidence'].append(f"Complete window: {m['full_window_fragment_names']} supporting fragment names; {m['expected_fragment_names']} expected / {m['opposing_fragment_names']} opposing orientation. Native review: {m['native_review_tier']}/{m['native_review_score']}.")
                region['notes'].extend(m['native_review_flags'])
            else:
                region['notes'].append('Other eligible splice window; not in the current plan')
        elif r['source'] == 'CTA/SELF':
            region['basis'] = 'Expressed gene · reference protein'
            region['evidence'] = ['Reference-protein CTA/self region. Normal-reference matches are expected and retained for this source.']
        else:
            region['basis'] = 'RNA fusion-derived region'
            region['evidence'] = ['One of the two recorded CDH1::NRG1 fusion-derived sequence regions. The windows are separate sequence hypotheses, not two independent fusion discoveries.']
            region['evidence'].append(f"Supplied gene values: CDH1 {expression['CDH1']:.5f}; NRG1 {expression['NRG1']:.5f}. A fusion-specific TPM is not supplied.")
            if i['affinity_gate_exception'] == 'True':
                region['notes'].append('Recorded class-I affinity exception: 508.95 nM')
        region['class_i'] = dict(binding=dict(value=number(i['best_affinity_nm']),peptide=i['best_affinity_peptide'],allele=i['best_affinity_allele']),
            presentation=dict(value=number(i['best_presentation_no_flanks']),peptide=i['best_presentation_no_flanks_peptide'],allele=i['best_presentation_no_flanks_allele']),
            source_context=dict(value=number(i['best_presentation_source_context']),peptide=i['best_presentation_source_context_peptide'],allele=i['best_presentation_source_context_allele'],full_flanks=i['best_presentation_source_full_flanks']=='True'),
            reference_affinity=number(i['reference_affinity_of_best_binder_nm']),
            mutant_over_reference=number(i['mutant_over_reference_of_best_binder']),
            splice_binding_rank=number(i.get('original_window_rank')),
            splice_presentation_rank=number(i.get('weighted_presentation_window_rank')))
        region['class_i']['by_allele'] = [by_allele[(r['annotation_id'], allele)] for allele in class_i_alleles]
        require(abs(min(p['affinity_nm'] for p in region['class_i']['by_allele']) - region['class_i']['binding']['value']) < 1e-7, 'Per-allele winner disagrees with region affinity')
        part = protein_parts.get(r['region_id'])
        region['protein_start'] = int(part['start_aa'])+1 if part else None
        region['protein_end'] = int(part['end_aa']) if part else None
        for target in req['targets']:
            values = by_target[r['annotation_id'],target['id']]
            region['class_ii'].append(dict(id=target['id'],primary=target['id'] in primary,**values))
        for name, primary_only in [('DR', True), ('DQ_DP', False)]:
            targets = [t for t in region['class_ii'] if t['primary'] == primary_only]
            region[name] = dict(binding=min((t['all_binding'] for t in targets), key=lambda p:(p['affinity'],p['allele'])),
                presentation=max((t['all_presentation'] for t in targets), key=lambda p:(p['presentation'],-p['affinity'])))
        if region['source'] != 'CTA/SELF' and region['DR']['binding']['self_match']:
            region['notes'].append('Best DR binding peptide has an exact searched self match')
        region['reference_comparison'] = reference_comparisons.get(region['region_id'])
        if region['reference_comparison'] is None and region['source'] == 'MUTATION':
            change = re.fullmatch(r'p\.([A-Z])\d+([A-Z])', region['change'])
            if change:
                lo, hi = region['target_start'], region['target_end']
                require(hi == lo+1 and region['sequence'][lo] == change[2], 'Substitution comparison position mismatch')
                reference = region['sequence'][:lo] + change[1] + region['sequence'][hi:]
                region['reference_comparison'] = dict(source='MUTATION', candidate=region['sequence'],
                    variants=[dict(reference=reference, effect='substitution', transcripts=[])])
        region['display_gene'] = region['gene'].replace('--', '::')
        region['display_change'] = region['change']
        if region['source'] == 'FUSION':
            original = recorded_selection.get(region['region_id'], {})
            region['display_change'] = 'Alternative fusion tail' if original.get('target_mode') == 'frameshift_tail' else 'In-frame junction'
        rows.append(region)
    order = {'MUTATION':0, 'SPLICE':1, 'CTA/SELF':2, 'FUSION':3}
    rows.sort(key=lambda r:(order[r['source']],r['gene'],r['change'],not r['planned'],r['id']))
    current = [r for r in rows if r['planned']]
    require(len(current)==30 and len(rows)==192, 'Incomplete review view')
    final, final_files = final_construct(current)
    sources = [DEFAULT/'verification.json',DEFAULT/'current-plan.tsv',DEFAULT/'region-target-annotations.tsv',
        ROOT/'data/hla/class-II.json',ROOT/'data/hla/class-I.tsv',ROOT/'data/expression/vendor-transformed-expression.tsv',
        ROOT/'data/shared-rna/splicing/junction-only/manual-evidence.json',
        ROOT/'designs/design-3/inputs/guardant/summary.tsv',ROOT/'designs/design-3/construct/protein/antigen-order.tsv',
        ROOT/'designs/design-3/selection/selected.json',ROOT/'designs/design-2/inputs/rna/primary-proteins.tsv'] + final_files
    payload = dict(title='Design #3 · antigen review',
        status='30-region construct assembled and replayed' if final else '30-region working plan; not assembled',
        final=final,
        current_plan_count=30,counts={k:sum(r['source']==k for r in current) for k in LABELS},regions=rows,
        class_i_alleles=[r['allele'] for r in table(ROOT/'data/hla/class-I.tsv')],
        typing=load(ROOT/'data/hla/class-II.json'),
        splice_comparison=table(DEFAULT.parent/'candidate-annotations/splice-gene-comparison.tsv'),
        source_sha256={str(p.relative_to(ROOT)):sha(p) for p in sources})
    md = ['# Design #3 — current antigen annotations','',
        ('**30 selected regions:** ' if final else '**30-region working plan:** ') +
        '18 mutation, 5 splice, 5 CTA/self and 2 fusion regions. ' +
        (f'The [assembled construct](../README.md) has {final["protein_length"]} amino acids; assembly, junction checks and RNA encoding passed offline replay. '
         if final else 'This plan has not been assembled. ') +
        'Design #3 is the lead candidate.','',
        '[Explore on the main site](https://rnafoundry.com/yash/#vaccine) · [Combined source table](class-II-annotations/current-plan.tsv) · '
        '[Class-I methods](../../../tools/candidate_annotations/README.md) · [Class-II methods](class-II-annotations/README.md)','',
        'Binding is predicted affinity in nM (lower is stronger). Presentation is a separate model score '
        '(higher is more ligand-like; not a probability of a response). Each column is the best within '
        'the region and can refer to a different peptide or allele. No new cutoff, combined score or selection change is applied.','',
        'DR summaries use DRB1*03:01/*08:02. Eight DQ/DP pairing hypotheses are kept separately in the interactive view; DRB3 is excluded. '
        'Both observed DRA peptide-binding domains are identical, so four possible DRA/DRB1 combinations map to two model targets.','']
    fmt = lambda v: '—' if v is None else f'{v:,.2f}'
    for source, label in LABELS.items():
        md.extend(['## '+label,'', '| Antigen | Supplied TPM | Class-I binding, nM | Class-I presentation | DR binding, nM | DR presentation | Notes |',
                   '|---|---:|---:|---:|---:|---:|---|'])
        for r in current:
            if r['source'] != source:
                continue
            note = '; '.join(r['notes']) or r['basis']
            md.append('| '+' | '.join([r['gene']+' '+r['change'],fmt(r['expression']),fmt(r['class_i']['binding']['value']),
                f"{r['class_i']['presentation']['value']:.3f}",fmt(r['DR']['binding']['affinity']),
                f"{r['DR']['presentation']['presentation']:.3f}",note])+' |')
        md.append('')
    md.extend(['## How to read these annotations','',
        '- Supplied TPM is the vendor-transformed gene value, shown as an annotation; it is not an input to either prediction model.',
        '- Class I uses MHCflurry 2.3.14. CapHLA separately scores 15–20-aa peptides for class II.',
        '- Reference comparisons restore only the nominated missense residue. Other background residues stay unchanged; these are not fully phased normal proteins.',
        '- Long-peptide exact matches are reported explicitly. CTA/self matches are intentional; for altered sources, a match prevents treating that particular peptide as unique.',
        '- Whole-peptide overlap with an alteration does not establish that the class-II binding core contains it. CapHLA does not assign that core here.',
        '- Personal-normal sequence models were built for ≤11-aa targets; their negative results are incomplete for class II. Reference searches do not measure presentation in healthy cells.','',
        'For allele-specific winners, counterpart predictions and the remaining splice pool, '
        'use the [interactive view on the main site](https://rnafoundry.com/yash/#vaccine).',''])
    payload.pop('source_sha256', None)
    other = [r for r in payload['regions'] if not r['planned']]
    payload['regions'] = [r for r in payload['regions'] if r['planned']]
    other_path = destination.with_name(destination.stem + '-other.json')
    payload['other_windows_url'] = other_path.name
    payload['other_window_count'] = len(other)
    compact = lambda data: (json.dumps(data, ensure_ascii=False, separators=(',', ':'))+'\n').encode()
    return {destination:compact(payload),other_path:compact(other),REVIEW:('\n'.join(md)).encode()}


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True,help='Website JSON asset, outside this repository');p.add_argument('--verify',action='store_true');a=p.parse_args()
    for path,body in build(a.output.resolve()).items():
        if a.verify:require(path.read_bytes()==body,'Review is stale: '+str(path))
        else:path.parent.mkdir(parents=True,exist_ok=True);path.write_bytes(body)
    print('Design #3 site data '+('verified' if a.verify else 'generated')+': 30 selected regions, 192 annotated windows, 10 class-II model targets; selection unchanged.')


if __name__=='__main__':main()
