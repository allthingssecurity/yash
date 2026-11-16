"""Reference-derived Guardant contexts, fresh HLA-I predictions and explicit additions.

The five plasma calls are user-nominated. This stage does not label their
sequences RNA-reconstructed or pretend they pass Design #2's DNA/RNA policy.
"""
from pathlib import Path
from collections import defaultdict
from dataclasses import asdict
import csv
import hashlib
import importlib.metadata
import json
import math
import os
import sys

os.environ['OMP_NUM_THREADS'] = '4'
os.environ['MKL_NUM_THREADS'] = '4'
PARENT = Path('/data/results/TN25-112103/matched-normal-20261107')
ROOT = Path(sys.argv[1])
MODELS = Path('/data/results/TN25-112103/openvax-mutations-20261018/references/mhcflurry/2.3.0/models_class1_presentation/models')
sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import sha, save, receipt


def table(path, rows):
    fields = list(dict.fromkeys(k for row in rows for k in row))
    with Path(path).open('w') as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, delimiter='\t', lineterminator='\n')
        writer.writeheader()
        for row in rows:
            writer.writerow({k: json.dumps(v) if isinstance(v, (dict, list, tuple)) else v for k, v in row.items()})


def prepare():
    from pyensembl import Genome
    genome = Genome(**json.loads((PARENT / 'analysis-genome-config.json').read_text()))
    report = json.loads((ROOT / 'inputs/guardant/results.json').read_text())
    config = json.loads((ROOT / 'policy.json').read_text())
    contexts = []
    for call in report:
        if call['scope'] != 'unmatched':
            continue
        assert call['candidate_count'] == 1
        mapping = call['mappings'][0]
        change = call['reported_protein_change']
        position, reference, alternate = int(change[1:-1]), change[0], change[-1]
        by_sequence = {}
        for model in mapping['compatible_transcripts']:
            transcript = genome.transcript_by_id(model['transcript_id'])
            protein = transcript.protein_sequence
            assert protein[position - 1] == reference
            lo, hi = max(0, position - 25), min(len(protein), position + 24)
            ref_context = protein[lo:hi]
            mutation = position - 1 - lo
            mutant = ref_context[:mutation] + alternate + ref_context[mutation + 1:]
            assert len(mutant) >= 25 and set(mutant) <= set('ACDEFGHIKLMNPQRSTVWY')
            key = mutant, mutation
            if key in by_sequence:
                by_sequence[key]['transcript_ids'].append(transcript.id)
                continue
            variant_id = 'guardant_' + '_'.join(map(str, mapping['allele']))
            by_sequence[key] = dict(context_id=variant_id + '_' + hashlib.sha256(mutant.encode()).hexdigest()[:12],
                variant_id=variant_id, gene=call['gene'], protein_change='p.' + change,
                normalized_allele=mapping['allele'], transcript_ids=[transcript.id],
                reference_sequence=ref_context, sequence=mutant,
                mutation_start=mutation, mutation_end=mutation + 1,
                reference_protein_start_0based=lo, guardant_cfDNA_percent=call['October_cfDNA_percent'],
                normal=mapping['normal'], tumour_DNA=mapping['tumour_DNA'], tumour_RNA=mapping['tumour_RNA'],
                sequence_basis='Ensembl-103 reference protein plus the inferred reported substitution',
                patient_RNA_sequence_confirmed=False, full_window_RNA_reconstruction_performed=False,
                Guardant_genomic_coordinates_supplied=False)
        contexts.extend(by_sequence.values())
    assert len({c['gene'] for c in contexts}) == config['additional_guardant_variants'] == 5
    out = ROOT / 'contexts'
    out.mkdir(exist_ok=True)
    save(out / 'contexts.json', contexts)
    table(out / 'contexts.tsv', contexts)
    (out / 'contexts.fasta').write_text(''.join('>' + r['context_id'] + ' ' + r['gene'] + '\n' + r['sequence'] + '\n' for r in contexts))
    for name in ['case.json', 'analysis-genome-config.json', 'affinity-model-manifest.json']:
        (ROOT / name).write_bytes((PARENT / name).read_bytes())
    # Reuse the already measured, explicitly scoped normal sequence library.
    if not (ROOT / 'normal-self').exists():
        (ROOT / 'normal-self').symlink_to(PARENT / 'normal-self', target_is_directory=True)
    receipt(out, [out / name for name in ['contexts.json', 'contexts.tsv', 'contexts.fasta']],
            variants=5, distinct_contexts=len(contexts), worker_sha256=sha(__file__),
            patient_RNA_reconstruction_claimed=False, parent_inputs_unchanged=True)


def predict_select():
    import numpy as np
    import torch
    from pyensembl import Genome
    from mhcflurry import Class1AffinityPredictor
    from mhctools.mhcflurry import MHCflurry_Affinity
    from mhctools.pred import Prediction, Kind
    from vaxrank.reference_proteome import self_reference_matches
    from vaxrank.candidate_epitope import CandidateEpitope
    from vaxrank.epitope_config import EpitopeConfig
    from vaxrank.epitope_dsl import attach_per_allele_scores, epitopes_for_ranking
    from vaxrank.epitope_logic import slice_epitopes
    from vaxrank.vaccine_antigen import VaccineAntigen, TargetableMask, AminoAcidInterval, TumorSpecificityAttestation
    from vaxrank.vaccine_peptide import VaccinePeptide
    from self_search import patient_self_matches
    sys.path.insert(0, str(ROOT / 'frozen/src'))
    from yash_vaccine.hydrophobicity import summarize
    torch.set_num_threads(4)
    torch.set_float32_matmul_precision('highest')
    torch.backends.cuda.matmul.allow_tf32 = False
    torch.backends.cudnn.allow_tf32 = False
    torch.manual_seed(0)
    np.random.seed(0)
    contexts = json.loads((ROOT / 'contexts/contexts.json').read_text())
    alleles = json.loads((ROOT / 'case.json').read_text())['hla_class_I']
    genome = Genome(**json.loads((ROOT / 'analysis-genome-config.json').read_text()))
    expected = json.loads((ROOT / 'affinity-model-manifest.json').read_text())
    checked = []
    for row in expected:
        if row['path'].startswith('affinity_predictor/'):
            path = MODELS / row['path']
            assert path.stat().st_size == row['bytes'] and sha(path) == row['sha256']
            checked.append(row)
    out = ROOT / 'rank'
    out.mkdir(exist_ok=True)
    assert not (out / 'complete.json').exists(), 'Completed inference must not be overwritten'
    occurrences = []
    for context in contexts:
        for length in [8, 9, 10, 11]:
            for start in range(len(context['sequence']) - length + 1):
                if start < context['mutation_end'] and start + length > context['mutation_start']:
                    occurrences.append(dict(context_id=context['context_id'], peptide=context['sequence'][start:start + length],
                                            start=start, end=start + length))
    peptides = sorted({r['peptide'] for r in occurrences})
    affinity, predictions = {}, []
    wrapper = MHCflurry_Affinity(alleles=alleles, default_peptide_lengths=[8, 9, 10, 11], models_path=str(MODELS / 'affinity_predictor'))
    for peptide, result in zip(peptides, wrapper.predict(peptides)):
        for pred in result.preds:
            nm, rank = float(pred.value), float(pred.percentile_rank)
            assert pred.kind == Kind.pMHC_affinity and pred.peptide == peptide
            assert math.isfinite(nm) and nm > 0 and math.isfinite(rank) and 0 <= rank <= 100
            key = peptide, pred.allele
            assert key not in affinity and pred.allele in alleles
            affinity[key] = nm, rank
            predictions.append(dict(peptide=peptide, allele=pred.allele, affinity_nm=nm, affinity_percentile=rank))
    assert set(affinity) == {(p, a) for p in peptides for a in alleles}
    print('NATIVE_MHCTOOLS_AFFINITY', len(affinity), 'pairs', flush=True)
    def antigen(context):
        return VaccineAntigen(kind='mutation', amino_acids=context['sequence'],
            targetable_mask=TargetableMask((AminoAcidInterval(context['mutation_start'], context['mutation_end']),)),
            tumor_specificity=TumorSpecificityAttestation(status='overridden', evidence_kind='user_selected_Guardant_VUS',
                evidence_source='Plasma report and matched-normal read/genotype review; reference-derived protein context',
                requires_review=True, override_reason='Explicit user-selected exploratory addition; not RNA-reconstructed and not established tumour specificity'))
    references = self_reference_matches(peptides, antigen(contexts[0]), genome)
    assert set(references) == set(peptides) and all(m.source_provenance_complete for m in references.values())
    normal_matches = patient_self_matches(ROOT, peptides)
    broad = json.loads((ROOT / 'inputs/sequence-review/reference-matches.json').read_text())
    assert set(peptides) <= set(broad), 'Broad-reference audit is missing prediction queries'
    assert all(broad[p]['extended'] is None and not broad[p]['personal']['occurs_in_searched_models'] for p in peptides), 'Revisit Guardant window scoring for new normal-self matches'
    save(out / 'reference-matches.json', {p: asdict(m) for p, m in references.items()})
    save(out / 'patient-normal-matches.json', normal_matches)
    table(out / 'peptide-HLA-affinity.tsv', predictions)
    table(out / 'peptide-occurrences.tsv', occurrences)
    expression = {}
    for line in Path('/data/inputs/TN25-112103/vendor-transformed-expression.tsv').read_text().splitlines():
        if line.strip() and not line.startswith('#'):
            gene, value = line.split('\t')
            expression[gene] = float(value)
    by_context = defaultdict(list)
    for row in occurrences:
        if not references[row['peptide']].occurs and not normal_matches[row['peptide']]['occurs_in_searched_models']:
            by_context[row['context_id']].append(row)
    config = EpitopeConfig()
    windows, ledger = [], []
    for context in contexts:
        eps = []
        for row in by_context[context['context_id']]:
            peptide = row['peptide']
            preds = tuple(Prediction(kind=Kind.pMHC_affinity, peptide=peptide, allele=a, score=affinity[peptide, a][0],
                value=affinity[peptide, a][0], predictor_name='mhcflurry', predictor_version='2.3.14') for a in alleles)
            eps.append(CandidateEpitope(sequence=peptide, source_sequence=context['sequence'], source_name=context['context_id'],
                offset=row['start'], predictions=preds, source_class='mutation', overlaps_mutation=True, overlaps_targetable=True,
                patient_alleles=tuple(alleles), prediction_id=context['context_id'] + ':' + str(row['start']) + ':' + str(row['end'])))
        scored = epitopes_for_ranking(attach_per_allele_scores(eps, config), config)
        for start in range(len(context['sequence']) - 24):
            end = start + 25
            if not start <= context['mutation_start'] < end:
                continue
            seq = context['sequence'][start:end]
            hydro = summarize(seq)
            inside = [o for o in by_context[context['context_id']] if start <= o['start'] and o['end'] <= end]
            vp = VaccinePeptide(antigen=antigen(context).sliced(start, end), epitopes=slice_epitopes(scored, start, end),
                num_target_epitopes_to_keep=0, combined_score_expr='target_epitope_score', ranking_rules=('target_epitope_score',))
            score = float(vp.combined_score)
            def contribution(nm):
                result = 0 if nm >= 5000 else (1 + math.exp(-350/150)) / (1 + math.exp((nm - 350)/150))
                return result if result >= 1e-5 else 0
            expected_score = sum(contribution(affinity[o['peptide'], a][0]) for o in inside for a in alleles)
            assert math.isclose(score, expected_score, abs_tol=1e-8, rel_tol=1e-8)
            best = min(((affinity[o['peptide'], a][0], o['peptide'], a) for o in inside for a in alleles), default=(None, '', ''))
            windows.append(dict(variant_id=context['variant_id'], context_id=context['context_id'], gene=context['gene'],
                protein_change=context['protein_change'], sequence=seq, reference_sequence=context['reference_sequence'][start:end],
                window_start=start, window_end=end, target_start_in_window=context['mutation_start'] - start,
                target_end_in_window=context['mutation_end'] - start, affinity_contribution_sum=score,
                best_target_affinity_nm=best[0], best_target_peptide=best[1], best_target_allele=best[2],
                passes_500nm_gate=best[0] is not None and best[0] <= 500, hydrophobicity=hydro,
                hydrophobicity_pass=hydro['status'] == 'PASS', distinct_target_peptides=len({o['peptide'] for o in inside}),
                centering_distance=abs((context['mutation_start'] - start) - 12)))
    chosen = []
    for variant_id in dict.fromkeys(c['variant_id'] for c in contexts):
        options = [r for r in windows if r['variant_id'] == variant_id and r['hydrophobicity_pass']]
        assert options, 'No hydropathy-compatible window: ' + variant_id
        best = min(options, key=lambda r: (-r['affinity_contribution_sum'], r['centering_distance'], r['context_id'], r['window_start']))
        context = next(c for c in contexts if c['context_id'] == best['context_id'])
        row = dict(best)
        row.update(source='MUTATION', source_subtype='GUARDANT_REFERENCE_DERIVED', event_id=variant_id,
            region_id='REG_' + hashlib.sha256((variant_id + ':' + row['sequence']).encode()).hexdigest()[:20],
            normalized_allele=context['normalized_allele'], transcript_ids=context['transcript_ids'],
            DNA_annotation_protein_change=context['protein_change'], DNA_annotation_transcript_id=context['transcript_ids'][0],
            native_affinity_sum=row['affinity_contribution_sum'], weighted_score=None, expression_proxy=None,
            supplied_gene_value=expression.get(row['gene']), rna_alternate_fraction=None,
            primary_DNA_eligible=False, primary_window_eligible=True, base_primary_window_eligible=False,
            manual_inclusion=True, research_only=True,
            admission_exception='GUARDANT_REPORTED_VUS_REFERENCE_DERIVED_NO_TUMOUR_RNA_RECONSTRUCTION',
            somatic_status='PLASMA_REPORTED_NORMAL_ALT_NOT_DETECTED', DNA_evidence_tier='GUARDANT_CALL_NORMAL_REFERENCE_SUPPORTED',
            patient_RNA_sequence_confirmed=False, sequence_basis=context['sequence_basis'],
            guardant_cfDNA_percent=context['guardant_cfDNA_percent'], normal_classification=context['normal']['classification'],
            normal_ref_fragments=context['normal']['read_evidence']['fragment_counts'].get('REF', 0),
            normal_alt_fragments=context['normal']['target_alternate_normal_fragments'],
            January_tumour_DNA=context['tumour_DNA']['fragment_counts'], January_tumour_RNA=context['tumour_RNA']['fragment_counts'])
        if row['normal_classification'] != 'NORMAL_REFERENCE_SUPPORTED':
            row['admission_exception'] += ';NORMAL_COVERAGE_BELOW_20_FRAGMENTS'
            row['DNA_evidence_tier'] = 'GUARDANT_CALL_NORMAL_COVERAGE_BORDERLINE'
        if not row['passes_500nm_gate']:
            row['admission_exception'] += ';AFFINITY_ABOVE_500NM_OR_NO_ELIGIBLE_NONSELF_PEPTIDE'
        chosen.append(row)
        print('SELECTED', row['gene'], row['sequence'], row['best_target_affinity_nm'], flush=True)
    assert len(chosen) == 5
    parent_selected = json.loads((ROOT / 'inputs/parent-selected.json').read_text())
    assert len(parent_selected) == 24
    policy = json.loads((ROOT / 'policy.json').read_text())
    audit = json.loads((ROOT / 'inputs/sequence-review/summary.json').read_text())
    excluded = set(policy['excluded_parent_region_ids'])
    assert {r['region_id'] for r in audit if r['review_status'] == 'NO_TARGET_REMAINS_AFTER_BROAD_SELF_EXCLUSION'} == excluded
    retained = [r for r in parent_selected if r['region_id'] not in excluded]
    assert len(retained) == policy['parent_regions_retained'] == 23
    selected = retained + [dict(r, antigen_id='AG' + str(i + 25).zfill(2)) for i, r in enumerate(chosen)]
    assert len({r['region_id'] for r in selected}) == 28
    integrate = ROOT / 'integrate'
    integrate.mkdir(exist_ok=True)
    save(integrate / 'selected.json', selected)
    table(integrate / 'selected.tsv', selected)
    (integrate / 'selected.fasta').write_text(''.join('>' + r['antigen_id'] + ' ' + r['gene'] + '\n' + r['sequence'] + '\n' for r in selected))
    table(out / 'all-windows.tsv', windows)
    table(out / 'selected-additions.tsv', chosen)
    save(out / 'selected-additions.json', chosen)
    save(out / 'inference.json', dict(status='COMPLETED', all_queries_freshly_predicted=True,
        peptides=len(peptides), peptide_HLA_pairs=len(predictions), alleles=alleles,
        gpu=torch.cuda.get_device_name() if torch.cuda.is_available() else None,
        device='CUDA' if torch.cuda.is_available() else 'CPU', class_II=False, presentation_model_used=False,
        exact_reference_matching_peptides=sum(x.occurs for x in references.values()),
        patient_self_matching_peptides=sum(x['occurs_in_searched_models'] for x in normal_matches.values()),
        native_Vaxrank_window_scores_verified=True, native_mhctools_affinity_pairs=len(affinity),
        affinity_adapter='mhctools.mhcflurry.MHCflurry_Affinity',
        versions={n: importlib.metadata.version(n) for n in ['vaxrank', 'mhctools', 'mhcflurry', 'torch', 'numpy']},
        parent_selected_input_sha256=sha(ROOT / 'inputs/parent-selected.json'), contexts_sha256=sha(ROOT / 'contexts/contexts.json'),
        model_manifest_sha256=sha(ROOT / 'affinity-model-manifest.json'), verified_model_files=checked,
        worker_sha256=sha(__file__)))
    receipt(out, [p for p in out.iterdir() if p.is_file() and p.name not in ['complete.json', 'status.json']],
            worker_sha256=sha(__file__), selected_additions=5, windows=len(windows))
    receipt(integrate, [integrate / name for name in ['selected.json', 'selected.tsv', 'selected.fasta']],
            retained_parent_regions=len(retained), excluded_parent_region_ids=sorted(excluded), explicit_Guardant_additions=len(chosen), total_regions=len(selected),
            source_inference_sha256=sha(out / 'complete.json'), worker_sha256=sha(__file__))


if __name__ == '__main__':
    {'prepare': prepare, 'predict-select': predict_select}[sys.argv[2]]()
