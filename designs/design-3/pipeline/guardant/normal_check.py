"""Map Guardant protein-only calls and measure the existing case BAMs read-only."""
from collections import Counter, defaultdict
from pathlib import Path
import csv
import datetime
import hashlib
import importlib.metadata
import itertools
import json
import sys

import pysam
import varcode
from pyensembl import Genome
from normal_evidence import NormalEvidence

ROOT = Path('/data/results/TN25-112103/matched-normal-20261107')
OUT = Path(sys.argv[1])
OUT.mkdir(parents=True, exist_ok=True)
QUERIES = [
    ('AR', 'P555S', 0.4, 'unmatched'),
    ('DHX36', 'R160K', 0.2, 'unmatched'),
    ('ERBB4', 'D113V', 0.3, 'unmatched'),
    ('KMT2D', 'R2789Q', 3.2, 'unmatched'),
    ('RAF1', 'G535S', 0.2, 'unmatched'),
    ('RXRA', 'A214S', 12.1, 'previously_assessed_control'),
    ('JUN', 'R276*', 14.8, 'previously_assessed_control'),
    ('MERTK', 'E192*', 19.4, 'previously_assessed_control'),
]
CODONS = dict(zip((''.join(x) for x in itertools.product('TCAG', repeat=3)),
                  'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'))
COMPLEMENT = str.maketrans('ACGT', 'TGCA')


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def save(name, value):
    path = OUT / name
    path.write_text(json.dumps(value, indent=2, default=str) + '\n')


def measured_counts(normal, variant, bam):
    original = normal.bam
    try:
        normal.bam = bam
        return normal.read_counts(variant)
    finally:
        normal.bam = original


def independent_bases(bam, chrom, position):
    """Independent direct alignment-position audit, without the Isovar collector."""
    fragments = defaultdict(set)
    orientations = defaultdict(set)
    for read in bam.fetch(chrom, position - 1, position):
        if (read.is_unmapped or read.is_secondary or read.is_supplementary or
                read.is_duplicate or read.is_qcfail or read.mapping_quality < 20):
            continue
        if read.has_tag('NH') and read.get_tag('NH') != 1:
            continue
        if read.query_sequence is None or read.query_qualities is None:
            continue
        for query_offset, reference_offset in read.get_aligned_pairs(matches_only=True):
            if reference_offset == position - 1 and read.query_qualities[query_offset] >= 20:
                fragments[read.query_name].add(read.query_sequence[query_offset].upper())
                orientations[read.query_name].add('reverse' if read.is_reverse else 'forward')
    bases = Counter(next(iter(values)) if len(values) == 1 else 'DISAGREEING_MATES'
                    for values in fragments.values())
    return dict(fragment_bases=dict(bases), qualifying_names=len(fragments),
                unit='Read names, not UMI-defined independent molecules',
                policy='MQ>=20, focal BQ>=20; no unmapped, secondary, supplementary, duplicate or QC-fail reads; NH=1 when present')


config_path = ROOT / 'analysis-genome-config.json'
config = json.loads(config_path.read_text())
genome = Genome(**config)
normal = NormalEvidence(ROOT)
normal_plan = json.loads((ROOT / 'normal-coding-plan.json').read_text())
normal_receipt = json.loads((ROOT / 'normal-coding/gather/complete.json').read_text())
for name in ['normal.coding.g.vcf.gz', 'normal.coding.g.vcf.gz.tbi']:
    assert sha(ROOT / 'normal-coding/gather' / name) == normal_receipt['outputs'][name]['sha256']
assert normal_plan['normal_bam_sha256'] == normal_receipt['normal_bam_sha256']
normal_bam = Path(normal_plan['normal_bam'])
assert normal_bam.stat().st_size == 143294724201
identity = json.loads((ROOT / 'identity/complete.json').read_text())
assert identity['identity'] == 'MATCH'
tumour = pysam.AlignmentFile('/data/inputs/TN25-112103/DNA_TN25-112103.bam', 'rb')
rna = pysam.AlignmentFile('/data/inputs/TN25-112103/RNA_TN25-112103.bam', 'rb')
tumour.check_index()
rna.check_index()
paired = pysam.VariantFile(ROOT / 'paired/filter/all-unique.vcf.gz')
for name, bam in [('normal', normal.bam), ('tumour', tumour), ('RNA', rna)]:
    assert 'chr1' in bam.references, (name, 'Unexpected reference names')

results, issues, audits = [], [], []
for gene_name, change, fraction, scope in QUERIES:
    print('MAPPING', gene_name, change, flush=True)
    ref_aa, alt_aa, aa_position = change[0], change[-1], int(change[1:-1])
    candidates = {}
    transcript_audit = []
    for gene in genome.genes_by_name(gene_name):
        for transcript in gene.transcripts:
            audit = dict(gene=gene_name, change=change, transcript_id=transcript.id)
            try:
                protein = transcript.protein_sequence
                if not protein or len(protein) < aa_position:
                    audit['status'] = 'NO_PROTEIN_OR_POSITION_OUT_OF_RANGE'
                    transcript_audit.append(audit)
                    continue
                if protein[aa_position - 1] != ref_aa:
                    audit.update(status='REFERENCE_AMINO_ACID_MISMATCH', observed=protein[aa_position - 1])
                    transcript_audit.append(audit)
                    continue
                cds = transcript.coding_sequence
                ranges = sorted(transcript.coding_sequence_position_ranges)
                positions = [pos for start, end in ranges for pos in range(start, end + 1)]
                if transcript.strand == '-':
                    positions.reverse()
                assert len(cds) == len(positions), (len(cds), len(positions))
                offset = 3 * (aa_position - 1)
                codon = cds[offset:offset + 3].upper()
                assert CODONS[codon] == ref_aa
                chrom = transcript.contig if transcript.contig.startswith('chr') else 'chr' + transcript.contig
                genomic_codon = ''.join(normal.fasta.fetch(chrom, pos - 1, pos).upper()
                                        for pos in positions[offset:offset + 3])
                if transcript.strand == '-':
                    genomic_codon = genomic_codon.translate(COMPLEMENT)
                assert genomic_codon == codon, (genomic_codon, codon)
                found = []
                for i, original in enumerate(codon):
                    for alternative in 'ACGT':
                        if alternative == original:
                            continue
                        mutant = codon[:i] + alternative + codon[i + 1:]
                        if CODONS[mutant] != alt_aa:
                            continue
                        position = positions[offset + i]
                        reference_base = original if transcript.strand == '+' else original.translate(COMPLEMENT)
                        alternate_base = alternative if transcript.strand == '+' else alternative.translate(COMPLEMENT)
                        key = chrom, position, reference_base, alternate_base
                        context = dict(transcript_id=transcript.id, protein_id=transcript.protein_id,
                                       strand=transcript.strand, reference_codon=codon,
                                       alternate_codon=mutant, codon_genomic_positions=positions[offset:offset + 3],
                                       local_reference_protein=protein[max(0, aa_position - 11):aa_position + 10])
                        candidates.setdefault(key, []).append(context)
                        found.append(list(key))
                audit.update(status='MAPPED' if found else 'NO_SINGLE_NUCLEOTIDE_SUBSTITUTION', alleles=found)
            except Exception as error:
                audit.update(status='TRANSCRIPT_MAPPING_REQUIRES_REVIEW', error=repr(error))
                issues.append(audit)
            transcript_audit.append(audit)
    mappings = []
    for key, transcripts in sorted(candidates.items()):
        chrom, position, reference_base, alternate_base = key
        variant = varcode.Variant(chrom, position, reference_base, alternate_base,
                                  genome=genome, convert_ucsc_contig_names=True)
        effects = [dict(gene=e.gene_name, transcript_id=getattr(e, 'transcript_id', None),
                        protein_change=getattr(e, 'short_description', None), consequence=type(e).__name__)
                   for e in variant.effects()]
        roundtrip = any(e['gene'] == gene_name and e['protein_change'] == 'p.' + change for e in effects)
        assert roundtrip, (key, gene_name, change, effects)
        native = []
        for record in paired.fetch(chrom, max(0, position - 2), position + 1):
            if record.pos == position and record.ref == reference_base and alternate_base in (record.alts or []):
                native.append(dict(position=record.pos, ref=record.ref, alts=record.alts,
                                   filters=list(record.filter), INFO=dict(record.info),
                                   samples={name: dict(record.samples[name]) for name in record.samples}))
        mapping = dict(allele=list(key), compatible_transcripts=transcripts,
                       protein_effect_roundtrip=True, effects=effects,
                       normal=normal.assess(list(key), variant),
                       normal_direct_pileup=independent_bases(normal.bam, chrom, position),
                       tumour_DNA=measured_counts(normal, variant, tumour),
                       tumour_DNA_direct_pileup=independent_bases(tumour, chrom, position),
                       tumour_RNA=measured_counts(normal, variant, rna),
                       paired_native_records=native)
        mappings.append(mapping)
        print('MEASURED', gene_name, change, list(key), mapping['normal']['classification'],
              mapping['normal']['read_evidence']['fragment_counts'], flush=True)
    results.append(dict(gene=gene_name, reported_protein_change=change, October_cfDNA_percent=fraction,
                        scope=scope, mapping_status='NO_MAPPING' if not mappings else
                        'UNIQUE_COMPATIBLE_ALLELE' if len(mappings) == 1 else 'MULTIPLE_COMPATIBLE_ALLELES',
                        candidate_count=len(mappings), mappings=mappings,
                        interpretation='Protein-only Guardant label; these are compatible Ensembl-103 interpretations, not verified Guardant genomic coordinates'))
    audits.extend(transcript_audit)
    save('results.json', results)
    save('transcript-mapping-audit.json', audits)
    save('mapping-issues.json', issues)

controls = {'RXRA': ('chr9', 134417268, 'G', 'T'),
            'JUN': ('chr1', 58782245, 'T', 'A'),
            'MERTK': ('chr2', 111945051, 'G', 'T')}
for gene, expected in controls.items():
    assert expected in {tuple(m['allele']) for r in results if r['gene'] == gene for m in r['mappings']}
flat = []
for row in results:
    for mapping in row['mappings']:
        normal_counts = mapping['normal']['read_evidence']['fragment_counts']
        tumour_counts = mapping['tumour_DNA']['fragment_counts']
        rna_counts = mapping['tumour_RNA']['fragment_counts']
        flat.append(dict(gene=row['gene'], reported_change=row['reported_protein_change'],
                         mapping_status=row['mapping_status'], allele=':'.join(map(str, mapping['allele'])),
                         normal_classification=mapping['normal']['classification'],
                         normal_REF=normal_counts.get('REF', 0), normal_ALT=normal_counts.get('ALT', 0),
                         normal_OTHER=normal_counts.get('OTHER', 0),
                         normal_genotype_status=mapping['normal']['genotype_evidence']['status'],
                         tumour_DNA_REF=tumour_counts.get('REF', 0), tumour_DNA_ALT=tumour_counts.get('ALT', 0),
                         tumour_RNA_REF=rna_counts.get('REF', 0), tumour_RNA_ALT=rna_counts.get('ALT', 0),
                         native_paired_record_count=len(mapping['paired_native_records'])))
with (OUT / 'summary.tsv').open('w') as stream:
    writer = csv.DictWriter(stream, fieldnames=list(flat[0]), delimiter='\t', lineterminator='\n')
    writer.writeheader()
    writer.writerows(flat)
save('complete.json', dict(status='COMPLETED', completed_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
                          inputs=dict(guardant_collection_date='2025-10-14', tumour_collection_date='2025-01-15',
                                      normal_collection_date='2026-10-30', normal_bam=str(normal_bam),
                                      normal_bam_sha256_from_recorded_receipt=normal_plan['normal_bam_sha256'],
                                      normal_bam_full_hash_recomputed=False, normal_bam_size_verified=normal_bam.stat().st_size,
                                      normal_GVCF_sha256_verified=normal_receipt['outputs']['normal.coding.g.vcf.gz']['sha256'],
                                      genome_config=config, genome_config_sha256=sha(config_path)),
                          versions={n: importlib.metadata.version(n) for n in ['pyensembl', 'varcode', 'isovar', 'pysam']},
                          worker_sha256=sha(__file__), helper_sha256=sha(Path(__file__).with_name('normal_evidence.py')),
                          queries=len(results), compatible_alleles=len(flat),
                          mapping_issues=len(issues), controls_roundtrip_passed=True,
                          normal_absence_is_not_proof_of_tumour_origin=True,
                          existing_case_outputs_modified=False))
print('COMPLETED', len(results), 'queries;', len(flat), 'compatible alleles', flush=True)
