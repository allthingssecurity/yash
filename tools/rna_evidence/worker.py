"""Read-only case RNA-sequence audit: normal DNA, broad references and personal-self models."""
from pathlib import Path
from collections import Counter, defaultdict
import argparse, bisect, csv, datetime, gzip, hashlib, json, sys
from Bio.Seq import Seq


def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda:f.read(8*1024*1024),b''):h.update(block)
    return h.hexdigest()


def save(path,value):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(value,indent=2,default=str)+'\n')


def translation_records(path):
    header=None;parts=[]
    with gzip.open(path,'rt') as f:
        for line in f:
            if line.startswith('>'):
                if header is not None:yield header,''.join(parts)
                header=line[1:].strip();parts=[]
            else:parts.append(line.strip())
    if header is not None:yield header,''.join(parts)


def primary_coding(header):
    import re
    m=re.search(r'chromosome:GRCh38:([^:]+):',header)
    return bool(m and m[1] in {str(i) for i in range(1,23)}|{'X','Y','MT'}
                and 'gene_biotype:protein_coding ' in header
                and 'transcript_biotype:protein_coding ' in header)


def library(records):
    headers=[];starts=[];parts=[];offset=0
    for header,sequence in records:
        headers.append(header);starts.append(offset);parts.append(sequence);offset+=len(sequence)+1
    return '!'.join(parts),starts,headers


def exact_hit(lib,peptide):
    sequence,starts,headers=lib
    offset=sequence.find(peptide)
    if offset<0:return None
    i=bisect.bisect_right(starts,offset)-1
    return dict(reference_header=headers[i],protein_offset_0based=offset-starts[i])


def direct_counts(bam,chrom,position):
    names=defaultdict(set)
    for read in bam.fetch(chrom,position-1,position):
        if (read.is_unmapped or read.is_secondary or read.is_supplementary or read.is_duplicate
                or read.is_qcfail or read.mapping_quality<20 or read.query_qualities is None):continue
        if read.has_tag('NH') and read.get_tag('NH')!=1:continue
        for q,r in read.get_aligned_pairs(matches_only=True):
            if r==position-1 and read.query_qualities[q]>=20:names[read.query_name].add(read.query_sequence[q].upper())
    counts=Counter(next(iter(v)) if len(v)==1 else 'DISAGREEING_MATES' for v in names.values())
    return dict(fragment_bases=dict(counts),qualifying_names=len(names),unit='Fragment names, not UMI-confirmed independent molecules',
                policy='MQ>=20, focal BQ>=20; exclude unmapped, secondary, supplementary, duplicate and QC-fail reads; NH=1 when present')


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root',type=Path,required=True)
    parser.add_argument('--request',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();root=args.root;out=args.output
    if (out/'complete.json').exists():raise FileExistsError('Use a new output directory; preserve prior measurements')
    out.mkdir(parents=True,exist_ok=True)
    request=json.loads(args.request.read_text())
    import pysam,varcode
    from pyensembl import Genome
    from normal_evidence import NormalEvidence
    from self_search import patient_self_matches
    normal=NormalEvidence(root)
    config=json.loads((root/'analysis-genome-config.json').read_text());genome=Genome(**config)
    identity=json.loads((root/'identity/complete.json').read_text());assert identity['identity']=='MATCH'
    plan=json.loads((root/'normal-coding-plan.json').read_text())
    assert Path(plan['normal_bam']).stat().st_size==143294724201
    normal_done=json.loads((root/'normal-coding/gather/complete.json').read_text())
    for name in ['normal.coding.g.vcf.gz','normal.coding.g.vcf.gz.tbi']:
        assert sha(root/'normal-coding/gather'/name)==normal_done['outputs'][name]['sha256']
    assert normal_done['normal_bam_sha256']==plan['normal_bam_sha256']
    tumour=pysam.AlignmentFile('/data/inputs/TN25-112103/DNA_TN25-112103.bam','rb')
    rna=pysam.AlignmentFile('/data/inputs/TN25-112103/RNA_TN25-112103.bam','rb')
    tumour.check_index();rna.check_index()
    candidates=[Path('/data/raft/references/homo_sapiens/pyensembl/GRCh38/ensembl103/Homo_sapiens.GRCh38.pep.all.fa.gz'),
                Path('/data/results/TN25-112103/mutation-evidence-v2-20261102/reference-cache/ensembl103/Homo_sapiens.GRCh38.pep.all.fa.gz')]
    fasta=next(p for p in candidates if p.is_file())
    fasta_sha=sha(fasta)
    assert fasta_sha=='f8df3d7cfc85718e9161579bd9168a3598547ee50022b5e1c272ed66293f0505'
    records=list(translation_records(fasta));primary_records=[r for r in records if primary_coding(r[0])]
    broad=library(records);primary=library(primary_records)
    windows=[];all_peptides=set();normal_rows=[]
    comp=str.maketrans('ACGTN','TGCAN')
    for w in request['windows']:
        assert str(Seq(w['coding_sequence']).translate())==w['sequence']
        chrom=w['coordinates'].split(':')[0];minus=w['coordinates'].endswith('(-)')
        assert None not in w['genomic_positions_0based']
        reference=''.join(normal.fasta.fetch(chrom,p,p+1).upper().translate(comp) if minus else normal.fasta.fetch(chrom,p,p+1).upper() for p in w['genomic_positions_0based'])
        aa=str(Seq(reference).translate());differences=[]
        for i,(a,b) in enumerate(zip(reference,w['coding_sequence'])):
            if a==b:continue
            pos=w['genomic_positions_0based'][i]+1;ref=a.translate(comp) if minus else a;alt=b.translate(comp) if minus else b
            allele=[chrom,pos,ref,alt]
            variant=varcode.Variant(chrom,pos,ref,alt,ensembl=genome)
            normal_result=normal.assess(allele,variant)
            counts={name:direct_counts(bam,chrom,pos) for name,bam in [('normal',normal.bam),('tumour_DNA',tumour),('tumour_RNA',rna)]}
            nc=counts['normal']['fragment_bases'];isovar=normal_result['read_evidence']['fragment_counts']
            assert nc.get(alt,0)==isovar.get('ALT',0),(w['gene'],allele,'Independent ALT counts disagree')
            codon=i//3
            d=dict(window_nt_position_1based=i+1,allele=allele,reference_codon=reference[3*codon:3*codon+3],rna_codon=w['coding_sequence'][3*codon:3*codon+3],
                   aa_position_1based=codon+1,reference_amino_acid=aa[codon],rna_inferred_amino_acid=w['sequence'][codon],
                   nonsynonymous=aa[codon]!=w['sequence'][codon],normal=normal_result,independent_counts=counts)
            differences.append(d)
            normal_rows.append(dict(gene=w['gene'],region_id=w['region_id'],**d))
            print('BASE',w['gene'],allele,'normal',nc,'native',normal_result['genotype_evidence']['status'],flush=True)
        lo,hi=w['target_start_in_window'],w['target_end_in_window']
        peptides=[]
        for length in range(8,12):
            for start in range(len(w['sequence'])-length+1):
                peptide=w['sequence'][start:start+length]
                peptides.append(dict(peptide=peptide,start_0based=start,end_0based_exclusive=start+length,
                  overlaps_recorded_target=start<hi and start+length>lo,
                  overlaps_nonsynonymous_RNA_difference=any(d['nonsynonymous'] and start<=d['aa_position_1based']-1<start+length for d in differences)))
                all_peptides.add(peptide)
        windows.append(dict(**w,reference_translation_on_recorded_path=aa,nucleotide_differences=differences,peptide_occurrences=peptides))
    all_peptides.update(request.get('additional_target_peptides',[]))
    for region in request.get('additional_regions',[]):
        for length in range(8,12):
            all_peptides.update(region['sequence'][i:i+length] for i in range(len(region['sequence'])-length+1))
    personal=patient_self_matches(root,all_peptides)
    search={p:dict(primary=exact_hit(primary,p),extended=exact_hit(broad,p),personal=personal[p]) for p in sorted(all_peptides)}
    summary=[]
    for w in windows:
        targets={r['peptide'] for r in w['peptide_occurrences'] if r['overlaps_recorded_target']}
        eligible={p for p in targets if search[p]['extended'] is None and not search[p]['personal']['occurs_in_searched_models']}
        full=exact_hit(broad,w['sequence'])
        w['broad_reference_full_window_hit']=full
        w['target_peptides_after_broad_reference_and_personal_self_exclusion']=sorted(eligible)
        summary.append(dict(gene=w['gene'],region_id=w['region_id'],sequence=w['sequence'],discovery_source='SPLICE',
          description='RNA-supported inferred amino-acid differences' if any(d['nonsynonymous'] for d in w['nucleotide_differences']) else 'RNA-derived junction-context hypothesis',
          full_window_fragment_names=w['fragment_count'],nonsynonymous_differences=[dict(position_1based=d['aa_position_1based'],reference=d['reference_amino_acid'],observed=d['rna_inferred_amino_acid']) for d in w['nucleotide_differences'] if d['nonsynonymous']],
          target_peptides=len(targets),primary_nonself_target_peptides=sum(search[p]['primary'] is None for p in targets),
          broad_nonself_target_peptides=sum(search[p]['extended'] is None for p in targets),
          broad_and_personal_nonself_target_peptides=len(eligible),
          full_window_extended_reference_match=full,
          review_status='NO_TARGET_REMAINS_AFTER_BROAD_SELF_EXCLUSION' if not eligible else 'RETAIN_FOR_REVIEW_WITH_SCOPE_LIMITATIONS',
          true_splice_origin_proven=False,normal_tissue_absence_proven=False,protein_translation_measured=False))
    for w in request.get('additional_regions',[]):
        lo,hi=w['target_start_in_window'],w['target_end_in_window']
        targets={w['sequence'][i:i+k] for k in range(8,12) for i in range(len(w['sequence'])-k+1) if i<hi and i+k>lo}
        eligible={p for p in targets if search[p]['extended'] is None and not search[p]['personal']['occurs_in_searched_models']}
        summary.append(dict(gene=w['gene'],region_id=w['region_id'],sequence=w['sequence'],discovery_source=w['source'],
            target_peptides=len(targets),primary_nonself_target_peptides=sum(search[p]['primary'] is None for p in targets),
            broad_nonself_target_peptides=sum(search[p]['extended'] is None for p in targets),
            broad_and_personal_nonself_target_peptides=len(eligible),full_window_extended_reference_match=exact_hit(broad,w['sequence']),
            review_status='NO_TARGET_REMAINS_AFTER_BROAD_SELF_EXCLUSION' if not eligible else 'RETAIN_FOR_REVIEW_WITH_SCOPE_LIMITATIONS',
            normal_DNA_difference_pileups_in_this_audit=False,normal_tissue_absence_proven=False,protein_translation_measured=False))
    save(out/'windows.json',windows);save(out/'normal-dna.json',normal_rows);save(out/'reference-matches.json',search);save(out/'summary.json',summary)
    save(out/'reference-scope.json',dict(reference=fasta.name,sha256=fasta_sha,primary_coding_translations=len(primary_records),all_provider_translations=len(records),
      extended_scope='All Ensembl 103 peptide translations, including alternative contigs and pseudogene annotations; these are predicted sequences, not demonstrated healthy-tissue proteins',
      patient_normal_scope='Recorded bounded normal-coding-genotype and typed-HLA sequence models; not exhaustive normal proteome or complete repeat haplotypes',
      normal_bam_sha256=plan['normal_bam_sha256'],normal_bam_bytes=143294724201,normal_coding_receipt_sha256=sha(root/'normal-coding/gather/complete.json'),
      reference_genome=plan['reference'],RNA_support_scope='Existing full-window reconstruction and BAM verification receipts; fresh direct pileups only at the sequence-difference positions',
      CTA_self_scope='Intentionally reference-derived CTA/self candidates are outside this altered-peptide novelty audit'))
    paths=sorted(p for p in out.iterdir() if p.is_file())
    save(out/'complete.json',dict(status='COMPLETED',completed_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
      request_sha256=sha(args.request),worker_sha256=sha(__file__),outputs={p.name:dict(bytes=p.stat().st_size,sha256=sha(p)) for p in paths}))
    print('COMPLETED',out,flush=True)


if __name__=='__main__':main()
