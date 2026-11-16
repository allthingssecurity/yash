"""Verify collected audit files, optionally recheck the provider FASTA, and render their findings."""
from pathlib import Path
import argparse,csv,gzip,hashlib,json,re

ROOT=Path(__file__).resolve().parents[2]
FOCAL={'F5':8,'MUC5B':12,'MUC3A':14,'MUC6':6}


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def rows(path):return json.loads(Path(path).read_text())


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--results',type=Path,default=Path('designs/design-3/inputs/sequence-review'))
    p.add_argument('--reference-proteins',type=Path,help='Optional Ensembl 103 provider FASTA for a new exact-match check; no model inference')
    a=p.parse_args();base=(ROOT/a.results).resolve()
    receipt=rows(base/'complete.json')
    assert receipt['status']=='COMPLETED' and receipt['request_sha256']==sha(base/'request.json')
    assert receipt['worker_sha256']==sha(ROOT/'tools/rna_evidence/worker.py')
    for name,r in receipt['outputs'].items():
        q=base/name;assert q.stat().st_size==r['bytes'] and sha(q)==r['sha256'],name
    req=rows(base/'request.json');assert req['selected_sha256']==sha(base/'candidate-pool.json')
    assert req['patient_window_sha256']==sha(ROOT/'data/shared-rna/splicing/reconstruction/patient-windows.jsonl.gz')
    assert req['existing_BAM_verification_sha256']==sha(ROOT/'data/shared-rna/splicing/read-verification/verification.json')
    search=rows(base/'reference-matches.json');scope=rows(base/'reference-scope.json')
    new_reference_check=False
    if a.reference_proteins:
        assert sha(a.reference_proteins)==scope['sha256']
        all_parts=[];primary_parts=[];header=None;parts=[]
        def add():
            if header is None:return
            sequence=''.join(parts);all_parts.append(sequence)
            m=re.search(r'chromosome:GRCh38:([^:]+):',header)
            if m and m[1] in {str(i) for i in range(1,23)}|{'X','Y','MT'} and 'gene_biotype:protein_coding ' in header and 'transcript_biotype:protein_coding ' in header:primary_parts.append(sequence)
        with gzip.open(a.reference_proteins,'rt') as f:
            for line in f:
                if line.startswith('>'):add();header=line[1:].strip();parts=[]
                else:parts.append(line.strip())
        add()
        assert len(all_parts)==scope['all_provider_translations'] and len(primary_parts)==scope['primary_coding_translations']
        broad='!'.join(all_parts);primary='!'.join(primary_parts)
        for peptide,r in search.items():
            assert (peptide in broad)==bool(r['extended']),peptide
            assert (peptide in primary)==bool(r['primary']),peptide
        new_reference_check=True
    summaries=rows(base/'summary.json');by_gene={r['gene']:r for r in summaries}
    windows=rows(base/'windows.json');normal=rows(base/'normal-dna.json')
    assert len(windows)==4 and len(summaries)==24
    assert {r['gene'] for r in summaries if r['review_status']=='NO_TARGET_REMAINS_AFTER_BROAD_SELF_EXCLUSION'}=={'MUC6'}
    normal_table=[]
    for gene,pos in FOCAL.items():
        r=next(x for x in normal if x['gene']==gene and x['aa_position_1based']==pos)
        ref,alt=r['allele'][2:];counts=r['independent_counts']['normal']['fragment_bases']
        assert counts.get(alt,0)==0 and counts.get(ref,0)>=20
        assert r['normal']['classification']=='NORMAL_REFERENCE_SUPPORTED'
        observations=r['normal']['genotype_evidence']['observations']
        gqs=[x['genotype'].get('GQ') for x in observations if x['high_quality']]
        normal_table.append(dict(gene=gene,reference_amino_acid=r['reference_amino_acid'],RNA_inferred_amino_acid=r['rna_inferred_amino_acid'],
            reference_RNA_codon=r['reference_codon'].replace('T','U'),observed_RNA_codon=r['rna_codon'].replace('T','U'),
            chromosome=r['allele'][0],position_1based=r['allele'][1],normal_reference_fragments=counts.get(ref,0),normal_alternate_fragments=counts.get(alt,0),
            native_genotype_status=r['normal']['genotype_evidence']['status'],minimum_native_GQ=min(gqs),
            primary_nonself_target_peptides=by_gene[gene]['primary_nonself_target_peptides'],broad_and_personal_nonself_target_peptides=by_gene[gene]['broad_and_personal_nonself_target_peptides'],
            decision='EXCLUDE_REFERENCE_MATCH' if gene=='MUC6' else 'RETAIN_RNA_SEQUENCE_HYPOTHESIS'))
    with (base/'summary.tsv').open('w') as f:
        w=csv.DictWriter(f,fieldnames=list(normal_table[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(normal_table)
    lines=['# RNA-supported amino-acid changes and normal-sequence screening','',
      'The RNA reads contain nucleotide sequences predicted to encode amino-acid changes relative to the nominated reference. Splice analysis discovered these candidates, but repeats leave their exact genomic placement and splice origin uncertain. The bases are useful sequence evidence; they are not direct protein measurements or confirmed RNA-editing events.','',
      '## Matched-normal DNA','',
      'Fresh focal-base pileups and native normal genotypes were checked at every nucleotide difference in the four selected RNA windows. Reference and alternate counts use fragment names, not UMI-confirmed independent molecules. Both the Isovar collector and a direct alignment-position implementation agree on alternate counts.','',
      '| Candidate | Reference → RNA-inferred amino acid | Reference → observed RNA codon | Normal ALT / REF fragments | Native normal GQ |','|---|---|---|---:|---:|']
    for r in normal_table:
        lines.append(f"| {r['gene']} | {r['reference_amino_acid']}→{r['RNA_inferred_amino_acid']} | `{r['reference_RNA_codon']}→{r['observed_RNA_codon']}` | {r['normal_alternate_fragments']} / {r['normal_reference_fragments']} | {r['minimum_native_GQ']} |")
    lines+=['','No alternate fragments were detected at these four positions. This is reassuring at the nominated loci; it does not establish the origin of a repetitive read, exclude every inherited repeat copy or prove tumour-specific translation. Normal DNA also cannot establish whether an unusual RNA join occurs in healthy tissue.','',
      'The full-window check additionally records MUC6 R→P at the first residue: 24 alternate normal fragments, but a reference native genotype, giving unresolved evidence rather than a confident inherited call. A third MUC6 nucleotide difference is synonymous and present in normal DNA. Neither changes the separate T→M result above.','',
      '## Broader reference and patient-self search','',
      'All 8–11-aa peptides in the 24 altered-source regions of the pre-filter pool were searched against all **115,262 Ensembl 103 peptide translations**, the **85,968 primary protein-coding subset**, and the recorded bounded patient-normal coding/HLA models. Target counts below use the existing target intervals. Intentional CTA/self regions are outside this novelty requirement.','',
      '| RNA candidate | Primary-nonself target peptides | Remain after broad-reference + personal-self exclusion | Decision |','|---|---:|---:|---|']
    for r in normal_table:lines.append(f"| {r['gene']} | {r['primary_nonself_target_peptides']} | {r['broad_and_personal_nonself_target_peptides']} | {'Exclude from altered-target selection' if r['gene']=='MUC6' else 'Retain as an RNA-sequence hypothesis'} |")
    lines+=['','MUC6’s entire selected 25-aa region matches **ENSP00000487059.1 / ENST00000627256.2**, on alternative contig **CHR_HSCHR11_2_CTG1**, annotated **polymorphic_pseudogene**. All 35 target-overlapping peptides match this broader reference; all 23 that were nonself against the primary coding subset therefore lose novelty under the expanded rule. Design #3 excludes this region.','',
      'This is a predicted alternative human translation, not a measurement that a healthy tissue expresses or displays it. The conservative policy excludes that exact sequence match even though the focal normal-DNA test supports the reference allele.','',
      'No additional target-peptide exclusions appeared in the other 23 altered-source regions relative to the primary-reference plus personal-self search. No negative search establishes absence from every normal protein: the reference release, predicted translations and bounded patient-normal models define the scope. No healthy-tissue immunopeptidome search or new HLA-binding inference was performed by this audit.','',
      '## Files and reruns','',
      '- [Summary table](summary.tsv) · [all altered-region decisions](summary.json).',
      '- [Normal-DNA genotypes and fresh normal/tumour/RNA base counts](normal-dna.json).',
      '- [RNA windows, codons and alignment paths](windows.json) · [reference witnesses](reference-matches.json).',
      '- [Reference/sample scope](reference-scope.json) · [audit request](request.json) · [candidate pool before exclusion](candidate-pool.json) · [execution receipt](complete.json).',
      '- [Runnable audit tools](../../../../tools/rna_evidence/README.md).','']
    if (base/'uniprot-reference-check.json').exists():
        from validate import verify_supplementary_reference
        u=verify_supplementary_reference(base)
        insert=lines.index('## Files and reruns')
        lines[insert:insert]=['## Second protein catalog: UniProt','',
          f"A supplementary exact search covered **{u['reference_sequences']:,} sequences** from UniProt's human reference proteome **UP000005640**, including canonical and additional/isoform FASTA files last modified 4 November 2026. All {u['reference_queries']:,} queried peptides were searched. It adds no target-peptide exclusions beyond the Ensembl/patient-normal result; the retained F5, MUC5B and MUC3A target sets are unchanged.",'',
          "MUC6's match is in the extended Ensembl inventory; neither its full window nor `SMAKTSTSL` matched this UniProt download. Its exclusion is a conservative response to a predicted alternative human translation, not a claim of measured healthy-tissue protein or a reviewed UniProt match.",'',
          '[Query-level decisions, provider URLs and file hashes](uniprot-reference-check.json) retain the complete search. Reference FASTA files come from [UniProt](https://www.uniprot.org/help/proteome), not from this repository.','']
    (base/'README.md').write_text('\n'.join(lines))
    result=dict(status='RNA_NORMAL_SELF_AUDIT_VERIFIED',focal_RNA_regions=4,altered_source_regions_searched=24,
        matched_normal_focal_changes_with_no_qualifying_ALT=4,excluded_region_ids=[by_gene['MUC6']['region_id']],
        fresh_local_reference_search_performed=new_reference_check,reference_queries=len(search),
        receipt_outputs_byte_verified=True,new_BAM_reading_in_this_verifier=False,normal_expression_or_presentation_proven=False)
    (base/'verification.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))


if __name__=='__main__':main()
