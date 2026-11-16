"""Supplement the frozen Ensembl/personal-self audit with a provider human reference proteome."""
from pathlib import Path
import argparse,hashlib,json
from worker import translation_records,library,exact_hit
ROOT=Path(__file__).resolve().parents[2]


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--reference-proteins',type=Path,required=True)
    p.add_argument('--download-receipt',type=Path,required=True)
    p.add_argument('--results',type=Path,default=Path('designs/design-3/inputs/sequence-review'))
    a=p.parse_args();base=(ROOT/a.results).resolve();meta=json.loads(a.download_receipt.read_text())
    meta=meta.get('reference_download',meta)
    assert hashlib.sha256(a.reference_proteins.read_bytes()).hexdigest()==meta['sha256']
    assert a.reference_proteins.stat().st_size==meta['bytes']
    records=list(translation_records(a.reference_proteins));lib=library(records)
    assert len(records)>20000
    original=json.loads((base/'reference-matches.json').read_text());matches={p:exact_hit(lib,p) for p in original}
    request=json.loads((base/'request.json').read_text())
    decisions=[]
    for r in request['windows']+request['additional_regions']:
        lo,hi=r['target_start_in_window'],r['target_end_in_window']
        target={r['sequence'][i:i+k] for k in range(8,12) for i in range(len(r['sequence'])-k+1) if i<hi and i+k>lo}
        previous={p for p in target if original[p]['extended'] is None and not original[p]['personal']['occurs_in_searched_models']}
        retained={p for p in previous if matches[p] is None}
        decisions.append(dict(gene=r['gene'],region_id=r['region_id'],previous_nonself_target_peptides=len(previous),
          combined_nonself_target_peptides=len(retained),new_reference_matching_target_peptides=sorted(previous-retained),
          full_window_reference_match=exact_hit(lib,r['sequence'])))
    headers={k.lower():v for k,v in meta['headers'].items()}
    result=dict(status='UNIPROT_REFERENCE_SEARCH_COMPLETED',provider='UniProt',reference_proteome='UP000005640',include_canonical_and_isoform_sequences=True,
      release=headers.get('x-uniprot-release'),release_date=headers.get('x-uniprot-release-date'),reference_download=meta,
      reference_sequences=len(records),reference_queries=len(matches),all_queries_exactly_searched=True,peptide_matches=matches,region_decisions=decisions,
      additional_excluded_region_ids=[r['region_id'] for r in decisions if r['previous_nonself_target_peptides'] and not r['combined_nonself_target_peptides']],
      source_audit_sha256=hashlib.sha256((base/'complete.json').read_bytes()).hexdigest(),
      scope='Exact sequence matches in the downloaded human reference proteome with isoforms; not exhaustive human variation or proof of healthy-tissue expression/presentation',
      provider_help='https://www.uniprot.org/help/proteome',new_model_inference=False)
    (base/'uniprot-reference-check.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ['peptide_matches','reference_download','region_decisions']},indent=2))
    for r in decisions:
        if r['new_reference_matching_target_peptides']:print('NEW_MATCHES',r['gene'],r['new_reference_matching_target_peptides'])


if __name__=='__main__':main()
