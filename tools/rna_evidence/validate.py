"""Fast integrity and exact-match-decision checks; no downloads or model inference."""
from pathlib import Path
import hashlib,json


def verify_supplementary_reference(directory):
    directory=Path(directory)
    result=json.loads((directory/'uniprot-reference-check.json').read_text())
    assert result['status']=='UNIPROT_REFERENCE_SEARCH_COMPLETED'
    assert result['source_audit_sha256']==hashlib.sha256((directory/'complete.json').read_bytes()).hexdigest()
    original=json.loads((directory/'reference-matches.json').read_text())
    assert set(original)==set(result['peptide_matches']) and len(original)==result['reference_queries']
    request=json.loads((directory/'request.json').read_text())
    by_id={r['region_id']:r for r in result['region_decisions']}
    assert len(by_id)==len(request['windows'])+len(request['additional_regions'])
    new_exclusions=[]
    for row in request['windows']+request['additional_regions']:
        lo,hi=row['target_start_in_window'],row['target_end_in_window']
        targets={row['sequence'][i:i+k] for k in range(8,12) for i in range(len(row['sequence'])-k+1) if i<hi and i+k>lo}
        previous={p for p in targets if original[p]['extended'] is None and not original[p]['personal']['occurs_in_searched_models']}
        remaining={p for p in previous if result['peptide_matches'][p] is None}
        decision=by_id[row['region_id']]
        assert decision['previous_nonself_target_peptides']==len(previous)
        assert decision['combined_nonself_target_peptides']==len(remaining)
        assert set(decision['new_reference_matching_target_peptides'])==previous-remaining
        if previous and not remaining:new_exclusions.append(row['region_id'])
    assert sorted(new_exclusions)==sorted(result['additional_excluded_region_ids'])
    # These case-specific selections have no extra matches; a changed reference
    # must trigger renewed scoring/selection rather than silently keep old scores.
    assert not any(r['new_reference_matching_target_peptides'] for r in result['region_decisions']), 'New normal-reference matches require renewed selection review'
    return dict(reference_sequences=result['reference_sequences'],reference_queries=result['reference_queries'],additional_excluded_regions=len(new_exclusions))
