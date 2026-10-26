"""Replay all algebraic coverage certificates and original-reply filters.

The verifier checks nilpotence by matrix powers rather than the generator's
characteristic polynomials. Both use SymPy exact arithmetic and the same
public incidence reduction; this is not an independent software backend.
"""
from pathlib import Path
from collections import Counter
import hashlib,json
from run_rur_study import verify_result
HERE=Path(__file__).resolve().parent
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def run():
    protocol=json.loads((HERE/'rur_protocol_v3.json').read_text())
    initial=json.loads((HERE/'rur_initial_v1/rur_protocol.json').read_text())
    assert [(x['id'],x['public_sha256']) for x in initial['cases']]==[(x['id'],x['public_sha256']) for x in protocol['cases']]
    for name,d in protocol['source_hashes'].items():assert sha(HERE/name)==d
    data=json.loads((HERE/'rur_results_v3.json').read_text())
    assert data['complete'] and data['protocol_sha256']==sha(HERE/'rur_protocol_v3.json')
    assert len(data['records'])==len(protocol['cases'])==37
    old={r['id']:r['result'] for r in json.loads((HERE/'theta_validation.json').read_text())['records']}
    counts=Counter();nonrational=0;arrays=0;degrees=[];tangent=None
    for case,row in zip(protocol['cases'],data['records']):
        assert case['id']==row['id'];assert row['status']=='complete',(row['id'],row['status'])
        path=HERE/'rur_public_v3'/(case['id']+'.json');assert sha(path)==case['public_sha256']
        public=json.loads(path.read_text());result=row['result'];verify_result(public,result,case['bounded'])
        if 'expected_count' in case:assert result['complete_count']==case['expected_count']
        if case['cohort']=='unchanged_theta27':
            assert {tuple(c['values']) for c in result['candidates']}=={tuple(c['values']) for c in old[case['id']]['candidates']}
        counts[result['status']]+=1;arrays+=result['complete_count']
        nonrational+=sum(not c.get('rational',True) for c in result.get('candidates',[]))
        if 'rur' in result and not result['rur']['empty']:
            degrees.append(result['rur']['quotient_dimension'])
        if case['id']=='theta3_c4_tangent':
            tangent={k:result['rur'][k] for k in ['quotient_dimension','distinct_complex_points']}
            assert tangent==dict(quotient_dimension=14,distinct_complex_points=13)
    assert arrays==41 and nonrational==10 and max(degrees)==14
    arithmetic=json.loads((HERE/'rur_arithmetic_validation.json').read_text())
    assert arithmetic['complete'] and len(arithmetic['controls'])==5 and len(arithmetic['rejected_corruptions'])==2
    assert arithmetic['solver_sha256']==sha(HERE/'exact_rur_decoder.py')
    result=dict(complete=True,jobs=37,distinct_public_inputs=34,old_theta_inputs=27,
                native_candidate_occurrences=arrays,nonrational_candidate_occurrences=nonrational,
                statuses=dict(counts),max_quotient_dimension=max(degrees),nonradical_native_control=tangent,
                arithmetic_controls=5,rejected_corruptions=2,initial_timeout_retained=True,
                public_hashes_unchanged_since_v1=True,protocol_sha256=sha(HERE/'rur_protocol_v3.json'),
                results_sha256=sha(HERE/'rur_results_v3.json'),verifier_sha256=sha(Path(__file__)),
                scope='same-machine exact certificate replay, including all real-root filters; no independent arithmetic backend or cross-solver timing claim')
    (HERE/'rur_validation.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result))
if __name__=='__main__':run()
