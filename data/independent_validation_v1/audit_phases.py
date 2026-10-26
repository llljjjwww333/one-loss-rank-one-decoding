"""Audit saved phase measurements and unchanged outcomes without repeating timings."""
from pathlib import Path
from collections import Counter
from math import isfinite,prod
import hashlib,json
HERE=Path(__file__).resolve().parent;ORIGINAL=HERE.parent/'energy_constraints_v1'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def run():
    protocol=json.loads((HERE/'phase_protocol.json').read_text());raw=json.loads((HERE/'phase_results.json').read_text())
    assert sha(HERE/'phase_study.py')==protocol['runner_sha256']
    assert sha(ORIGINAL/'coordinate_decoder.py')==protocol['decoder_sha256']
    assert sha(ORIGINAL/'prospective_inputs.jsonl')==protocol['inputs_sha256']
    assert sha(ORIGINAL/'prospective_outcomes.jsonl')==protocol['reference_sha256']
    assert raw['complete'] and len(raw['records'])==288 and raw['protocol_sha256']==sha(HERE/'phase_protocol.json')
    refs={r['id']:r for r in map(json.loads,(ORIGINAL/'prospective_outcomes.jsonl').read_text().splitlines())};assert len(refs)==96
    seen=set();counts=Counter()
    for r in raw['records']:
        key=r['id'],r['repeat'];assert key not in seen and r['repeat'] in [1,2,3];seen.add(key)
        expected=refs[r['id']];assert r['r']==expected['r'] and r['N']==expected['N']
        assert r['outcome']=='complete' and r['result']==expected['result']
        assert isfinite(r['seconds']) and r['seconds']>=sum(r['phases'].values())
        assert all(isfinite(x) and x>=0 for x in r['phases'].values())
        assert r['calls']['reduction']==r['calls']['coordinates']==r['calls']['propagate']==1
        assert r['calls']['real_roots']==r['r']
        d=r['result'];assert not d['truth_passed_to_solver'] and d['cartesian_tuples_enumerated']==0
        assert d['cartesian_before']==prod(len(c['roots']) for c in d['cycles'])
        assert d['cartesian_after']==prod(map(len,d['domains'])) and d['deletions']<=3*r['r']
        counts[d['status']]+=1
    assert counts=={'CertifiedUnique':252,'Unresolved':36}
    out=dict(complete=True,source_sha256=sha(Path(__file__)),records=288,distinct_inputs=96,
      outcomes=dict(counts),phase_results_sha256=sha(HERE/'phase_results.json'),timings_repeated=False,
      scope='Saved timing records and exact scientific outcomes audited; new timing values are machine dependent.')
    (HERE/'phase_validation.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out))
if __name__=='__main__':run()
