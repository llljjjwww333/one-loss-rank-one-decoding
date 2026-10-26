"""Frozen public-input study of exact algebraic reply decoding.

Per-job 90-second budgets; no solver speed comparison. Old 27-input
cohort unchanged; new native inputs are constructed from rational margins
and norms, never by passing a hidden target to the decoder.
"""
from pathlib import Path
from fractions import Fraction as Q
import hashlib,json,subprocess,sys,time
from exact_rur_decoder import decode,verify_rur,public_system,polynomial,Z,reply_identities
from exact_cactus_decoder import actual_budget
import sympy as s
HERE=Path(__file__).resolve().parent
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def put(p,x):p.write_text(json.dumps(x,indent=2)+'\n')
def public(m,n,edges,values,extra_norm=Q(0)):
    u=[Q(1,m)]*m;v=[Q(1,n)]*n;gu=[Q(0)]*m;gv=[Q(0)]*n;loss=extra_norm/2
    for (i,j),value in zip(edges,values):
        residual=u[i]*v[j]-value;gu[i]+=residual*v[j];gv[j]+=residual*u[i];loss+=residual**2/2
    return dict(m=m,n=n,edges=edges,u=list(map(str,u)),v=list(map(str,v)),
                gu=list(map(str,gu)),gv=list(map(str,gv)),loss=str(loss))
def source_hashes():
    return {name:sha(HERE/name) for name in ['exact_rur_decoder.py','normalized_block_decoder.py','run_rur_study.py']}

def prepare():
    destination=HERE/'rur_public_v3';destination.mkdir(exist_ok=False);cases=[]
    old=json.loads((HERE/'theta_frozen_inputs.json').read_text())['cases']
    for case in old:
        source=HERE/'theta_public'/(case['id']+'.json');target=destination/(case['id']+'.json')
        assert sha(source)==case['public_sha256'];target.write_bytes(source.read_bytes())
        cases.append(dict(id=case['id'],bounded=False,public_sha256=sha(target),cohort='unchanged_theta27'))
    c4=[(0,0),(0,1),(1,1),(1,0)]
    theta=[edge for i in range(1,4) for edge in [(0,i),(i,i),(i,0)]]
    attached=[(0,4),(0,5),(4,5),(4,4)];scale=Q(1,64)
    new=[('c4_irrational',public(2,2,c4,[Q(0)]*4,Q(1,32)),2),
         ('c4_tangent',public(2,2,c4,[Q(0)]*4),1),
         ('theta3_c4_irrational',public(5,6,theta+attached,[value*scale for _ in range(3) for value in [3,1,3]]+[Q(0)]*4,2*scale**2),2),
         ('theta3_c4_tangent',public(5,6,theta+attached,[value*scale for _ in range(3) for value in [3,1,3]]+[Q(0)]*4),1),
         ('c4_irrational_outside_budget',public(2,2,c4,[Q(0)]*4,Q(8)),2),
         ('nonactual_forest',public(2,2,[(0,0),(0,1),(1,1)],[Q(1),Q(0),Q(1)]),0)]
    inconsistent=public(2,2,[(0,0)],[Q(0)]);inconsistent['gu'][1]='1'
    new.append(('inconsistent_margins',inconsistent,0))
    for name,p,expected in new:
        target=destination/(name+'.json');put(target,p)
        cases.append(dict(id=name,bounded=False,public_sha256=sha(target),cohort='new_native_algebraic',expected_count=expected))
    for base,expected in [('c4_irrational',2),('c4_irrational_outside_budget',0),('theta3_c4_irrational',2)]:
        target=destination/(base+'_bounded.json');target.write_bytes((destination/(base+'.json')).read_bytes())
        cases.append(dict(id=base+'_bounded',bounded=True,public_sha256=sha(target),cohort='optional_budget',expected_count=expected))
    assert len(cases)==37
    put(HERE/'rur_protocol_v3.json',dict(version=3,complete=False,per_worker_seconds=90,max_quotient_dimension=64,
        platform='Windows; exact reference certificates, no timing ranking',cases=cases,source_hashes=source_hashes(),
        truth_passed_to_solver=False,old_native_timings_repeated=False))
    print('Frozen 37 jobs on 34 distinct public inputs, including 27 unchanged theta inputs')

def verify_result(p,result,bounded):
    assert result['complete'];system=public_system(p)
    if 'rur' not in result:
        check=decode(p,bounded,64);assert check==result;return
    certificate=result['rur'];equations=system['equations']+[system['energy']]
    verify_rur(equations,system['variables'],certificate,64)
    if certificate['empty']:assert result['complete_count']==0;return
    P=polynomial(certificate['polynomial']);coords=list(map(polynomial,certificate['coordinates']))
    substitution={v:w.as_expr() for v,w in zip(system['variables'],coords)}
    edges=[s.Poly(f.subs(substitution,simultaneous=True),Z,domain=s.QQ).rem(P) for f in system['values']]
    assert [list(map(str,p.all_coeffs())) for p in edges]==result['edge_polynomials']
    reply_identities(system,edges,P);expected=[]
    for index,root in enumerate(P.real_roots(radicals=False)):
        values=[v.as_expr().subs(Z,root) for v in edges]
        accepted,_=actual_budget(p['m'],p['n'],p['edges'],values,bounded)
        if accepted:expected.append((index,list(map(str,values))))
    assert [(x['real_root_ordinal'],x['values']) for x in result['candidates']]==expected
    assert result['complete_count']==len(expected)

def worker(case_id):
    protocol=json.loads((HERE/'rur_protocol_v3.json').read_text());case=next(x for x in protocol['cases'] if x['id']==case_id)
    path=HERE/'rur_public_v3'/(case_id+'.json');assert sha(path)==case['public_sha256']
    p=json.loads(path.read_text());result=decode(p,case['bounded'],protocol['max_quotient_dimension'])
    if result['complete']:verify_result(p,result,case['bounded'])
    put(HERE/'rur_outputs_v3'/(case_id+'.json'),result)

def run():
    protocol=json.loads((HERE/'rur_protocol_v3.json').read_text());assert protocol['source_hashes']==source_hashes()
    (HERE/'rur_outputs_v3').mkdir(exist_ok=False);records=[]
    for case in protocol['cases']:
        began=time.perf_counter();path=HERE/'rur_public_v3'/(case['id']+'.json');assert sha(path)==case['public_sha256']
        try:
            r=subprocess.run([sys.executable,'-X','utf8',str(Path(__file__)),'worker',case['id']],capture_output=True,text=True,encoding='utf-8',timeout=protocol['per_worker_seconds'])
            if r.returncode:row=dict(id=case['id'],status='error',stdout=r.stdout,stderr=r.stderr)
            else:
                result=json.loads((HERE/'rur_outputs_v3'/(case['id']+'.json')).read_text())
                row=dict(id=case['id'],status='complete' if result['complete'] else 'unresolved',result=result)
        except subprocess.TimeoutExpired:row=dict(id=case['id'],status='timeout',budget_seconds=protocol['per_worker_seconds'])
        row['wall_seconds_including_startup_and_verifier']=time.perf_counter()-began;records.append(row)
        print(case['id'],row['status'],row.get('result',{}).get('complete_count'),flush=True)
        put(HERE/'rur_results_v3.json',dict(complete=len(records)==len(protocol['cases']),protocol_sha256=sha(HERE/'rur_protocol_v3.json'),records=records))

if __name__=='__main__':
    if sys.argv[1]=='prepare':prepare()
    elif sys.argv[1]=='worker':worker(sys.argv[2])
    else:run()
