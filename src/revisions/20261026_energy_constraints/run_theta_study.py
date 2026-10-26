"""Frozen public-input theta study; retain timeouts and unresolved replies."""
from pathlib import Path
from fractions import Fraction as Q
from collections import Counter
import hashlib,json,random,subprocess,sys,time
HERE=Path(__file__).resolve().parent


def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def mask(blocks):
    edges=[];m,n=1,0
    for p in blocks:
        if p=='C6':
            edges.extend([(0,n),(m,n),(m,n+1),(m+1,n+1),(m+1,n+2),(0,n+2)])
            m+=2;n+=3
        else:
            for i in range(p):edges.extend([(0,n+1+i),(m+i,n+1+i),(m+i,n)])
            m+=p;n+=p+1
    return m,n,edges


def data(blocks,kind,N,rng):
    m,n,edges=mask(blocks)
    if kind=='positive':a=[1+i%3 for i in range(m)];b=[1+j%4 for j in range(n)]
    elif kind=='mixed':a=[(-1)**i*(1+i%2) for i in range(m)];b=[(-1)**(j//2)*(1+j%3) for j in range(n)]
    else:a=[(-1)**(i//2)*(i*i+2) for i in range(m)];b=[(-1)**j*(j+2) for j in range(n)]
    a=[Q(z,4*sum(map(abs,a))) for z in a];b=[Q(z,4*sum(map(abs,b))) for z in b]
    u=[Q(1,m*rng.randint(1,N)) for _ in range(m)];v=[Q(1,n*rng.randint(1,N)) for _ in range(n)]
    B=[a[i]*b[j] for i,j in edges];gu=[Q(0)]*m;gv=[Q(0)]*n;loss=Q(0)
    for (i,j),value in zip(edges,B):
        z=u[i]*v[j]-value;gu[i]+=z*v[j];gv[j]+=z*u[i];loss+=z*z/2
    public=dict(m=m,n=n,edges=edges,u=list(map(str,u)),v=list(map(str,v)),
                gu=list(map(str,gu)),gv=list(map(str,gv)),loss=str(loss))
    return public,B


def run():
    protocol_path=HERE/'theta_protocol.json';protocol=json.loads(protocol_path.read_text())
    pubdir=HERE/'theta_public';pubdir.mkdir(exist_ok=True)
    dest=HERE/('theta_replay_outputs' if '--replay' in sys.argv else 'theta_checkpoint_outputs')
    dest.mkdir(exist_ok=False)
    rng=random.Random(protocol['seed']);cases=[]
    for p in [3,4]:
        for kind in ['positive','mixed','asymmetric']:
            for N in protocol['query_grids']:cases.append(dict(id=f'theta{p}_{kind}_N{N}',blocks=[p],kind=kind,N=N))
    for blocks in [[3,'C6'],[3,3],[3,4]]:
        cases.append(dict(id='glued_'+'_'.join(map(str,blocks)),blocks=blocks,kind='asymmetric',N=8))
    assert len(cases)==27
    inputs=[]
    for case in cases:
        public,B=data(case['blocks'],case['kind'],case['N'],rng)
        path=pubdir/(case['id']+'.json');path.write_text(json.dumps(public,indent=2)+'\n')
        inputs.append((case,path,B))
    frozen=dict(protocol_sha256=sha(protocol_path),cases=[dict(case,public_sha256=sha(path)) for case,path,_ in inputs],
                source_hashes={p.name:sha(p) for p in [Path(__file__),HERE/'normalized_block_decoder.py']},
                truth_passed_to_solver=False)
    (HERE/'theta_frozen_inputs.json').write_text(json.dumps(frozen,indent=2)+'\n')
    records=[]
    def verify_result(result,B,case,path):
        assert result['public_sha256']==sha(path)
        assert result['status']!='Empty',case['id']
        if result['complete']:
            assert any(list(map(Q,c['values']))==B and c['unit_factor_budget_passed'] for c in result['candidates']),case['id']
    for case,path,B in inputs:
        output=dest/(case['id']+'.json');start=time.perf_counter()
        try:
            child=subprocess.run([sys.executable,'-X','utf8',str(Path(__file__)),'--worker',str(path),str(output),
                                  str(int(len(case['blocks'])==1))],capture_output=True,text=True,encoding='utf-8',
                                 timeout=protocol['timeout_seconds'])
            record=dict(case,exit_code=child.returncode)
            if child.returncode:
                record.update(status='error',stderr=child.stderr,stdout=child.stdout)
            else:
                result=json.loads(output.read_text());verify_result(result,B,case,path)
                record.update(status='complete' if result['complete'] else 'unresolved',result=result,
                              hidden_truth_comparison=result['complete'])
        except subprocess.TimeoutExpired as exc:
            record=dict(case,status='timeout',timeout_seconds=protocol['timeout_seconds'])
            (dest/(case['id']+'.timeout.txt')).write_text('Prespecified 30-second timeout; retained.\n')
            if output.exists():
                result=json.loads(output.read_text());verify_result(result,B,case,path)
                record.update(result=result,norm_certificate_recovered=result['complete'],
                              hidden_truth_comparison=result['complete'],gradient_audit_timeout=True)
        record['elapsed_seconds']=time.perf_counter()-start;records.append(record)
        (HERE/'theta_progress.json').write_text(json.dumps(dict(completed=len(records),total=len(inputs),
            statuses=dict(Counter(x['status'] for x in records))),indent=2)+'\n')
        print(case['id'],record['status'],record.get('result',{}).get('gradient',{}).get('real_actual_candidates'),flush=True)
    result=dict(complete=True,scope='all 27 prescribed exact public inputs, with all unresolved cases and timeouts retained',
                cases=len(records),statuses=dict(Counter(x['status'] for x in records)),records=records,
                norm_statuses=dict(Counter(x.get('result',{}).get('status','unresolved') for x in records)),
                protocol_sha256=sha(protocol_path),frozen_inputs_sha256=sha(HERE/'theta_frozen_inputs.json'),
                source_hashes=frozen['source_hashes'],truth_passed_to_solver=False)
    (HERE/'theta_validation.json').write_text(json.dumps(result,indent=2)+'\n')
    print(result['statuses'],flush=True)


if __name__=='__main__':
    if len(sys.argv)>1 and sys.argv[1]=='--worker':
        from normalized_block_decoder import decode,gradient_count,public_system
        public=json.loads(Path(sys.argv[2]).read_text())
        audit=bool(int(sys.argv[4]));result=decode(public,with_gradient_count=False)
        result.update(public_sha256=sha(Path(sys.argv[2])),gradient_phase='pending' if audit else 'not_requested')
        Path(sys.argv[3]).write_text(json.dumps(result,indent=2)+'\n')
        if audit:
            result['gradient']=gradient_count(public_system(public));result['gradient_phase']='complete'
            Path(sys.argv[3]).write_text(json.dumps(result,indent=2)+'\n')
    else:run()
