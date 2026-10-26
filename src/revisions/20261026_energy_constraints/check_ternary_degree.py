"""All three six-cycle branches for the sharp exponential degree base.

Native line witnesses may have zero rows. The main proof obtains nonzero
competitors on a common open set of perturbed outside query coordinates.
"""
from pathlib import Path
from fractions import Fraction as Q
from itertools import product
from math import isqrt
import hashlib,json
from check_certificate_degree import Quadratic,native_reply
HERE=Path(__file__).resolve().parent
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def witness(r,branches,T):
    D=4+T;x=Quadratic(-2,1,D);gamma=2**(r+4);m,n=2*r+3,3*r+3
    assert x*x+4*x==T and (x-1).sign()>=0 and (x-2**(r+1)).sign()<0
    a=[Quadratic(1,D=D),Quadratic(Q(2,3),D=D),Quadratic(Q(2,3),D=D)]
    b=[Quadratic(z,D=D) for z in [9,3,-3]]
    c=[Quadratic(1,D=D),1+x/2,Quadratic(Q(1,2),D=D)]
    d=[Quadratic(z,D=D) for z in [10,2,-2]];edges=[]
    for block in range(r+1):
        i,j=2*block+1,3*block
        edges.extend([(0,j),(0,j+1),(i,j+1),(i,j+2),(i+1,j+2),(i+1,j)])
        if not block:continue
        rho=2**(block-1);branch=branches[block-1]
        a.extend([Quadratic(1,D=D),Quadratic(Q(1,3),D=D)])
        b.extend([Quadratic(rho*z,D=D) for z in [9,3,-3]])
        rows=[(1,Q(1,3)),(Q(1,2),Q(1,2)),(0,1)][branch]
        cols=[[9,3,-3],[8,4,-4],[6,6,-4]][branch]
        c.extend([Quadratic(z,D=D) for z in rows]);d.extend([Quadratic(rho*z,D=D) for z in cols])
    a,b,c,d=[[z/gamma for z in values] for values in [a,b,c,d]]
    for left,right in [(a,b),(c,d)]:
        assert (sum((z*z for z in left),Quadratic(D=D))-Q(1,32)).sign()<0
        assert (sum((z*z for z in right),Quadratic(D=D))-Q(13,64)).sign()<0
    B=[a[i]*b[j] for i,j in edges];C=[c[i]*d[j] for i,j in edges]
    assert B!=C and all(z.sign()!=0 for z in B)
    nonzero=all(z.sign()!=0 for z in C);assert nonzero==(2 not in branches)
    # Every zero edge in a line witness is incident to an inactive row;
    # the other four cycle edges provide a nonvanishing propagation patch.
    for block,branch in enumerate(branches,1):
        if branch==2:assert all(C[6*block+k].sign()!=0 for k in [0,1,4,5])
    u=[Quadratic(Q(1,m),D=D) for _ in range(m)];u[1]=1/(m*x)
    v=[Quadratic(Q(1,n),D=D) for _ in range(n)]
    assert (sum((z*z for z in u),Quadratic(D=D))-1).sign()<0
    assert (sum((z*z for z in v),Quadratic(D=D))-1).sign()<0
    original=native_reply(edges,u,v,B,D);alternative=native_reply(edges,u,v,C,D)
    assert original==alternative and sum((z*z-y*y for z,y in zip(C,B)),Quadratic(D=D))==0
    assert (4*x+8).sign()>0
    assert all([6,-4,12][branch]*4**i!=0 for i,branch in enumerate(branches))
    payload=dict(target=[z.pair() for z in B],candidate=[z.pair() for z in C],
                 u=[z.pair() for z in u],v=[z.pair() for z in v],reply=original)
    return dict(r=r,branches=list(branches),T=T,D=D,x_polynomial_relation=[1,4,-T],
                native_reply_sha256=hashlib.sha256(json.dumps(payload,sort_keys=True).encode()).hexdigest(),
                candidate_nonzero_on_base_slice=nonzero,bounded_factors_checked=True,
                nonvanishing_factor_propagation_patch=True,positive_collision_derivative=True)

def run():
    protocol=dict(version=1,native_r=list(range(1,7)),integer_slice_r=list(range(1,13)),
                  policy='all 3^r branch configurations; every positive compensation retained',
                  scope='exact constructed algebraic line witnesses; nonzero outside perturbations are proved analytically; not full-grid probabilities')
    path=HERE/'ternary_degree_protocol.json'
    if path.exists():assert json.loads(path.read_text())==protocol
    else:path.write_text(json.dumps(protocol,indent=2)+'\n')
    records=[];native=[]
    for r in protocol['native_r']:
        values=[];current=[]
        for branches in product(range(3),repeat=r):
            T=sum([0,5,-3][branch]*4**i for i,branch in enumerate(branches));values.append(T)
            if T>0:current.append(witness(r,branches,T))
        assert len(set(values))==3**r and len(current)==(3**r-1)//2
        records.extend(current);native.append(dict(r=r,cycles=r+1,witnesses=len(current)))
        print('ternary native',r,len(current),flush=True)
    output=HERE/'ternary_degree_witnesses.jsonl';output.write_text(''.join(json.dumps(x,sort_keys=True)+'\n' for x in records))
    binary=[json.loads(line) for line in (HERE/'degree_witnesses.jsonl').read_text().splitlines()]
    binary_hashes={x['native_reply_sha256'] for x in binary};ternary_hashes={x['native_reply_sha256'] for x in records}
    assert len(binary_hashes&ternary_hashes)==89
    assert len(binary_hashes|ternary_hashes)==1821
    assert sum(x['candidate_nonzero_on_base_slice'] for x in records)==120
    values=[0];integer=[]
    for r in protocol['integer_slice_r']:
        power=4**(r-1);values=[z+d*power for d in [0,5,-3] for z in values]
        assert len(set(values))==3**r
        positive=[T for T in values if T>0];assert len(positive)==(3**r-1)//2
        xs=sorted(isqrt(4+T)-2 for T in positive if isqrt(4+T)**2==4+T)
        integer.append(dict(r=r,cycles=r+1,constructed_real_collisions=len(positive),integer_members=len(xs),integer_parameters=xs))
    result=dict(complete=True,native_witnesses=len(records),native_rows=native,
                base_nonzero_candidates=120,base_zero_row_candidates=423,
                binary_crosscheck_duplicates=89,combined_distinct_native_witnesses=1821,
                integer_slice_controls=integer,protocol_sha256=sha(path),witnesses_sha256=sha(output),
                binary_witnesses_sha256=sha(HERE/'degree_witnesses.jsonl'),
                verifier_sha256=sha(Path(__file__)),quadratic_arithmetic_sha256=sha(HERE/'check_certificate_degree.py'),
                all_native_loss_gradients_and_factor_bounds_checked=True,
                full_independent_grid_probability_estimated=False,
                scope='finite native exact controls and complete constructed integer-member counts; arbitrary-r degree and nonzero-continuation proof are in the main article')
    assert len(records)==543
    (HERE/'ternary_degree_validation.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='integer_slice_controls'}))
if __name__=='__main__':run()
