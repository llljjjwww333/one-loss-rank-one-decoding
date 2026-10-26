"""Exact rational reflection witnesses, plus exhaustive finite-grid counts."""
from fractions import Fraction as F
from itertools import product
from pathlib import Path
import hashlib,json
HERE=Path(__file__).resolve().parent
def reflect(a,q):
    c=2*sum(x*y for x,y in zip(a,q))/sum(x*x for x in q)
    return [c*y-x for x,y in zip(a,q)]
def check(r,N,coords):
    x=coords[:r+1];y=coords[r+1:]
    u=[F(1,(r+1)*z) for z in x];v=[F(1,2*r*z) for z in y]
    lam=F(1,2*r+2);a=[lam]+[-lam]*r;b=[lam,-lam]*r
    edges=[(i,j) for k in range(r) for i in [0,k+1] for j in [2*k,2*k+1]]
    B=[a[i]*b[j] for i,j in edges]
    rowevents=[k for k in range(r) if x[k+1]==x[0]]
    colevents=[k for k in range(r) if y[2*k]==y[2*k+1]]
    if not rowevents and not colevents:
        for k in range(r):
            # The four-cycle polynomial has degree one and nonzero slope.
            slope=lam**2*(u[0]-u[k+1])*(v[2*k]-v[2*k+1])
            assert slope
        return 0,0
    aa=a[:];bb=b[:]
    if rowevents:
        k=rowevents[0];bb[2*k:2*k+2]=reflect(b[2*k:2*k+2],v[2*k:2*k+2])
        assert sum(z*z for z in aa)<1 and sum(z*z for z in bb)<1
    else:
        k=colevents[0];ar=reflect([lam,-lam],[u[0],u[k+1]])
        assert ar[0]!=0
        aa[k+1]=lam*ar[1]/ar[0]
        bb[2*k:2*k+2]=[z*ar[0]/lam for z in b[2*k:2*k+2]]
    C=[aa[i]*bb[j] for i,j in edges];assert C!=B
    assert sum(z*z for z in C)==sum(z*z for z in B)
    for i in range(r+1):
        assert sum((c-z)*v[j] for (ii,j),c,z in zip(edges,C,B) if ii==i)==0
    for j in range(2*r):
        assert sum((c-z)*u[i] for (i,jj),c,z in zip(edges,C,B) if jj==j)==0
    return 1,int(bool(rowevents))
def run():
    protocol={'version':1,'exhaustive':[[1,n] for n in range(1,7)]+[[2,n] for n in range(1,4)]+[[3,2]],
              'analytic_r':[1,2,4,8,12,32,64],'analytic_N':[1,2,4,8,16,32,64,256,1024],
              'test':'rational reflection witnesses; nonzero linear cycle equations for all other grid points'}
    (HERE/'grid_protocol.json').write_text(json.dumps(protocol,indent=2)+'\n')
    rows=[]
    for r,N in protocol['exhaustive']:
        collisions=bounded=0
        for coords in product(range(1,N+1),repeat=3*r+1):
            c,b=check(r,N,coords);collisions+=c;bounded+=b
        total=N**(3*r+1)
        assert collisions==total-N**(r+1)*(N-1)**(2*r)
        assert F(bounded,total)==1-(1-F(1,N))**r
        rows.append(dict(r=r,N=N,queries=total,collisions=collisions,bounded_witnesses=bounded))
    out=dict(complete=True,exhaustive=rows,total_queries=sum(x['queries'] for x in rows),
             collision_witnesses=sum(x['collisions'] for x in rows),
             source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
             analytic=[dict(r=r,N=N,probability=str(1-(1-F(1,N))**(2*r)),bound=str(min(F(1),F(2*r,N))))
                       for r in protocol['analytic_r'] for N in protocol['analytic_N']])
    (HERE/'grid_validation.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({k:v for k,v in out.items() if k not in ['analytic','exhaustive']}))
if __name__=='__main__':run()
