"""Public-input rational six-cycle preprocessing and exact capped list counts.

This experimental backend reports unsupported irrational roots explicitly.
The general paper algorithm has a separate algebraic-number implementation.
"""
from fractions import Fraction as F
from functools import reduce
from itertools import product
from math import gcd,lcm
import sympy as s

def preprocess(public):
    m,n=public['m'],public['n'];r=n//3
    assert m==2*r+1 and n==3*r
    u=list(map(F,public['u']));v=list(map(F,public['v']))
    gu=list(map(F,public['gu']));gv=list(map(F,public['gv']))
    edges=[tuple(e) for e in public['edges']]
    assert len(edges)==6*r and len(set(u))==len(set(v))==1
    rm=[(sum(u[i]*v[j]**2 for ii,j in edges if ii==i)-gu[i])/v[0] for i in range(m)]
    cm=[(sum(v[j]*u[i]**2 for i,jj in edges if jj==j)-gv[j])/u[0] for j in range(n)]
    t=s.Symbol('t');energies=[];branches=[];roots_total=rejects=0
    for k in range(r):
        b20=cm[3*k];b22=rm[2*k+2]-b20;b12=cm[3*k+2]-b22
        b11=rm[2*k+1]-b12;b01=cm[3*k+1]-b11
        base=[F(0),b01,b11,b12,b22,b20]
        forms=[s.Rational(z.numerator,z.denominator)+(-1)**i*t for i,z in enumerate(base)]
        poly=s.Poly(s.prod(forms[::2])-s.prod(forms[1::2]),t)
        roots=s.polys.polytools.ground_roots(poly);assert sum(roots.values())==poly.degree(), 'unsupported nonrational roots'
        local=[];ee=[]
        for root in sorted(roots):
            vals=[F(z) for z in [f.subs(t,root) for f in forms]];roots_total+=1
            local_edges=[(0,0),(0,1),(1,1),(1,2),(2,2),(2,0)]
            active_r={i for (i,j),z in zip(local_edges,vals) if z};active_c={j for (i,j),z in zip(local_edges,vals) if z}
            if any(not z and i in active_r and j in active_c for (i,j),z in zip(local_edges,vals)):
                rejects+=1;continue
            assert all(vals), 'unsupported zero branch after local filter'
            local.append(vals);ee.append(sum(z*z for z in vals))
        assert local;branches.append(local);energies.append(ee)
    norm=2*F(public['loss'])+sum(u[i]**2*v[j]**2 for i,j in edges)-sum(x*y for x,y in zip(u,gu))-sum(x*y for x,y in zip(v,gv))
    Q=reduce(lcm,[z.denominator for row in energies for z in row]+[norm.denominator],1)
    minima=[min(row) for row in energies]
    costs=[[int(Q*(z-low)) for z in row] for row,low in zip(energies,minima)]
    U=int(Q*(norm-sum(minima)));g=reduce(gcd,[z for row in costs for z in row],0) or 1
    if U%g: return dict(costs=costs,U=-1,branches=branches,roots=roots_total,rejected=rejects,lattice_obstruction=True)
    return dict(costs=[[z//g for z in row] for row in costs],U=U//g,branches=branches,roots=roots_total,rejected=rejects,lattice_obstruction=False)

def enumeration(costs,U):
    count=checks=0;witness=None
    for path in product(*(range(len(row)) for row in costs)):
        checks+=1
        if sum(costs[i][j] for i,j in enumerate(path))==U:
            count=min(2,count+1);witness=path
    return dict(count=count,witness=list(witness) if count==1 else None,operations=checks,storage=len(costs))

def partial(rows):
    ans={0:(1,[])};ops=0
    for row in rows:
        new={}
        for total,(count,path) in ans.items():
            for j,c in enumerate(row):
                ops+=1;z=total+c;prev=new.get(z,(0,None))
                new[z]=(min(2,prev[0]+count),path+[j])
        ans=new
    return ans,ops

def mitm(costs,U):
    h=len(costs)//2;left,ops1=partial(costs[:h]);right,ops2=partial(costs[h:])
    count=0;witness=None
    for val,(a,p) in left.items():
        if U-val in right:
            b,q=right[U-val];count=min(2,count+a*b);witness=p+q
    return dict(count=count,witness=witness if count==1 else None,
                operations=ops1+ops2+len(left),storage=len(left)+len(right))

def dynamic(costs,U):
    if U<0:return dict(count=0,witness=None,operations=0,storage=0)
    counts=[0]*(U+1);counts[0]=1;pred=[];ops=0
    for row in costs:
        new=[0]*(U+1);back=[None]*(U+1)
        for total in range(U+1):
            for j,c in enumerate(row):
                ops+=1
                if c<=total and counts[total-c]:
                    new[total]=min(2,new[total]+counts[total-c]);back[total]=(total-c,j)
        counts=new;pred.append(back)
    witness=None
    if counts[U]==1:
        witness=[];total=U
        for back in reversed(pred):total,j=back[total];witness.append(j)
        witness.reverse();assert total==0
    return dict(count=counts[U],witness=witness,operations=ops,storage=len(costs)*(U+1))

def verify_witness(public,branches,path):
    m,n=public['m'],public['n'];u=list(map(F,public['u']));v=list(map(F,public['v']))
    edges=[tuple(e) for e in public['edges']];C=[z for row,j in zip(branches,path) for z in row[j]]
    assert len(C)==len(edges)
    assert all(sum((u[i]*v[j]-z)*v[j] for (ii,j),z in zip(edges,C) if ii==i)==F(public['gu'][i]) for i in range(m))
    assert all(sum((u[i]*v[j]-z)*u[i] for (i,jj),z in zip(edges,C) if jj==j)==F(public['gv'][j]) for j in range(n))
    assert sum((u[i]*v[j]-z)**2 for (i,j),z in zip(edges,C))/2==F(public['loss'])
    aa=[F(1)];bb=[]
    for vals in [row[j] for row,j in zip(branches,path)]:
        z0,z1,z2,z3,z4,z5=vals;b0,b1=z0,z1;a1=z2/b1;b2=z3/a1;a2=z4/b2
        assert a2*b0==z5;aa.extend([a1,a2]);bb.extend([b0,b1,b2])
    assert sum(z*z for z in aa)*sum(z*z for z in bb)<1
    return True
