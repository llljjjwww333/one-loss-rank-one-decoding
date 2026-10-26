"""Exact parsimonious SUBSET-SUM instances for the native query interface."""
from fractions import Fraction as Q
from itertools import product


def cycle_parameters(item):
    w=2*int(item)
    assert w>=2
    r,s=Q(w+1,2),Q(w-1,2)
    return r,s,r*s-r-s,r*s+r+s


def cycle_values(item,high=False):
    r,s,d,e=cycle_parameters(item);t=r if high else s
    return [t,1-t,1+t,-t,d+t,e-t]


def factors(items,bits):
    a=[Q(1)];b=[]
    for item,bit in zip(items,bits):
        r,s,_,_=cycle_parameters(item);t=r if bit else s;other=s if bit else r
        a.extend([(1+t)/(1-t),other*(1+t)/t])
        b.extend([t,1-t,t*(t-1)/(t+1)])
    scale=Q(1,10*(1+2*sum(items)))
    return [x*scale for x in a],[x*scale for x in b]


def construction(items,target):
    items=list(map(int,items));assert items and min(items)>0 and target>=0
    count=len(items);m,n=2*count+1,3*count
    edges=[];low=[];scale=Q(1,10*(1+2*sum(items)))
    for k,item in enumerate(items):
        rows=[0,2*k+1,2*k+2];cols=[3*k,3*k+1,3*k+2]
        edges.extend((rows[i],cols[j]) for i,j in [(0,0),(0,1),(1,1),(1,2),(2,2),(2,0)])
        low.extend(x*scale**2 for x in cycle_values(item))
    u=[Q(1,m)]*m;v=[Q(1,n)]*n
    gu=[Q(0)]*m;gv=[Q(0)]*n
    for (i,j),x in zip(edges,low):
        residual=u[i]*v[j]-x;gu[i]+=residual*v[j];gv[j]+=residual*u[i]
    baseline=sum((u[i]*v[j]-x)**2 for (i,j),x in zip(edges,low))/2
    # Replacing s by r changes squared array norm by 4*item*scale**4.
    loss=baseline+2*target*scale**4
    public=dict(m=m,n=n,edges=edges,u=list(map(str,u)),v=list(map(str,v)),
        gu=list(map(str,gu)),gv=list(map(str,gv)),loss=str(loss))
    return public,dict(scale=str(scale),baseline_norm=str(sum(x*x for x in low)),
        required_norm=str(sum(x*x for x in low)+4*target*scale**4),
        item_weights=items,target=int(target))


def solution_bits(items,target):
    return [bits for bits in product([0,1],repeat=len(items)) if sum(a*b for a,b in zip(items,bits))==target]
