"""Exact public-input counts for at most two cactus cycles of degree <=3.

QQ factorization and quadratic radicals suffice for the prespecified rational
targets. The true target is never supplied to this solver. Unsupported factors
and zero cycle polynomials are returned as unresolved, never counted as unique.
There is no factor-budget restriction: counts are among all actual real arrays.
"""
from fractions import Fraction as Q
from itertools import product
from math import isqrt
import sympy as s
from exact_cactus_decoder import canonical_public, cycle_basis, forest_base


def sign(x): return (x > 0)-(x < 0)


class Quad:
    """a+b sqrt(d), with positive nonsquare rational d when b != 0."""
    def __init__(self,a=0,b=0,d=0):
        self.a,self.b,self.d=map(Q,(a,b,d))
        if self.b:
            assert self.d>0
            p,q=isqrt(self.d.numerator),isqrt(self.d.denominator)
            if p*p==self.d.numerator and q*q==self.d.denominator:
                self.a+=self.b*Q(p,q); self.b=Q(0)
        if not self.b:self.d=Q(0)

    def __add__(self,other):
        other=other if isinstance(other,Quad) else Quad(other)
        assert not self.b or not other.b or self.d==other.d
        return Quad(self.a+other.a,self.b+other.b,self.d or other.d)
    __radd__=__add__
    def __mul__(self,other):
        other=other if isinstance(other,Quad) else Quad(other)
        assert not self.b or not other.b or self.d==other.d
        d=self.d or other.d
        return Quad(self.a*other.a+self.b*other.b*d,self.a*other.b+self.b*other.a,d)
    __rmul__=__mul__
    def zero(self):return not self.a and not self.b
    def encode(self):return [str(self.a),str(self.b),str(self.d)]
    def expr(self):return s.Rational(self.a)+s.Rational(self.b)*s.sqrt(s.Rational(self.d))


def radical_sign(a,b,d):
    if not b:return sign(a)
    if not a:return sign(b)
    if sign(a)==sign(b):return sign(a)
    return sign(a)*sign(a*a-b*b*d)


def energy_zero(parts,target):
    """Exact zero test for a sum of at most two quadratic-field elements."""
    c=sum((p.a for p in parts),Q(0))-target
    radicals=[(p.b,p.d) for p in parts if p.b]
    if not radicals:return not c
    if len(radicals)==1:return not c and not radicals[0][0]
    (b,d),(e,f)=radicals
    if not c:return sign(b)==-sign(e) and b*b*d==e*e*f
    u=c*c+b*b*d-e*e*f
    return (u*u==4*c*c*b*b*d and sign(-u)==sign(c*b)
            and radical_sign(c,b,d)==-sign(e))


def real_roots(poly):
    roots=[]
    for factor,multiplicity in s.factor_list(poly)[1]:
        c=list(map(Q,factor.all_coeffs()))
        if factor.degree()==1: roots.append(Quad(-c[1]/c[0]))
        elif factor.degree()==2:
            a,b,z=c; delta=b*b-4*a*z
            if delta<0:continue
            if delta==0:roots.append(Quad(-b/(2*a)))
            else:roots.extend([Quad(-b/(2*a),sgn/(2*a),delta) for sgn in (-1,1)])
        else:raise ValueError('unsupported irreducible factor degree')
    assert len({tuple(root.encode()) for root in roots})==len(roots)
    return roots


def decode(public):
    p=canonical_public(public); m,n,edges=p['m'],p['n'],p['edges']
    u,v,gu,gv=[p[key] for key in ('u','v','gu','gv')]
    cycles=cycle_basis(m,n,edges)
    assert len(cycles)<=2
    weights=[u[i]*v[j] for i,j in edges]
    pu=[sum(u[i]*v[j]**2 for a,j in edges if a==i) for i in range(m)]
    pv=[sum(v[j]*u[i]**2 for i,b in edges if b==j) for j in range(n)]
    rhs=[u[i]*(pu[i]-gu[i]) for i in range(m)]+[v[j]*(pv[j]-gv[j]) for j in range(n)]
    base=forest_base(m,n,edges,cycles,rhs)
    if base is None:return dict(complete=False,reason='inconsistent margins')
    total=Q(2*p['loss']+sum(w*w for w in weights)-sum(a*b for a,b in zip(u,gu))-sum(a*b for a,b in zip(v,gv)))
    t=s.Symbol('t'); branch_sets=[]; records=[]
    for cycle in cycles:
        z=[base[e]+(-1)**j*t for j,e in enumerate(cycle)]
        poly=s.Poly(s.prod(z[::2])-s.prod(z[1::2]),t,domain=s.QQ)
        if poly.is_zero:return dict(complete=False,reason='identically zero cycle')
        if poly.degree()>3:return dict(complete=False,reason='cycle degree above three')
        try:roots=real_roots(poly)
        except ValueError as exc:return dict(complete=False,reason=str(exc))
        branches=[]
        for root in roots:
            # Exact substitution checks every enumerated polynomial root.
            val=Quad()
            for coef in poly.all_coeffs():val=val*root+Q(coef)
            assert val.zero()
            values=[Quad(Q(base[e]/weights[e]))+root*Q((-1)**j/weights[e]) for j,e in enumerate(cycle)]
            energy=sum((c*c for c in values),Quad())
            branches.append(dict(root=root,values=values,energy=energy))
        branch_sets.append(branches)
        records.append(dict(edges=cycle,coefficients=list(map(str,poly.all_coeffs())),
                            roots=[root.encode() for root in roots]))
    cycle_edges={e for cycle in cycles for e in cycle}
    fixed={e:Quad(Q(base[e]/weights[e])) for e in range(len(edges)) if e not in cycle_edges}
    fixed_energy=sum((x.a*x.a for x in fixed.values()),Q(0))
    decisions=[];full=[];gradient=[]
    for choice in product(*(range(len(branches)) for branches in branch_sets)):
        values=dict(fixed);parts=[]
        for cycle,branches,k in zip(cycles,branch_sets,choice):
            branch=branches[k]; parts.append(branch['energy'])
            values.update(zip(cycle,branch['values']))
        active_vertices=set()
        for e,(i,j) in enumerate(edges):
            if not values[e].zero():active_vertices.update((i,m+j))
        actual=not any(values[e].zero() and i in active_vertices and m+j in active_vertices for e,(i,j) in enumerate(edges))
        # On a cactus, the root equations check every all-nonzero cycle;
        # deleting any zero edge breaks that cycle. This support test is
        # precisely the remaining finite-factor condition of Lemma 3.1.
        match=energy_zero(parts,total-fixed_energy)
        decisions.append(dict(choice=choice,actual=actual,energy_match=match))
        if actual:
            candidate=[values[e].encode() for e in range(len(edges))]
            gradient.append(candidate)
            if match:full.append(candidate)
    return dict(complete=True,gradient_count=len(gradient),complete_count=len(full),
                cycles=records,decisions=decisions,gradient_candidates=gradient,complete_candidates=full)
