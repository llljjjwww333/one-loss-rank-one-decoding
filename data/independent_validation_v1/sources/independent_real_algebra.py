"""Rational polynomial, Sturm and real-field predicates without a CAS import.

Univariate factors are supplied by Singular. Every sign uses an exact
zero test plus rational interval refinement; budgets use an integer norm
gap, including equality. No floating-point tolerance is used.
"""
from fractions import Fraction as Q
from math import gcd,isqrt,lcm


def trim(p):
    p=list(map(Q,p))
    while p and not p[-1]: p.pop()
    return p


def add(p,q):
    return trim([(p[i] if i<len(p) else Q(0))+(q[i] if i<len(q) else Q(0)) for i in range(max(len(p),len(q)))])


def neg(p): return [-x for x in p]


def mul(p,q):
    if not p or not q: return []
    out=[Q(0)]*(len(p)+len(q)-1)
    for i,a in enumerate(p):
        for j,b in enumerate(q): out[i+j]+=a*b
    return trim(out)


def divrem(p,q):
    p=trim(p); q=trim(q); assert q
    out=[Q(0)]*max(0,len(p)-len(q)+1)
    while p and len(p)>=len(q):
        k=len(p)-len(q); a=p[-1]/q[-1]; out[k]=a
        for j,b in enumerate(q): p[k+j]-=a*b
        p=trim(p)
    return trim(out),p


def monic(p):
    p=trim(p)
    return [x/p[-1] for x in p] if p else []


def derivative(p): return trim([i*p[i] for i in range(1,len(p))])


def pgcd(p,q):
    while q: p,q=q,divrem(p,q)[1]
    return monic(p)


def inverse(p,f):
    r0,r1=f,divrem(p,f)[1]; s0,s1=[],[Q(1)]
    while r1:
        quotient,remainder=divrem(r0,r1)
        r0,r1=r1,remainder; s0,s1=s1,add(s0,neg(mul(quotient,s1)))
    assert len(r0)==1, 'A supplied factor did not define a field for this inverse'
    return divrem([x/r0[0] for x in s0],f)[1]


def evaluate(p,x):
    value=Q(0)
    for coefficient in reversed(p): value=value*x+coefficient
    return value


def integer_primitive(p):
    p=trim(p); denominator=lcm(*(x.denominator for x in p))
    coefficients=[int(x*denominator) for x in p]
    content=gcd(*coefficients)
    coefficients=[x//content for x in coefficients]
    return coefficients if coefficients[-1]>0 else [-x for x in coefficients]


def sturm(p):
    p=trim(p); result=[p,derivative(p)]
    while result[-1]:
        remainder=neg(divrem(result[-2],result[-1])[1])
        if not remainder: break
        result.append(remainder)
    return [p for p in result if p]


def variations(sequence,x):
    signs=[1 if value>0 else -1 for p in sequence if (value:=evaluate(p,x))]
    return sum(a!=b for a,b in zip(signs,signs[1:]))


def count(sequence,a,b): return variations(sequence,a)-variations(sequence,b)


def product_interval(a,b):
    values=[x*y for x in a for y in b]
    return min(values),max(values)


def interval(p,a,b):
    low=high=Q(0)
    for coefficient in reversed(p):
        low,high=product_interval((low,high),(a,b))
        low+=coefficient; high+=coefficient
    return low,high


def cauchy(p):
    ratio=max([abs(x/p[-1]) for x in p[:-1]]+[Q(0)])
    return 2+(ratio.numerator//ratio.denominator)


class RealRoot:
    def __init__(self,f,a,b,index):
        self.f=monic(f); self.a=Q(a); self.b=Q(b); self.index=index
        self.sequence=sturm(self.f)
        self.integer_polynomial=integer_primitive(self.f)
        self.degree=len(self.f)-1
        self.complex_bound=cauchy(self.f)
        assert self.a==self.b and not evaluate(self.f,self.a) or count(self.sequence,self.a,self.b)==1

    def refine(self,bits):
        while self.b-self.a>Q(1,1<<bits):
            mid=(self.a+self.b)/2; value=evaluate(self.f,mid)
            if not value: self.a=self.b=mid; break
            if count(self.sequence,self.a,mid): self.b=mid
            else: self.a=mid

    def sign(self,p):
        p=divrem(p,self.f)[1]
        if not p: return 0
        common=pgcd(self.f,p)
        if len(common)>1:
            if self.a==self.b and not evaluate(common,self.a): return 0
            if self.a!=self.b and count(sturm(common),self.a,self.b): return 0
        bits=8
        while True:
            self.refine(bits); low,high=interval(p,self.a,self.b)
            if low>0: return 1
            if high<0: return -1
            bits*=2

    def value(self,p): return Element(self,p)

    def json(self):
        return dict(integer_polynomial=self.integer_polynomial,real_root_index=self.index,
                    isolating_interval=[str(self.a),str(self.b)])


def isolate(f):
    f=monic(f); assert len(f)>1 and len(pgcd(f,derivative(f)))==1
    sequence=sturm(f); bound=cauchy(f)
    queue=[(Q(-bound),Q(bound),count(sequence,Q(-bound),Q(bound)))]; boxes=[]
    while queue:
        a,b,number=queue.pop()
        if not number: continue
        if number==1: boxes.append((a,b)); continue
        mid=(a+b)/2
        if not evaluate(f,mid):
            mid=next((k*a+b)/(k+1) for k in range(2,len(f)+3) if evaluate(f,(k*a+b)/(k+1)))
        left=count(sequence,a,mid)
        queue.extend([(a,mid,left),(mid,b,number-left)])
    boxes.sort()
    return [RealRoot(f,a,b,index) for index,(a,b) in enumerate(boxes)]


class Element:
    def __init__(self,root,p):
        self.root=root; self.p=divrem(trim(p),root.f)[1]

    def coerce(self,x):
        if isinstance(x,Element): assert x.root is self.root; return x
        return Element(self.root,[Q(x)])

    def __add__(self,x): return Element(self.root,add(self.p,self.coerce(x).p))
    __radd__=__add__
    def __neg__(self): return Element(self.root,neg(self.p))
    def __sub__(self,x): return self+-self.coerce(x)
    def __rsub__(self,x): return self.coerce(x)+-self
    def __mul__(self,x): return Element(self.root,mul(self.p,self.coerce(x).p))
    __rmul__=__mul__
    def __truediv__(self,x):
        other=self.coerce(x); assert other.sign()!=0
        return Element(self.root,mul(self.p,inverse(other.p,self.root.f)))
    def __rtruediv__(self,x): return self.coerce(x)/self
    def __pow__(self,n):
        assert n>=0
        result=self.coerce(1); base=self
        while n:
            if n&1: result=result*base
            base=base*base; n//=2
        return result
    def sign(self): return self.root.sign(self.p)
    def __bool__(self): return self.sign()!=0
    def __eq__(self,x): return (self-self.coerce(x)).sign()==0
    def bounds(self): return interval(self.p,self.root.a,self.root.b)
    def json(self): return list(map(str,self.p))


def sqrt_bound(x,bits,upper):
    assert x>=0
    scale=1<<bits; k=isqrt((x.numerator*scale*scale)//x.denominator)
    if upper and Q(k,scale)**2!=x: k+=1
    return Q(k,scale)


def budget(radicands):
    if not radicands: return True,dict(sign=-1,components=0,kind='empty')
    root=radicands[0].root; assert all(value.root is root and value.sign()>0 for value in radicands)
    denominator=1
    for value in radicands:
        for coefficient in value.p: denominator*=coefficient.denominator
    leading=abs(root.integer_polynomial[-1])
    q=denominator*leading**(root.degree-1)
    conjugate=max([sum(abs(c)*root.complex_bound**i for i,c in enumerate(value.p)) for value in radicands]+[Q(1)])
    crude=q*(1+len(radicands)*max(Q(1),conjugate))
    bound=(crude.numerator+crude.denominator-1)//crude.denominator
    precision=q.bit_length()+(root.degree*(1<<len(radicands))-1)*bound.bit_length()
    gap=Q(1,1<<precision); bits=8
    while True:
        root.refine(bits); low=high=Q(-1)
        for value in radicands:
            a,b=value.bounds(); assert b>=0
            low+=sqrt_bound(max(Q(0),a),bits,False)
            high+=sqrt_bound(b,bits,True)
        if high<0: sign=-1
        elif low>0: sign=1
        elif high-low<gap: sign=0
        else: bits*=2; continue
        return sign<=0,dict(kind='algebraic_radical_integer_norm_gap',sign=sign,
            components=len(radicands),precision_bits=bits,nonzero_gap_bits=precision,
            integer_scale_bits=q.bit_length(),conjugate_bound_bits=bound.bit_length(),
            interval=[str(low),str(high)],radicands=[x.json() for x in radicands])
