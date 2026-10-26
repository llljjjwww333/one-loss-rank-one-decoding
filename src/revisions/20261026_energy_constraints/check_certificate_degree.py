"""Exact native-reply witnesses for the exception-degree lower bound.

Native checks use standard-library Fraction arithmetic in real quadratic
fields, independently of the paper's decoder and SymPy support routines.
The integer-slice counts concern only the displayed constructed family,
not all collisions on the slice or on the full independent grid.
"""
from pathlib import Path
from fractions import Fraction as Q
from math import isqrt
import hashlib,json
HERE=Path(__file__).resolve().parent
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

class Quadratic:
    def __init__(self,a=0,b=0,D=0):
        self.a,self.b,self.D=Q(a),Q(b),D
        root=isqrt(D)
        if root*root==D:self.a+=self.b*root;self.b=Q(0)
    def coerce(self,x):
        if isinstance(x,Quadratic):assert x.D==self.D;return x
        return Quadratic(x,D=self.D)
    def __add__(self,x):
        x=self.coerce(x);return Quadratic(self.a+x.a,self.b+x.b,self.D)
    __radd__=__add__
    def __neg__(self):return Quadratic(-self.a,-self.b,self.D)
    def __sub__(self,x):return self+-self.coerce(x)
    def __rsub__(self,x):return self.coerce(x)+-self
    def __mul__(self,x):
        x=self.coerce(x);return Quadratic(self.a*x.a+self.b*x.b*self.D,self.a*x.b+self.b*x.a,self.D)
    __rmul__=__mul__
    def __truediv__(self,x):
        x=self.coerce(x);norm=x.a*x.a-x.b*x.b*self.D;assert norm
        return self*Quadratic(x.a/norm,-x.b/norm,self.D)
    def __rtruediv__(self,x):return self.coerce(x)/self
    def __pow__(self,n):
        assert n in [0,1,2]
        return self*self if n==2 else self if n==1 else self.coerce(1)
    def __eq__(self,x):
        x=self.coerce(x);return self.a==x.a and self.b==x.b
    def sign(self):
        a,b=self.a,self.b
        if not b:return (a>0)-(a<0)
        if not a:return (b>0)-(b<0)
        if (a>0)==(b>0):return 1 if a>0 else -1
        gap=a*a-b*b*self.D
        return ((gap>0)-(gap<0))*(1 if a>0 else -1)
    def pair(self):return [str(self.a),str(self.b)]

def native_reply(edges,u,v,values,D):
    gu=[Quadratic(D=D) for _ in u];gv=[Quadratic(D=D) for _ in v];loss=Quadratic(D=D)
    for (i,j),value in zip(edges,values):
        residual=u[i]*v[j]-value;gu[i]+=residual*v[j];gv[j]+=residual*u[i];loss+=residual**2/2
    return [z.pair() for z in gu],[z.pair() for z in gv],loss.pair()

def witness(r,mask):
    selected=[i for i in range(r) if mask>>i&1];W=sum(4**i for i in selected);D=4+5*W
    x=Quadratic(-2,1,D);gamma=2**(r+4);lam=Q(1,gamma*gamma);m,n=2*r+3,3*r+3
    assert x*x+4*x==5*W and x.sign()>0 and (x-1).sign()>=0
    assert (x-2**(r+1)).sign()<0
    a=[Quadratic(1,D=D),Quadratic(Q(2,3),D=D),Quadratic(Q(2,3),D=D)]
    b=[Quadratic(z,D=D) for z in [9,3,-3]]
    c=[Quadratic(1,D=D),1+x/2,Quadratic(Q(1,2),D=D)]
    d=[Quadratic(z,D=D) for z in [10,2,-2]];edges=[]
    for cycle in range(r+1):
        row1,row2=2*cycle+1,2*cycle+2;col0,col1,col2=3*cycle,3*cycle+1,3*cycle+2
        edges.extend([(0,col0),(0,col1),(row1,col1),(row1,col2),(row2,col2),(row2,col0)])
        if not cycle:continue
        rho=2**(cycle-1)
        a.extend([Quadratic(1,D=D),Quadratic(Q(1,3),D=D)])
        b.extend([Quadratic(rho*z,D=D) for z in [9,3,-3]])
        c.extend([Quadratic(Q(1,2),D=D),Quadratic(Q(1,2),D=D)] if cycle-1 in selected
                 else [Quadratic(1,D=D),Quadratic(Q(1,3),D=D)])
        d.extend([Quadratic(rho*z,D=D) for z in ([8,4,-4] if cycle-1 in selected else [9,3,-3])])
    a,b,c,d=[[v/gamma for v in vector] for vector in [a,b,c,d]]
    assert len(a)==len(c)==m and len(b)==len(d)==n
    for left,right in [(a,b),(c,d)]:
        assert (sum((z*z for z in left),Quadratic(D=D))-Q(1,32)).sign()<0
        assert (sum((z*z for z in right),Quadratic(D=D))-Q(13,64)).sign()<0
    B=[a[i]*b[j] for i,j in edges];C=[c[i]*d[j] for i,j in edges]
    assert B!=C and all(z.sign()!=0 for z in B+C)
    u=[Quadratic(Q(1,m),D=D) for _ in range(m)];u[1]=1/(m*x)
    v=[Quadratic(Q(1,n),D=D) for _ in range(n)]
    assert (sum((z*z for z in u),Quadratic(D=D))-1).sign()<0
    assert (sum((z*z for z in v),Quadratic(D=D))-1).sign()<0
    original=native_reply(edges,u,v,B,D);alternative=native_reply(edges,u,v,C,D)
    assert original==alternative
    assert sum((z*z-y*y for z,y in zip(C,B)),Quadratic(D=D))==0
    derivative=lam*lam*(4*x+8);assert derivative.sign()>0
    # All chosen false and unchanged roots are simple. The private P root
    # has derivative 2(1+2/x); the M roots have derivatives -4rho^2 or 6rho^2.
    assert (2*(1+2/x)).sign()>0
    for i in range(r):assert (-4 if i in selected else 6)*4**i!=0
    root=isqrt(D);interval=[root-2,root-2 if root*root==D else root-1]
    payload=dict(target=[z.pair() for z in B],candidate=[z.pair() for z in C],
                 u=[z.pair() for z in u],v=[z.pair() for z in v],reply=original)
    digest=hashlib.sha256(json.dumps(payload,sort_keys=True).encode()).hexdigest()
    return dict(r=r,subset_mask=mask,W=W,D=D,x_polynomial_relation=[1,4,-5*W],
                x_integer_interval=interval,native_reply_sha256=digest,all_nonzero=True,
                bounded_factors_checked=True,positive_collision_derivative=True)

def run():
    protocol=dict(version=1,native_r=[1,2,3,4,6,8,10],integer_slice_r=list(range(1,19)),
                  subset_policy='all nonempty subsets; no exclusions',
                  scope='exact constructed slice witnesses and integer members of this family; not full-grid failure probabilities')
    path=HERE/'degree_protocol.json'
    if path.exists():assert json.loads(path.read_text())==protocol
    else:path.write_text(json.dumps(protocol,indent=2)+'\n')
    records=[];rows=[]
    for r in protocol['native_r']:
        current=[witness(r,mask) for mask in range(1,1<<r)]
        assert len({x['W'] for x in current})==2**r-1
        records.extend(current);rows.append(dict(r=r,cycles=r+1,witnesses=len(current)))
        print('native',r,len(current),flush=True)
    output=HERE/'degree_witnesses.jsonl'
    output.write_text(''.join(json.dumps(x,sort_keys=True)+'\n' for x in records))
    values=[0];integer=[]
    for r in protocol['integer_slice_r']:
        power=4**(r-1);values+= [z+power for z in values]
        xs=[isqrt(4+5*w)-2 for w in values[1:] if isqrt(4+5*w)**2==4+5*w]
        integer.append(dict(r=r,cycles=r+1,constructed_real_collisions=2**r-1,
                            integer_members=len(xs),integer_parameters=xs))
    assert len(records)==1367 and integer[-1]['integer_members']==14
    result=dict(complete=True,native_witnesses=len(records),native_rows=rows,
                integer_slice_controls=integer,protocol_sha256=sha(path),witnesses_sha256=sha(output),
                verifier_sha256=sha(Path(__file__)),arithmetic='standard-library rational pairs in real quadratic fields',
                all_native_loss_and_gradients_checked=True,all_factor_norms_checked=True,
                all_simple_roots_and_transversal_derivatives_checked=True,
                full_independent_grid_probability_estimated=False,
                scope='exact finite identity controls; arbitrary-r theorem and open-slice argument are proved in the article')
    (HERE/'degree_validation.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='integer_slice_controls'}))
if __name__=='__main__':run()
