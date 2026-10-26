"""Exact normalization witnesses for cycles and length-three parallel paths.

Counts concern the unit specialization, not generic ML-degree claims.
This proof-development control does not import the public decoder.
"""
from pathlib import Path
from itertools import combinations
import hashlib,json
import sympy as s
HERE=Path(__file__).resolve().parent


def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def paths(p):
    x=s.symbols('t1:'+str(p))
    t=(-sum(x),)+x
    pairs=[s.expand((1+t[i])**2*(1-t[j])-(1+t[j])**2*(1-t[i]))
           for i,j in combinations(range(p),2)]
    return x,t,pairs


def run():
    out=dict(complete=False,cycle_controls=[],parallel_path_controls=[],boundary_controls=[])
    t,z=s.symbols('t z')
    for k in range(2,17):
        F=s.Poly((1+t)**k-(1-t)**k,t)
        G=F.exquo(s.Poly(t,t));d=k if k%2 else k-1
        R=s.Poly(s.expand(G.as_expr().subs(t,s.I*z)),z)
        assert all(not c.has(s.I) for c in R.all_coeffs())
        assert F.degree()==d and R.degree()==d-1
        assert R.eval(0)!=0 and len(R.intervals(eps=s.Rational(1,10**6)))==d-1
        assert s.gcd(R,R.diff()).degree()==0
        out['cycle_controls'].append(dict(half_length=k,degree=d,nonzero_roots=d-1,
                                         all_nonzero_roots_pure_imaginary=True))
    for p in [2,3,4]:
        x,flows,F=paths(p)
        highest=[sum(c*s.prod(v**a for v,a in zip(x,mon))
                     for mon,c in s.Poly(f,*x).terms() if sum(mon)==3) for f in F]
        patches=[]
        for variable in x:
            ideal=s.groebner(highest+[variable-1],*x,order='lex')
            assert list(ideal)==[1]
            patches.append(str(variable))
        roots=[(s.S.Zero,)*p];patterns=[]
        for j in range(1,p):
            c=s.Rational(p-j,j);P=c*z*z+(1-c)*z+3
            assert s.discriminant(P,z)<0
            beta_re=(c-1)/(2*c)
            beta2_re=s.simplify((c*c-8*c+1)/(2*c*c))
            changes=[]
            for multiplier in [s.S.One,-c]:
                for sign in [-1,1]:
                    gap=s.simplify(multiplier**2*beta2_re+2*sign*multiplier*beta_re)
                    assert gap<0,(p,j,multiplier,sign,gap)
                    changes.append(str(gap))
            patterns.append(dict(multiplicity=j,ratio=str(c),edge_real_energy_changes=changes))
            for beta in s.solve(P,z):
                for selected in combinations(range(p),j):
                    root=tuple(s.simplify(-c*beta if i in selected else beta) for i in range(p))
                    assert s.simplify(sum(root))==0
                    for f in F:assert s.simplify(f.subs(dict(zip(x,root[1:]))))==0
                    if not any(all(s.simplify(a-b)==0 for a,b in zip(root,other)) for other in roots):roots.append(root)
        assert len(roots)==2**p-1
        jac=s.Matrix(F).jacobian(x)
        assert jac.subs(dict.fromkeys(x,0)).rank()==p-1
        zz=s.Symbol('projection')
        for offset in range(1,10):
            linear=sum((offset**i)*v for i,v in enumerate(x))
            gb=s.groebner(F+[zz-linear],*x,zz,order='lex',domain=s.QQ)
            uni=[g for g in gb if g.free_symbols<={zz}]
            assert len(uni)==1
            if s.degree(uni[0],zz)==len(roots):break
        else:raise AssertionError('no separating projection')
        assert len(gb)==len(x)+1
        for variable in x:
            relations=[g for g in gb if g.has(variable) and
                       not any(g.has(other) for other in x if other!=variable)]
            assert len(relations)==1 and s.diff(relations[0],variable)==1
        assert s.gcd(s.Poly(uni[0],zz),s.Poly(uni[0],zz).diff()).degree()==0
        out['parallel_path_controls'].append(dict(paths=p,vertices=2*p+2,edges=3*p,
            affine_points=len(roots),groebner_projection_degree=s.degree(uni[0],zz),
            projection_offset=offset,projection_polynomial=str(uni[0]),
            no_roots_at_infinity_patches=patches,true_root_simple=True,patterns=patterns))
    c=s.Integer(4);beta_re=(c-1)/(2*c);beta2_re=(c*c-8*c+1)/(2*c*c)
    gap=s.simplify(beta2_re+2*beta_re);assert gap==s.Rational(9,32)
    out['boundary_controls'].append(dict(paths=5,one_outer_edge_change=str(gap),
        interpretation='edgewise contraction proof fails; this is not an identifiability counterexample'))
    out.update(complete=True,verifier_sha256=sha(Path(__file__)),
               scope='exact cycle root checks, all 7/15 unit theta points, projective boundary and contraction identities')
    (HERE/'unit_block_validation.json').write_text(json.dumps(out,indent=2,default=str)+'\n')
    print(json.dumps(out,indent=2,default=str),flush=True)


if __name__=='__main__':run()
