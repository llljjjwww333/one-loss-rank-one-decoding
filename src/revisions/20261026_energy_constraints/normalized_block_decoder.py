"""Target-blind exact reduction for small cycle/parallel-path block masks.

Groebner bases are classical. Complete rational singletons are certified
from the full norm-augmented ideal. Shape failures remain unresolved.
"""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from itertools import product
import sympy as s
from exact_cactus_decoder import canonical_public,forest_base,actual_budget


def simple_cycles(m,n,edges):
    adj=[[] for _ in range(m+n)]
    for e,(i,j) in enumerate(edges):
        adj[i].append((m+j,e));adj[m+j].append((i,e))
    cycles={}
    def search(start,node,visited,path):
        for other,e in sorted(adj[node]):
            if other==start and len(path)>=3:
                cycle=path+[e];key=tuple(sorted(cycle))
                cycles.setdefault(key,cycle)
            elif other>start and other not in visited:
                search(start,other,visited|{other},path+[e])
    for node in range(m+n):search(node,node,{node},[])
    return [cycles[key] for key in sorted(cycles)]


def public_system(public):
    p=canonical_public(public);m,n,edges=(p[k] for k in ['m','n','edges'])
    u,v,gu,gv=(p[k] for k in ['u','v','gu','gv'])
    A=s.zeros(m+n,len(edges));parent=list(range(m+n));chords=[]
    def find(x):
        while parent[x]!=x:
            parent[x]=parent[parent[x]];x=parent[x]
        return x
    for e,(i,j) in enumerate(edges):
        A[i,e]=A[m+j,e]=1;a,b=find(i),find(m+j)
        if a==b:chords.append([e])
        else:parent[a]=b
    weights=[u[i]*v[j] for i,j in edges]
    rhs=[-u[i]*gu[i] for i in range(m)]+[-v[j]*gv[j] for j in range(n)]
    for e,(i,j) in enumerate(edges):rhs[i]+=weights[e]**2;rhs[m+j]+=weights[e]**2
    base=forest_base(m,n,edges,chords,rhs)
    if base is None:return None
    basis=A.nullspace();variables=s.symbols('t0:'+str(len(basis)))
    z=[base[e]+sum(t*b[e] for t,b in zip(variables,basis)) for e in range(len(edges))]
    equations=[]
    cycles=simple_cycles(m,n,edges)
    for cycle in cycles:
        f=s.expand(s.prod(z[e] for e in cycle[::2])-s.prod(z[e] for e in cycle[1::2]))
        if f!=0:equations.append(f)
    values=[s.expand(x/w) for x,w in zip(z,weights)]
    norm=2*p['loss']+sum(w*w for w in weights)-sum(x*y for x,y in zip(u,gu))-sum(x*y for x,y in zip(v,gv))
    energy=s.expand(sum(x*x for x in values)-norm)
    return dict(public=p,variables=variables,values=values,equations=equations,energy=energy,cycles=cycles,norm=norm)


def exact_reply(public,values):
    p=canonical_public(public);gu=[s.S.Zero]*p['m'];gv=[s.S.Zero]*p['n'];loss=s.S.Zero
    for (i,j),value in zip(p['edges'],values):
        z=p['u'][i]*p['v'][j]-value
        gu[i]+=z*p['v'][j];gv[j]+=z*p['u'][i];loss+=z*z/2
    return gu==p['gu'] and gv==p['gv'] and loss==p['loss']


def gradient_count(system):
    variables=system['variables'];F=system['equations'];z=s.Symbol('projection')
    if not variables:
        return dict(complete=True,real_actual_candidates=1,complex_degree=1)
    for offset in range(1,7):
        linear=sum(offset**i*t for i,t in enumerate(variables))
        G=s.groebner(F+[z-linear],*variables,z,order='lex',domain=s.QQ)
        if list(G)==[1]:return dict(complete=True,real_actual_candidates=0,complex_degree=0)
        univariate=[g for g in G if g.free_symbols<={z}]
        if len(univariate)!=1 or len(G)!=len(variables)+1:continue
        P=s.Poly(univariate[0],z,domain=s.QQ)
        if s.gcd(P,P.diff()).degree()!=0:continue
        substitution={}
        for variable in variables:
            relations=[g for g in G if g.has(variable) and not any(g.has(t) for t in variables if t!=variable)]
            if len(relations)!=1 or s.diff(relations[0],variable)!=1:break
            substitution[variable]=s.expand(variable-relations[0])
        if len(substitution)!=len(variables):continue
        # Every root is a nonzero toric point if no edge vanishes on P.
        edge_polynomials=[s.Poly(x.subs(substitution,simultaneous=True),z,domain=s.QQ).rem(P) for x in system['values']]
        if any(s.gcd(P,x).degree()!=0 for x in edge_polynomials):
            count=0;roots=P.real_roots(radicals=False)
            for root in roots:
                values=[x.as_expr().subs(z,root) for x in edge_polynomials]
                actual,_=actual_budget(system['public']['m'],system['public']['n'],system['public']['edges'],values,False)
                count+=int(actual)
            return dict(complete=True,real_actual_candidates=count,complex_degree=P.degree(),
                        real_algebraic_closure_points=len(roots),nonactual_real_points=len(roots)-count,
                        projection_offset=offset,zero_supports_checked=True)
        for f in F:
            assert s.Poly(f.subs(substitution,simultaneous=True),z,domain=s.QQ).rem(P).is_zero
        return dict(complete=True,real_actual_candidates=int(P.count_roots(-s.oo,s.oo)),
                    complex_degree=P.degree(),projection_offset=offset,
                    separating_polynomial=list(map(str,P.all_coeffs())),all_branches_nonzero=True)
    return dict(complete=False,reason='separating shape not obtained in six fixed trials')


def decode(public,with_gradient_count=True):
    system=public_system(public)
    if system is None:return dict(complete=True,status='Empty',reason='inconsistent margins')
    p=system['public'];variables=system['variables']
    if len(variables)>6:return dict(complete=False,status='Unresolved',reason='prespecified six-variable reference budget')
    if not variables:
        points=[()];G=[]
    else:
        G=s.groebner(system['equations']+[system['energy']],*variables,order='lex',domain=s.QQ)
        if list(G)==[1]:return dict(complete=True,status='Empty',norm_ideal=list(map(str,G)))
        if len(G)==len(variables) and all(s.Poly(g,*variables).total_degree()==1 for g in G):
            points=list(s.linsolve(list(G),variables))
        else:
            try:points=s.solve_poly_system(list(G),*variables,strict=True)
            except (NotImplementedError,s.polys.polyerrors.UnsolvableFactorError):points=None
        if points is None or any(not x.is_Rational for point in points for x in point):
            return dict(complete=False,status='Unresolved',reason='norm ideal not resolved into rational points',norm_ideal=list(map(str,G)))
    candidates=[]
    for point in points:
        subs=dict(zip(variables,point));values=[s.expand(x.subs(subs)) for x in system['values']]
        assert s.expand(sum(x*x for x in values)-system['norm'])==0
        unbounded_actual,_=actual_budget(p['m'],p['n'],p['edges'],values,False)
        if not unbounded_actual:continue
        actual,budget=actual_budget(p['m'],p['n'],p['edges'],values,True)
        assert exact_reply(public,values)
        candidates.append(dict(values=list(map(str,values)),unit_factor_budget_passed=actual,budget=str(budget)))
    result=dict(complete=True,status='Unique' if len(candidates)==1 else 'Ambiguous' if candidates else 'Empty',
                variables=len(variables),simple_cycles=len(system['cycles']),candidates=candidates,
                complete_count=len(candidates),norm_ideal=list(map(str,G)),original_reply_checked=True,
                truth_passed_to_solver=False,candidate_class='all actual real arrays without factor budget')
    if len(candidates)==1:result.update(candidates[0])
    if with_gradient_count:result['gradient']=gradient_count(system)
    return result
