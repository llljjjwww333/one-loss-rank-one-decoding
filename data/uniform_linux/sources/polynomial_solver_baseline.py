"""Exact generic polynomial baselines from PUBLIC replies only.

Factor mode uses the ORIGINAL gradient and loss equations, with support
strata and component gauges. Reduced mode uses affine cycle coordinates.
Both call SymPy's lex Groebner/solve_triangulated(extension=True) backend.
No claim of representing the fastest computer algebra solver is made.
"""
import argparse
from itertools import product
import json
from pathlib import Path
from time import perf_counter
import sympy as s
from exact_cactus_decoder import (canonical_public,cycle_basis,forest_base,
                                  actual_budget,exact_sign)


def expressions(values):
    result=[]
    for value in values:
        value=s.sympify(value)
        while value.has(s.AlgebraicNumber):
            value=value.xreplace({a:a.as_expr() for a in value.atoms(s.AlgebraicNumber)})
        result.append(s.cancel(value))
    return result


def is_real(value):
    if value.is_real is not None:
        return bool(value.is_real)
    return all(c==0 for c in s.to_number_field(value-s.conjugate(value)).coeffs())


def equal(a,b):
    return all(exact_sign(s.expand(x-y))==0 for x,y in zip(a,b))


def new_stats():
    return dict(groebner_calls=0,groebner_seconds=0.0,triangular_seconds=0.0,
                filter_seconds=0.0,max_variables=0,max_equation_degree=0,
                complex_tuples=0,real_tuples=0,unit_ideals=0)


class PositiveDimension(Exception):
    pass


def solve_exact(equations,variables,stats):
    equations=[s.expand(e) for e in equations if e!=0]
    if not variables:
        return [()] if all(e==0 for e in equations) else []
    primitive=[]
    for expr in equations:
        poly=s.Poly(expr,*variables,domain=s.QQ)
        primitive.append(poly.clear_denoms(convert=True)[1].primitive()[1].as_expr())
        stats['max_equation_degree']=max(stats['max_equation_degree'],poly.total_degree())
    stats['max_variables']=max(stats['max_variables'],len(variables))
    start=perf_counter()
    gb=s.groebner(primitive,*variables,domain=s.QQ,order='lex')
    stats['groebner_seconds']+=perf_counter()-start
    stats['groebner_calls']+=1
    if any(poly.total_degree()==0 and not poly.is_zero for poly in gb.polys):
        stats['unit_ideals']+=1
        return []
    if not gb.is_zero_dimensional:
        raise PositiveDimension('Gauge-fixed stratum has a positive-dimensional complex solution set')
    start=perf_counter()
    # This API recomputes the basis from the supplied basis. That additional
    # work is included in triangular_seconds, rather than hidden in preprocessing.
    solutions=s.solve_triangulated([p.as_expr() for p in gb.polys],*variables,extension=True)
    stats['triangular_seconds']+=perf_counter()-start
    stats['complex_tuples']+=len(solutions)
    real=[]
    for values in solutions:
        values=expressions(values)
        if not all(is_real(value) for value in values):
            continue
        substitution=dict(zip(variables,values))
        # Filter backend solutions against ALL original stratum equations.
        assert all(exact_sign(s.expand(e.subs(substitution)))==0 for e in equations)
        real.append(values)
    stats['real_tuples']+=len(real)
    return real


def weighted_margins(p):
    pu,pv=[s.S.Zero]*p['m'],[s.S.Zero]*p['n']
    for i,j in p['edges']:
        pu[i]+=p['u'][i]*p['v'][j]**2
        pv[j]+=p['v'][j]*p['u'][i]**2
    return [p['u'][i]*(pu[i]-p['gu'][i]) for i in range(p['m'])]+[
            p['v'][j]*(pv[j]-p['gv'][j]) for j in range(p['n'])]


def accept(p,values,bounded,arrays,stats):
    values=expressions(values)
    start=perf_counter()
    # solve_triangulated can represent extension roots in several forms.
    # Convert pure radicals if needed for the shared public decoder's filter.
    if any(not x.is_Rational and not x.has(s.CRootOf) for x in values):
        values=[s.to_number_field(x).to_root(radicals=False) if not x.is_Rational else x for x in values]
    admissible,_=actual_budget(p['m'],p['n'],p['edges'],values,bounded)
    if admissible and not any(equal(values,previous) for previous in arrays):
        arrays.append(values)
    stats['filter_seconds']+=perf_counter()-start


def factor_decode(p,bounded,stats):
    m,n,edges=p['m'],p['n'],p['edges']
    factors=list(s.symbols('a:'+str(m)))+list(s.symbols('b:'+str(n)))
    margins=weighted_margins(p)
    incident={vertex for i,j in edges for vertex in [i,m+j]}
    forced={i for i,value in enumerate(margins) if value!=0}
    if not forced.issubset(incident):
        return []
    optional=sorted(incident-forced)
    stats.update(optional_zero_margin_vertices=len(optional),support_strata_planned=2**len(optional),
                 support_strata_solved=0,support_strata_skipped=0)
    residuals=[p['u'][i]*p['v'][j]-factors[i]*factors[m+j] for i,j in edges]
    original=[sum((residuals[e]*p['v'][j] for e,(ii,j) in enumerate(edges) if ii==i),s.S.Zero)-p['gu'][i] for i in range(m)]
    original +=[sum((residuals[e]*p['u'][i] for e,(i,jj) in enumerate(edges) if jj==j),s.S.Zero)-p['gv'][j] for j in range(n)]
    original +=[sum((r*r for r in residuals),s.S.Zero)-2*p['loss']]
    arrays=[]
    for on in product([False,True],repeat=len(optional)):
        active=forced|{vertex for vertex,flag in zip(optional,on) if flag}
        adjacency={i:[] for i in active}
        for i,j in edges:
            if i in active and m+j in active:
                adjacency[i].append(m+j)
                adjacency[m+j].append(i)
        if any(not neighbors for neighbors in adjacency.values()):
            # Such a nonzero factor has no nonzero observed edge; setting it
            # to zero gives the same array in an enumerated smaller support.
            stats['support_strata_skipped']+=1
            continue
        fixed={factors[i]:s.S.Zero for i in range(m+n) if i not in active}
        visited=set()
        for first in sorted(active):
            if first in visited:
                continue
            fixed[factors[first]]=s.S.One
            visited.add(first)
            stack=[first]
            while stack:
                node=stack.pop()
                for neighbor in adjacency[node]:
                    if neighbor not in visited:
                        visited.add(neighbor)
                        stack.append(neighbor)
        variables=[x for x in factors if x not in fixed]
        equations=[e.subs(fixed) for e in original]
        # A nonzero public margin already forces its vertex factor nonzero.
        # Saturate only the remaining declared-active free factors.
        nonzero=[factors[i] for i in sorted(active-forced) if factors[i] not in fixed]
        if nonzero:
            auxiliary=s.Symbol('inverse_support')
            equations.append(auxiliary*s.prod(nonzero)-1)
            variables=[auxiliary]+variables
        roots=solve_exact(equations,variables,stats)
        stats['support_strata_solved']+=1
        for values in roots:
            solved=fixed|dict(zip(variables,values))
            candidate=[solved[factors[i]]*solved[factors[m+j]] for i,j in edges]
            accept(p,candidate,bounded,arrays,stats)
    return arrays


def reduced_decode(p,bounded,stats):
    cycles=cycle_basis(p['m'],p['n'],p['edges'])
    base=forest_base(p['m'],p['n'],p['edges'],cycles,weighted_margins(p))
    if base is None:
        return []
    variables=list(s.symbols('t:'+str(len(cycles))))
    z=list(base)
    equations=[]
    for t,cycle in zip(variables,cycles):
        for offset,e in enumerate(cycle):
            z[e]+=(-1)**offset*t
        equations.append(s.prod(z[e] for e in cycle[::2])-s.prod(z[e] for e in cycle[1::2]))
    weights=[p['u'][i]*p['v'][j] for i,j in p['edges']]
    norm=2*p['loss']+sum(w*w for w in weights)-sum(p['u'][i]*p['gu'][i] for i in range(p['m']))-sum(p['v'][j]*p['gv'][j] for j in range(p['n']))
    candidates=[entry/w for entry,w in zip(z,weights)]
    equations.append(sum(c*c for c in candidates)-norm)
    roots=solve_exact(equations,variables,stats)
    arrays=[]
    for values in roots:
        substitution=dict(zip(variables,values))
        accept(p,[c.subs(substitution) for c in candidates],bounded,arrays,stats)
    return arrays


def decode(public,method,bounded=False):
    p=canonical_public(public)
    stats=new_stats()
    try:
        arrays=(factor_decode if method=='factor_groebner' else reduced_decode)(p,bounded,stats)
    except PositiveDimension as error:
        return dict(complete=False,status='positive_dimensional',reason=str(error),stats=stats,arrays=[])
    return dict(complete=True,status='Empty' if not arrays else 'Unique' if len(arrays)==1 else 'Ambiguous',
                arrays=[[str(x) for x in values] for values in arrays],stats=stats)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('public',type=Path)
    parser.add_argument('output',type=Path)
    parser.add_argument('--method',choices=['factor_groebner','reduced_groebner','cycle_enumerate','cycle_mitm'],required=True)
    parser.add_argument('--bounded',action='store_true')
    args=parser.parse_args()
    start=perf_counter()
    public=json.loads(args.public.read_text(encoding='utf-8'))
    if args.method.startswith('cycle_'):
        from exact_cactus_decoder import decode as cycle_decode
        cert=cycle_decode(public,bounded=args.bounded,search='mitm' if args.method=='cycle_mitm' else 'enumerate')
        result=dict(complete=cert['complete'],status=cert['status'],arrays=[c['values'] for c in cert['candidates']],
                    stats=dict(tuple_bound=cert.get('tuple_bound'),energy_search=cert.get('energy_search')))
    else:
        result=decode(public,args.method,bounded=args.bounded)
    result.update(method=args.method,bounded=args.bounded,solver_core_seconds=perf_counter()-start,
                  sympy=s.__version__,truth_passed_to_solver=False)
    args.output.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(result['status'],len(result['arrays']))


if __name__=='__main__':
    main()
