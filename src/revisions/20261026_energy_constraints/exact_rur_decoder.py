"""Classical trace RUR plus a replayable nilpotence coverage certificate.

Rational public inputs; arbitrary real algebraic candidates. Nonradical
zero-dimensional ideals are allowed. No radical formula is required.
Quotient-size limits are reported as Unresolved, never as emptiness.
"""
from collections import deque
from pathlib import Path
import hashlib,json,sys
from math import gcd,lcm
import sympy as s
from normalized_block_decoder import public_system as raw_public_system
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from exact_cactus_decoder import actual_budget,canonical_public

Z=s.Symbol('_rur_z')

class Unresolved(Exception):pass

def public_system(public):
    system=raw_public_system(public)
    if system is None:return None
    variables=system['variables'];p=system['public'];zero=dict.fromkeys(variables,0)
    constants=[s.Rational(f.subs(zero)*p['u'][i]*p['v'][j])
               for f,(i,j) in zip(system['values'],p['edges'])]
    denominator=lcm(*(int(x.q) for x in constants));numerator=0
    for x in constants:numerator=gcd(numerator,abs(int(x*denominator)))
    scale=s.Rational(numerator,denominator) if numerator else s.S.One
    system['rational_flow_scale']=str(scale)
    if variables:
        substitution={v:scale*v for v in variables}
        def primitive(f):
            return s.Poly(f.subs(substitution,simultaneous=True),*variables,domain=s.QQ).primitive()[1].as_expr()
        system['equations']=[primitive(f) for f in system['equations']]
        system['energy']=primitive(system['energy'])
        system['values']=[s.expand(f.subs(substitution,simultaneous=True)) for f in system['values']]
    return system

def coefficients(p):return list(map(str,p.all_coeffs()))
def polynomial(c):return s.Poly.from_list(list(map(s.Rational,c)),Z,domain=s.QQ)

def multiplication_algebra(equations,variables,max_dimension=None):
    G=s.groebner(equations,*variables,order='grevlex',domain=s.QQ)
    if list(G)==[1]:return G,[],[]
    if not G.is_zero_dimensional:raise Unresolved('positive-dimensional norm ideal')
    leads=[p.LM(order=G.order).exponents for p in G.polys]
    zero=(0,)*len(variables);seen={zero};queue=deque([zero]);basis=[]
    while queue:
        mon=queue.popleft()
        if any(all(a>=b for a,b in zip(mon,lm)) for lm in leads):continue
        basis.append(mon)
        if max_dimension is not None and len(basis)>max_dimension:
            raise Unresolved('quotient dimension exceeds supplied budget')
        for i in range(len(variables)):
            child=tuple(a+int(j==i) for j,a in enumerate(mon))
            if child not in seen:seen.add(child);queue.append(child)
    basis.sort(key=lambda m:(sum(m),m));index={m:i for i,m in enumerate(basis)}
    monomials=[s.prod(v**k for v,k in zip(variables,m)) for m in basis]
    matrices=[]
    for variable in variables:
        M=s.zeros(len(basis))
        for col,monomial in enumerate(monomials):
            remainder=s.Poly(G.reduce(variable*monomial)[1],*variables,domain=s.QQ)
            for mon,value in remainder.terms():M[index[mon],col]=value
        matrices.append(M)
    return G,basis,matrices

def evaluate_matrix(poly,M):
    identity=s.eye(M.rows);out=s.zeros(M.rows)
    for a in poly.all_coeffs():out=out*M+a*identity
    return out

def nilpotent(M,dimension,method='charpoly'):
    if M.is_zero_matrix:return True
    if method=='power':return (M**dimension).is_zero_matrix is True
    return M.charpoly(Z).all_coeffs()==[s.S.One]+[s.S.Zero]*dimension

def rational_univariate(equations,variables,max_dimension=None):
    G,basis,matrices=multiplication_algebra(equations,variables,max_dimension)
    if not basis:return dict(empty=True,groebner=list(map(str,G)))
    dimension=len(basis);identity=s.eye(dimension);trials=[]
    bound=(len(variables)-1)*dimension*(dimension-1)//2+1
    for offset in range(1,bound+1):
        weights=[s.Integer(offset)**i for i in range(len(variables))]
        M=sum((a*m for a,m in zip(weights,matrices)),s.zeros(dimension))
        characteristic=M.charpoly(Z).as_poly();P=characteristic.sqf_part().monic();degree=P.degree()
        powers=[identity]
        for _ in range(1,degree):powers.append(powers[-1]*M)
        def trace_numerator(operator):
            traces=[s.trace(operator*m) for m in powers]
            return s.Poly(sum(P.nth(k)*traces[j]*Z**(k-1-j)
                              for k in range(1,degree+1) for j in range(k)),Z,domain=s.QQ)
        denominator=trace_numerator(identity)
        assert s.gcd(denominator,P).degree()==0
        inverse=s.invert(denominator,P)
        coordinates=[(trace_numerator(m)*inverse).rem(P) for m in matrices]
        tests=[nilpotent(m-evaluate_matrix(v,M),dimension) for m,v in zip(matrices,coordinates)]
        trials.append(dict(offset=offset,distinct_projected_values=degree,separating=all(tests)))
        if not all(tests):continue
        substitution={v:w.as_expr() for v,w in zip(variables,coordinates)}
        for f in equations:
            assert s.Poly(f.subs(substitution,simultaneous=True),Z,domain=s.QQ).rem(P).is_zero
        assert s.Poly(sum(a*v.as_expr() for a,v in zip(weights,coordinates))-Z,Z).rem(P).is_zero
        return dict(empty=False,quotient_dimension=dimension,distinct_complex_points=degree,
                    projection_weights=list(map(str,weights)),projection_bound=bound,trials=trials,
                    polynomial=coefficients(P),coordinates=[coefficients(v) for v in coordinates],
                    characteristic=coefficients(characteristic),standard_monomials=[list(m) for m in basis],
                    nilpotence_exponent=dimension,groebner=list(map(str,G)),
                    method='classical multiplicity-weighted trace RUR; nilpotence-certified coverage')
    raise AssertionError('deterministic separating-form bound exhausted')

def verify_rur(equations,variables,certificate,max_dimension=None):
    """Check coverage with matrix powers, without recomputing trace formulas."""
    G,basis,matrices=multiplication_algebra(equations,variables,max_dimension)
    if certificate['empty']:
        assert list(G)==[1];return
    dimension=len(basis);assert dimension==certificate['quotient_dimension']
    assert [list(m) for m in basis]==certificate['standard_monomials']
    assert certificate['nilpotence_exponent']==dimension
    P=polynomial(certificate['polynomial']);coords=list(map(polynomial,certificate['coordinates']))
    assert P.LC()==1 and s.gcd(P,P.diff()).degree()==0 and len(coords)==len(variables)
    assert P.degree()==certificate['distinct_complex_points']
    assert all(v.degree()<P.degree() for v in coords)
    weights=list(map(s.Rational,certificate['projection_weights']))
    assert len(weights)==len(variables)
    M=sum((a*m for a,m in zip(weights,matrices)),s.zeros(dimension))
    assert nilpotent(evaluate_matrix(P,M),dimension,'power')
    for variable,m,v in zip(variables,matrices,coords):
        assert nilpotent(m-evaluate_matrix(v,M),dimension,'power'),str(variable)
    substitution={v:w.as_expr() for v,w in zip(variables,coords)}
    for f in equations:assert s.Poly(f.subs(substitution,simultaneous=True),Z).rem(P).is_zero
    assert s.Poly(sum(a*v.as_expr() for a,v in zip(weights,coords))-Z,Z).rem(P).is_zero

def reply_identities(system,values,P):
    p=system['public'];gu=[s.S.Zero]*p['m'];gv=[s.S.Zero]*p['n'];loss=s.S.Zero
    for (i,j),v in zip(p['edges'],values):
        residual=p['u'][i]*p['v'][j]-v.as_expr()
        gu[i]+=residual*p['v'][j];gv[j]+=residual*p['u'][i];loss+=residual**2/2
    for value,expected in list(zip(gu,p['gu']))+list(zip(gv,p['gv']))+[(loss,p['loss'])]:
        assert s.Poly(s.expand(value-expected),Z,domain=s.QQ).rem(P).is_zero

def decode(public,bounded=False,max_dimension=None):
    system=public_system(public)
    if system is None:return dict(complete=True,status='Empty',complete_count=0,reason='inconsistent margins')
    variables=system['variables'];p=system['public']
    if not variables:
        values=system['values'];ok=system['energy']==0 and all(f==0 for f in system['equations'])
        actual,budget=actual_budget(p['m'],p['n'],p['edges'],values,bounded) if ok else (False,None)
        candidates=[dict(values=list(map(str,values)),budget=str(budget))] if actual else []
        return dict(complete=True,status='Unique' if candidates else 'Empty',complete_count=len(candidates),candidates=candidates)
    equations=system['equations']+[system['energy']]
    try:certificate=rational_univariate(equations,variables,max_dimension)
    except Unresolved as error:return dict(complete=False,status='Unresolved',reason=str(error))
    if certificate['empty']:return dict(complete=True,status='Empty',complete_count=0,candidates=[],rur=certificate)
    P=polynomial(certificate['polynomial']);coordinates=list(map(polynomial,certificate['coordinates']))
    substitution={v:w.as_expr() for v,w in zip(variables,coordinates)}
    edges=[s.Poly(f.subs(substitution,simultaneous=True),Z,domain=s.QQ).rem(P) for f in system['values']]
    reply_identities(system,edges,P)
    roots=P.real_roots(radicals=False);intervals=P.intervals(eps=s.Rational(1,65536))
    assert len(roots)==len(intervals)
    candidates=[];decisions=[]
    for ordinal,(root,(interval,multiplicity)) in enumerate(zip(roots,intervals)):
        assert multiplicity==1
        values=[v.as_expr().subs(Z,root) for v in edges]
        actual,budget=actual_budget(p['m'],p['n'],p['edges'],values,bounded)
        decision=dict(real_root_ordinal=ordinal,isolating_interval=list(map(str,interval)),actual_and_budget_passed=actual)
        decisions.append(decision)
        if actual:candidates.append(dict(**decision,values=list(map(str,values)),budget=str(budget),
                                         rational=all(v.is_Rational for v in values)))
    return dict(complete=True,status='Unique' if len(candidates)==1 else 'Ambiguous' if candidates else 'Empty',
                complete_count=len(candidates),candidates=candidates,rur=certificate,
                edge_polynomials=[coefficients(v) for v in edges],real_root_decisions=decisions,
                real_closure_points=len(roots),original_reply_checked=True,truth_passed_to_solver=False,
                rational_flow_scale=system['rational_flow_scale'],
                candidate_class='actual real arrays'+(' with unit factor budget' if bounded else ' without factor budget'))

if __name__=='__main__':
    import argparse
    parser=argparse.ArgumentParser();parser.add_argument('public');parser.add_argument('--output',required=True)
    parser.add_argument('--bounded',action='store_true');parser.add_argument('--max-dimension',type=int)
    args=parser.parse_args();public=json.loads(Path(args.public).read_text())
    result=decode(public,args.bounded,args.max_dimension)
    result['public_sha256']=hashlib.sha256(Path(args.public).read_bytes()).hexdigest()
    Path(args.output).write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ['rur','candidates','edge_polynomials','real_root_decisions']}))
