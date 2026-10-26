"""Exact, parameter-uniform checks of the scaled theta proof witness.

This is proof development. Hidden target factors are never decoder inputs.
The coefficient certificates cover every real c >= 1, epsilon >= c+1;
the finite path counts are additional controls, not the proof of that range.
"""
from pathlib import Path
from fractions import Fraction as Q
from itertools import combinations
import hashlib,json
import sympy as s

HERE=Path(__file__).resolve().parent

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def glue_control():
    # Every new block attaches at one vertex, including endpoint/internal
    # attachments on both sides of the bipartition and a bridge.
    assigned={};blocks=[];next_id=0
    for index,kind in enumerate([('theta',5),('cycle',4),('theta',15),
                                ('bridge',1),('theta',32),('cycle',3)]):
        name,p=kind
        if name=='theta':
            row=['U']+['r'+str(i) for i in range(p)]
            col=['V']+['c'+str(i) for i in range(p)]
            factors={v:Q(p if v in ['U','V'] else 1) for v in row+col}
            edges=[e for i in range(p) for e in [('U','c'+str(i)),
                   ('r'+str(i),'c'+str(i)),('r'+str(i),'V')]]
        elif name=='cycle':
            row=['r'+str(i) for i in range(p)];col=['c'+str(i) for i in range(p)]
            factors={v:Q(1) for v in row+col}
            edges=[e for i in range(p) for e in [(row[i],col[i]),(row[(i+1)%p],col[i])]]
        else:
            row=['r0'];col=['c0'];factors={'r0':Q(1),'c0':Q(1)};edges=[('r0','c0')]
        mapping={}
        if assigned:
            # Alternate row/column cut vertices and force internal local
            # attachments rather than testing shared endpoints alone.
            side='r' if index%2 else 'c'
            cut=next(v for v in reversed(list(assigned)) if v[0]==side)
            local=(row[-1] if side=='r' else col[-1]);mapping[local]=cut
            scale=assigned[cut]/factors[local]
            gauge=scale if side=='r' else 1/scale
        else:gauge=Q(1)
        for local in row+col:
            if local not in mapping:
                mapping[local]=( 'r' if local in row else 'c', next_id);next_id+=1
            value=factors[local]*(gauge if local in row else 1/gauge)
            if mapping[local] in assigned:assert assigned[mapping[local]]==value
            else:assigned[mapping[local]]=value
        for a,b in edges:assert assigned[mapping[a]]*assigned[mapping[b]]==factors[a]*factors[b]
        blocks.append(dict(kind=name,parameter=p,edges=len(edges),products_preserved=True))
    return dict(vertices=len(assigned),blocks=blocks,all_positive=all(v>0 for v in assigned.values()),
                shared_row_and_column_internal_vertices_checked=True)

def run():
    c,e,z,A,B,x,d=s.symbols('c epsilon z alpha beta x d')
    K=e*e+2*e;P=c*z*z+(1-c)*z+K
    F=(e+A)**2*(1-B)-(e+B)**2*(1-A)
    assert s.expand(F-(A-B)*(A+B-A*B+K))==0
    assert s.expand((A+B-A*B+K).subs({A:-c*z,B:z})-P)==0
    re=(c-1)/(2*c);re2=((c-1)**2-2*c*K)/(2*c*c)
    gaps=dict(beta_outer=s.factor(re2+2*e*re),
              alpha_outer=s.factor(c*c*re2-2*e*c*re),
              beta_middle=s.factor(re2-2*re),
              alpha_middle=s.factor(c*c*re2+2*c*re))
    certificates={}
    # Positive-coefficient expansions prove strict signs for all x,d>=0,
    # hence c=1+x>=1 and epsilon=c+1+d>=c+1, without a finite grid argument.
    expressions={'negative_discriminant':4*c*K-(c-1)**2,
                 'negative_real_beta_square':2*c*K-(c-1)**2}
    for name,gap in gaps.items():expressions['negative_'+name]=-s.fraction(s.cancel(gap))[0]
    for name,expr in expressions.items():
        poly=s.Poly(s.expand(expr.subs({c:1+x,e:2+x+d})),x,d)
        assert all(a>=0 for a in poly.coeffs()) and poly.eval({x:0,d:0})>0,(name,poly)
        certificates[name]=str(poly.as_expr())
    derivative=s.diff(F,A).subs({A:0,B:0})
    assert s.expand(derivative-K)==0
    cases=[]
    for p in [2,3,4,5,8,14,15,32,64,1000]:
        maximum=None;patterns=0
        for j in range(1,p//2+1):
            ratio=Q(p-j,j);eps=Q(p);k=eps*eps+2*eps
            assert (ratio-1)**2-4*ratio*k<0
            r=(ratio-1)/(2*ratio);r2=((ratio-1)**2-2*ratio*k)/(2*ratio**2)
            physical=[(r2+2*eps*r)/eps**2,(ratio**2*r2-2*eps*ratio*r)/eps**2,
                      r2-2*r,ratio**2*r2+2*ratio*r]
            assert max(physical)<0
            maximum=max(physical) if maximum is None else max(maximum,*physical)
            patterns+=1
        cases.append(dict(paths=p,epsilon=p,multiplicity_patterns=patterns,
                          largest_real_relative_edge_energy_change=str(maximum)))
    # The former five-path obstruction disappears under this scaling.
    assert cases[3]['largest_real_relative_edge_energy_change']!='9/32'
    # Small complete affine fibers are computed independently, including
    # every complex branch, with a generic separating linear projection.
    spectra=[]
    for p in [2,3,4]:
        variables=s.symbols('t1:'+str(p));t=(-sum(variables),)+variables
        equations=[s.expand((p+t[i])**2*(1-t[j])-(p+t[j])**2*(1-t[i]))
                   for i,j in combinations(range(p),2)]
        jac=s.Matrix(equations).jacobian(variables)
        assert jac.subs(dict.fromkeys(variables,0)).rank()==p-1
        projection=s.Symbol('projection')
        for offset in range(1,12):
            linear=sum(offset**i*v for i,v in enumerate(variables))
            gb=s.groebner(equations+[projection-linear],*variables,projection,order='lex',domain=s.QQ)
            univariate=[f for f in gb if f.free_symbols<={projection}]
            if len(univariate)==1 and s.degree(univariate[0],projection)==2**p-1:break
        else:raise AssertionError('no separating projection')
        assert len(gb)==p
        for v in variables:
            relations=[f for f in gb if f.has(v) and not any(f.has(w) for w in variables if w!=v)]
            assert len(relations)==1 and s.diff(relations[0],v)==1
        poly=s.Poly(univariate[0],projection)
        assert s.gcd(poly,poly.diff()).degree()==0
        spectra.append(dict(paths=p,epsilon=p,affine_points=poly.degree(),
                            projection_offset=offset,projection_polynomial=str(poly.as_expr())))
    result=dict(complete=True,universal_coefficient_certificates=certificates,
                coefficient_domain='c=1+x, epsilon=c+1+d, x>=0, d>=0',
                cycle_identity_verified=True,true_jacobian_factor=str(K),finite_path_controls=cases,
                scaled_affine_spectra=spectra,block_cut_gluing=glue_control(),
                previous_five_path_obstruction_removed=True,verifier_sha256=sha(Path(__file__)),
                scope='exact symbolic inequalities for all path counts plus finite spectra and rational block gluing; no observer query uses hidden factors')
    (HERE/'scaled_block_validation.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2),flush=True)

if __name__=='__main__':run()
