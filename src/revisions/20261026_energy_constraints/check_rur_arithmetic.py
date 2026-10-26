"""Exact RUR arithmetic controls and deliberate certificate corruptions."""
from copy import deepcopy
from pathlib import Path
import hashlib,json
import sympy as s
from exact_rur_decoder import rational_univariate,verify_rur,polynomial
HERE=Path(__file__).resolve().parent
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def run():
    x,y=s.symbols('x y');cases=[
        ('colliding_projection',[y+x,x**3-x],(x,y),3,3,3),
        ('nonradical_algebraic',[(x*x-2)**2,y-x],(x,y),4,2,2),
        ('nongorenstein_local',[x*x,x*y,y*y],(x,y),3,1,1),
        ('rational_nonradical',[(x-s.Rational(1,3))**2,y-x],(x,y),2,1,1),
        ('nonsplit_quintic',[x**5-x-1],(x,),5,5,1)]
    records=[]
    for name,equations,variables,dimension,degree,real in cases:
        cert=rational_univariate(equations,variables,64);verify_rur(equations,variables,cert,64)
        assert cert['quotient_dimension']==dimension and cert['distinct_complex_points']==degree
        assert len(polynomial(cert['polynomial']).real_roots(radicals=False))==real
        if name=='colliding_projection':
            assert len(cert['trials'])==2 and cert['trials'][0]['separating'] is False
        records.append(dict(id=name,quotient_dimension=dimension,complex_points=degree,real_points=real,
                            projection_trials=cert['trials'],certificate=cert))
    # A zero mean can itself be a solution of a collapsed fiber. Checking
    # equations at the mean is insufficient; the nilpotence gate rejects it.
    collapsed=dict(empty=False,quotient_dimension=3,distinct_complex_points=1,
                   standard_monomials=[[0,0],[0,1],[0,2]],nilpotence_exponent=3,
                   polynomial=['1','0'],coordinates=[['0'],['0']],projection_weights=['1','1'])
    collapsed['standard_monomials']=records[0]['certificate']['standard_monomials']
    rejected=[]
    for name,eqs,vars_,cert in [
        ('collapsed_projection_mean',[y+x,x**3-x],(x,y),collapsed),
        ('deleted_real_root',[x**3-x],(x,),dict(rational_univariate([x**3-x],(x,),64),
                                              polynomial=['1','0','-1'],distinct_complex_points=2))]:
        try:verify_rur(eqs,vars_,cert,64)
        except AssertionError:rejected.append(name)
        else:raise AssertionError('corrupted certificate accepted: '+name)
    assert len(rejected)==2
    out=dict(complete=True,controls=records,rejected_corruptions=rejected,
             verifier_sha256=sha(Path(__file__)),solver_sha256=sha(HERE/'exact_rur_decoder.py'),
             scope='exact synthetic quotient-algebra controls; separate from native loss/gradient experiments')
    (HERE/'rur_arithmetic_validation.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({k:v for k,v in out.items() if k!='controls'}))
if __name__=='__main__':run()
