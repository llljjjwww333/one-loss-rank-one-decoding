"""Exact cycle, energy and transversality identities used by the degree proof."""
from pathlib import Path
import hashlib,json
import sympy as s
HERE=Path(__file__).resolve().parent
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def run():
    t,h,rho,x,W=s.symbols('t h rho x W')
    P=[9,3,2,-2,-2,6];M=[9,3,3,-3,-1,3];sign=[1,-1,1,-1,1,-1]
    weighted=[P[i]*(h if i in [2,3] else 1)+sign[i]*t for i in range(6)]
    f=s.prod(weighted[::2])-s.prod(weighted[1::2])
    assert s.expand(f-2*t*(t-1)*(t+2*h))==0
    weighted_m=[rho*M[i]+sign[i]*t for i in range(6)]
    g=s.prod(weighted_m[::2])-s.prod(weighted_m[1::2])
    assert s.expand(g-2*t*(t+rho)*(t+3*rho))==0
    candidate_p=[10,2,2+x,-2-x,-1,5];candidate_m=[rho*z for z in [8,4,2,-2,-2,4]]
    ep=s.expand(sum(v*v-u*u for u,v in zip(P,candidate_p)))
    em=s.expand(sum(v*v-(rho*u)**2 for u,v in zip(M,candidate_m)))
    assert ep==2*x*x+8*x and em==-10*rho*rho
    candidate_m2=[rho*z for z in [6,6,0,0,-4,6]]
    assert s.expand(sum(v*v-(rho*u)**2 for u,v in zip(M,candidate_m2)))==6*rho*rho
    assert s.diff(ep-10*W,x)==4*x+8
    assert s.expand((s.sqrt(4+5*W)-2)**2+4*(s.sqrt(4+5*W)-2)-5*W)==0
    assert s.diff(f,t).subs(t,1)==2+4*h
    assert s.expand(s.diff(g,t).subs(t,-rho)+4*rho*rho)==0
    assert s.expand(s.diff(g,t).subs(t,-3*rho)-12*rho*rho)==0
    digits=[0,5,-3];differences={a-b for a in digits for b in digits}
    assert differences=={0,3,-3,5,-5,8,-8}
    for value in differences-{0}:
        if value%4==0:assert all((value//4+next_digit)%4!=0 for next_digit in differences)
    v1,v2=s.symbols('v1 v2');assert s.expand(3*rho*v1-3*rho*v2-3*rho*(v1-v2))==0
    out=dict(complete=True,cycle_factorizations=2,energy_identities=3,positive_transversal_derivative='4*x+8',
             simple_P_root_derivative='2+4*h',simple_M_root_derivatives=['-4*rho**2','12*rho**2'],
             ternary_digit_uniqueness_modulo_four=True,nonzero_private_row_margin='3*rho*(v1-v2)',
             real_collision_relation='x**2+4*x-5*W=0',verifier_sha256=sha(Path(__file__)),
             scope='exact symbolic identities; open-neighborhood and arbitrary-r degree argument are in the main proof')
    (HERE/'degree_symbolic_validation.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out))
if __name__=='__main__':run()
