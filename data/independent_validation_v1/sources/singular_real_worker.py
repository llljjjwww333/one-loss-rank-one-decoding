"""Independent zero-dimensional real-fiber decoding via a separating form.

Singular supplies radical lex bases and rational univariate factors.
Standard-library code isolates every real root, checks actualness and
decides the algebraic radical-sum budget with an integer-norm gap.
"""
import argparse
from collections import deque
from fractions import Fraction as Q
import json
import os
from pathlib import Path
import subprocess
from time import perf_counter
from singular_fiber_worker import equations,parse_basis,sha,direct_reply
from independent_real_algebra import (trim,add,mul,divrem,monic,derivative,pgcd,
                                      isolate,budget)


def source(public,gamma):
    equations_list,cycles=equations(public); e=max(1,len(public['edges']))
    result='LIB "primdec.lib";\n'
    result+=f'ring r=0,(c(1..{e})),dp;\noption(redSB);\n'
    result+='ideal I='+',\n'.join(equations_list)+';\nideal B=std(I);\n'
    result+='print("RAW_DIM|"+string(dim(B)));\nprint("RAW_VDIM|"+string(vdim(B)));\n'
    result+='ideal R=std(radical(I));\nint pointCount=vdim(R);\n'
    result+=f'ring s=0,(c(1..{e}),zz),lp;\nideal J=imap(r,R);\n'
    linear='+'.join(f'({gamma**i})*c({i+1})' for i in range(e))
    result+=f'J[size(J)+1]=zz-({linear});\nideal G=std(J);\n'
    result+='''print("POINT_COUNT|"+string(pointCount));
print("BEGIN_BASIS");
print("DIM|"+string(dim(G)));
print("VDIM|"+string(vdim(G)));
int i; poly p;
for(i=1;i<=size(G);i++){
 print("POLY"); p=G[i];
 while(p!=0){ print("TERM|"+string(leadcoef(p))+"|"+string(leadexp(p))); p=p-lead(p); }
}
print("END_BASIS");
list FF=factorize(G[1]);
ideal factors=FF[1]; intvec powers=FF[2];
print("BEGIN_FACTORS");
for(i=1;i<=size(factors);i++){
 if(deg(factors[i])>0){
  print("FACTOR|"+string(powers[i])); p=factors[i];
  while(p!=0){ print("TERM|"+string(leadcoef(p))+"|"+string(leadexp(p))); p=p-lead(p); }
 }
}
print("END_FACTORS");
quit;
'''
    return result,cycles


def univariate(poly,e):
    out=[]
    for coefficient,powers in poly:
        if any(powers[:e]): return None
        while len(out)<=powers[-1]: out.append(Q(0))
        out[powers[-1]]+=Q(coefficient)
    return trim(out)


def shape(basis,e,points):
    pure=[p for poly in basis if (p:=univariate(poly,e)) is not None and len(p)>1]
    if len(pure)!=1: return None
    g=monic(pure[0])
    if len(g)-1!=points or len(pgcd(g,derivative(g)))!=1: return None
    coordinates=[None]*e
    for poly in basis:
        polynomial=[]; linear=[]
        for coefficient,powers in poly:
            if any(powers[:e]):
                if sum(powers[:e])!=1 or powers[-1]: break
                linear.append((powers[:e].index(1),Q(coefficient)))
            else:
                while len(polynomial)<=powers[-1]: polynomial.append(Q(0))
                polynomial[powers[-1]]+=Q(coefficient)
        else:
            if len(linear)==1:
                index,coefficient=linear[0]
                assert coordinates[index] is None
                coordinates[index]=divrem([-x/coefficient for x in polynomial],g)[1]
    if any(p is None for p in coordinates): return None
    return g,coordinates


def factors_from_output(output,e,g):
    data=output.split('BEGIN_FACTORS\n',1)[1].split('END_FACTORS',1)[0]
    polys=[]
    for line in data.strip().splitlines():
        if line.startswith('FACTOR|'):
            assert int(line.split('|')[1])==1
            polys.append([])
        else:
            tag,coefficient,powers=line.split('|'); assert tag=='TERM'
            vector=list(map(int,powers.split(','))); assert len(vector)==e+1 and not any(vector[:e])
            while len(polys[-1])<=vector[-1]: polys[-1].append(Q(0))
            polys[-1][vector[-1]]+=Q(coefficient)
    product=[Q(1)]
    for p in polys: product=mul(product,monic(p))
    assert monic(product)==g
    return [monic(p) for p in polys]


def actual_filter(public,values,bounded):
    m=public['m']; edges=public['edges']; adjacency={}
    for (i,j),value in zip(edges,values):
        if value:
            adjacency.setdefault(i,[]).append((m+j,value))
            adjacency.setdefault(m+j,[]).append((i,value))
    for (i,j),value in zip(edges,values):
        if not value and i in adjacency and m+j in adjacency:
            return False,dict(reason='zero_edge_active_endpoints')
    factors={}; radicands=[]
    for first in sorted(adjacency):
        if first in factors: continue
        factors[first]=values[0].coerce(1); queue=deque([first]); vertices=[]
        while queue:
            v=queue.popleft(); vertices.append(v)
            for w,value in adjacency[v]:
                predicted=value/factors[v]
                if w in factors:
                    if factors[w]!=predicted: return False,dict(reason='inconsistent_factor_propagation')
                else: factors[w]=predicted; queue.append(w)
        a=sum(factors[v]**2 for v in vertices if v<m)
        b=sum(factors[v]**2 for v in vertices if v>=m)
        radicands.append(a*b)
    if not bounded: return True,dict(kind='unbounded_support',components=len(radicands))
    return budget(radicands)


def basis_at_point(basis,values,parameter):
    all_values=values+[parameter]
    for poly in basis:
        total=parameter.coerce(0)
        for coefficient,powers in poly:
            term=parameter.coerce(Q(coefficient))
            for value,power in zip(all_values,powers): term=term*value**power
            total=total+term
        assert total==0


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('public',type=Path); parser.add_argument('output',type=Path)
    parser.add_argument('--runtime',type=Path,required=True); parser.add_argument('--timeout',type=int,default=45)
    parser.add_argument('--bounded',action='store_true'); args=parser.parse_args()
    public=json.loads(args.public.read_text()); runtime=args.runtime
    binary=runtime/'runtime/usr/bin/Singular'
    manifest=json.loads((runtime/'runtime-manifest.json').read_text()); assert sha(binary)==manifest['binary_sha256']
    module=runtime/'runtime/usr/libexec/x86_64-linux-gnu/singular/MOD'
    env=os.environ.copy(); env['LD_LIBRARY_PATH']=str(runtime/'runtime/usr/lib/x86_64-linux-gnu')
    env['SINGULARPATH']=str(runtime/'runtime/usr/share/singular/LIB')+':'+str(module)
    env['SINGULAR_PROCS_DIR']=str(module)
    result=dict(complete=False,uses_sympy=False,public_sha256=sha(args.public),bounded=args.bounded,
                backend_binary_sha256=sha(binary),attempts=[])
    start=perf_counter(); cas_seconds=0
    try:
        gamma=2; attempts_bound=None
        while attempts_bound is None or gamma<2+attempts_bound:
            text,cycles=source(public,gamma)
            script=args.output.with_name(args.output.stem+f'--gamma{gamma}.sing'); script.write_text(text)
            before=perf_counter()
            process=subprocess.run([str(binary),'-q','--no-tty',str(script)],env=env,
                                   capture_output=True,text=True,timeout=max(.001,args.timeout-cas_seconds))
            cas_seconds+=perf_counter()-before
            script.with_suffix('.stdout.txt').write_text(process.stdout); script.with_suffix('.stderr.txt').write_text(process.stderr)
            assert process.returncode==0,process.stderr
            e=max(1,len(public['edges'])); dim,vdim,basis=parse_basis(process.stdout,e+1)
            points=int(next(line.split('|')[1] for line in process.stdout.splitlines() if line.startswith('POINT_COUNT|')))
            if dim==-1 and points==0:
                result.update(complete=True,outcome='complete_list',status='Empty',complex_points=0,
                              candidates=[],decisions=[],wall_seconds=perf_counter()-start,singular_seconds=cas_seconds)
                args.output.write_text(json.dumps(result,indent=2)+'\n')
                print('complete_list',0,0,flush=True)
                raise SystemExit(0)
            assert dim==0 and points>0
            attempts_bound=(e-1)*points*(points-1)//2+1
            representation=shape(basis,e,points)
            result['attempts'].append(dict(gamma=gamma,script=script.name,script_sha256=sha(script),
                 stdout_sha256=sha(script.with_suffix('.stdout.txt')),points=points,dimension=dim,
                 quotient_dimension=vdim,shape_verified=representation is not None))
            if representation is not None: break
            gamma+=1
        assert representation is not None, 'No separating form within the deterministic bound'
        g,coordinates=representation; factor_polys=factors_from_output(process.stdout,e,g)
        result.update(gamma=gamma,separation_trials_bound=attempts_bound,complex_points=points,
            eliminator=list(map(str,g)),coordinate_polynomials=[list(map(str,p)) for p in coordinates],
            factor_polynomials=[list(map(str,p)) for p in factor_polys],basis=basis,decisions=[],candidates=[])
        for factor_index,f in enumerate(factor_polys):
            for root in isolate(f):
                parameter=root.value([0,1]); all_values=[root.value(p) for p in coordinates]
                values=all_values[:len(public['edges'])]
                assert sum((gamma**i)*v for i,v in enumerate(all_values))==parameter
                basis_at_point(basis,all_values,parameter); direct_reply(public,values)
                accepted,certificate=actual_filter(public,values,args.bounded)
                decision=dict(factor_index=factor_index,root=root.json(),accepted=accepted,
                              filter_certificate=certificate,values=[v.json() for v in values])
                result['decisions'].append(decision)
                if accepted: result['candidates'].append(decision)
        result.update(complete=True,outcome='complete_list',
                      status='Empty' if not result['candidates'] else 'Unique' if len(result['candidates'])==1 else 'Ambiguous')
    except subprocess.TimeoutExpired as error:
        cas_seconds+=perf_counter()-before
        for stream in ['stdout','stderr']:
            value=getattr(error,stream) or b''
            if isinstance(value,bytes): value=value.decode(errors='replace')
            args.output.with_suffix('.timeout.'+stream+'.txt').write_text(value)
        result.update(outcome='timeout')
    except Exception as error: result.update(outcome='error',error=repr(error))
    result.update(wall_seconds=perf_counter()-start,singular_seconds=cas_seconds)
    args.output.write_text(json.dumps(result,indent=2)+'\n')
    print(result['outcome'],len(result.get('decisions',[])),len(result.get('candidates',[])),flush=True)
    raise SystemExit(result['outcome']=='error')


if __name__=='__main__': main()
