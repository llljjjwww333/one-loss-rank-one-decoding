"""Independent observed-entry polynomial fiber over QQ; no SymPy imports.

Singular computes a reduced lex basis. Unit ideals and rational singleton
ideals give complete actual lists after independent support/budget tests.
Other bases are retained as algebraic-fiber checks, not full real decoding.
"""
import argparse
from collections import deque
from fractions import Fraction as Q
import hashlib
import json
from math import isqrt
import os
from pathlib import Path
import subprocess
from time import perf_counter


def sha(path): return hashlib.sha256(path.read_bytes()).hexdigest()


def cycle_edges(m,n,edges):
    # All simple even cycles; no shared incidence elimination or SymPy.
    adj={v:[] for v in range(m+n)}
    for k,(i,j) in enumerate(edges):adj[i].append((m+j,k));adj[m+j].append((i,k))
    cycles={}
    for start in sorted(adj):
        def dfs(vertex,seen,path):
            for other,k in adj[vertex]:
                if other==start and len(path)>=3:
                    c=path+[k];key=tuple(sorted(c))
                    assert len(c)%2==0
                    cycles.setdefault(key,c)
                elif other>start and other not in seen:
                    dfs(other,seen|{other},path+[k])
        dfs(start,{start},[])
    return [cycles[k] for k in sorted(cycles)]


def equations(public):
    m,n=public['m'],public['n']; edges=public['edges']; e=len(edges)
    u,v,gu,gv=([Q(x) for x in public[k]] for k in ['u','v','gu','gv'])
    loss=Q(public['loss']); assert len(u)==len(gu)==m and len(v)==len(gv)==n
    weights=[u[i]*v[j] for i,j in edges]
    row=[u[i]*(sum(u[i]*v[j]**2 for ii,j in edges if ii==i)-gu[i]) for i in range(m)]
    col=[v[j]*(sum(v[j]*u[i]**2 for i,jj in edges if jj==j)-gv[j]) for j in range(n)]
    norm=2*loss+sum(w*w for w in weights)-sum(x*y for x,y in zip(u,gu))-sum(x*y for x,y in zip(v,gv))
    result=[]
    for i in range(m):
        result.append('+'.join(f'({weights[k]})*c({k+1})' for k,(ii,j) in enumerate(edges) if ii==i) or '0')
        result[-1]+=f'-({row[i]})'
    for j in range(n):
        result.append('+'.join(f'({weights[k]})*c({k+1})' for k,(i,jj) in enumerate(edges) if jj==j) or '0')
        result[-1]+=f'-({col[j]})'
    cycles=cycle_edges(m,n,edges)
    for cycle in cycles:
        result.append('*'.join(f'c({k+1})' for k in cycle[::2])+'-'+ '*'.join(f'c({k+1})' for k in cycle[1::2]))
    result.append(('+'.join(f'c({k+1})^2' for k in range(e)) or '0')+f'-({norm})')
    if not e: result.append('c(1)')
    return result,cycles


def singular_source(public):
    eqs,cycles=equations(public)
    source=f'ring r=0,(c(1..{max(1,len(public["edges"]))})),lp;\noption(redSB);\n'
    source+='ideal I='+',\n'.join(eqs)+';\nideal G=std(I);\n'
    source+='''print("BEGIN_BASIS");
print("DIM|"+string(dim(G)));
print("VDIM|"+string(vdim(G)));
int i; poly p;
for(i=1;i<=size(G);i++){
 print("POLY"); p=G[i];
 while(p!=0){
  print("TERM|"+string(leadcoef(p))+"|"+string(leadexp(p)));
  p=p-lead(p);
 }
}
print("END_BASIS");
quit;
'''
    return source,cycles


def parse_basis(output,variables):
    assert 'BEGIN_BASIS' in output and 'END_BASIS' in output
    assert '?' not in output and 'Could not find' not in output,output[:1000]
    data=output.split('BEGIN_BASIS\n',1)[1].split('END_BASIS',1)[0]
    basis=[]; dimension=vdim=None
    for line in data.strip().splitlines():
        if line.startswith('DIM|'): dimension=int(line.split('|')[1])
        elif line.startswith('VDIM|'): vdim=int(line.split('|')[1])
        elif line=='POLY': basis.append([])
        else:
            tag,coef,powers=line.split('|'); assert tag=='TERM'
            exponents=list(map(int,powers.split(','))); assert len(exponents)==variables
            basis[-1].append([str(Q(coef)),exponents])
    assert dimension is not None and vdim is not None
    return dimension,vdim,basis


def linear_point(basis,n):
    rows=[]
    for poly in basis:
        row=[Q(0)]*(n+1)
        for coef,powers in poly:
            degree=sum(powers); assert degree<=1
            if degree: row[powers.index(1)]+=Q(coef)
            else: row[-1]-=Q(coef)
        rows.append(row)
    rank=0
    for j in range(n):
        found=next((i for i in range(rank,len(rows)) if rows[i][j]),None)
        assert found is not None
        rows[rank],rows[found]=rows[found],rows[rank]
        pivot=rows[rank][j]; rows[rank]=[x/pivot for x in rows[rank]]
        for i in range(len(rows)):
            if i!=rank:
                a=rows[i][j]; rows[i]=[x-a*y for x,y in zip(rows[i],rows[rank])]
        rank+=1
    assert all(any(row[:-1]) or row[-1]==0 for row in rows)
    return [rows[i][-1] for i in range(n)]


def radical_budget(radicands):
    """Exact sign of sum(sqrt(nonnegative rationals))-1 via integer norm.

    With T the product of radicand denominators, T*S is integral of degree
    <=2^h. The explicit norm gap bounds every possible nonzero S, so a
    sufficiently narrow interval containing zero certifies equality.
    """
    if not radicands: return True,dict(kind='empty',sign=-1)
    assert all(x>0 for x in radicands)
    scale=1
    for x in radicands: scale*=x.denominator
    def sqrt_bounds(x,bits):
        mult=1<<bits; k=isqrt((x.numerator*mult*mult)//x.denominator)
        low=Q(k,mult)
        return low,low if low*low==x else Q(k+1,mult)
    bound=scale*(1+sum(sqrt_bounds(x,0)[1] for x in radicands))
    assert bound.denominator==1
    precision=scale.bit_length()+((1<<len(radicands))-1)*int(bound).bit_length()
    gap=Q(1,1<<precision)
    bits=8
    while True:
        intervals=[sqrt_bounds(x,bits) for x in radicands]
        low=sum(a for a,b in intervals)-1; high=sum(b for a,b in intervals)-1
        if high<0: sign=-1
        elif low>0: sign=1
        elif high-low<gap: sign=0
        else: bits*=2; continue
        return sign<=0,dict(kind='rational_radical_integer_norm_gap',sign=sign,
                           components=len(radicands),precision_bits=bits,nonzero_gap_bits=precision,
                           interval=[str(low),str(high)])


def actual_and_budget(public,values,bounded):
    m,n=public['m'],public['n']; edges=public['edges']; adjacency={}
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
        factors[first]=Q(1); queue=deque([first]); vertices=[]
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
    accepted,certificate=radical_budget(radicands)
    certificate['radicands']=list(map(str,radicands))
    return accepted,certificate


def direct_reply(public,values):
    m,n=public['m'],public['n']; u,v,gu,gv=([Q(x) for x in public[k]] for k in ['u','v','gu','gv'])
    a=[Q(0)]*m; b=[Q(0)]*n; loss=Q(0)
    for (i,j),value in zip(public['edges'],values):
        residual=u[i]*v[j]-value
        a[i]+=residual*v[j]; b[j]+=residual*u[i]; loss+=residual*residual/2
    assert a==gu and b==gv and loss==Q(public['loss'])


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('public',type=Path); parser.add_argument('output',type=Path)
    parser.add_argument('--runtime',required=True,type=Path)
    parser.add_argument('--timeout',type=int,default=45)
    parser.add_argument('--bounded',action='store_true')
    args=parser.parse_args(); start=perf_counter()
    public=json.loads(args.public.read_text()); source,cycles=singular_source(public)
    script=args.output.with_suffix('.sing'); script.write_text(source)
    runtime=args.runtime; binary=runtime/'runtime/usr/bin/Singular'
    manifest=json.loads((runtime/'runtime-manifest.json').read_text())
    assert sha(binary)==manifest['binary_sha256']
    module=runtime/'runtime/usr/libexec/x86_64-linux-gnu/singular/MOD'
    env=os.environ.copy(); env['LD_LIBRARY_PATH']=str(runtime/'runtime/usr/lib/x86_64-linux-gnu')
    env['SINGULARPATH']=str(runtime/'runtime/usr/share/singular/LIB')+':'+str(module)
    env['SINGULAR_PROCS_DIR']=str(module)
    base=dict(public_sha256=sha(args.public),uses_sympy=False,bounded=args.bounded,
              model='All observed entries; independent weighted margins, cycle binomials and norm',
              variables=max(1,len(public['edges'])),cycle_lengths=list(map(len,cycles)),
              script_sha256=sha(script),backend_binary_sha256=sha(binary))
    try:
        process=subprocess.run([str(binary),'-q','--no-tty',str(script)],env=env,capture_output=True,text=True,timeout=args.timeout)
        args.output.with_suffix('.stdout.txt').write_text(process.stdout)
        args.output.with_suffix('.stderr.txt').write_text(process.stderr)
        assert process.returncode==0,process.stderr
        dim,vdim,basis=parse_basis(process.stdout,base['variables'])
        base.update(dimension=dim,quotient_dimension=vdim,basis=basis,complete_actual_list=False)
        if dim<0:
            base.update(outcome='complete_list',complete_actual_list=True,status='Empty',arrays=[])
        elif dim==0 and vdim==1:
            values=linear_point(basis,base['variables'])[:len(public['edges'])]
            direct_reply(public,values)
            accepted,certificate=actual_and_budget(public,values,args.bounded)
            base.update(outcome='complete_list',complete_actual_list=True,status='Unique' if accepted else 'Empty',
                        rational_fiber_point=list(map(str,values)),filter_certificate=certificate,
                        arrays=[list(map(str,values))] if accepted else [])
        else:
            base.update(outcome='algebraic_fiber_only' if dim==0 else 'positive_dimensional',
                        status='No independent full real-list claim',arrays=None)
    except subprocess.TimeoutExpired as error:
        for stream in ['stdout','stderr']:
            value=getattr(error,stream) or ''
            if isinstance(value,bytes): value=value.decode(errors='replace')
            args.output.with_suffix('.'+stream+'.txt').write_text(value)
        base.update(outcome='timeout',complete_actual_list=False)
    except Exception as error:
        base.update(outcome='error',complete_actual_list=False,error=repr(error))
    base['wall_seconds']=perf_counter()-start
    args.output.write_text(json.dumps(base,indent=2)+'\n')
    print(base['outcome'],round(base['wall_seconds'],3),flush=True)
    raise SystemExit(base['outcome']=='error')


if __name__=='__main__': main()
