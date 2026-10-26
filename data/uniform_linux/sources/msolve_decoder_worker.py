"""Exact public-input msolve rational-parametrization route, with QQ replay."""
from pathlib import Path
from time import perf_counter
import argparse,ast,json,subprocess
import sympy as s
import exact_cactus_decoder as cd
import polynomial_solver_baseline as pb
import msolve_runtime as runtime


def parse_output(raw):
    tree=ast.parse(raw.strip().rstrip(':').replace('^','**'),mode='eval')
    def visit(node):
        if isinstance(node,ast.Expression):return visit(node.body)
        if isinstance(node,ast.List):return [visit(x) for x in node.elts]
        if isinstance(node,ast.Constant) and type(node.value) in (int,str):return node.value
        if isinstance(node,ast.UnaryOp) and isinstance(node.op,(ast.USub,ast.UAdd)):
            return -visit(node.operand) if isinstance(node.op,ast.USub) else visit(node.operand)
        if isinstance(node,ast.BinOp) and isinstance(node.op,(ast.Div,ast.Pow)):
            a,b=visit(node.left),visit(node.right)
            return s.Rational(a,b) if isinstance(node.op,ast.Div) else a**b
        raise ValueError('Unsupported msolve syntax: '+ast.dump(node))
    return visit(tree)


def roots_from_output(raw,variables,equations):
    data=parse_output(raw)
    if data==[-1]:return [],dict(dimension=-1)
    # With -P 2, msolve 0.6.5 emits [0, [dimension,...]].
    assert data[0]==0 and len(data)==2,data[:2]
    dim,nvars,degree,names,form,block=data[1]
    if dim!=0:return None,dict(dimension=dim)
    assert block[0]==1 and len(block)==2
    wdata,ddata,params=block[1]
    t=s.Symbol('_T')
    def polynomial(encoded):
        deg,coeff=encoded;assert deg+1==len(coeff)
        return s.Poly(sum(s.Rational(c)*t**j for j,c in enumerate(coeff)),t,domain=s.QQ)
    w,den=polynomial(wdata),polynomial(ddata)
    assert w.degree()==degree and s.gcd(w,den).degree()==0
    assert len(params)==len(names)-1 and nvars==len(names)
    functions=[-polynomial(p).as_expr()/s.Rational(c)/den.as_expr() for p,c in params]+[t]
    mapping={name:expr for name,expr in zip(names,functions)}
    assert {str(v) for v in variables}.issubset(mapping)
    coord=[mapping[str(v)] for v in variables]
    # Replay every input equation identically in QQ[T]/(w), before any
    # root substitution. No approximate interval endpoint is a candidate.
    for eq in equations:
        numerator=s.cancel(eq.subs(dict(zip(variables,coord)),simultaneous=True)).as_numer_denom()[0]
        assert s.rem(s.Poly(numerator,t),w).is_zero
    roots=[]
    for root,mult in w.real_roots(multiple=False,radicals=False):
        roots.append([s.cancel(f.subs(t,root)) for f in coord])
    return roots,dict(dimension=dim,parametrization_degree=degree,real_parameter_roots=len(roots),equations_replayed=len(equations))


def decode(public,output):
    p=cd.canonical_public(public)
    cycles=cd.cycle_basis(p['m'],p['n'],p['edges'])
    base=cd.forest_base(p['m'],p['n'],p['edges'],cycles,pb.weighted_margins(p))
    assert base is not None and cycles
    variables=list(s.symbols('t:'+str(len(cycles))));z=list(base);equations=[]
    for t,cycle in zip(variables,cycles):
        for offset,e in enumerate(cycle):z[e]+=(-1)**offset*t
        equations.append(s.prod(z[e] for e in cycle[::2])-s.prod(z[e] for e in cycle[1::2]))
    weights=[p['u'][i]*p['v'][j] for i,j in p['edges']]
    norm=2*p['loss']+sum(w*w for w in weights)-sum(x*y for x,y in zip(p['u'],p['gu']))-sum(x*y for x,y in zip(p['v'],p['gv']))
    arrays_expr=[entry/w for entry,w in zip(z,weights)]
    equations.append(sum(c*c for c in arrays_expr)-norm)
    primitive=[s.Poly(e,*variables,domain=s.QQ).clear_denoms(convert=True)[1].primitive()[1].as_expr() for e in equations if e!=0]
    source=','.join(map(str,variables))+'\n0\n'+',\n'.join(str(e).replace('**','^') for e in primitive)+'\n'
    input_path=output.with_suffix('.msolve.in');native_path=output.with_suffix('.msolve.out')
    input_path.write_text(source)
    start=perf_counter()
    proc=subprocess.run([str(runtime.binary()),'-f',str(input_path),'-o',str(native_path),'-P','2','-t','1'],env=runtime.environment(),capture_output=True,text=True)
    elapsed=perf_counter()-start
    output.with_suffix('.cas-stdout.txt').write_text(proc.stdout);output.with_suffix('.cas-stderr.txt').write_text(proc.stderr)
    assert proc.returncode==0,(proc.returncode,proc.stderr)
    roots,info=roots_from_output(native_path.read_text(),variables,equations)
    if roots is None:return dict(complete=False,outcome='positive_dimension',arrays=[],cas_seconds=elapsed,**info)
    arrays=[];stats=pb.new_stats()
    for values in roots:
        sub=dict(zip(variables,values))
        pb.accept(p,[c.subs(sub) for c in arrays_expr],False,arrays,stats)
    return dict(complete=True,outcome='complete',arrays=[[str(v) for v in a] for a in arrays],
                cas_seconds=elapsed,root_extraction='msolve exact rational parametrization; SymPy real-root isolation and actual-support filtering',**info)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('public',type=Path);parser.add_argument('output',type=Path)
    args=parser.parse_args();s.core.cache.clear_cache();print('READY',flush=True);start=perf_counter()
    try:result=decode(json.loads(args.public.read_text()),args.output)
    except Exception as exc:result=dict(complete=False,outcome='error',error=repr(exc),arrays=[])
    result.update(core_seconds=perf_counter()-start,truth_passed_to_solver=False,bounded=False,method='msolve_0.6.5')
    args.output.write_text(json.dumps(result,indent=2)+'\n');print(result['outcome'],flush=True)
