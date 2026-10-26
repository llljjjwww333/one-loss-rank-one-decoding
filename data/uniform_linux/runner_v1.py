"""One Linux protocol for exact cycle MITM and msolve reduced systems.

Imported modules are shared through fork; every measured job is a fresh
child, serially scheduled and confined to the same available CPU.
"""
from pathlib import Path
import hashlib,json,os,platform,shutil,signal,sys,time,multiprocessing as mp
HERE=Path(__file__).resolve().parent;PAPER=HERE.parents[1];RUN=HERE/'uniform_linux'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def dump(p,x):p.write_text(json.dumps(x,indent=2)+'\n')
def freeze():
    assert not (RUN/'protocol.json').exists()
    for name in ['public','sources','outputs']:(RUN/name).mkdir(parents=True,exist_ok=True)
    old=PAPER/'review3_msolve_v2';p=json.loads((old/'protocol.json').read_text())
    cases=p['cases']
    for case in cases:shutil.copy2(old/'public'/(case['id']+'.json'),RUN/'public')
    for name in ['exact_cactus_decoder.py','polynomial_solver_baseline.py','msolve_decoder_worker.py','msolve_runtime.py']:
        shutil.copy2(PAPER/name,RUN/'sources'/name)
    protocol=dict(version=1,cases=[dict(id=c['id'],repeats=3) for c in cases],jobs=6*len(cases),
        methods=['cycle_mitm','msolve_reduced'],seconds_after_ready=30,execution='serial fresh fork children after common imports',
        includes='public IO, construction, exact root isolation, energy matching or native msolve startup, parametrization replay, support checks, serialization',
        excludes='one shared import stage, recorded separately; no Windows timing pooled',single_available_CPU=True,
        source_hashes={p.name:sha(p) for p in (RUN/'sources').glob('*.py')},
        public_hashes={p.name:sha(p) for p in (RUN/'public').glob('*.json')},runner_sha256=sha(Path(__file__)))
    dump(RUN/'protocol.json',protocol);print('frozen',protocol['jobs'],flush=True)
def child(conn,method,public_path,output):
    os.setsid();import sympy as s
    s.core.cache.clear_cache();conn.send('READY');started=time.perf_counter()
    try:
        public=json.loads(Path(public_path).read_text())
        if method=='cycle_mitm':
            import exact_cactus_decoder as cd
            out=cd.decode(public,bounded=False,search='mitm')
        else:
            import msolve_decoder_worker as mw
            out=mw.decode(public,Path(output))
        dump(Path(output),out)
        conn.send(dict(outcome='complete' if out.get('complete') else out.get('outcome','unresolved'),
                       seconds=time.perf_counter()-started,arrays=len(out.get('arrays',[])),output_sha256=sha(Path(output))))
    except BaseException as ex:
        conn.send(dict(outcome='error',seconds=time.perf_counter()-started,error=repr(ex)))
    finally:conn.close()
def run():
    assert os.name=='posix'
    protocol=json.loads((RUN/'protocol.json').read_text());assert not (RUN/'results.json').exists()
    for name,d in protocol['source_hashes'].items():assert sha(RUN/'sources'/name)==d
    for name,d in protocol['public_hashes'].items():assert sha(RUN/'public'/name)==d
    os.environ.update(OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1')
    allowed=sorted(os.sched_getaffinity(0));os.sched_setaffinity(0,{allowed[0]})
    sys.path.insert(0,str(RUN/'sources'));start=time.perf_counter()
    import sympy as s,exact_cactus_decoder,msolve_decoder_worker,msolve_runtime
    imports=time.perf_counter()-start
    runtime=dict(python=sys.version,sympy=s.__version__,platform=platform.platform(),cpu=next(x for x in Path('/proc/cpuinfo').read_text().splitlines() if x.startswith('model name')),
                 affinity=sorted(os.sched_getaffinity(0)),imports_seconds=imports,msolve_sha256=sha(msolve_runtime.binary()),
                 native_manifest_sha256=sha(msolve_runtime.LOCAL/'runtime-manifest.json'),parent_pid=os.getpid())
    dump(RUN/'runtime.json',runtime);records=[];ctx=mp.get_context('fork')
    for rep in range(1,4):
        for case in protocol['cases']:
            # Alternating method order balances serial order across repetitions.
            for method in (protocol['methods'] if rep%2 else list(reversed(protocol['methods']))):
                job=f"{case['id']}--{method}--{rep}";output=RUN/'outputs'/(job+'.json')
                parent,conn=ctx.Pipe();proc=ctx.Process(target=child,args=(conn,method,RUN/'public'/(case['id']+'.json'),output))
                load=Path('/proc/loadavg').read_text().strip();start=time.perf_counter();proc.start();conn.close()
                record=dict(id=case['id'],method=method,repeat=rep,load_before=load)
                if parent.poll(60) and parent.recv()=='READY':
                    record['child_startup_seconds']=time.perf_counter()-start
                    if parent.poll(protocol['seconds_after_ready']):record.update(parent.recv())
                    else:record.update(outcome='timeout',censor_seconds=protocol['seconds_after_ready'])
                else:record.update(outcome='startup_error')
                if proc.is_alive() and record['outcome'] in ['timeout','startup_error']:os.killpg(proc.pid,signal.SIGKILL)
                proc.join(5);parent.close();records.append(record)
                dump(RUN/'results.json',dict(complete=len(records)==protocol['jobs'],protocol_sha256=sha(RUN/'protocol.json'),records=records))
                print(json.dumps(record),flush=True)
if __name__=='__main__':
    {'freeze':freeze,'run':run}[sys.argv[1]]()
