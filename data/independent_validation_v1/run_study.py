"""One serial Linux protocol; native-edge Singular has no SymPy dependency."""
from pathlib import Path
from collections import Counter
import hashlib,json,os,sys,shutil,subprocess,time,signal,platform,multiprocessing as mp
HERE=Path(__file__).resolve().parent;PAPER=HERE.parent;ROOT=PAPER.parents[1]
PREV=PAPER/'revisions/20261026_selfcheck';SOURCES=HERE/'sources'
SINGULAR=ROOT/'tools/third_party/jsc_singular_4.3.2_v2'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def put(p,x):p.write_text(json.dumps(x,indent=2)+'\n')
def freeze():
    assert not (HERE/'protocol.json').exists()
    old=json.loads((PREV/'uniform_linux/protocol.json').read_text())
    rur=json.loads((PREV/'rur_protocol_v3.json').read_text())
    controls=[dict(c,filename=c['id']+'.json') for c in rur['cases']]
    # The previous bounded filenames encode the same rational public input.
    for c in controls:assert sha(HERE/'controls'/c['filename'])==c['public_sha256']
    p=dict(version=1,methods=['cycle_mitm','msolve_reduced','singular_observed'],cases=old['cases'],
      repetitions=3,comparison_jobs=126,controls=controls,control_jobs=len(controls),seconds_after_ready=30,
      execution='serial fresh process groups, one common parent import, fixed CPU affinity, rotated method order',
      timing='from READY through public IO, equation construction, root work, native executable and Python subworker startup, real support tests and serialization',
      exclusions='shared parent import only; child startup is recorded separately; no pooling of earlier timings',
      runner_sha256=sha(Path(__file__)),source_hashes={p.name:sha(p) for p in sorted(SOURCES.glob('*.py'))},
      public_hashes={p.name:sha(p) for p in sorted((HERE/'public').glob('*.json'))},
      control_hashes={p.name:sha(p) for p in sorted((HERE/'controls').glob('*.json'))},
      singular_manifest_sha256=sha(SINGULAR/'runtime-manifest.json'))
    put(HERE/'protocol.json',p);print('FROZEN',p['comparison_jobs'],p['control_jobs'],flush=True)
def native(public,output,bounded):
    command=[sys.executable,str(SOURCES/'singular_real_worker.py'),str(public),str(output),'--runtime',str(SINGULAR),'--timeout','30']
    if bounded:command.append('--bounded')
    done=subprocess.run(command,capture_output=True,text=True)
    output.with_suffix('.worker.stdout.txt').write_text(done.stdout);output.with_suffix('.worker.stderr.txt').write_text(done.stderr)
    if not output.exists():raise RuntimeError(('missing native outcome',done.returncode,done.stderr[-1000:]))
    result=json.loads(output.read_text())
    if done.returncode or result.get('outcome')=='error':raise RuntimeError(result.get('error',done.stderr))
    return result
def child(conn,method,public,output,bounded=False):
    os.setsid();conn.send('READY');start=time.perf_counter()
    try:
        if method=='singular_observed':
            x=native(public,output,bounded);arrays=len(x.get('candidates',[]));status='complete' if x['complete'] else x['outcome']
        else:
            data=json.loads(public.read_text())
            if method=='cycle_mitm':
                import exact_cactus_decoder as cd
                x=cd.decode(data,bounded=False,search='mitm');arrays=len(x.get('candidates',[]))
            else:
                import msolve_decoder_worker as mw
                x=mw.decode(data,output);arrays=len(x.get('arrays',[]))
            put(output,x);status='complete' if x.get('complete') else x.get('outcome','unresolved')
        conn.send(dict(outcome=status,seconds=time.perf_counter()-start,arrays=arrays,output_sha256=sha(output)))
    except BaseException as e:conn.send(dict(outcome='error',seconds=time.perf_counter()-start,error=repr(e)))
    finally:conn.close()
def run():
    assert os.name=='posix' and not (HERE/'results.json').exists()
    protocol=json.loads((HERE/'protocol.json').read_text());assert sha(Path(__file__))==protocol['runner_sha256']
    for group,folder in [('source_hashes',SOURCES),('public_hashes',HERE/'public'),('control_hashes',HERE/'controls')]:
        for name,d in protocol[group].items():assert sha(folder/name)==d
    os.environ.update(OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1')
    allowed=sorted(os.sched_getaffinity(0));os.sched_setaffinity(0,{allowed[0]})
    sys.path.insert(0,str(SOURCES));start=time.perf_counter()
    import sympy,exact_cactus_decoder,msolve_decoder_worker,msolve_runtime
    put(HERE/'runtime.json',dict(python=sys.version,sympy=sympy.__version__,platform=platform.platform(),
      cpu=next(x for x in Path('/proc/cpuinfo').read_text().splitlines() if x.startswith('model name')),
      affinity=sorted(os.sched_getaffinity(0)),parent_import_seconds=time.perf_counter()-start,
      singular_sha256=sha(SINGULAR/'runtime/usr/bin/Singular'),msolve_sha256=sha(msolve_runtime.binary())))
    records=[];ctx=mp.get_context('fork')
    def job(case,method,rep,group,bounded=False):
        identifier=case['id'];public=HERE/('public' if group=='comparison' else 'controls')/(identifier+'.json')
        output=HERE/'outputs'/f'{group}--{identifier}--{method}--{rep}.json'
        parent,other=ctx.Pipe();proc=ctx.Process(target=child,args=(other,method,public,output,bounded))
        start=time.perf_counter();proc.start();other.close();row=dict(group=group,id=identifier,method=method,repeat=rep,bounded=bounded,output=output.name,load_before=Path('/proc/loadavg').read_text().strip())
        if parent.poll(60) and parent.recv()=='READY':
            row['child_startup_seconds']=time.perf_counter()-start
            if parent.poll(30):row.update(parent.recv())
            else:row.update(outcome='timeout',censor_seconds=30)
        else:row.update(outcome='startup_error')
        if proc.is_alive() and row['outcome'] in ['timeout','startup_error']:os.killpg(proc.pid,signal.SIGKILL)
        proc.join(5);parent.close()
        if row['outcome']=='timeout':row['censor_seconds']=30
        records.append(row);put(HERE/'results.json',dict(complete=len(records)==126+len(protocol['controls']),protocol_sha256=sha(HERE/'protocol.json'),records=records))
        print(json.dumps({k:v for k,v in row.items() if k not in ['load_before','output_sha256']}),flush=True)
    # All independent coverage controls run first; timings are a separate task.
    for c in protocol['controls']:job(c,'singular_observed',1,'coverage',c['bounded'])
    for rep in range(1,4):
        methods=protocol['methods'][rep-1:]+protocol['methods'][:rep-1]
        for c in protocol['cases']:
            for method in methods:job(c,method,rep,'comparison')
    print('DONE',dict(Counter((x['group'],x['method'],x['outcome']) for x in records)),flush=True)
if __name__=='__main__':{'freeze':freeze,'run':run}[sys.argv[1]]()
