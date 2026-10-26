"""Instrument the unchanged coordinate decoder on all 96 prospective inputs."""
from pathlib import Path
from collections import Counter
import hashlib,json,sys,os,time,multiprocessing as mp,signal
HERE=Path(__file__).resolve().parent;ORIGINAL=HERE.parent/'energy_constraints_v1'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def put(p,x):p.write_text(json.dumps(x,indent=2)+'\n')
def freeze():
    assert not (HERE/'phase_protocol.json').exists()
    put(HERE/'phase_protocol.json',dict(inputs=96,repetitions=3,jobs=288,seconds_after_ready=30,
      phases=['reduction','real_roots','coordinates','propagate','rational_candidate'],
      timing='warm shared imports, fresh fork children, clear SymPy cache; time decode only, public parse and output-file writing excluded',
      affinity='same first available CPU as the solver protocol; execute only after that protocol finishes',
      decoder_sha256=sha(ORIGINAL/'coordinate_decoder.py'),inputs_sha256=sha(ORIGINAL/'prospective_inputs.jsonl'),
      reference_sha256=sha(ORIGINAL/'prospective_outcomes.jsonl'),runner_sha256=sha(Path(__file__))))
def child(conn,public):
    os.setsid();import coordinate_decoder as d, sympy
    sympy.core.cache.clear_cache();phases={};calls={}
    for name in ['reduction','real_roots','coordinates','propagate','rational_candidate']:
        fn=getattr(d,name)
        def wrapper(*args,fn=fn,name=name,**kwargs):
            start=time.perf_counter()
            try:return fn(*args,**kwargs)
            finally:phases[name]=phases.get(name,0)+time.perf_counter()-start;calls[name]=calls.get(name,0)+1
        setattr(d,name,wrapper)
    conn.send('READY');start=time.perf_counter()
    try:
        result=d.decode(public);elapsed=time.perf_counter()-start
        conn.send(dict(outcome='complete',seconds=elapsed,phases=phases,calls=calls,result=result))
    except BaseException as e:conn.send(dict(outcome='error',error=repr(e)))
    finally:conn.close()
def run():
    assert os.name=='posix' and not (HERE/'phase_results.json').exists()
    protocol=json.loads((HERE/'phase_protocol.json').read_text());assert sha(Path(__file__))==protocol['runner_sha256']
    assert json.loads((HERE/'results.json').read_text())['complete']
    assert sha(ORIGINAL/'coordinate_decoder.py')==protocol['decoder_sha256']
    inputs=[json.loads(s) for s in (ORIGINAL/'prospective_inputs.jsonl').read_text().splitlines()]
    references={r['id']:r['result'] for r in map(json.loads,(ORIGINAL/'prospective_outcomes.jsonl').read_text().splitlines())}
    assert len(inputs)==96
    os.environ.update(OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1')
    os.sched_setaffinity(0,{min(os.sched_getaffinity(0))});sys.path.insert(0,str(ORIGINAL));import coordinate_decoder
    ctx=mp.get_context('fork');records=[]
    for rep in range(1,4):
        for c in inputs:
            parent,other=ctx.Pipe();proc=ctx.Process(target=child,args=(other,c['public']));start=time.perf_counter();proc.start();other.close()
            row=dict(id=c['id'],r=c['r'],N=c['N'],repeat=rep)
            if parent.poll(30) and parent.recv()=='READY':
                row['startup_seconds']=time.perf_counter()-start
                if parent.poll(30):row.update(parent.recv())
                else:row.update(outcome='timeout',censor_seconds=30)
            else:row.update(outcome='startup_error')
            if proc.is_alive() and row['outcome'] in ['timeout','startup_error']:os.killpg(proc.pid,signal.SIGKILL)
            proc.join(5);parent.close()
            if row['outcome']=='complete':assert row['result']==references[c['id']],c['id']
            records.append(row);put(HERE/'phase_results.json',dict(complete=len(records)==288,protocol_sha256=sha(HERE/'phase_protocol.json'),records=records))
            if len(records)%16==0:print('PHASE',len(records),'of 288',flush=True)
    assert all(r['outcome']=='complete' for r in records)
    print('PHASE DONE 288; unchanged outcomes throughout',flush=True)
if __name__=='__main__':{'freeze':freeze,'run':run}[sys.argv[1]]()
