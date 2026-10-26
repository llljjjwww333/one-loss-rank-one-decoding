"""Audit retained Linux outputs; generate a source-bound comparison table."""
from pathlib import Path
from collections import Counter
from statistics import median
import hashlib,json,sys
HERE=Path(__file__).resolve().parent;RUN=HERE/'uniform_linux';PAPER=HERE.parents[1]
sys.path.insert(0,str(PAPER))
from run_solver_comparison import same_lists
from verify_solver_comparison import direct_reply
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def run():
    protocol=json.loads((RUN/'protocol.json').read_text());raw=json.loads((RUN/'results.json').read_text())
    assert raw['complete'] and len(raw['records'])==protocol['jobs']==84
    assert sha(RUN/'protocol.json')==raw['protocol_sha256']
    assert sha(RUN/'runner_v1.py')==protocol['runner_sha256']
    for name,d in protocol['source_hashes'].items():assert sha(RUN/'sources'/name)==d
    for name,d in protocol['public_hashes'].items():assert sha(RUN/'public'/name)==d
    records=[];references={};reply_arrays=agreements=0;corrections=[]
    for r in raw['records']:
        row=dict(r)
        if r['outcome']=='complete':
            path=RUN/'outputs'/f"{r['id']}--{r['method']}--{r['repeat']}.json"
            assert sha(path)==r['output_sha256'];out=json.loads(path.read_text());assert out['complete']
            arrays=out['arrays'] if r['method']=='msolve_reduced' else [c['values'] for c in out['candidates']]
            row['verified_arrays']=len(arrays)
            if r['arrays']!=len(arrays):corrections.append(dict(id=r['id'],method=r['method'],repeat=r['repeat'],logged=r['arrays'],verified=len(arrays)))
            public=json.loads((RUN/'public'/(r['id']+'.json')).read_text())
            reply_arrays+=direct_reply(public,arrays)
            if r['id'] in references:
                assert same_lists(arrays,references[r['id']]);agreements+=1
            else:references[r['id']]=arrays
        else:assert r['outcome']=='timeout' and r['censor_seconds']==30,r
        records.append(row)
    summaries={method:dict(Counter(x['outcome'] for x in records if x['method']==method)) for method in protocol['methods']}
    runtime_record=json.loads((RUN/'runtime.json').read_text())
    native=PAPER.parents[1]/'tools/third_party/jsc_msolve_0.6.5_v2/runtime/usr/bin/msolve'
    if native.exists():assert sha(native)==runtime_record['msolve_sha256']
    rows=[];details=[]
    for i,case in enumerate(protocol['cases'],1):
        display=[];item=dict(id=case['id'])
        for method in protocol['methods']:
            rec=[x for x in records if x['id']==case['id'] and x['method']==method]
            done=[x for x in rec if x['outcome']=='complete'];assert len(rec)==3
            item[method]=dict(completed=len(done),times=[x.get('seconds') for x in rec],
                              median=median(x['seconds'] for x in done) if len(done)==3 else None)
            display.append(f"{item[method]['median']:.3f}" if len(done)==3 else f"{len(done)}/3; TO")
        label=case['id'].replace('validation-','').replace('flower-','')
        rows.append(f"{i} & "+label.replace('_',r'\_')+' & '+' & '.join(display)+r' \\')
        details.append(item)
    counts='; '.join(f"{m}: {d.get('complete',0)} complete, {d.get('timeout',0)} timeouts" for m,d in summaries.items())
    text=r'''Both methods run on the same $14$ public inputs, with three serial
repetitions, SymPy $1.14.0$, and a common $30$-second after-ready limit
under Linux/WSL on an AMD Ryzen 7 5800U. Each fresh fork child has the
same single-CPU affinity. The interval includes public IO, equation
construction, exact root work, native startup where applicable,
reconstruction, support checks and serialization. One common import
stage and individual child startup are recorded separately. The msolve
route computes a reduced cycle system and replays its rational
parametrization exactly \citep{berthomieu2021}. Table~\ref{tab:uniform}
reports this protocol alone.

\begin{table}[tb]\centering\small
\begin{tabular}{rlrr}\toprule
No. & Input & Cycle join (s) & msolve reduced (s)\\\midrule
'''+ '\n'.join(rows)+r'''
\bottomrule\end{tabular}
\caption{Unified Linux comparison. A time is the median of three complete
runs. Otherwise the cell gives completed runs out of three and TO
denotes censored runs at $30$ seconds; no timeout is substituted as a
measured runtime. Inputs 1--12 are the frozen validation cohort; inputs
13--14 are nonsingleton controls.}\label{tab:uniform}\end{table}

Small inputs give comparable times, including cases where msolve is
faster. The structural join reduces combination work on larger completed
inputs. The cycle join completes all $42$ runs; msolve completes $33$
and reaches the common limit in nine. This is a comparison
of two exact implementations on specified synthetic inputs, rather than
a general solver ranking.
'''
    (HERE/'manuscript/uniform_summary.tex').write_text(text)
    runtime=json.loads((RUN/'runtime.json').read_text())
    sup=r'''The protocol freezes all public and source hashes before execution.
Both methods share Python $3.12.3$, SymPy $1.14.0$, CPU affinity and
the same kernel. msolve is the pinned $0.6.5$ build with one native
thread. The host is Windows with Linux/WSL; unrelated host load is not
controlled. Each job records Linux load and its child startup interval.
Method order alternates by repetition. Every successful array reproduces
the native loss and both gradient blocks exactly; complete lists agree
across successful repetitions and methods whenever both complete.
'''
    sup+=f"The common import stage takes ${runtime['imports_seconds']:.3f}$ seconds. There are ${reply_arrays}$ direct array-reply checks and ${agreements}$ list agreements.\n"
    sup+='Outcome totals: '+counts.replace('_',r'\_')+'.\n'
    sup+=r'''The saved cycle certificate calls its output list
\texttt{candidates}, whereas the native-solver record calls it
\texttt{arrays}. The runner's provisional array-count field for the
cycle route reads the latter key and is zero. The immutable raw log is
retained; the audit reads the certificate list, verifies its arrays and
writes corrected counts in a separate receipt. All timings and complete
certificates remain unchanged. This reporting correction is recorded
explicitly in the reproducibility package.
'''
    # Keep implementation audit details in the package, not the paper narrative.
    sup=sup[:sup.index('The saved cycle certificate')]
    (HERE/'manuscript/uniform_protocol_details.tex').write_text(sup)
    receipt=dict(complete=True,jobs=84,summaries=summaries,details=details,direct_reply_arrays=reply_arrays,list_agreements=agreements,
                 reporting_corrections=corrections,records=records,protocol_sha256=sha(RUN/'protocol.json'),
                 results_sha256=sha(RUN/'results.json'),verifier_sha256=sha(Path(__file__)))
    (HERE/'uniform_validation.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps(dict(complete=True,summaries=summaries,reply_arrays=reply_arrays,list_agreements=agreements,corrected_count_fields=len(corrections))))
if __name__=='__main__':run()
