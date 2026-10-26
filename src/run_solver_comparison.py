"""Freeze a serial four-method comparison; never discard failures/timeouts.

The cohort is the already-used exact development/control cohort, NOT a new
holdout. Every method starts in a fresh process with the same wall limit.
"""
import argparse
import hashlib
import json
from pathlib import Path
import platform
import subprocess
import sys
from time import perf_counter
import sympy as s
from exact_cactus_decoder import exact_sign

HERE=Path(__file__).resolve().parent
METHODS=['cycle_enumerate','cycle_mitm','reduced_groebner','factor_groebner']


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read(path):
    return json.loads(path.read_text(encoding='utf-8'))


def write(path,data):
    path.write_text(json.dumps(data,indent=2)+'\n',encoding='utf-8')


def same_lists(first,second):
    if len(first)!=len(second):
        return False
    remaining=[list(map(s.sympify,row)) for row in second]
    for row in first:
        values=list(map(s.sympify,row))
        for i,other in enumerate(remaining):
            if len(values)==len(other) and all(exact_sign(s.expand(a-b))==0 for a,b in zip(values,other)):
                remaining.pop(i)
                break
        else:
            return False
    return not remaining


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('folder',type=Path)
    parser.add_argument('--timeout',type=int,default=45)
    args=parser.parse_args()
    folder=args.folder.resolve()
    if folder.exists():
        raise ValueError('Existing run folder is immutable; choose a fresh name')
    for name in ['public','sources','outputs','logs']:
        (folder/name).mkdir(parents=True,exist_ok=True)
    previous=HERE/'study_exact_v6'
    old_protocol=read(previous/'protocol.json')
    cases=old_protocol['cases']
    for case in cases:
        source=previous/'public'/(case['id']+'.json')
        assert sha(source)==old_protocol['input_hashes'][source.name]
        (folder/'public'/source.name).write_bytes(source.read_bytes())
    references={}
    for case in cases:
        cert=read(previous/'certificates'/(case['id']+'.json'))
        references[case['id']]=dict(complete=cert['complete'],status=cert['status'],
                                  arrays=[c['values'] for c in cert['candidates']])
    write(folder/'references.json',references)
    names=['polynomial_solver_baseline.py','exact_cactus_decoder.py','run_solver_comparison.py']
    for name in names:
        (folder/'sources'/name).write_bytes((HERE/name).read_bytes())
    protocol=dict(format_version=1,cases=cases,methods=METHODS,wall_timeout_seconds=args.timeout,
                  timing='Fresh subprocess wall time includes startup/imports/JSON IO; core time begins after imports and includes parsing, equations, solving and final filters',
                  process_policy='Serial execution, one fresh worker per method and case; no repeat/best-of timing',
                  source_hashes={name:sha(folder/'sources'/name) for name in names},
                  input_hashes={p.name:sha(p) for p in sorted((folder/'public').glob('*.json'))},
                  references_sha256=sha(folder/'references.json'),
                  reference_protocol_sha256=sha(previous/'protocol.json'),
                  cohort='The same 33 previously used exact development/control cases; NOT independent holdout data',
                  solver_backend='SymPy 1.14.0 lex Groebner and solve_triangulated(extension=True); no Singular/msolve performance claim',
                  truth_passed_to_solver=False,python=platform.python_version(),sympy=s.__version__,platform=platform.platform())
    write(folder/'protocol.json',protocol)
    records=[]
    for case in cases:
        # Fixed rotating order avoids giving one method a fixed within-case position.
        offset=len(records)//len(METHODS)%len(METHODS)
        order=METHODS[offset:]+METHODS[:offset]
        for method in order:
            job=case['id']+'--'+method
            output=folder/'outputs'/(job+'.json')
            command=[sys.executable,str(folder/'sources/polynomial_solver_baseline.py'),
                     str(folder/'public'/(case['id']+'.json')),str(output),'--method',method]
            if case['bounded']:
                command.append('--bounded')
            record=dict(id=case['id'],method=method,bounded=case['bounded'])
            started=perf_counter()
            try:
                completed=subprocess.run(command,capture_output=True,text=True,timeout=args.timeout)
                record.update(wall_seconds=perf_counter()-started,exit_code=completed.returncode)
                (folder/'logs'/(job+'.stdout.txt')).write_text(completed.stdout,encoding='utf-8')
                (folder/'logs'/(job+'.stderr.txt')).write_text(completed.stderr,encoding='utf-8')
                if completed.returncode:
                    record['outcome']='error'
                else:
                    result=read(output)
                    assert not result['truth_passed_to_solver']
                    record.update(outcome='complete' if result['complete'] else 'incomplete',status=result['status'],
                                  retained=len(result['arrays']),solver_core_seconds=result['solver_core_seconds'],
                                  stats=result['stats'],output_sha256=sha(output))
                    validation_start=perf_counter()
                    if result['complete'] and references[case['id']]['complete']:
                        record['same_complete_list']=same_lists(result['arrays'],references[case['id']]['arrays'])
                        if not record['same_complete_list']:
                            record['outcome']='list_mismatch'
                    else:
                        record['same_complete_list']=None
                    record['validation_seconds']=perf_counter()-validation_start
            except subprocess.TimeoutExpired as error:
                record.update(outcome='timeout',wall_seconds=perf_counter()-started)
                for stream in ['stdout','stderr']:
                    value=getattr(error,stream) or b''
                    if isinstance(value,bytes): value=value.decode('utf-8',errors='replace')
                    (folder/'logs'/(job+'.'+stream+'.txt')).write_text(value,encoding='utf-8')
            except Exception as error:
                record.update(outcome='validation_error',wall_seconds=perf_counter()-started,error=repr(error))
            records.append(record)
            summary={method:{state:sum(r['method']==method and r['outcome']==state for r in records)
                             for state in ['complete','incomplete','timeout','error','list_mismatch','validation_error']}
                     for method in METHODS}
            report=dict(experiment_finished=len(records)==len(cases)*len(METHODS),
                        all_solvers_completed=all(r['outcome']=='complete' for r in records),
                        scope='Recorded development benchmark, including noncompletion; no fastest-solver or population performance claim',
                        protocol_sha256=sha(folder/'protocol.json'),records=records,summary=summary,
                        issues=[dict(id=r['id'],method=r['method'],outcome=r['outcome']) for r in records
                                if r['outcome'] in ['error','list_mismatch','validation_error']])
            write(folder/'results.json',report)
            print(case['id'],method,record['outcome'],round(record['wall_seconds'],3),flush=True)
    raise SystemExit(bool(report['issues']))


if __name__=='__main__':
    main()
