"""Validate source bindings and complete-list/original-reply comparisons.

Does not rerun timed solvers or replace the algebra backend's correctness.
"""
import hashlib
import json
from pathlib import Path
import sympy as s
from exact_cactus_decoder import canonical_public,exact_sign
from run_solver_comparison import same_lists,METHODS

HERE=Path(__file__).resolve().parent
RUN=HERE/'solver_comparison_v1'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read(path):
    return json.loads(path.read_text(encoding='utf-8'))


def direct_reply(public,arrays):
    p=canonical_public(public)
    for row in arrays:
        assert len(row)==len(p['edges'])
        values=list(map(s.sympify,row))
        gu,gv,loss=[s.S.Zero]*p['m'],[s.S.Zero]*p['n'],s.S.Zero
        for (i,j),c in zip(p['edges'],values):
            residual=p['u'][i]*p['v'][j]-c
            gu[i]+=residual*p['v'][j]
            gv[j]+=residual*p['u'][i]
            loss+=residual**2/2
        assert all(exact_sign(s.expand(a-b))==0 for a,b in zip(gu+gv+[loss],p['gu']+p['gv']+[p['loss']]))
    return len(arrays)


def main():
    protocol,result=read(RUN/'protocol.json'),read(RUN/'results.json')
    assert result['experiment_finished'] and not result['issues']
    assert result['protocol_sha256']==sha(RUN/'protocol.json')
    assert not protocol['truth_passed_to_solver']
    assert protocol['references_sha256']==sha(RUN/'references.json')
    for name,digest in protocol['source_hashes'].items():
        assert sha(RUN/'sources'/name)==digest and sha(HERE/name)==digest,name
    for name,digest in protocol['input_hashes'].items():
        assert sha(RUN/'public'/name)==digest,name
    references=read(RUN/'references.json')
    records={(r['id'],r['method']):r for r in result['records']}
    assert len(records)==len(protocol['cases'])*len(METHODS)
    complete_pairs,reply_checks,cross_method_pairs=0,0,0
    for case in protocol['cases']:
        public=read(RUN/'public'/(case['id']+'.json'))
        outputs=[]
        for method in METHODS:
            record=records[case['id'],method]
            assert record['bounded']==case['bounded']
            if record['outcome']=='timeout':
                assert record['wall_seconds']>=protocol['wall_timeout_seconds']
                continue
            assert record['outcome'] in ['complete','incomplete']
            path=RUN/'outputs'/(case['id']+'--'+method+'.json')
            assert sha(path)==record['output_sha256']
            output=read(path)
            assert output['method']==method and output['bounded']==case['bounded']
            assert output['truth_passed_to_solver'] is False
            assert 0<=output['solver_core_seconds']<=record['wall_seconds']
            if not output['complete']:
                assert output['status'] in ['Degenerate','positive_dimensional']
                continue
            reply_checks+=direct_reply(public,output['arrays'])
            outputs.append(output)
            if references[case['id']]['complete']:
                assert same_lists(output['arrays'],references[case['id']]['arrays'])
                complete_pairs+=1
        if outputs:
            for other in outputs[1:]:
                assert same_lists(outputs[0]['arrays'],other['arrays'])
                cross_method_pairs+=1
    # The all-zero four-cycle exposes a useful conservative exit: the
    # cycle routine is Degenerate, but norm zero alone forces every REAL
    # observed entry to zero. Both polynomial baselines should recover it.
    for method in ['factor_groebner','reduced_groebner']:
        path=RUN/'outputs'/('zero-polynomial--'+method+'.json')
        if path.exists() and read(path)['complete']:
            assert read(path)['status']=='Unique'
            assert read(path)['arrays']==[['0']*4]
    table=['\\begin{tabular}{rrrrrrr}','\\toprule',
           '$r$ & $k$ & raw bits & enum. & join & reduced GB & factor GB\\\\',
           '\\midrule']
    for case in protocol['cases']:
        if not case['id'].startswith('grid-'):
            continue
        _,r,k,b=case['id'].split('-')
        values=[]
        for method in METHODS:
            rec=records[case['id'],method]
            values.append(f"{rec['wall_seconds']:.2f}" if rec['outcome']=='complete' else
                          '$\\geq '+str(protocol['wall_timeout_seconds'])+'$' if rec['outcome']=='timeout' else rec['outcome'])
        table.append(' & '.join([r[1:],k[1:],b[1:]]+values)+'\\\\')
    table+=['\\bottomrule','\\end{tabular}']
    target=HERE/'solver_comparison_summary.tex'
    target.write_text('\n'.join(table)+'\n',encoding='utf-8')
    audit=dict(complete=True,scope='Source-bound benchmark audit; same SymPy backend, not a fastest-solver result',
               jobs=len(records),cases=len(protocol['cases']),summary=result['summary'],
               complete_reference_list_checks=complete_pairs,original_candidate_reply_checks=reply_checks,
               complete_cross_method_list_checks=cross_method_pairs,
               protocol_sha256=sha(RUN/'protocol.json'),results_sha256=sha(RUN/'results.json'),
               audit_source_sha256=sha(Path(__file__)),table_sha256=sha(target),
               source_hashes=protocol['source_hashes'])
    (HERE/'checks/solver-comparison-validation.json').write_text(json.dumps(audit,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(audit,indent=2))


if __name__=='__main__':
    main()
