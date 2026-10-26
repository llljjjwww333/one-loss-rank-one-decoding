"""Whole-cohort native generation and separate exact certificate replay."""
from pathlib import Path
from collections import Counter, defaultdict
from fractions import Fraction as Q
from itertools import product
import copy, hashlib, json, random, sys
from coordinate_decoder import decode, propagate, coordinates, Quad, reduction
from verify_coordinate import verify, real_root_count
import sympy as s

HERE = Path(__file__).resolve().parent
PAPER = HERE.parent
PRIOR = PAPER/'revisions/20261026_mainline'


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def controls():
    rng = random.Random(2026102606)
    arithmetic = []
    # Exhaustive oracles test every feasible tuple, not only a chosen witness.
    for index in range(500):
        dims = rng.randrange(1, 5)
        vectors = [[[Q(rng.randrange(-5, 6), rng.randrange(1, 4)) for _ in range(dims)]
                    for _ in range(rng.randrange(1, 4))] for _ in range(rng.randrange(0, 6))]
        target = [Q(rng.randrange(-8, 9), rng.randrange(1, 4)) for _ in range(dims)]
        if vectors and index % 2 == 0:
            path = [rng.randrange(len(row)) for row in vectors]
            target = [sum(row[j][k] for row, j in zip(vectors, path)) for k in range(dims)]
        feasible = [path for path in product(*(range(len(row)) for row in vectors))
                    if all(sum(vectors[i][j][k] for i, j in enumerate(path)) == target[k] for k in range(dims))]
        domains, trace = propagate(vectors, target)
        assert all(all(j in domains[i] for i, j in enumerate(path)) for path in feasible)
        if not all(domains):
            assert not feasible
        if all(len(dom) == 1 for dom in domains):
            assert len(feasible) <= 1
        arithmetic.append(dict(index=index, cycles=len(vectors), feasible=len(feasible),
                               retained=list(map(len, domains)), rounds=len(trace)))
    radical_cases = [
        [[Quad(2), Quad(2, 1, 2)], [Quad(3), Quad(3, -1, 8)]],
        [[Quad(1), Quad(0, 1, 2)], [Quad(2), Quad(0, 1, 3)], [Quad(3), Quad(0, 1, 6)]],
        [[Quad(1), Quad(1, 1, 2), Quad(1, -1, 2)],
         [Quad(1), Quad(1, 2, 2), Quad(1, -2, 2)]],
    ]
    radicals = []
    for groups in radical_cases:
        reps, vectors = coordinates(groups)
        domains, trace = propagate(vectors, [sum(g[0].a for g in groups)]+[Q(0)]*len(reps))
        radicals.append(dict(representatives=list(map(str, reps)), domains=domains, rounds=len(trace)))
    t = s.Symbol('t')
    sturm = []
    for expr, expected in [(t-1,1),(t*t+1,0),(t*t-2,2),((t-1)**3,1),((t-1)**2*(t+2),2),(t**3-2,1)]:
        count = real_root_count(s.Poly(expr, t, domain=s.QQ))
        assert count == expected
        sturm.append(dict(polynomial=str(expr), distinct_real_roots=count))
    native = []
    empty = dict(m=0,n=0,edges=[],u=[],v=[],gu=[],gv=[],loss='0')
    for loss, want in [('0','CertifiedUnique'),('1','CertifiedEmpty')]:
        public = dict(empty, loss=loss)
        result = decode(public)
        assert result['status'] == want and verify(public, result)
        native.append(dict(name='empty_loss_'+loss, result=result))
    return dict(seed=2026102606, vector_oracles=arithmetic, radical_controls=radicals,
                sturm_controls=sturm, native_controls=native)


def independent_native_check(public, result):
    if result['status'] != 'CertifiedUnique':
        return
    vals = list(map(Q, result['values']))
    m, n = public['m'], public['n']
    u, v = list(map(Q, public['u'])), list(map(Q, public['v']))
    gu, gv = [Q(0)]*m, [Q(0)]*n
    active_rows, active_cols = set(), set()
    loss = Q(0)
    for (i,j), val in zip(public['edges'], vals):
        residual = u[i]*v[j]-val
        gu[i] += residual*v[j]
        gv[j] += residual*u[i]
        loss += residual*residual/2
        if val:
            active_rows.add(i)
            active_cols.add(j)
    assert gu == list(map(Q,public['gu'])) and gv == list(map(Q,public['gv'])) and loss == Q(public['loss'])
    assert not any(not val and i in active_rows and j in active_cols for (i,j),val in zip(public['edges'],vals))


def small_full_oracle(public, result):
    p,cycles,weights,base,_,_=reduction(public)
    vectors=[[[Q(z) for z in vec] for vec in row] for row in result['vectors']]
    target=list(map(Q,result['target']));roots=[[Quad(*z) for z in row['roots']] for row in result['cycles']]
    count=0
    for path in product(*(range(len(row)) for row in roots)):
        if any(sum(vectors[i][j][k] for i,j in enumerate(path))!=target[k] for k in range(len(target))):
            continue
        values=[Quad(Q(z/w)) for z,w in zip(base,weights)]
        for cycle,row,j in zip(cycles,roots,path):
            for k,e in enumerate(cycle):values[e]=values[e]+row[j]*Q((-1)**k/weights[e])
        active_r={i for (i,j),val in zip(p['edges'],values) if not val.zero()}
        active_c={j for (i,j),val in zip(p['edges'],values) if not val.zero()}
        if not any(val.zero() and i in active_r and j in active_c for (i,j),val in zip(p['edges'],values)):
            count+=1
    return count


def run():
    protocol = json.loads((HERE/'protocol.json').read_text())
    assert not (HERE/'validation.json').exists(), 'use a fresh study namespace'
    for name, digest in protocol['source_hashes'].items():
        assert sha(HERE/name) == digest
    queries = PAPER/'review4_probability_v1/queries.jsonl'
    reference = PAPER/'review4_probability_v1/outcomes.jsonl'
    assert sha(queries) == protocol['queries_sha256'] and sha(reference) == protocol['reference_sha256']
    prior = {r['id']:r for r in map(json.loads,(PRIOR/'squareclass_outcomes.jsonl').open())}
    refs = {r['id']:r['result'] for r in map(json.loads,reference.open())}
    rows, summary = [], defaultdict(Counter)
    with (HERE/'outcomes.jsonl').open('x',encoding='utf-8') as output:
        for index, query in enumerate(map(json.loads,queries.open()),1):
            result = decode(query['public'])
            assert verify(query['public'],result)
            independent_native_check(query['public'],result)
            # These annotations are read only after generation and replay.
            ref = refs[query['id']]
            assert ref['complete']
            if result['status']=='CertifiedUnique':
                assert ref['complete_count']==1
            if result['status']=='CertifiedEmpty':
                assert ref['complete_count']==0
            if prior[query['id']]['result']['status']=='CertifiedUnique':
                assert result['status']=='CertifiedUnique'
            row = dict(id=query['id'],target=query['target'],N=query['N'],reference_count=ref['complete_count'],
                       prior_status=prior[query['id']]['result']['status'],result=result)
            output.write(json.dumps(row,sort_keys=True)+'\n')
            rows.append(row)
            summary[f"{query['target']}_N{query['N']}"][result['status']]+=1
            if index % 100==0:
                output.flush()
                print(index,dict(Counter(r['result']['status'] for r in rows)),flush=True)
    assert len(rows)==len(refs)==1605
    additional=[]
    for cohort,directory in [('uniform',PRIOR/'uniform_linux/public'),('lattice',PRIOR/'lattice_public')]:
        for path in sorted(directory.glob('*.json')):
            public=json.loads(path.read_text());result=decode(public)
            assert verify(public,result)
            independent_native_check(public,result)
            additional.append(dict(cohort=cohort,id=path.stem,public_sha256=sha(path),result=result))
    assert sha(HERE/'prospective_inputs.jsonl')==protocol['prospective_inputs_sha256']
    prospective=[]
    with (HERE/'prospective_outcomes.jsonl').open('x',encoding='utf-8') as output:
        for query in map(json.loads,(HERE/'prospective_inputs.jsonl').open()):
            result=decode(query['public']);assert verify(query['public'],result)
            independent_native_check(query['public'],result)
            assert result['status']!='CertifiedEmpty'
            if result['status']=='CertifiedUnique':assert result['values']==query['truth_array']
            row=dict(id=query['id'],r=query['r'],N=query['N'],result=result)
            if query['r']<=4:
                count=small_full_oracle(query['public'],result);assert count>=1
                if result['status']=='CertifiedUnique':assert count==1
                row['complete_actual_count']=count
            prospective.append(row);output.write(json.dumps(row,sort_keys=True)+'\n')
        print('prospective',len(prospective),dict(Counter(r['result']['status'] for r in prospective)),flush=True)
    chosen=next(r for r in rows if r['result']['status']=='CertifiedUnique' and r['result']['trace'] and r['result']['representatives'])
    public=next(r['public'] for r in map(json.loads,queries.open()) if r['id']==chosen['id'])
    mutations=[]
    for label in ['omitted_root','false_bound','wrong_candidate','merged_classes']:
        bad=copy.deepcopy(chosen['result'])
        if label=='omitted_root':
            bad['cycles'][0]['roots'].pop()
        elif label=='false_bound':
            bad['trace'][0][0]['need']=str(Q(bad['trace'][0][0]['need'])+1)
        elif label=='wrong_candidate':
            bad['values'][0]=str(Q(bad['values'][0])+1)
        else:
            bad['representatives'].append(bad['representatives'][0])
        try:
            verify(public,bad)
        except AssertionError:
            mutations.append(dict(mutation=label,rejected=True))
        else:
            raise AssertionError('corrupted certificate accepted: '+label)
    report=dict(complete=True,protocol_sha256=sha(HERE/'protocol.json'),source_hashes=protocol['source_hashes'],
                queries=1605,no_exclusions=True,truth_passed_to_solver=False,
                statuses=dict(Counter(r['result']['status'] for r in rows)),
                prior_certified=708,newly_certified=sum(r['result']['status']=='CertifiedUnique' and r['prior_status']=='Unresolved' for r in rows),
                unresolved_reference_counts=dict(Counter('ambiguous' if r['reference_count']>1 else 'singleton' for r in rows if r['result']['status']=='Unresolved')),
                max_rounds=max(r['result'].get('rounds',0) for r in rows),
                max_deletions=max(r['result'].get('deletions',0) for r in rows),
                summary={k:dict(v) for k,v in sorted(summary.items())},additional=additional,
                prospective=dict(inputs=len(prospective),statuses=dict(Counter(r['result']['status'] for r in prospective)),
                                 small_complete_oracles=sum('complete_actual_count' in r for r in prospective),
                                 max_root_product=max(r['result']['cartesian_before'] for r in prospective),
                                 outcomes_sha256=sha(HERE/'prospective_outcomes.jsonl')),
                arithmetic_controls=controls(),corrupted_certificates=mutations,outcomes_sha256=sha(HERE/'outcomes.jsonl'),
                scope='exact local public-input certificate study; classical propagation; no new CAS timing claim or average-query complexity theorem')
    (HERE/'validation.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k not in ['summary','additional','arithmetic_controls']},indent=2))


if __name__=='__main__':
    run()
