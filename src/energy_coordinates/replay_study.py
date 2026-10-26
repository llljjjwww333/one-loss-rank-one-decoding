"""Read-only replay of all native certificates and exhaustive controls."""
from pathlib import Path
from collections import Counter
from fractions import Fraction as Q
import copy,hashlib,json
from verify_coordinate import verify
from run_study import controls, independent_native_check, small_full_oracle

HERE=Path(__file__).resolve().parent;PAPER=HERE.parent;PRIOR=PAPER/'revisions/20261026_mainline'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def records(path):return list(map(json.loads,path.open(encoding='utf-8')))

def run():
    protocol=json.loads((HERE/'protocol.json').read_text());report=json.loads((HERE/'validation.json').read_text())
    assert report['complete'] and report['protocol_sha256']==sha(HERE/'protocol.json')
    for name,digest in protocol['source_hashes'].items():assert sha(HERE/name)==digest
    for name,digest in protocol['shared_library_hashes'].items():assert sha(PAPER/name)==digest
    assert sha(HERE/'exploratory_cached.json')==protocol['exploratory_sha256']
    assert sha(HERE/'outcomes.jsonl')==report['outcomes_sha256']
    queries=records(PAPER/'review4_probability_v1/queries.jsonl');rows=records(HERE/'outcomes.jsonl')
    assert sha(PAPER/'review4_probability_v1/queries.jsonl')==protocol['queries_sha256']
    assert sha(PAPER/'review4_probability_v1/outcomes.jsonl')==protocol['reference_sha256']
    assert sha(PRIOR/'squareclass_outcomes.jsonl')==protocol['prior_outcomes_sha256']
    old={r['id']:r for r in records(PRIOR/'squareclass_outcomes.jsonl')}
    refs={r['id']:r['result'] for r in records(PAPER/'review4_probability_v1/outcomes.jsonl')}
    assert len(queries)==len(rows)==1605 and [r['id'] for r in queries]==[r['id'] for r in rows]
    already=nontrivial=max_gradient=0
    for query,row in zip(queries,rows):
        result=row['result'];assert verify(query['public'],result)
        independent_native_check(query['public'],result)
        reference=refs[row['id']];assert row['reference_count']==reference['complete_count']
        assert row['prior_status']==old[row['id']]['result']['status']
        if row['prior_status']=='CertifiedUnique':assert result['status']=='CertifiedUnique'
        if result['status']=='CertifiedUnique':
            assert reference['complete_count']==1
            already+=reference['gradient_count']==1
            nontrivial+=reference['gradient_count']>1
            max_gradient=max(max_gradient,reference['gradient_count'])
        if result['status']=='CertifiedEmpty':assert reference['complete_count']==0
    counts=dict(Counter(r['result']['status'] for r in rows))
    assert counts==report['statuses']=={'CertifiedUnique':1486,'Unresolved':119}
    assert (already,nontrivial,max_gradient)==(125,1361,531441)
    unresolved=dict(Counter('ambiguous' if r['reference_count']>1 else 'singleton' for r in rows if r['result']['status']=='Unresolved'))
    assert unresolved==report['unresolved_reference_counts']=={'ambiguous':65,'singleton':54}
    summary={}
    for r in [1,2,4,8,12]:
        summary[str(r)]=[sum(row['result']['status']=='CertifiedUnique' for row in rows if row['target']==f'PM_r{r}' and row['N']==N) for N in [1,2,4,8,16,32]]
    assert summary=={'1':[1,59,64,64,64,64],'2':[0,60,64,64,64,64],'4':[0,54,64,64,64,64],
                     '8':[0,33,57,63,64,64],'12':[0,19,52,64,64,64]}
    lattice={r['id']:r for r in json.loads((PRIOR/'lattice_outcomes.json').read_text())['records']}
    additional=report['additional'];assert len(additional)==46
    extra_summary={}
    for row in additional:
        subdir='uniform_linux/public' if row['cohort']=='uniform' else 'lattice_public'
        path=PRIOR/subdir/(row['id']+'.json');assert sha(path)==row['public_sha256']
        public=json.loads(path.read_text());assert verify(public,row['result'])
        independent_native_check(public,row['result'])
        if row['cohort']=='lattice' and row['result']['status']!='Unresolved':
            want='1' if row['result']['status']=='CertifiedUnique' else '0'
            assert lattice[row['id']]['reference_count']==want
    for cohort in ['uniform','lattice']:
        extra_summary[cohort]=dict(Counter(r['result']['status'] for r in additional if r['cohort']==cohort))
    fresh=records(HERE/'prospective_inputs.jsonl');outcomes=records(HERE/'prospective_outcomes.jsonl')
    assert sha(HERE/'prospective_inputs.jsonl')==protocol['prospective_inputs_sha256']
    assert sha(HERE/'prospective_outcomes.jsonl')==report['prospective']['outcomes_sha256']
    assert len(fresh)==len(outcomes)==96 and [r['id'] for r in fresh]==[r['id'] for r in outcomes]
    small=0
    for query,row in zip(fresh,outcomes):
        result=row['result'];assert verify(query['public'],result)
        independent_native_check(query['public'],result)
        assert result['status']!='CertifiedEmpty'
        if result['status']=='CertifiedUnique':assert result['values']==query['truth_array']
        if query['r']<=4:
            assert small_full_oracle(query['public'],result)==row['complete_actual_count']
            small+=1
    assert small==48
    assert dict(Counter(r['result']['status'] for r in outcomes))=={'CertifiedUnique':84,'Unresolved':12}
    largest_cert=max(r['result']['cartesian_before'] for r in outcomes if r['result']['status']=='CertifiedUnique')
    assert largest_cert==31381059609
    assert controls()==report['arithmetic_controls']
    chosen=next(r for r in rows if r['result']['status']=='CertifiedUnique' and r['result']['trace'] and r['result']['representatives'])
    public=next(r['public'] for r in queries if r['id']==chosen['id'])
    rejected=[]
    for label in ['omitted_root','false_bound','wrong_candidate','merged_classes']:
        bad=copy.deepcopy(chosen['result'])
        if label=='omitted_root':bad['cycles'][0]['roots'].pop()
        elif label=='false_bound':bad['trace'][0][0]['need']=str(Q(bad['trace'][0][0]['need'])+1)
        elif label=='wrong_candidate':bad['values'][0]=str(Q(bad['values'][0])+1)
        else:bad['representatives'].append(bad['representatives'][0])
        try:verify(public,bad)
        except AssertionError:rejected.append(label)
        else:raise AssertionError('accepted corrupted certificate')
    result=dict(complete=True,frozen_queries=1605,certified_unique=1486,prior_certified=708,newly_certified=778,
                nontrivial_pruning=nontrivial,already_singleton_gradients=already,max_actual_gradient_candidates=max_gradient,
                unresolved_singletons=54,retained_ambiguities=65,prospective_inputs=96,prospective_certified=84,
                prospective_unresolved=12,small_complete_oracles=48,largest_certified_real_root_product=largest_cert,
                additional=extra_summary,vector_oracles=500,sturm_controls=6,radical_controls=3,empty_masks=2,
                rejected_corrupted_certificates=rejected,table_counts=summary,
                verifier_sha256=sha(Path(__file__)),protocol_sha256=sha(HERE/'protocol.json'),
                source_hashes=protocol['source_hashes'],outcomes_sha256=sha(HERE/'outcomes.jsonl'),
                prospective_outcomes_sha256=sha(HERE/'prospective_outcomes.jsonl'),
                scope='local exact certificate replay; shares rational arithmetic and incidence reduction; not a new priority or random-query complexity result')
    (HERE/'replay_validation.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))
if __name__=='__main__':run()
