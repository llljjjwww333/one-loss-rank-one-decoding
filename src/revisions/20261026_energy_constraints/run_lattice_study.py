"""Freeze and run a public-input, exact uniqueness study (no CAS speed claim)."""
from pathlib import Path
from time import perf_counter
from math import comb,prod
import hashlib,json,platform,sys
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'20261025_structural_depth'))
from hardness_family import construction
from lattice_decoder import preprocess,enumeration,mitm,dynamic,verify_witness
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def dump(p,x):p.write_text(json.dumps(x,indent=2)+'\n')
def exact_reference(items,T,kind):
    if kind=='binary':return 1
    if kind=='equal':return comb(len(items),T)
    if kind=='empty':return 0
    # Uncapped independent integer generating-function coefficients.
    coeff=[1]
    for a in items:
        new=coeff+[0]*a
        for j,c in enumerate(coeff):new[j+a]+=c
        coeff=new
    return coeff[T] if T<len(coeff) else 0
def run():
    cases=[]
    for r in [1,2,4,8,12,16,24,32]:
        for kind in ['dense','binary','empty','equal']:
            if kind=='dense':items=[1+(i*7%9) for i in range(r)];T=sum(items)//2
            elif kind=='binary':items=[2**i for i in range(r)];T=sum(a for i,a in enumerate(items) if i%3==0)
            elif kind=='empty':items=[4]+[6]*(r-1);T=2
            else:items=[1]*r;T=r//2
            cases.append(dict(id=f'{kind}_r{r}',r=r,kind=kind,items=items,target=T))
    protocol=dict(version=1,cases=cases,enumeration_tuple_cap=65536,mitm_half_tuple_cap=131072,
                  dp_cell_cap=200000,timing='same Python process after imports; preprocessing and combination separately; operation counts primary',
                  methods='exact capped-count decision and singleton recovery, not complete list printing',
                  source_hashes={p.name:sha(p) for p in [Path(__file__),HERE/'lattice_decoder.py',HERE.parent/'20261025_structural_depth/hardness_family.py']})
    dump(HERE/'lattice_protocol.json',protocol)
    pubdir=HERE/'lattice_public';pubdir.mkdir(exist_ok=True);records=[]
    for case in cases:
        public,_=construction(case['items'],case['target']);p=pubdir/(case['id']+'.json');dump(p,public)
        start=perf_counter();data=preprocess(public);prep=perf_counter()-start
        costs,U=data['costs'],data['U'];r=case['r'];D=prod(map(len,costs))
        record=dict(id=case['id'],r=r,public_sha256=sha(p),root_count=data['roots'],locally_rejected=data['rejected'],
                    actual_gradient_candidates=D,integer_target=U,preprocessing_seconds=prep,
                    reference_count=str(exact_reference(case['items'],case['target'],case['kind'])),methods={})
        for name,method,within in [('enumeration',enumeration,D<=65536),('mitm',mitm,2**((r+1)//2)<=131072),
                                   ('dynamic',dynamic,(U+1)*r<=200000)]:
            if not within:
                record['methods'][name]=dict(outcome='prespecified_budget_skip');continue
            start=perf_counter();out=method(costs,U);elapsed=perf_counter()-start
            assert out['count']==min(2,int(record['reference_count'])),(case,out)
            if out['count']==1:assert verify_witness(public,data['branches'],out['witness'])
            out.update(outcome='complete',combination_seconds=elapsed);record['methods'][name]=out
        records.append(record);print(case['id'],flush=True)
    dump(HERE/'lattice_outcomes.json',dict(complete=True,protocol_sha256=sha(HERE/'lattice_protocol.json'),
        python=sys.version,platform=platform.platform(),records=records))
    dump(HERE/'lattice_validation.json',dict(complete=True,cases=len(records),
        method_completions=sum(o['outcome']=='complete' for x in records for o in x['methods'].values()),
        retained_budget_skips=sum(o['outcome']!='complete' for x in records for o in x['methods'].values()),
        singleton_reply_replays=sum(o.get('count')==1 for x in records for o in x['methods'].values()),
        public_only_preprocessing=True,outcomes_sha256=sha(HERE/'lattice_outcomes.json')))
if __name__=='__main__':run()
