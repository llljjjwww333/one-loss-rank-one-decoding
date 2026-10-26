"""Replay native certificates and compare RURs using rational polynomial identities."""
from pathlib import Path
from fractions import Fraction as Q
from collections import Counter
import hashlib,json,sys
HERE=Path(__file__).resolve().parent;sys.path.insert(0,str(HERE/'sources'))
from independent_real_algebra import (trim,add,mul,divrem,monic,pgcd,derivative,isolate,RealRoot,sturm,count)
from singular_real_worker import actual_filter,basis_at_point,direct_reply
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def put(p,x):p.write_text(json.dumps(x,indent=2)+'\n')
def compose(p,q,f):
    z=[]
    for c in reversed(p):z=divrem(add(mul(z,q),[Q(c)]),f)[1]
    return z
def minus(a,b):return add(a,[-x for x in b])
def reference():
    assert not (HERE/'reference_snapshot.json').exists()
    p=HERE.parent/'revisions/20261026_selfcheck'
    old=json.loads((p/'rur_results_v3.json').read_text())
    lists={}
    for f in sorted((p/'uniform_linux/outputs').glob('*--cycle_mitm--1.json')):
        x=json.loads(f.read_text());lists[f.name.split('--')[0]]=dict(arrays=[c['values'] for c in x['candidates']],source_sha256=sha(f))
    assert len(lists)==14
    put(HERE/'reference_snapshot.json',dict(rur_source_sha256=sha(p/'rur_results_v3.json'),rur_records=old['records'],baseline_lists=lists))
def native_replay(public,out,bounded):
    assert out['complete'] and out['uses_sympy'] is False
    if out.get('complex_points')==0:
        assert out['status']=='Empty' and not out['candidates'];return dict(points=0,real_points=0,actual=0)
    f=monic(out['eliminator']);assert len(f)-1==out['complex_points'] and len(pgcd(f,derivative(f)))==1
    product=[Q(1)]
    for factor in out['factor_polynomials']:product=mul(product,monic(factor))
    assert monic(product)==f
    expected=sum(len(isolate(q)) for q in out['factor_polynomials'])
    assert len(out['decisions'])==expected
    seen=set();accepted=[]
    for d in out['decisions']:
        root_data=d['root'];factor=out['factor_polynomials'][d['factor_index']]
        assert monic(factor)==monic(root_data['integer_polynomial'])
        key=(d['factor_index'],root_data['real_root_index']);assert key not in seen;seen.add(key)
        canonical=isolate(factor);idx=root_data['real_root_index'];assert 0<=idx<len(canonical)
        ca,cb=canonical[idx].a,canonical[idx].b;a,b=map(Q,root_data['isolating_interval'])
        if a==b:assert ca<=a<=cb
        else:assert max(ca,a)<min(cb,b) and count(sturm(factor),max(ca,a),min(cb,b))==1
        root=RealRoot(factor,*root_data['isolating_interval'],root_data['real_root_index'])
        assert all(divrem(p,factor)[1]==trim(v) for p,v in zip(out['coordinate_polynomials'],d['values']))
        values=[root.value(v) for v in d['values']];direct_reply(public,values)
        basis_at_point(out['basis'],values,root.value([0,1]))
        actual,certificate=actual_filter(public,values,bounded);assert actual==d['accepted']
        if actual:accepted.append(d)
    assert accepted==out['candidates']
    return dict(points=out['complex_points'],real_points=expected,actual=len(accepted))
def compare_rurs(ref,out):
    if not ref.get('rur') or ref['rur'].get('empty'):
        assert not out['candidates'] and ref['complete_count']==0
        return dict(empty_actual_list=True,native_closure_points=out['complex_points'])
    f=list(map(Q,reversed(ref['rur']['polynomial'])))
    edges=[list(map(Q,reversed(p))) for p in ref['edge_polynomials']]
    assert len(f)-1==out['complex_points']==ref['rur']['distinct_complex_points']
    projection=[]
    for i,e in enumerate(edges):projection=add(projection,[c*out['gamma']**i for c in e])
    assert not compose(list(map(Q,out['eliminator'])),projection,f)
    for a,b in zip(out['coordinate_polynomials'],edges):assert not minus(compose(list(map(Q,a)),projection,f),divrem(b,f)[1])
    matched=set()
    for d in ref['real_root_decisions']:
        root=RealRoot(f,*d['isolating_interval'],d['real_root_ordinal']);z=root.value(projection);hits=[]
        for j,native in enumerate(out['decisions']):
            data=native['root']
            if root.sign(compose(list(map(Q,data['integer_polynomial'])),projection,f)):continue
            a,b=map(Q,data['isolating_interval'])
            if a==b:
                if (z-a).sign():continue
            else:
                for bits in [16,32,64,128,256,512,1024,2048]:
                    root.refine(bits);lo,hi=z.bounds()
                    if hi<a or lo>b:break
                    if a<lo and hi<b:break
                if not(a<lo and hi<b):continue
            hits.append(j)
        assert len(hits)==1,('real embedding not identified',d,hits)
        j=hits[0];assert j not in matched;matched.add(j)
        assert bool(d['actual_and_budget_passed'])==out['decisions'][j]['accepted']
    assert len(matched)==len(out['decisions'])
    assert len(out['candidates'])==ref['complete_count']
    return dict(complex_points=out['complex_points'],real_embeddings=len(matched),polynomial_identities=len(edges)+1,actual=len(out['candidates']))
def run(require_complete=True):
    p=json.loads((HERE/'protocol.json').read_text());raw=json.loads((HERE/'results.json').read_text())
    assert raw['protocol_sha256']==sha(HERE/'protocol.json') and sha(HERE/'run_study.py')==p['runner_sha256']
    for group,folder in [('source_hashes',HERE/'sources'),('public_hashes',HERE/'public'),('control_hashes',HERE/'controls')]:
        for name,d in p[group].items():assert sha(folder/name)==d
    references=json.loads((HERE/'reference_snapshot.json').read_text());refs={r['id']:r['result'] for r in references['rur_records']}
    coverage=[];comparisons=[];native_arrays=0
    for r in raw['records']:
        if r['outcome']=='complete':
            output=HERE/'outputs'/r['output'];assert sha(output)==r['output_sha256'];out=json.loads(output.read_text())
            public=json.loads((HERE/('controls' if r['group']=='coverage' else 'public')/(r['id']+'.json')).read_text())
            if r['method']=='singular_observed':
                checks=native_replay(public,out,r['bounded']);native_arrays+=checks['actual']
                if r['group']=='coverage':
                    cross=compare_rurs(refs[r['id']],out);coverage.append(dict(id=r['id'],**checks,cross_rur=cross));continue
                arrays=[]
                for d in out['candidates']:
                    assert all(len(v)<=1 for v in d['values']),'nonrational baseline array needs field comparison'
                    arrays.append([v[0] if v else '0' for v in d['values']])
            else:arrays=out['arrays'] if r['method']=='msolve_reduced' else [c['values'] for c in out['candidates']]
            expected=references['baseline_lists'][r['id']]['arrays']
            canonical=lambda values:sorted(tuple(Q(v) for v in a) for a in values)
            assert canonical(arrays)==canonical(expected)
            for a in arrays:
                # Recheck every rational reply independently of primary arithmetic.
                root=RealRoot([Q(0),Q(1)],Q(0),Q(0),0)
                direct_reply(public,[root.value([Q(v)]) for v in a])
            comparisons.append(dict(id=r['id'],method=r['method'],repeat=r['repeat'],arrays=len(arrays)))
        else:assert r['outcome']=='timeout' and r['censor_seconds']==30,r
    assert len(coverage)==37
    if require_complete:assert raw['complete'] and len(raw['records'])==163
    controls=[]
    good=json.loads((HERE/'outputs/coverage--c4_irrational--singular_observed--1.json').read_text())
    public=json.loads((HERE/'controls/c4_irrational.json').read_text())
    for kind in ['missing_real_root','changed_coordinate','changed_support_decision']:
        bad=json.loads(json.dumps(good))
        if kind=='missing_real_root':bad['decisions'].pop()
        if kind=='changed_coordinate':bad['decisions'][0]['values'][0]=['1']
        if kind=='changed_support_decision':bad['decisions'][0]['accepted']=not bad['decisions'][0]['accepted']
        try:native_replay(public,bad,False)
        except AssertionError:controls.append(kind)
        else:raise AssertionError(('corrupt certificate accepted',kind))
    result=dict(complete=raw['complete'],source_sha256=sha(Path(__file__)),results_sha256=sha(HERE/'results.json'),reference_sha256=sha(HERE/'reference_snapshot.json'),
      coverage=coverage,comparison_list_agreements=comparisons,native_arrays_rechecked=native_arrays,
      outcomes={group:{method:dict(Counter(r['outcome'] for r in raw['records'] if r['group']==group and r['method']==method)) for method in p['methods']} for group in ['coverage','comparison']},
      corruption_controls_rejected=controls,
      scope='Native equations, rational/Sturm predicates and cross-RUR modular identities are independent of SymPy and the primary incidence reduction; native Groebner/radical computations trust Singular.')
    put(HERE/('validation.json' if require_complete else 'coverage_validation.json'),result)
    print(json.dumps({k:v for k,v in result.items() if k not in ['coverage','comparison_list_agreements']}))
if __name__=='__main__':reference() if sys.argv[1]=='reference' else run(sys.argv[1]!='coverage')
