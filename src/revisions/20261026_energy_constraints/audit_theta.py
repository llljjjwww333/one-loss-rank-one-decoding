"""Replay exact norm ideals and independently check native replies/supports.

Saved gradient-audit timeouts are retained; expensive timed audits are
not repeated by this integrity replay.
"""
from pathlib import Path
from fractions import Fraction as Q
from math import isqrt
from collections import Counter
import hashlib,json
import sympy as s
from normalized_block_decoder import public_system
HERE=Path(__file__).resolve().parent


def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def rational_check(public,values):
    m,n,edges=public['m'],public['n'],public['edges']
    u,v=[list(map(Q,public[k])) for k in ['u','v']]
    C=list(map(Q,values));gu=[Q(0)]*m;gv=[Q(0)]*n;loss=Q(0)
    adj=[[] for _ in range(m+n)]
    for e,((i,j),value) in enumerate(zip(edges,C)):
        residual=u[i]*v[j]-value;gu[i]+=residual*v[j];gv[j]+=residual*u[i];loss+=residual*residual/2
        if value:adj[i].append((m+j,e));adj[m+j].append((i,e))
    assert gu==list(map(Q,public['gu'])) and gv==list(map(Q,public['gv'])) and loss==Q(public['loss'])
    assert not any(value==0 and adj[i] and adj[m+j] for (i,j),value in zip(edges,C))
    factors={};radicands=[]
    for first in range(m+n):
        if first in factors or not adj[first]:continue
        factors[first]=Q(1);stack=[first];component=[]
        while stack:
            node=stack.pop();component.append(node)
            for other,e in adj[node]:
                if other not in factors:factors[other]=C[e]/factors[node];stack.append(other)
                else:assert factors[node]*factors[other]==C[e]
        left=sum((factors[x]**2 for x in component if x<m),Q(0))
        right=sum((factors[x]**2 for x in component if x>=m),Q(0));radicands.append(left*right)
    bounded=None
    for bits in [16,32,64,128,256]:
        scale=1<<bits;lo=hi=Q(0)
        for d in radicands:
            k=isqrt(d.numerator*scale*scale//d.denominator)
            lo+=Q(k,scale);hi+=Q(k if k*k*d.denominator==d.numerator*scale*scale else k+1,scale)
        if hi<=1:bounded=True;break
        if lo>1:bounded=False;break
    assert bounded is not None,'budget comparison unresolved in audit'
    return bounded


def run():
    data=json.loads((HERE/'theta_validation.json').read_text())
    frozen=json.loads((HERE/'theta_frozen_inputs.json').read_text())
    original=json.loads((HERE/'theta_initial_v1/theta_frozen_inputs.json').read_text())
    assert [(x['id'],x['public_sha256']) for x in frozen['cases']]==[(x['id'],x['public_sha256']) for x in original['cases']]
    assert data['protocol_sha256']==sha(HERE/'theta_protocol.json')
    for name,digest in data['source_hashes'].items():assert sha(HERE/name)==digest
    assert len(data['records'])==len(frozen['cases'])==27
    counts=[];checked=0;gradient=[]
    for row,case in zip(data['records'],frozen['cases']):
        assert row['id']==case['id'];path=HERE/'theta_public'/(row['id']+'.json')
        assert sha(path)==case['public_sha256']
        public=json.loads(path.read_text());cert=row['result']
        assert cert['complete'] and cert['public_sha256']==sha(path)
        system=public_system(public);variables=system['variables']
        G=s.groebner(system['equations']+[system['energy']],*variables,order='lex',domain=s.QQ)
        assert list(map(str,G))==cert['norm_ideal']
        points=s.solve_poly_system(list(G),*variables,strict=True)
        expected={tuple(str(s.expand(value.subs(dict(zip(variables,point))))) for value in system['values']) for point in points}
        # Every rational ideal point in this fixed cohort is actual.
        actual={tuple(c['values']) for c in cert['candidates']}
        assert actual==expected and len(actual)==cert['complete_count']
        for candidate in cert['candidates']:
            assert rational_check(public,candidate['values'])==candidate['unit_factor_budget_passed'];checked+=1
        counts.append(cert['complete_count'])
        if row['status']=='timeout':
            assert row['gradient_audit_timeout'] and row['norm_certificate_recovered']
        if 'gradient' in cert and cert['gradient']['complete']:
            gradient.append(dict(id=row['id'],count=cert['gradient']['real_actual_candidates']))
    assert Counter(counts)=={1:25,2:2} and checked==29
    assert data['statuses']=={'complete':25,'timeout':2} and len(gradient)==22
    receipt=dict(complete=True,cases=27,exact_arrays_checked=29,norm_singletons=25,norm_ambiguous=2,
                 retained_worker_timeouts=2,completed_gradient_audits=22,public_inputs_unchanged=True,
                 scope='fresh exact norm-ideal replay and independent rational native reply/support/budget checks; timed gradient audits not repeated',
                 data_sha256=sha(HERE/'theta_validation.json'),inputs_sha256=sha(HERE/'theta_frozen_inputs.json'),
                 protocol_sha256=sha(HERE/'theta_protocol.json'),
                 verifier_sha256=sha(Path(__file__)),gradient_counts=gradient)
    (HERE/'theta_audit_validation.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps({k:v for k,v in receipt.items() if k!='gradient_counts'}))


if __name__=='__main__':run()
