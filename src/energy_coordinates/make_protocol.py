"""Freeze all new native inputs and sources before the native study."""
from pathlib import Path
from fractions import Fraction as Q
import hashlib,json,random
HERE=Path(__file__).resolve().parent;PAPER=HERE.parent;PRIOR=PAPER/'revisions/20261026_mainline'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
assert not (HERE/'protocol.json').exists()
rng=random.Random(2026102607)
rows=[]
for r in [1,2,4,8,12,24]:
    for N in [2,4,8,32]:
        for repetition in range(4):
            m,n=2*r+1,3*r
            edges=[]
            for k in range(r):edges.extend([(0,3*k),(0,3*k+1),(2*k+1,3*k+1),(2*k+1,3*k+2),(2*k+2,3*k+2),(2*k+2,3*k)])
            a=[Q(rng.choice([-5,-3,-2,-1,1,2,3,5]),8*m) for _ in range(m)]
            b=[Q(rng.choice([-7,-5,-3,-2,-1,1,2,3,5,7]),8*n) for _ in range(n)]
            assert sum(z*z for z in a)<1 and sum(z*z for z in b)<1
            u=[Q(1,m*rng.randrange(1,N+1)) for _ in range(m)]
            v=[Q(1,n*rng.randrange(1,N+1)) for _ in range(n)]
            values=[a[i]*b[j] for i,j in edges];gu=[Q(0)]*m;gv=[Q(0)]*n;loss=Q(0)
            for (i,j),value in zip(edges,values):
                residual=u[i]*v[j]-value;gu[i]+=residual*v[j];gv[j]+=residual*u[i];loss+=residual*residual/2
            p=dict(m=m,n=n,edges=edges,u=list(map(str,u)),v=list(map(str,v)),gu=list(map(str,gu)),gv=list(map(str,gv)),loss=str(loss))
            rows.append(dict(id=f'fresh_r{r}_N{N}_{repetition}',r=r,N=N,public=p,truth_array=list(map(str,values))))
with (HERE/'prospective_inputs.jsonl').open('x',encoding='utf-8') as output:
    for row in rows:output.write(json.dumps(row,sort_keys=True)+'\n')
sources=['coordinate_decoder.py','verify_coordinate.py','run_study.py','make_protocol.py']
protocol=dict(version='coordinate-propagation-2026.10.26.1',status='frozen_before_native_run',
              exploratory_cache_disclosed=True,exploratory_sha256=sha(HERE/'exploratory_cached.json'),
              development_cohort_previously_observed=True,no_claim_of_heldout_old_cohort=True,
              queries_sha256=sha(PAPER/'review4_probability_v1/queries.jsonl'),
              reference_sha256=sha(PAPER/'review4_probability_v1/outcomes.jsonl'),
              prior_outcomes_sha256=sha(PRIOR/'squareclass_outcomes.jsonl'),
              source_hashes={name:sha(HERE/name) for name in sources},
              shared_library_hashes={name:sha(PAPER/name) for name in ['exact_cactus_decoder.py','cubic_probability_decoder.py']},
              frozen_queries=1605,additional_public_inputs=46,prospective_inputs=len(rows),prospective_seed=2026102607,
              prospective_inputs_sha256=sha(HERE/'prospective_inputs.jsonl'),
              vector_oracles=500,sturm_controls=6,radical_controls=3,empty_mask_controls=2,corrupted_certificates=4,
              no_input_exclusions=True,no_cartesian_join_by_decoder=True,no_new_comparative_timing_study=True,
              candidate_class='all actual real targets; no factor budget',
              prior_art=[dict(title='Finite Domain Bounds Consistency Revisited',authors='Chiu Wo Choi; Warwick Harvey; Jimmy Ho-Man Lee; Peter J. Stuckey',year=2004,url='https://arxiv.org/abs/cs/0412021',scope='classical bounds propagation; no priority claim for that operation')])
(HERE/'protocol.json').write_text(json.dumps(protocol,indent=2)+'\n')
print(json.dumps(protocol,indent=2))
