"""Separate nontrivial pruning from already-singleton gradient fibers."""
from pathlib import Path
import hashlib, json
HERE=Path(__file__).resolve().parent


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def run():
    path=HERE/'squareclass_outcomes.jsonl'
    previous=HERE.parents[1]/'review4_probability_v1/outcomes.jsonl'
    old={x['id']:x['result'] for x in map(json.loads,previous.open())}
    rows=list(map(json.loads,path.open()))
    successful=[x for x in rows if x['result']['status']=='CertifiedUnique']
    nontrivial=[x for x in successful if old[x['id']]['gradient_count']>1]
    assert len(rows)==len(old)==1605 and len(successful)==708 and len(nontrivial)==583
    assert all(x['result']['cartesian_after']==1 for x in successful)
    assert sum(x['result']['cartesian_before']>1 for x in successful)==583
    largest=max(old[x['id']]['gradient_count'] for x in nontrivial)
    assert largest==177147
    result=dict(complete=True,queries=1605,certified_singletons=708,
                already_singleton_gradient_fibers=125,multiple_gradient_candidates_pruned=583,
                largest_gradient_fiber_certified=largest,
                largest_examples=[x['id'] for x in nontrivial if old[x['id']]['gradient_count']==largest],
                results_sha256=sha(path),reference_sha256=sha(previous),verifier_sha256=sha(Path(__file__)),
                scope='whole frozen cohort; distinguishes certification frequency from nontrivial pruning; no timing claim')
    (HERE/'squareclass_effect_validation.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result))


if __name__=='__main__':
    run()
