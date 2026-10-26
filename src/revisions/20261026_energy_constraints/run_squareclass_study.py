"""Replay a frozen whole cohort for an exact sufficient uniqueness certificate."""
from pathlib import Path
from collections import Counter, defaultdict
from fractions import Fraction as Q
import hashlib, json, sys
from squareclass_decoder import decode, prune, Quad

HERE = Path(__file__).resolve().parent
PAPER = HERE.parents[1]


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def controls():
    cases = []
    def check(name, groups, expected):
        keep, basis = prune(groups)
        assert keep == expected, (name, keep, expected)
        cases.append(dict(name=name, retained=keep, basis=basis))
    check('same square class, different representatives',
          [[Quad(0, 1, 2)], [Quad(0, -1, 8)]], [[0], [0]])
    check('distinct but multiplicatively dependent classes',
          [[Quad(1), Quad(0, 1, 2)], [Quad(2), Quad(0, 1, 3)],
           [Quad(3), Quad(0, 1, 6)]], [[0], [0], [0]])
    check('rational energy needs no radical label',
          [[Quad(3, 0, 2), Quad(0, 1, 2)]], [[0]])
    check('shared-class cancellation must remain possible',
          [[Quad(2), Quad(2, 1, 2)], [Quad(3), Quad(3, -1, 2)]],
          [[0, 1], [0, 1]])
    check('private class leaves no branch', [[Quad(0, 1, 2)]], [[]])
    empty = dict(m=0,n=0,edges=[],u=[],v=[],gu=[],gv=[],loss='0')
    assert decode(empty)['status'] == 'CertifiedUnique'
    empty['loss'] = '1'
    assert decode(empty)['status'] == 'CertifiedEmpty'
    return dict(squareclass_cases=cases, empty_mask_cases=2)


def run():
    protocol = json.loads((HERE/'squareclass_protocol.json').read_text())
    queries = PAPER/'review4_probability_v1/queries.jsonl'
    outcomes = PAPER/'review4_probability_v1/outcomes.jsonl'
    old = {x['id']:x['result'] for x in map(json.loads, outcomes.open())}
    records, summaries = [], defaultdict(Counter)
    report = dict(complete=False, protocol_sha256=sha(HERE/'squareclass_protocol.json'),
                  source_hashes={p.name:sha(p) for p in
                                 [HERE/'squareclass_decoder.py',Path(__file__)]},
                  library_hashes={p.name:sha(p) for p in
                                  [PAPER/'cubic_probability_decoder.py',PAPER/'exact_cactus_decoder.py']},
                  input_sha256=sha(queries),reference_sha256=sha(outcomes),
                  arithmetic_controls=controls())
    dest = HERE/'squareclass_outcomes.jsonl'
    with dest.open('w', encoding='utf-8') as stream:
        for index, query in enumerate(map(json.loads, queries.open()), 1):
            result = decode(query['public'])
            reference = old[query['id']]
            assert reference['complete']
            if result['status'] == 'CertifiedUnique':
                assert reference['complete_count'] == 1, query['id']
            if result['status'] == 'CertifiedEmpty':
                assert reference['complete_count'] == 0, query['id']
            row = dict(id=query['id'], target=query['target'], N=query['N'],
                       reference_count=reference['complete_count'], result=result)
            stream.write(json.dumps(row, sort_keys=True)+'\n')
            records.append(row)
            summaries[f"{query['target']}_N{query['N']}"][result['status']] += 1
            if index % 100 == 0:
                stream.flush()
                print('queries', index, dict(Counter(x['result']['status'] for x in records)), flush=True)
    assert len(records) == len(old) == 1605
    negative = []
    for label, directory in [('uniform',HERE/'uniform_linux/public'), ('lattice',HERE/'lattice_public')]:
        for path in sorted(directory.glob('*.json')):
            result = decode(json.loads(path.read_text()))
            negative.append(dict(cohort=label,id=path.stem,public_sha256=sha(path),result=result))
    assert Counter(x['cohort'] for x in negative) == {'uniform':14,'lattice':32}
    report.update(complete=True, queries=len(records),
                  summary={k:dict(v) for k,v in sorted(summaries.items())},
                  statuses=dict(Counter(x['result']['status'] for x in records)),
                  unresolved_reference_counts=dict(Counter(
                      'singleton' if x['reference_count']==1 else 'ambiguous' if x['reference_count']>1 else 'empty'
                      for x in records if x['result']['status']=='Unresolved')),
                  negative_controls=negative, outcomes_sha256=sha(dest),
                  no_exclusions=True, truth_passed_to_solver=False)
    (HERE/'squareclass_validation.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k not in
                      ['negative_controls','summary','arithmetic_controls']},indent=2),flush=True)


if __name__ == '__main__':
    run()
