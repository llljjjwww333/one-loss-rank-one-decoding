"""Replay the v1.2.0 supplementary checks in an isolated temporary workspace."""
from pathlib import Path
import argparse
import hashlib
import json
import subprocess
import sys
import tempfile


ROOT = Path(__file__).resolve().parent


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read(path):
    return json.loads(path.read_text(encoding='utf-8'))


def main():
    import zipfile
    import sympy
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=ROOT / 'reproduction-output-v120')
    args = parser.parse_args()
    assert sympy.__version__ == '1.14.0', 'Use the pinned SymPy 1.14.0.'
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    report = {'complete': False, 'version': 'v1.2.0', 'sympy': sympy.__version__,
              'timings_remeasured': False, 'archives': {}, 'checks': []}
    receipt = output / 'validation.json'

    def save():
        receipt.write_text(json.dumps(report, indent=2) + '\n', encoding='utf-8')

    def run(name, command, cwd):
        result = subprocess.run([sys.executable, '-X', 'utf8', *command], cwd=cwd,
                                capture_output=True, encoding='utf-8', timeout=600)
        (output / (name + '.log')).write_text(result.stdout + result.stderr, encoding='utf-8')
        assert result.returncode == 0, (name, result.stdout, result.stderr)
        return result.stdout

    save()
    with tempfile.TemporaryDirectory(prefix='rankone_v120_') as temporary:
        workspace = Path(temporary).resolve()
        for name, expected in read(ROOT / 'supplementary/manifest.json').items():
            source = ROOT / 'supplementary' / name
            assert sha(source) == expected['sha256'], name
            assert source.stat().st_size == expected['bytes'], name
            destination = workspace / source.stem
            destination.mkdir()
            with zipfile.ZipFile(source) as archive:
                assert len(archive.infolist()) == expected['members'], name
                assert archive.testzip() is None, name
                for member in archive.infolist():
                    target = (destination / member.filename).resolve()
                    assert target.is_relative_to(destination), member.filename
                    assert not ((member.external_attr >> 16) & 0o170000) == 0o120000
                archive.extractall(destination)
            report['archives'][name] = expected
        save()
        verification = workspace / 'Supplementary_verification'
        expected_theta = read(verification / 'theta_specialization_checks.json')
        run('theta-specialization', ['check_theta_specialization.py'], verification)
        actual_theta = read(verification / 'theta_specialization_checks.json')
        assert len(actual_theta['rows']) == len(expected_theta['rows']) == 4
        # Numeric residual order can differ; require a small residual and exact
        # agreement of every algebraic count, rank, identity and energy bound.
        for actual, expected in zip(actual_theta['rows'], expected_theta['rows']):
            assert sympy.Float(actual.pop('numeric_residual_max'), 50) < sympy.Float('1e-40')
            expected.pop('numeric_residual_max')
            assert actual == expected
        assert actual_theta == expected_theta
        report['checks'].append({'name': 'theta-specialization', 'passed': True,
                                 'false_points': sum(row['false_points'] for row in actual_theta['rows'])})
        save()
        print('theta-specialization PASS', flush=True)
        expected_minimal = read(verification / 'minimal_two_cycles_result.json')
        run('minimal-two-cycles', ['minimal_two_cycles.py'], verification)
        assert read(verification / 'minimal_two_cycles_result.json') == expected_minimal
        report['checks'].append({'name': 'minimal-two-cycles', 'passed': True,
                                 'gradient_candidates': expected_minimal['gradient_candidates']})
        save()
        print('minimal-two-cycles PASS', flush=True)
        ablation = workspace / 'Ablation_reproducibility'
        ablation_report = json.loads(run('ablation', ['tools/replay_ablation.py', '.'], ablation))
        assert ablation_report == {'status': 'PASS', 'timed_runs_retained': 192,
                                   'complete_lists_checked': 180,
                                   'direct_reply_and_support_arrays': 180,
                                   'pruning_traces_replayed': 72,
                                   'timings_remeasured': False}
        report['checks'].append({'name': 'ablation', 'passed': True, **ablation_report})
        report['complete'] = True
        report['scope'] = 'Exact supplementary reruns and saved-output replay; no benchmark timing repetition or independent human proof review.'
        save()
        print('ablation PASS; report:', receipt, flush=True)


if __name__ == '__main__':
    main()
