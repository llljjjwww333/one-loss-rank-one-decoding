"""Reproduce exact checks in an isolated workspace without altering snapshots."""
from pathlib import Path
import argparse,hashlib,json,subprocess,sys,tempfile,time
ROOT=Path(__file__).resolve().parent
CURRENT='paper/cactus_global_loss/revisions/20261026_energy_constraints'
TASKS=[
 ('grid','check_grid_law.py','grid_validation.json'),
 ('lattice','run_lattice_study.py','lattice_outcomes.json'),
 ('saved-linux','analyze_uniform.py','uniform_validation.json'),
 ('private-classes','run_squareclass_study.py','squareclass_validation.json'),
 ('private-effect','analyze_squareclass.py','squareclass_effect_validation.json'),
 ('unit-blocks','check_unit_blocks.py','unit_block_validation.json'),
 ('scaled-blocks','check_scaled_blocks.py','scaled_block_validation.json'),
 ('theta-lists','audit_theta.py','theta_audit_validation.json'),
 ('rur-arithmetic','check_rur_arithmetic.py','rur_arithmetic_validation.json'),
 ('rur-coverage','audit_rur.py','rur_validation.json'),
 ('degree-symbolics','check_degree_symbolics.py','degree_symbolic_validation.json'),
 ('degree-binary','check_certificate_degree.py','degree_validation.json'),
 ('degree-ternary','check_ternary_degree.py','ternary_degree_validation.json'),
 ('energy-coordinates','../../energy_constraints_v1/replay_study.py','../../energy_constraints_v1/replay_validation.json'),
]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def semantic(name,obj):
 if name=='lattice':
  obj=json.loads(json.dumps(obj));obj.pop('python');obj.pop('platform')
  for row in obj['records']:
   row.pop('preprocessing_seconds')
   for method in row['methods'].values():method.pop('combination_seconds',None)
 return obj
def main():
 parser=argparse.ArgumentParser();parser.add_argument('--full',action='store_true');parser.add_argument('--output',default='reproduction-output');args=parser.parse_args()
 output=Path(args.output).resolve();output.mkdir(parents=True,exist_ok=True)
 manifest=json.loads((ROOT/'snapshot_manifest.json').read_text())
 work=Path(tempfile.mkdtemp(prefix='rank_one_replay_'));seen=set()
 for entry in manifest['files']:
  source=ROOT/entry['path'];assert sha(source)==entry['sha256'],entry['path']
  target=(work/entry['target']).resolve();assert target.is_relative_to(work.resolve())
  assert str(target) not in seen;seen.add(str(target));target.parent.mkdir(parents=True,exist_ok=True);target.write_bytes(source.read_bytes())
 current=work/CURRENT;(current/'manuscript').mkdir(exist_ok=True)
 selected=TASKS if args.full else [t for t in TASKS if t[0] in ['grid','rur-arithmetic','energy-coordinates']]
 report=dict(complete=False,mode='full' if args.full else 'quick',manifest_sha256=sha(ROOT/'snapshot_manifest.json'),snapshot_files_verified=len(seen),tasks=[],native_solver_timings_repeated=False)
 receipt=output/'validation.json'
 receipt.write_text(json.dumps(report,indent=2)+'\n')
 for name,script,result in selected:
  result_path=current/result;before=json.loads(result_path.read_text());start=time.monotonic()
  run=subprocess.run([sys.executable,'-X','utf8',str(current/script)],capture_output=True,text=True,encoding='utf-8',timeout=600)
  (output/(name+'.log')).write_text(run.stdout+run.stderr,encoding='utf-8')
  assert run.returncode==0,(name,run.stdout,run.stderr)
  after=json.loads(result_path.read_text());assert semantic(name,before)==semantic(name,after),name
  row=dict(name=name,passed=True,scientific_results_equal=True,elapsed_seconds=time.monotonic()-start)
  report['tasks'].append(row);receipt.write_text(json.dumps(report,indent=2)+'\n');print(name,'PASS',flush=True)
 report['complete']=True;report['scope']='exact local replay using shared SymPy and incidence reduction; saved Linux outcomes audited, timed native solvers not rerun'
 receipt.write_text(json.dumps(report,indent=2)+'\n');print('All',len(selected),'tasks passed; report:',receipt)
if __name__=='__main__':main()
