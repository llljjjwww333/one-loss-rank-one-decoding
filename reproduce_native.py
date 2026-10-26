"""Reconstruct a fresh Ubuntu workspace and repeat the native protocols."""
from pathlib import Path
import argparse,hashlib,json,os,shutil,subprocess,sys,urllib.request
REPO=Path(__file__).resolve().parent
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def run():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--workspace',type=Path,required=True)
    ap.add_argument('--package-cache',type=Path,action='append',default=[])
    ap.add_argument('--prepare-only',action='store_true');args=ap.parse_args()
    assert sys.platform.startswith('linux'),'Native repetition requires Linux/Ubuntu; saved-certificate replay works on Windows.'
    work=args.workspace.resolve();assert not work.exists(),'Use a fresh workspace; prior experiments are preserved.';work.mkdir(parents=True)
    manifest=json.loads((REPO/'snapshot_manifest.json').read_text())
    study=Path('paper/cactus_global_loss/independent_validation_v1')
    excluded={'results.json','phase_results.json','validation.json','phase_validation.json','coverage_validation.json'}
    for e in manifest['files']:
        src=REPO/e['path'];assert sha(src)==e['sha256']
        target=Path(e['target'])
        if target.is_relative_to(study) and ('outputs' in target.parts or target.name in excluded):continue
        dst=(work/target).resolve();assert dst.is_relative_to(work);dst.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(src,dst)
    caches={}
    for folder in args.package_cache:
        for p in folder.resolve().rglob('*.deb'):caches[sha(p)]=p
    provision=[]
    for name in ['jsc_singular_4.3.2_v2','jsc_msolve_0.6.5_v2']:
        local=work/'tools/third_party'/name;runtime=local/'runtime';runtime.mkdir(parents=True)
        original=json.loads((REPO/'native'/name/'runtime-manifest.json').read_text())
        for package,meta in original['packages'].items():
            src=caches.get(meta['sha256'])
            if src is None:
                src=local/Path(meta['repository_file']).name
                urllib.request.urlretrieve('https://archive.ubuntu.com/ubuntu/'+meta['repository_file'],src)
            assert sha(src)==meta['sha256'] and src.stat().st_size==meta['bytes'],package
            subprocess.run(['dpkg-deb','-x',str(src),str(runtime)],check=True)
        binary=runtime/'usr/bin'/('Singular' if 'singular' in name else 'msolve')
        assert sha(binary)==original['binary_sha256']
        shutil.copy2(REPO/'native'/name/'runtime-manifest.json',local/'runtime-manifest.json')
        provision.append(dict(runtime=name,binary_sha256=sha(binary),packages=len(original['packages'])))
    report=dict(complete=True,scope='Runtime provisioning only; no scientific experiment has run yet.',runtimes=provision)
    (work/'runtime_provision_validation.json').write_text(json.dumps(report,indent=2)+'\n')
    if args.prepare_only:print(json.dumps(report));return
    env=os.environ.copy();env.update(OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1')
    source=work/study
    for script,mode in [('run_study.py','run'),('verify_study.py','final'),('phase_study.py','run'),('audit_phases.py',None)]:
        command=[sys.executable,'-u',str(source/script)]+([mode] if mode else [])
        subprocess.run(command,cwd=work,env=env,check=True)
    print('Native comparison, independent coverage and phase repetition complete:',work)
if __name__=='__main__':run()
