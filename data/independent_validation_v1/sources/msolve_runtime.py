"""Unpack a hash-checked Ubuntu msolve runtime under this workspace only."""
from pathlib import Path
import hashlib,json,os,subprocess,sys
ROOT=next(p for p in Path(__file__).resolve().parents if (p/'tools/third_party').is_dir())
LOCAL=ROOT/'tools/third_party/jsc_msolve_0.6.5_v2'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def environment():
    env=os.environ.copy()
    env['LD_LIBRARY_PATH']=str(LOCAL/'runtime/usr/lib/x86_64-linux-gnu')+':'+str(ROOT/'tools/third_party/jsc_singular_4.3.2_v2/runtime/usr/lib/x86_64-linux-gnu')
    env['OMP_NUM_THREADS']='1'
    return env
def binary():return LOCAL/'runtime/usr/bin/msolve'
def setup():
    assert os.name=='posix' and str(ROOT)=='/mnt/f/autoresearch'
    assert not LOCAL.exists()
    (LOCAL/'packages').mkdir(parents=True);(LOCAL/'runtime').mkdir()
    metadata={}
    for name,version in [('msolve','0.6.5-1build2'),('libmsolve-0.6.5','0.6.5-1build2'),('libgomp1','14-20240412-0ubuntu1')]:
        raw=subprocess.check_output(['apt-cache','show',name+'='+version],text=True)
        info=dict(line.split(': ',1) for line in raw.splitlines() if ': ' in line and not line.startswith(' '))
        metadata[name]=dict(version=version,sha256=info['SHA256'],bytes=int(info['Size']),repository_file=info['Filename'])
    (LOCAL/'package-plan.json').write_text(json.dumps(metadata,indent=2)+'\n')
    out=subprocess.run(['apt-get','download']+[name+'='+info['version'] for name,info in metadata.items()],cwd=LOCAL/'packages',capture_output=True,text=True)
    (LOCAL/'download.stdout.txt').write_text(out.stdout);(LOCAL/'download.stderr.txt').write_text(out.stderr)
    assert out.returncode==0,out.stderr
    for p in (LOCAL/'packages').glob('*.deb'):
        name=subprocess.check_output(['dpkg-deb','-f',str(p),'Package'],text=True).strip()
        assert sha(p)==metadata[name]['sha256'] and p.stat().st_size==metadata[name]['bytes']
        subprocess.run(['dpkg-deb','-x',str(p),str(LOCAL/'runtime')],check=True)
    finalize(metadata)

def finalize(metadata):
    out=subprocess.run([str(binary()),'-h'],env=environment(),capture_output=True,text=True)
    (LOCAL/'help.stdout.txt').write_text(out.stdout);(LOCAL/'help.stderr.txt').write_text(out.stderr)
    assert out.returncode in (0,1) and 'msolve library' in out.stdout and not out.stderr,(out.returncode,out.stderr)
    info=dict(complete=True,packages=metadata,binary_sha256=sha(binary()),scope='Local extraction only; no system package installation',
              source_sha256=sha(Path(__file__)),linked_libraries=subprocess.check_output(['ldd',str(binary())],env=environment(),text=True))
    (LOCAL/'runtime-manifest.json').write_text(json.dumps(info,indent=2)+'\n')
    print(json.dumps(info,indent=2))
if __name__=='__main__':
    if '--finalize' in sys.argv:
        metadata=json.loads((LOCAL/'package-plan.json').read_text())
        for p in (LOCAL/'packages').glob('*.deb'):
            name=subprocess.check_output(['dpkg-deb','-f',str(p),'Package'],text=True).strip()
            assert sha(p)==metadata[name]['sha256']
        finalize(metadata)
    else:setup()
