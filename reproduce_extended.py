"""Extract and verify the complete v1.3.0 extended report without changing snapshots."""
from pathlib import Path
import argparse, hashlib, json, zipfile
ROOT=Path(__file__).resolve().parent
EXPECTED='380a4173755443d1daaffa75107fdca4e28af7368e76b406068c9dacba3148d1'
def main():
    parser=argparse.ArgumentParser(); parser.add_argument('--output',type=Path,default=ROOT/'reproduction-output-v130')
    args=parser.parse_args(); archive=ROOT/'supplementary/Extended_materials_v1.3.0.zip'
    assert hashlib.sha256(archive.read_bytes()).hexdigest()==EXPECTED
    out=args.output.resolve(); assert not out.exists(),'Use a fresh output folder.'; out.mkdir(parents=True)
    with zipfile.ZipFile(archive) as z:
        assert z.testzip() is None
        for member in z.infolist(): assert (out/member.filename).resolve().is_relative_to(out)
        z.extractall(out)
    files=json.loads((out/'manifest.json').read_text())
    for name,digest in files.items():
        assert hashlib.sha256((out/name).read_bytes()).hexdigest()==digest,name
    report=dict(status='PASS',verified_members=len(files),archive_sha256=EXPECTED,source_build_command='pwsh -File source/build_report.ps1',scope='Archive extraction and byte integrity; saved scientific results and timings are unchanged.')
    (out/'extraction_validation.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(report,indent=2))
if __name__=='__main__':main()
