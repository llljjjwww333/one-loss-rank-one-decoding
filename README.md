# One Loss Resolves Signed Rank-One Gradient Fibers

Exact reconstruction software and synthetic computational data accompanying
*One Loss Resolves Signed Rank-One Gradient Fibers: Exact Decoding on Cycle and Theta Blocks*.
This release corresponds to manuscript version `jsc-certified-2026.10.29.17`.

Software release **v1.3.0** adds the indexed extended report and all of its
LaTeX sources while retaining every earlier scientific snapshot unchanged.
The version-specific v1.3.0 DOI is **https://doi.org/10.5281/zenodo.23754950**.
The concept DOI is **https://doi.org/10.5281/zenodo.23660117**.
The [v1.2.0 archive](https://doi.org/10.5281/zenodo.23753925),
[v1.1.0 archive](https://doi.org/10.5281/zenodo.23664952) and
[v1.0.0 archive](https://doi.org/10.5281/zenodo.23660118) remain unchanged.

## Extended report and archive map (v1.3.0)

```bash
python reproduce_extended.py
```

This verifies the ZIP hash, safely extracts the 15-page extended report and
its source, and checks every member against its embedded SHA-256 manifest.
Use `--output FRESH_FOLDER` to select an extraction directory. The first
page of the report maps each section to the data and reproduction command.
Its LaTeX source rebuilds independently with `source/build_report.ps1`.

`supplementary/Extended_materials_v1.3.0.zip` contains expanded matched
examples, the integer-lattice uniqueness refinement, separation-certificate
coefficient heights, all per-input comparison and ablation tables, exact-check
and finite-grid outcomes, and the coordinate-certificate phase figure.
The journal retains the core mathematical proofs in its essential supplement;
the main article and essential supplement are not part of this software release.
The report source was independently extracted and rebuilt: all 15 pages had
equal text, dimensions and rendered pixels. This verifies artifact consistency,
not independent mathematical proof or performance claims.

The complete earlier code/data snapshot remains byte-preserved, including
all timeout, fallback, empty and multiple controls. The original 192 ablation
runs, 180 completed lists, 12 timeouts and 72 traces are unchanged.

## Added checks and same-input ablation (v1.2.0)

```bash
python reproduce_supplementary.py
```

Use Python 3.11+ and SymPy 1.14.0. This entry point verifies the two frozen
ZIPs, extracts them into an isolated workspace, reruns both standalone exact
checkers, and replays all saved ablation outputs. The report and logs are
written to `reproduction-output-v120/`. It checks 112 false theta points
for p=3,4,5,6, a two-six-cycle base point with nine gradient candidates
and complete-reply layer sizes 1/4/4, and the saved 192-run ablation.

The ablation holds 16 public inputs fixed and varies propagation on/off
and enumeration/meet-in-the-middle, with three repetitions each. Completion
counts are 42/48, 48/48, 42/48 and 48/48; the remaining 12 runs are censored
at 30 seconds. All completed outputs, empty and multiple lists, propagation
fallbacks and timeouts are retained. Saved-output replay checks 180 complete
lists, 180 returned arrays and 72 pruning traces. This cohort shows no
additional completed cases from propagation over MITM alone.

- `supplementary/Supplementary_verification.zip`: standalone checkers and
  exact receipts for the theta specialization and minimal two-cycle example.
- `supplementary/Ablation_reproducibility.zip`: frozen source, inputs,
  references, all timed records, completed outputs and pruning traces.
- `supplementary/manifest.json`: attachment sizes and SHA-256 hashes.

To repeat the ablation timings on Linux/WSL, extract the ablation ZIP into
a fresh folder and follow its README. The saved-output replay does not
repeat benchmark timings. The exact checkers and ablation use SymPy;
their checks do not constitute independent human proof review or a new
independent algebra-system comparison.

## Reproduce

The byte-preserving replay has been validated on Windows with Python 3.11.
Other platforms have not yet been validated; some original outcome hashes
include platform-dependent text line endings. Use the pinned packages:

```bash
python -m pip install -r requirements.txt
python reproduce.py
python reproduce.py --full
```

The quick run checks the four-cycle grid law, RUR arithmetic controls and
energy-coordinate certificates. The full run executes 16 groups of checks,
including all original input cases, failed sufficient conditions, saved
solver timeouts, independent native certificates and saved phase records. Each run checks SHA-256 for every snapshot file, creates
an isolated temporary workspace, and compares scientific results with the
saved outcomes. It writes logs and a validation report to
`reproduction-output/`. The original snapshots remain unchanged.

Lattice runtime measurements and machine-description strings are recorded
separately from scientific decisions. The full run audits saved Linux
solver outcomes; it does not rerun or reinterpret native solver timings.

## Contents

- `src/`: original implementation bytes, including exact cactus and block
  decoders, real-algebraic coverage checks and energy-coordinate propagation.
- `data/six_cycle_grids/`: all 1,605 public queries and complete-fiber references.
- `data/energy_coordinates/`: propagation outcomes and 96 additional public
  inputs, fixed before their run.
- `data/theta_public/`, `data/rur_public_v3/`: public inputs for block and
  algebraic decoding studies.
- `data/uniform_linux/`: complete saved comparison results and protocols,
  including nine msolve timeouts.
- `data/lattice_public/`, `data/results/`: lattice cases, certificate-degree
  witnesses, protocols and complete outcomes.
- `snapshot_manifest.json`: file hashes and compatibility paths. These
  paths preserve the original programs' import and input conventions;
  they are recreated only inside the temporary workspace.

The repository contains software, synthetic data and the extended computational report.
The main journal article, essential journal supplement, private reviews, account
credentials and author-approval records are not distributed here.

## Main computational results

Energy-coordinate propagation certifies 1,486 singleton replies among all
1,605 six-cycle inputs. Its 119 unresolved replies include 54 singleton
references and all 65 ambiguous references. On 96 additional inputs,
84 are certified and 12 remain unresolved. These are certificate counts,
not estimates of a general uniqueness probability or average complexity.

Other checks cover the exact four-cycle collision law, 1,821 distinct
algebraic collision witnesses, block normalization and real-algebraic
list coverage, including nonradical finite ideals. These are exact local
computational checks, not independent human proof verification.

The primary generators and verifiers share SymPy rational arithmetic and
the incidence reduction. The additional native path independently constructs
observed-edge equations, uses Singular for Groebner/radical computations,
and uses separate standard-library rational and Sturm code for exact checks.
It does not import SymPy or the primary incidence reduction. Classical root isolation, rational univariate representations,
meet-in-the-middle and constraint propagation are used; the repository
does not claim priority for those general algorithms.

## AI assistance

OpenAI Codex assisted with implementation, mathematical exploration,
computational checks and documentation. The provided exact checks and
their limitations are described above. Publishing this software does not
represent final author approval of a journal manuscript.

## Citation and license

Use this repository's version and, once available, the Zenodo software DOI.
The software DOI is separate from any future journal-article DOI.
Metadata currently records the confirmed creator 楼俊文, without inventing
an English name, affiliation, ORCID or a journal DOI.

Original software and accompanying synthetic data are released under MIT.
Third-party dependencies retain their own licenses. No rights in the
separate manuscript are granted by this software license.

## Independent validation and timing protocols (v1.1.0)

`data/independent_validation_v1/` contains:

- All 37 native coverage jobs on 34 distinct public inputs, including
  irrational points, nonradical ideals and three optional budget jobs.
- Exact rational modular identities comparing 35 RUR representations,
  real-embedding/support checks and two direct empty-list controls.
- All 126 Linux comparison runs: cycle join 42 complete; msolve 33
  complete and nine timeouts; observed-edge Singular 27 complete and
  fifteen timeouts. All completed lists agree with their references.
- All 288 phase measurements (96 inputs, three repetitions) and their
  unchanged exact results: 84 singleton certificates and 12 unresolved
  replies in each repetition.
- Pinned input/source hashes, raw native systems and outputs, protocol
  definitions, recorded machine/affinity information and three rejected
  corrupted certificates.

The full replay audits these saved native certificates and timing records
without requiring Linux, Singular or msolve. It does not repeat timings.
To repeat the timed native study, use Ubuntu 24.04, Python 3.12.3 and
SymPy 1.14.0, provision the pinned Singular 4.3.2 and msolve 0.6.5
runtimes described by the archived native manifests, and reconstruct the
compatibility paths in `snapshot_manifest.json` in a fresh workspace.
Run `run_study.py run` before `phase_study.py run`, with results absent;
the recorded after-ready and phase timing scopes differ and are not pooled.
The primary library folder must be on PYTHONPATH. Runtimes are external
dependencies; executable binaries are not redistributed here.

All inputs are synthetic. Certificate frequency, exact list coverage,
single-host timings and mathematical correctness are distinct claims.

An executable native repetition entry point is included:

```bash
python reproduce_native.py --workspace /path/to/fresh-workspace
```

It reconstructs the snapshot paths, downloads the pinned Ubuntu packages,
checks package and executable hashes, extracts them locally without system
installation, and runs the comparison, exact verification and phase study
in order. Use `--package-cache /path/to/debs` for an offline cache, or
`--prepare-only` to test provisioning without claiming an experiment.
Use Ubuntu 24.04 amd64 with the Python dependencies installed. Native
runs can take several minutes because all thirty-second timeouts remain.

The certificate and phase plot is reproducible from the saved profiles:

```bash
python -m pip install -r requirements-plot.txt
python plot_profiles.py --output-dir plots
```
