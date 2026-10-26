# One Loss Resolves Signed Rank-One Gradient Fibers

Exact reconstruction software and synthetic computational data accompanying
*One Loss Resolves Signed Rank-One Gradient Fibers: Exact Decoding on Cycle and Theta Blocks*.
This release corresponds to manuscript version `jsc-coordinate-2026.10.26.6`.

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
energy-coordinate certificates. The full run executes 14 groups of checks,
including all original input cases, failed sufficient conditions and saved
solver timeouts. Each run checks SHA-256 for every snapshot file, creates
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

The repository contains software and synthetic data. The unpublished
manuscript, private review reports, account credentials and author-approval
records are not distributed here.

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

Generators and verifiers share SymPy rational arithmetic and the incidence
reduction. Classical root isolation, rational univariate representations,
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
