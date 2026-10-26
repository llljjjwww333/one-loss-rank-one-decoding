"""Replay exact certificates without running the generator's propagation.

Root coverage is checked by a rational Sturm chain, not by refactorization.
The checker shares SymPy rational arithmetic and the incidence reduction.
"""
from fractions import Fraction as Q
from math import prod
import sympy as s
from coordinate_decoder import reduction, Quad, qsqrt, rational_candidate


def real_root_count(poly):
    if poly.degree() <= 0:
        return 0
    f = poly.exquo(poly.gcd(poly.diff()))
    chain = [f, f.diff()]
    while not chain[-1].is_zero:
        remainder = -chain[-2].rem(chain[-1])
        if remainder.is_zero:
            break
        chain.append(remainder)
    def variations(at_minus_infinity):
        signs = [s.sign(f.LC())*((-1)**f.degree() if at_minus_infinity else 1) for f in chain]
        return sum(a != b for a, b in zip(signs, signs[1:]))
    return variations(True)-variations(False)


def same_root(a, b):
    if a.a != b.a:
        return False
    if not a.b or not b.b:
        return a.b == b.b
    return (a.b > 0) == (b.b > 0) and a.b*a.b*a.d == b.b*b.b*b.d


def verify(public, certificate):
    p, cycles, weights, base, polynomials, required = reduction(public)
    status, reason = certificate['status'], certificate.get('reason')
    assert status in ['CertifiedEmpty', 'CertifiedUnique', 'Unresolved']
    assert certificate['cartesian_tuples_enumerated'] == 0
    assert certificate['candidate_class'] == 'all actual real targets; no factor budget'
    assert certificate['truth_passed_to_solver'] is False
    if base is None:
        assert status == 'CertifiedEmpty' and reason == 'inconsistent margins'
        return True
    if any(poly.degree() == 0 and not poly.is_zero for poly in polynomials):
        assert status == 'CertifiedEmpty' and reason == 'nonzero constant cycle'
        return True
    if 'cycles' not in certificate:
        assert status == 'Unresolved'
        assert (any(poly.is_zero for poly in polynomials) and reason == 'zero cycle polynomial'
                or any(poly.degree() > 3 for poly in polynomials) and reason == 'cycle degree exceeds three'
                or reason == 'irreducible cubic factor' and any(
                    factor.degree() > 2 for poly in polynomials for factor, _ in s.factor_list(poly)[1]))
        return True
    records = certificate['cycles']
    assert len(records) == len(cycles)
    reps = list(map(Q, certificate['representatives']))
    for k, rep in enumerate(reps):
        assert rep > 0 and qsqrt(rep) is None
        assert all(qsqrt(rep/other) is None for other in reps[:k])
    roots, vectors = [], []
    for cycle, poly, record in zip(cycles, polynomials, records):
        assert record['edges'] == cycle and record['coefficients'] == list(map(str, poly.all_coeffs()))
        local = [Quad(*encoded) for encoded in record['roots']]
        assert [root.encode() for root in local] == record['roots']
        assert len(local) == real_root_count(poly)
        row, encoded_energies = [], []
        for j, root in enumerate(local):
            assert all(not same_root(root, other) for other in local[:j])
            residual = Quad()
            for coef in poly.all_coeffs():
                residual = residual*root+Q(coef)
            assert residual.zero()
            values = [Quad(Q(base[e]/weights[e]))+root*Q((-1)**k/weights[e])
                      for k, e in enumerate(cycle)]
            energy = sum((z*z for z in values), Quad())
            encoded_energies.append(energy.encode())
            vec = [energy.a]+[Q(0)]*len(reps)
            if energy.b:
                matching = [(k, qsqrt(energy.d/rep)) for k, rep in enumerate(reps)
                            if qsqrt(energy.d/rep) is not None]
                assert len(matching) == 1
                k, ratio = matching[0]
                vec[k+1] = energy.b*ratio
            row.append(vec)
        assert record['energies'] == encoded_energies
        roots.append(local)
        vectors.append(row)
    assert certificate['vectors'] == [[[str(z) for z in vec] for vec in row] for row in vectors]
    target = [required]+[Q(0)]*len(reps)
    assert certificate['target'] == list(map(str, target))
    domains = [list(range(len(row))) for row in vectors]
    deletions = 0
    for cuts in certificate['trace']:
        assert cuts and all(domains)
        removed = set()
        for cut in cuts:
            i, j, k = cut['cycle'], cut['branch'], cut['coordinate']
            assert 0 <= i < len(domains) and j in domains[i] and 0 <= k < len(target)
            assert (i, j) not in removed
            # Recompute each bound directly from all OTHER cycle domains.
            lo = sum(min(vectors[a][b][k] for b in dom) for a, dom in enumerate(domains) if a != i)
            hi = sum(max(vectors[a][b][k] for b in dom) for a, dom in enumerate(domains) if a != i)
            need = target[k]-vectors[i][j][k]
            assert list(map(str, (need, lo, hi))) == [cut['need'], cut['lower'], cut['upper']]
            assert need < lo or need > hi
            removed.add((i, j))
        domains = [[j for j in dom if (i, j) not in removed] for i, dom in enumerate(domains)]
        deletions += len(removed)
    assert certificate['domains'] == domains
    assert certificate['cartesian_before'] == prod(map(len, roots))
    assert certificate['cartesian_after'] == prod(map(len, domains))
    assert certificate['rounds'] == len(certificate['trace']) and certificate['deletions'] == deletions
    if any(not dom for dom in domains):
        assert status == 'CertifiedEmpty' and reason == 'exhausted cycle domain'
        return True
    if any(len(dom) != 1 for dom in domains):
        assert status == 'Unresolved' and reason == 'multiple retained branches'
        return True
    values, valid = rational_candidate(p, cycles, weights, base, roots, domains)
    if values is None:
        assert status == 'Unresolved' and reason == 'retained nonrational root'
    elif not valid:
        assert status == 'CertifiedEmpty' and reason == 'final reply/support test failed'
    else:
        assert status == 'CertifiedUnique' and certificate['values'] == list(map(str, values))
        assert certificate['reply_recomputed'] and certificate['actual_support_checked']
    return True
