"""Exact finite-domain propagation in rational energy coordinates.

Public input only. This is a sufficient empty/singleton certificate, with
Unresolved for every unsupported or unfinished case. It uses classical
finite-domain interval propagation; no Cartesian-product tuples are joined.
"""
from pathlib import Path
from fractions import Fraction as Q
from math import isqrt, prod
import sys
import sympy as s

PAPER = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PAPER))
from cubic_probability_decoder import Quad, real_roots
from exact_cactus_decoder import canonical_public, cycle_basis, forest_base, actual_budget


def qsqrt(value):
    if value < 0:
        return None
    a, b = isqrt(value.numerator), isqrt(value.denominator)
    return Q(a, b) if a*a == value.numerator and b*b == value.denominator else None


def reduction(public):
    p = canonical_public(public)
    m, n, edges = (p[key] for key in ('m', 'n', 'edges'))
    u, v, gu, gv = (p[key] for key in ('u', 'v', 'gu', 'gv'))
    cycles = cycle_basis(m, n, edges)
    weights = [u[i]*v[j] for i, j in edges]
    rhs = [-u[i]*gu[i] for i in range(m)]+[-v[j]*gv[j] for j in range(n)]
    for k, (i, j) in enumerate(edges):
        rhs[i] += weights[k]**2
        rhs[m+j] += weights[k]**2
    base = forest_base(m, n, edges, cycles, rhs)
    if base is None:
        return p, cycles, weights, None, [], None
    t = s.Symbol('t')
    polynomials = []
    for cycle in cycles:
        z = [base[e]+(-1)**j*t for j, e in enumerate(cycle)]
        polynomials.append(s.Poly(s.prod(z[::2])-s.prod(z[1::2]), t, domain=s.QQ))
    covered = {e for cycle in cycles for e in cycle}
    forest_energy = sum((Q(base[e]/weights[e])**2 for e in range(len(edges)) if e not in covered), Q(0))
    total = Q(2*p['loss']+sum(w*w for w in weights)
              -sum(a*b for a, b in zip(u, gu))-sum(a*b for a, b in zip(v, gv)))
    return p, cycles, weights, base, polynomials, total-forest_energy


def coordinates(groups):
    reps, labels = [], []
    for group in groups:
        row = []
        for energy in group:
            if not energy.b:
                row.append((None, Q(0)))
                continue
            for k, representative in enumerate(reps):
                ratio = qsqrt(energy.d/representative)
                if ratio is not None:
                    break
            else:
                k, ratio = len(reps), Q(1)
                reps.append(energy.d)
            row.append((k, energy.b*ratio))
        labels.append(row)
    vectors = []
    for group, row in zip(groups, labels):
        vector_row = []
        for energy, (k, b) in zip(group, row):
            vec = [energy.a]+[Q(0)]*len(reps)
            if k is not None:
                vec[k+1] = b
            vector_row.append(vec)
        vectors.append(vector_row)
    return reps, vectors


def propagate(vectors, target):
    domains = [list(range(len(row))) for row in vectors]
    trace = []
    while all(domains):
        lows = [[min(vectors[i][j][k] for j in dom) for k in range(len(target))]
                for i, dom in enumerate(domains)]
        highs = [[max(vectors[i][j][k] for j in dom) for k in range(len(target))]
                 for i, dom in enumerate(domains)]
        totals_low = [sum(row[k] for row in lows) for k in range(len(target))]
        totals_high = [sum(row[k] for row in highs) for k in range(len(target))]
        cuts, updated = [], []
        for i, dom in enumerate(domains):
            keep = []
            for j in dom:
                for k, required in enumerate(target):
                    need = required-vectors[i][j][k]
                    lo, hi = totals_low[k]-lows[i][k], totals_high[k]-highs[i][k]
                    if need < lo or need > hi:
                        cuts.append(dict(cycle=i, branch=j, coordinate=k,
                                         need=str(need), lower=str(lo), upper=str(hi)))
                        break
                else:
                    keep.append(j)
            updated.append(keep)
        if not cuts:
            break
        trace.append(cuts)
        domains = updated
    return domains, trace


def rational_candidate(p, cycles, weights, base, roots, domains):
    chosen = [row[dom[0]] for row, dom in zip(roots, domains)]
    if any(root.b for root in chosen):
        return None, None
    values = [Q(z/w) for z, w in zip(base, weights)]
    for cycle, root in zip(cycles, chosen):
        for j, e in enumerate(cycle):
            values[e] += root.a*Q((-1)**j/weights[e])
    gu, gv = [Q(0)]*p['m'], [Q(0)]*p['n']
    loss = Q(0)
    for (i, j), w, value in zip(p['edges'], weights, values):
        residual = Q(w)-value
        gu[i] += residual*Q(p['v'][j])
        gv[j] += residual*Q(p['u'][i])
        loss += residual**2/2
    reply_ok = gu == list(map(Q, p['gu'])) and gv == list(map(Q, p['gv'])) and loss == Q(p['loss'])
    actual, _ = actual_budget(p['m'], p['n'], p['edges'], list(map(s.Rational, values)), False)
    return values, bool(reply_ok and actual)


def decode(public):
    p, cycles, weights, base, polynomials, required = reduction(public)
    out = dict(candidate_class='all actual real targets; no factor budget',
               truth_passed_to_solver=False, cartesian_tuples_enumerated=0)
    if base is None:
        return dict(out, status='CertifiedEmpty', reason='inconsistent margins')
    if any(poly.degree() == 0 and not poly.is_zero for poly in polynomials):
        return dict(out, status='CertifiedEmpty', reason='nonzero constant cycle')
    if any(poly.is_zero for poly in polynomials):
        return dict(out, status='Unresolved', reason='zero cycle polynomial')
    if any(poly.degree() > 3 for poly in polynomials):
        return dict(out, status='Unresolved', reason='cycle degree exceeds three')
    roots, energies = [], []
    for cycle, poly in zip(cycles, polynomials):
        try:
            local = real_roots(poly)
        except ValueError:
            return dict(out, status='Unresolved', reason='irreducible cubic factor')
        roots.append(local)
        energies.append([])
        for root in local:
            values = [Quad(Q(base[e]/weights[e]))+root*Q((-1)**j/weights[e])
                      for j, e in enumerate(cycle)]
            energies[-1].append(sum((z*z for z in values), Quad()))
    reps, vectors = coordinates(energies)
    target = [required]+[Q(0)]*len(reps)
    domains, trace = propagate(vectors, target)
    out.update(cycles=[dict(edges=cycle, coefficients=list(map(str, poly.all_coeffs())),
                            roots=[z.encode() for z in local], energies=[z.encode() for z in group])
                       for cycle, poly, local, group in zip(cycles, polynomials, roots, energies)],
               representatives=list(map(str, reps)), vectors=[[[str(z) for z in vec] for vec in row] for row in vectors],
               target=list(map(str, target)), trace=trace, domains=domains,
               cartesian_before=prod(map(len, roots)), cartesian_after=prod(map(len, domains)),
               rounds=len(trace), deletions=sum(map(len, trace)))
    if any(not dom for dom in domains):
        return dict(out, status='CertifiedEmpty', reason='exhausted cycle domain')
    if not all(len(dom) == 1 for dom in domains):
        return dict(out, status='Unresolved', reason='multiple retained branches')
    values, valid = rational_candidate(p, cycles, weights, base, roots, domains)
    if values is None:
        return dict(out, status='Unresolved', reason='retained nonrational root')
    if not valid:
        return dict(out, status='CertifiedEmpty', reason='final reply/support test failed')
    return dict(out, status='CertifiedUnique', values=list(map(str, values)),
                reply_recomputed=True, actual_support_checked=True)
