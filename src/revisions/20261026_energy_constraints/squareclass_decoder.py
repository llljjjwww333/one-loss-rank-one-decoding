"""Public-input quadratic-energy pruning; unresolved is never a uniqueness claim.

The certificate concerns all actual real targets, without a factor budget.
No target, floating-point comparison or integer factorization is used.
"""
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from fractions import Fraction as Q
from math import isqrt, prod
import sympy as s
from cubic_probability_decoder import Quad, real_roots
from exact_cactus_decoder import canonical_public, cycle_basis, forest_base, actual_budget


def rational_sqrt(x):
    if x < 0:
        return None
    a, b = isqrt(x.numerator), isqrt(x.denominator)
    return Q(a, b) if a*a == x.numerator and b*b == x.denominator else None


def prune(groups):
    """Group *energy* classes, including rational-energy irrational roots."""
    reps, owners, labels = [], [], []
    for j, group in enumerate(groups):
        row = []
        for energy in group:
            if not energy.b:
                row.append(None)
                continue
            for k, representative in enumerate(reps):
                if rational_sqrt(energy.d/representative) is not None:
                    break
            else:
                k = len(reps)
                reps.append(energy.d)
                owners.append(set())
            owners[k].add(j)
            row.append(k)
        labels.append(row)
    keep = [[i for i, k in enumerate(row) if k is None or len(owners[k]) > 1]
            for row in labels]
    return keep, dict(representatives=list(map(str, reps)),
                      owners=[sorted(group) for group in owners], labels=labels)


def decode(public):
    p = canonical_public(public)
    m, n, edges = (p[key] for key in ('m', 'n', 'edges'))
    u, v, gu, gv = (p[key] for key in ('u', 'v', 'gu', 'gv'))
    cycles = cycle_basis(m, n, edges)
    weights = [u[i]*v[j] for i, j in edges]
    rhs = [-u[i]*gu[i] for i in range(m)]+[-v[j]*gv[j] for j in range(n)]
    for e, (i, j) in enumerate(edges):
        rhs[i] += weights[e]**2
        rhs[m+j] += weights[e]**2
    base = forest_base(m, n, edges, cycles, rhs)
    out = dict(candidate_class='all actual real targets; no factor budget',
               truth_passed_to_solver=False)
    if base is None:
        return dict(out, status='CertifiedEmpty', reason='inconsistent margins')
    t = s.Symbol('t')
    polynomials = []
    for cycle in cycles:
        z = [base[e]+(-1)**k*t for k, e in enumerate(cycle)]
        polynomials.append(s.Poly(s.prod(z[::2])-s.prod(z[1::2]), t, domain=s.QQ))
    if any(poly.degree() == 0 and not poly.is_zero for poly in polynomials):
        return dict(out, status='CertifiedEmpty', reason='nonzero constant cycle')
    if any(poly.is_zero for poly in polynomials):
        return dict(out, status='Unresolved', reason='zero cycle polynomial')
    if any(poly.degree() > 3 for poly in polynomials):
        return dict(out, status='Unresolved', reason='cycle degree exceeds three')
    roots, groups = [], []
    for cycle, poly in zip(cycles, polynomials):
        try:
            local_roots = real_roots(poly)
        except ValueError:
            return dict(out, status='Unresolved', reason='irreducible cubic factor')
        roots.append(local_roots)
        energies = []
        for root in local_roots:
            residual = Quad()
            for coefficient in poly.all_coeffs():
                residual = residual*root+Q(coefficient)
            assert residual.zero()
            values = [Quad(Q(base[e]/weights[e]))+root*Q((-1)**k/weights[e])
                      for k, e in enumerate(cycle)]
            energies.append(sum((x*x for x in values), Quad()))
        groups.append(energies)
    keep, basis = prune(groups)
    out.update(square_classes=basis, before=[len(g) for g in groups],
               after=list(map(len, keep)), cartesian_before=prod(map(len, groups)),
               cartesian_after=prod(map(len, keep)),
               cycles=[dict(edges=c, coefficients=list(map(str, poly.all_coeffs())),
                            roots=[x.encode() for x in local], energies=[x.encode() for x in g],
                            retained=k)
                       for c, poly, local, g, k in zip(cycles, polynomials, roots, groups, keep)])
    if any(not group for group in keep):
        return dict(out, status='CertifiedEmpty', reason='all branches excluded on a cycle')
    if any(len(group) != 1 for group in keep):
        return dict(out, status='Unresolved', reason='multiple retained branches')
    chosen = [group[k[0]] for group, k in zip(roots, keep)]
    if any(root.b for root in chosen):
        return dict(out, status='Unresolved', reason='retained nonrational root')
    values = [Q(z/w) for z, w in zip(base, weights)]
    for cycle, root in zip(cycles, chosen):
        for k, e in enumerate(cycle):
            values[e] += root.a*Q((-1)**k/weights[e])
    norm_target = Q(2*p['loss']+sum(w*w for w in weights)
                    -sum(x*y for x, y in zip(u, gu))-sum(x*y for x, y in zip(v, gv)))
    norm_ok = sum((x*x for x in values), Q(0)) == norm_target
    actual, _ = actual_budget(m, n, edges, list(map(s.Rational, values)), False)
    if not norm_ok or not actual:
        return dict(out, status='CertifiedEmpty', reason='rational final norm/support test failed')
    # Check the supplied reply independently of the extracted norm identity.
    check_gu, check_gv = [s.S.Zero]*m, [s.S.Zero]*n
    check_loss = s.S.Zero
    for (i, j), weight, value in zip(edges, weights, values):
        residual = weight-s.Rational(value)
        check_gu[i] += residual*v[j]
        check_gv[j] += residual*u[i]
        check_loss += residual**2/2
    assert check_gu == gu and check_gv == gv and check_loss == p['loss']
    return dict(out, status='CertifiedUnique', values=list(map(str, values)),
                reply_recomputed=True, actual_support_checked=True)
