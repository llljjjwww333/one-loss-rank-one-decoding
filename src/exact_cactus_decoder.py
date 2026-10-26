"""Exact QQ public-reply decoder, including nonsplit real cycle polynomials.

No target, numerical tolerance, or caller-supplied cycle basis is used.
The JSON certificate is replayable, not a formally verified proof object.
SymPy supplies rational root isolation and primitive-element zero tests.
"""
from functools import cmp_to_key, lru_cache
from itertools import product
from math import isqrt, prod
import argparse
import hashlib
import json
from pathlib import Path
import sympy as s


def rational(value):
    value = s.sympify(value)
    if value.is_Rational is not True:
        raise ValueError('Public coordinates must be exact rationals (use strings for fractions)')
    return value


def cycle_basis(m, n, edges):
    """Fundamental cycles; reject overlapping cycles (noncactus input)."""
    adjacency = [[] for _ in range(m+n)]
    cycles, used = [], set()
    for index, (row, col) in enumerate(edges):
        a, b = row, m+col
        parent = {a: (None, None)}
        stack = [a]
        while stack and b not in parent:
            node = stack.pop()
            for neighbor, edge in adjacency[node]:
                if neighbor not in parent:
                    parent[neighbor] = (node, edge)
                    stack.append(neighbor)
        if b in parent:
            path, node = [], b
            while node != a:
                node, edge = parent[node]
                path.append(edge)
            cycle = list(reversed(path))+[index]
            if used.intersection(cycle):
                raise ValueError('Mask is not a cactus')
            used.update(cycle)
            cycles.append(cycle)
        else:
            adjacency[a].append((b, index))
            adjacency[b].append((a, index))
    return cycles


def forest_base(m, n, edges, cycles, rhs):
    """Set chord entries to zero and solve the forest margins by leaf removal."""
    chords = {cycle[-1] for cycle in cycles}
    adjacency = [set() for _ in range(m+n)]
    for e, (i, j) in enumerate(edges):
        if e not in chords:
            adjacency[i].add(e)
            adjacency[m+j].add(e)
    remaining = list(rhs)
    base = [s.S.Zero]*len(edges)
    leaves = [node for node, incident in enumerate(adjacency) if len(incident) == 1]
    while leaves:
        node = leaves.pop()
        if not adjacency[node]:
            continue
        e = adjacency[node].pop()
        i, j = edges[e]
        neighbor = m+j if node == i else i
        base[e] = remaining[node]
        remaining[node] = s.S.Zero
        remaining[neighbor] -= base[e]
        adjacency[neighbor].remove(e)
        if len(adjacency[neighbor]) == 1:
            leaves.append(neighbor)
    if any(remaining):
        return None
    return base


def _mul(a, b):
    values = [x*y for x in a for y in b]
    return min(values), max(values)


def _sqrt_bound(q, bits, upper):
    q = s.Rational(q)
    scale = 1 << bits
    k = isqrt(int(q.p)*scale*scale//int(q.q))
    if upper and s.Rational(k*k, scale*scale) != q:
        k += 1
    return s.Rational(k, scale)


@lru_cache(maxsize=32768)
def enclosure(expr, bits):
    """Outward rational interval arithmetic; no floating-point decisions."""
    if expr.is_Rational:
        return expr, expr
    if isinstance(expr, s.CRootOf):
        interval = expr._get_interval().refine_size(s.Rational(1, 1 << bits))
        return s.Rational(interval.a), s.Rational(interval.b)
    if expr.is_Add:
        pairs = [enclosure(arg, bits) for arg in expr.args]
        return sum(a for a, _ in pairs), sum(b for _, b in pairs)
    if expr.is_Mul:
        result = (s.S.One, s.S.One)
        for arg in expr.args:
            result = _mul(result, enclosure(arg, bits))
        return result
    if expr.func == s.Abs:
        lo, hi = enclosure(expr.args[0], bits)
        return (s.S.Zero if lo <= 0 <= hi else min(abs(lo), abs(hi)), max(abs(lo), abs(hi)))
    if expr.is_Pow:
        base, exponent = expr.args
        lo, hi = enclosure(base, bits)
        if exponent.is_Rational and exponent.q == 2:
            # Budget radicands are provably nonnegative sums of squares.
            if hi < 0:
                raise ValueError('Negative square-root radicand')
            lo, hi = _sqrt_bound(max(0, lo), bits, False), _sqrt_bound(max(0, hi), bits, True)
            exponent = s.Integer(exponent.p)
        if exponent.is_Integer:
            k = int(exponent)
            if k < 0:
                if lo <= 0 <= hi:
                    raise ZeroDivisionError('Refine denominator interval')
                lo, hi, k = 1/hi, 1/lo, -k
            if k == 0:
                return s.S.One, s.S.One
            if k % 2:
                return lo**k, hi**k
            return (s.S.Zero if lo <= 0 <= hi else min(lo**k, hi**k), max(lo**k, hi**k))
    raise NotImplementedError('Interval expression: '+str(expr))


def exact_sign(expr):
    expr = s.sympify(expr)
    if expr.is_Rational:
        return int(s.sign(expr))
    bits, zero_tested = 16, False
    while True:
        try:
            lo, hi = enclosure(expr, bits)
            if lo > 0:
                return 1
            if hi < 0:
                return -1
            if lo == hi == 0:
                return 0
        except ZeroDivisionError:
            pass
        if not zero_tested:
            # AlgebraicField.is_positive is NOT an ordering for this embedding.
            if all(c == 0 for c in s.to_number_field(expr).coeffs()):
                return 0
            zero_tested = True
        bits *= 2


def actual_budget(m, n, edges, values, bounded):
    # Reduce propagation in one exact field. Repeated primitive-element
    # tests on unreduced rational expressions can otherwise be very costly.
    generators = sorted(set().union(*(c.atoms(s.CRootOf) for c in values)), key=str)
    field = s.QQ.algebraic_field(*generators) if generators else s.QQ
    converted = [field.from_sympy(c) for c in values]
    active = [c != field.zero for c in converted]
    adjacency = [[] for _ in range(m+n)]
    for edge, ((i, j), on) in enumerate(zip(edges, active)):
        if on:
            adjacency[i].append((m+j, edge))
            adjacency[m+j].append((i, edge))
    if any(not on and adjacency[i] and adjacency[m+j]
           for (i, j), on in zip(edges, active)):
        return False, None
    factors, terms = {}, []
    for first in range(m+n):
        if first in factors or not adjacency[first]:
            continue
        factors[first], stack, component = field.one, [first], []
        while stack:
            node = stack.pop()
            component.append(node)
            for neighbor, edge in adjacency[node]:
                if neighbor not in factors:
                    factors[neighbor] = converted[edge]/factors[node]
                    stack.append(neighbor)
                elif factors[node]*factors[neighbor] != converted[edge]:
                    return False, None
        if bounded:
            left = sum((factors[i]**2 for i in component if i < m), field.zero)
            right = sum((factors[i]**2 for i in component if i >= m), field.zero)
            terms.append(s.sqrt(field.to_sympy(left*right)))
    budget = sum(terms, s.S.Zero) if bounded else None
    return (not bounded or exact_sign(budget-1) <= 0), budget


def canonical_public(public):
    m, n = public['m'], public['n']
    if type(m) is not int or type(n) is not int or min(m, n) < 0:
        raise ValueError('Invalid dimensions')
    edges = [tuple(edge) for edge in public['edges']]
    if len(edges) != len(set(edges)) or any(
        len(edge) != 2 or any(type(k) is not int for k in edge) or
        not (0 <= edge[0] < m and 0 <= edge[1] < n) for edge in edges):
        raise ValueError('Invalid mask')
    result = dict(m=m, n=n, edges=edges)
    for key, size in [('u', m), ('v', n), ('gu', m), ('gv', n)]:
        if len(public[key]) != size:
            raise ValueError('Reply/query dimension mismatch')
        result[key] = list(map(rational, public[key]))
    if any(q <= 0 for q in result['u']+result['v']):
        raise ValueError('Query must be positive')
    result['loss'] = rational(public['loss'])
    return result


def json_public(public):
    p = canonical_public(public)
    return {key: [str(c) for c in value] if key in {'u','v','gu','gv'}
            else str(value) if key == 'loss' else value for key, value in p.items()}


def digest_public(public):
    return hashlib.sha256(json.dumps(json_public(public), sort_keys=True).encode()).hexdigest()


def energy_choices(branches, target):
    """Exact meet-in-the-middle energy join, including equal-energy blocks."""
    groups, sizes, comparisons = [[], []], [1, 1], 0
    for cycle in sorted(range(len(branches)), key=lambda j: -len(branches[j])):
        side = 0 if sizes[0] <= sizes[1] else 1
        groups[side].append(cycle)
        sizes[side] *= len(branches[cycle])

    def compare(a, b):
        nonlocal comparisons
        comparisons += 1
        return exact_sign(a[0]-b[0])

    lists = []
    for side, group in enumerate(groups):
        partial = []
        for choice in product(*(range(len(branches[j])) for j in group)):
            total = sum((branches[j][k] for j, k in zip(group, choice)), s.S.Zero)
            partial.append((total if side == 0 else target-total, choice))
        partial.sort(key=cmp_to_key(compare))
        lists.append(partial)
    left, right = lists
    i, j, matches = 0, 0, []
    while i < len(left) and j < len(right):
        order = compare(left[i], right[j])
        if order < 0:
            i += 1
        elif order > 0:
            j += 1
        else:
            ii, jj = i+1, j+1
            while ii < len(left) and compare(left[i], left[ii]) == 0:
                ii += 1
            while jj < len(right) and compare(right[j], right[jj]) == 0:
                jj += 1
            for a in left[i:ii]:
                for b in right[j:jj]:
                    choice = [0]*len(branches)
                    for group, item in zip(groups, (a,b)):
                        for cycle, k in zip(group, item[1]):
                            choice[cycle] = k
                    matches.append(tuple(choice))
            i, j = ii, jj
    return sorted(matches), dict(groups=groups, partial_sizes=sizes,
                                  exact_comparisons=comparisons, energy_matches=len(matches))


def decode(public, bounded=False, search='enumerate'):
    if search not in {'enumerate','mitm'}:
        raise ValueError('Unknown energy search')
    p = canonical_public(public)
    m, n, edges = p['m'], p['n'], p['edges']
    u, v, gu, gv = [p[k] for k in ('u', 'v', 'gu', 'gv')]
    cycles = cycle_basis(m, n, edges)
    certificate = dict(format_version=1, public_sha256=digest_public(p), bounded=bounded, search=search,
                       arithmetic='QQ real algebraic', cycles=[], decisions=[], candidates=[])

    def finish(status, reason=None):
        certificate.update(status=status, complete=status != 'Degenerate')
        if reason:
            certificate['reason'] = reason
        return certificate

    weights = []
    pu, pv = [s.S.Zero]*m, [s.S.Zero]*n
    for edge, (i, j) in enumerate(edges):
        pu[i] += u[i]*v[j]**2
        pv[j] += v[j]*u[i]**2
        weights.append(u[i]*v[j])
    rhs = [u[i]*(pu[i]-gu[i]) for i in range(m)]+[v[j]*(pv[j]-gv[j]) for j in range(n)]
    base = forest_base(m, n, edges, cycles, rhs)
    if base is None:
        return finish('Empty', 'inconsistent_margins')
    norm = 2*p['loss']+sum(w*w for w in weights)-sum(u[i]*gu[i] for i in range(m))-sum(v[j]*gv[j] for j in range(n))
    certificate.update(norm=str(norm), base=list(map(str, base)))
    t, root_sets, polynomials = s.Symbol('t'), [], []
    for cycle in cycles:
        z = [base[edge]+(-1)**offset*t for offset, edge in enumerate(cycle)]
        poly = s.Poly(s.prod(z[::2])-s.prod(z[1::2]), t, domain=s.QQ)
        polynomials.append(poly)
        record = dict(edges=cycle, coefficients=list(map(str, poly.all_coeffs())))
        certificate['cycles'].append(record)
    if any(not poly.is_zero and poly.degree() == 0 for poly in polynomials):
        return finish('Empty', 'nonzero_constant_cycle')
    if any(poly.is_zero for poly in polynomials):
        return finish('Degenerate', 'identically_zero_cycle')
    for poly, record in zip(polynomials, certificate['cycles']):
        roots = poly.real_roots(multiple=False, radicals=False)
        intervals = poly.intervals(eps=s.Rational(1, 65536))
        assert len(roots) == len(intervals)
        record['degree'] = int(poly.degree())
        record['roots'] = [dict(expression=str(root), multiplicity=int(mult),
                                interval=list(map(str, interval)))
                           for (root, mult), (interval, count) in zip(roots, intervals)]
        assert all(mult == count for (_, mult), (_, count) in zip(roots, intervals))
        root_sets.append([root for root, _ in roots])
    certificate['tuple_bound'] = prod(map(len, root_sets))
    cycle_edges = {edge for cycle in cycles for edge in cycle}
    fixed_energy = sum((base[e]/weights[e])**2 for e in range(len(edges)) if e not in cycle_edges)
    branch_energies = []
    for cycle, roots in zip(cycles, root_sets):
        energy = s.Poly(sum(((base[e]+(-1)**j*t)/weights[e])**2 for j, e in enumerate(cycle)), t)
        branch_energies.append([energy.eval(root) if root.is_Rational else
                                energy.rem(s.Poly(root.poly.as_expr().subs(root.poly.gen, t), t)).eval(root)
                                for root in roots])
    if search == 'mitm':
        choices, stats = energy_choices(branch_energies, norm-fixed_energy)
        certificate['energy_search'] = stats
    else:
        choices = product(*(range(len(roots)) for roots in root_sets))
    for choice in choices:
        energy = fixed_energy+sum(branch_energies[j][k] for j, k in enumerate(choice))
        decision = dict(roots=list(choice), energy_sign=exact_sign(energy-norm))
        certificate['decisions'].append(decision)
        if decision['energy_sign']:
            continue
        values = [base[e]/weights[e] for e in range(len(edges))]
        for cycle, roots, k in zip(cycles, root_sets, choice):
            for j, e in enumerate(cycle):
                values[e] += (-1)**j*roots[k]/weights[e]
        admissible, budget = actual_budget(m, n, edges, values, bounded)
        decision['admissible'] = admissible
        if budget is not None:
            decision['factor_budget'] = str(budget)
        if admissible:
            certificate['candidates'].append(dict(roots=list(choice), values=list(map(str, values))))
    return finish('Empty' if not certificate['candidates'] else
                  'Unique' if len(certificate['candidates']) == 1 else 'Ambiguous')


def verify_certificate(public, certificate):
    """Recompute every root and decision; no trust in the supplied candidate list.

    Reuses the decoder and SymPy backend: replay, not independent proof checking.
    """
    return certificate == decode(public, bounded=certificate['bounded'],search=certificate['search'])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('public', type=Path)
    parser.add_argument('certificate', type=Path)
    parser.add_argument('--bounded', action='store_true')
    parser.add_argument('--search', choices=['enumerate','mitm'],default='enumerate')
    parser.add_argument('--verify', action='store_true')
    args = parser.parse_args()
    public = json.loads(args.public.read_text(encoding='utf-8'))
    if args.verify:
        if not verify_certificate(public, json.loads(args.certificate.read_text(encoding='utf-8'))):
            raise SystemExit('Certificate replay FAILED')
        print('Certificate replay PASS')
    else:
        result = decode(public, bounded=args.bounded,search=args.search)
        args.certificate.write_text(json.dumps(result, indent=2)+'\n', encoding='utf-8')
        print(result['status'], len(result['candidates']))


if __name__ == '__main__':
    main()
