#!/usr/bin/env python3
"""Independent, stdlib-only exact checker. Run: python3 check.py [certificate].

Only the raw bitmasks, ground size, and port labels feed the computation.
Certificate matrix/rate fields are compared afterwards, never used as evidence.
The default input and REPORT_openai_v12.md are resolved relative to this script.
"""
import argparse
from collections import Counter
from fractions import Fraction
from hashlib import sha256
from itertools import combinations
import json
from pathlib import Path


def validate(sets, n):
    if type(n) is not int or n < 2:
        raise ValueError('invalid ground size')
    if not sets or any(type(s) is not int or not 0 <= s < (1 << n) for s in sets):
        raise ValueError('empty family or invalid bitmask')
    if len(set(sets)) != len(sets):
        raise ValueError('duplicate sets')


def containment(sets):
    """Visit every unordered pair, testing both strict inclusion directions."""
    pairs = 0
    for i, a in enumerate(sets):
        for j in range(i + 1, len(sets)):
            b = sets[j]
            pairs += 1
            if a != b and (a & b == a or a & b == b):
                return (a, b), pairs
    return None, pairs


def sunflower(sets, n):
    """Exact exhaustive pair/complement search; no size heuristics or sampling.

    For any fixed A,B, all possible third sets are precisely
    (A & B) | S for S subset of the n-bit complement of (A | B).
    Searching even candidates below B avoids relying on an ordering shortcut.
    """
    masks = sorted(sets)
    lookup = set(masks)
    full = (1 << n) - 1
    pairs = candidates = 0
    for i, a in enumerate(masks):
        for j in range(i + 1, len(masks)):
            b = masks[j]
            core = a & b
            free = full ^ (a | b)
            s = free
            pairs += 1
            while True:
                c = core | s
                candidates += 1
                if c != a and c != b and c in lookup:
                    return (a, b, c), pairs, candidates
                if s == 0:
                    break
                s = (s - 1) & free
    return None, pairs, candidates


def matrix_from_sets(sets, u, v):
    counts = Counter(((s >> u) & 1, (s >> v) & 1) for s in sets)
    # Column state is (old-port present, old-port absent).
    # Row state is (new-port present, new-port absent).
    return counts, [[counts[0, 1], counts[1, 1]],
                    [counts[0, 0], counts[1, 0]]]


def polynomial_data(m):
    a, b = m[0]
    c, d = m[1]
    return a + d, a * d - b * c


def evaluate(x, trace, determinant):
    return x * x - trace * x + determinant


def compare_rate(rate, m, exponent):
    """Return exact relation of rate**exponent to the larger eigenvalue.

    Both off-diagonal entries must be positive. The larger root is strictly
    above trace/2. On that half-line the characteristic polynomial increases.
    """
    if rate <= 0 or exponent <= 0 or min(m[0][1], m[1][0]) <= 0:
        raise ValueError('positive rate/exponent and off-diagonal entries required')
    trace, det = polynomial_data(m)
    x = rate ** exponent
    p = evaluate(x, trace, det)
    if x <= Fraction(trace, 2):
        relation = -1
    else:
        relation = (p > 0) - (p < 0)
    return relation, x, p


def controls(sets, n, m, exponent):
    # Exercise the same algorithms on all 256 families on a three-point set,
    # against direct definitions (all pairs/all triples), including empty ones.
    for family_bits in range(1 << 8):
        f = [s for s in range(8) if family_bits >> s & 1]
        brute_c = any(a != b and (a & b == a or a & b == b)
                      for a, b in combinations(f, 2))
        brute_s = any(a & b == a & c == b & c
                      for a, b, c in combinations(f, 3))
        if (containment(f)[0] is not None) != brute_c:
            raise AssertionError('containment disagrees with direct definition')
        if (sunflower(f, 3)[0] is not None) != brute_s:
            raise AssertionError('sunflower disagrees with direct definition')
    planted_s = sorted(set(sets) | {1, 2, 4})
    sw = sunflower(planted_s, n)[0]
    if sw is None or not (len(set(sw)) == 3 and sw[0] & sw[1] == sw[0] & sw[2] == sw[1] & sw[2]):
        raise AssertionError('planted sunflower not rejected')
    planted_c = sorted(set(sets) | {0})
    cw = containment(planted_c)[0]
    if cw is None:
        raise AssertionError('planted containment not rejected')
    inflated = Fraction('1.564773096083')
    relation, x, p = compare_rate(inflated, m, exponent)
    if relation != 1:
        raise AssertionError('specified inflated claim is not strictly above the rate')
    # Duplicate and out-of-ground negative controls use the real validator.
    for f in (sets + [sets[0]], sets + [1 << n]):
        try:
            validate(f, n)
        except ValueError:
            pass
        else:
            raise AssertionError('invalid input not rejected')
    return sw, cw, inflated, x, p


def fmt(q):
    return f'{q.numerator}/{q.denominator}'


def main():
    here = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('certificate', nargs='?', type=Path,
                        default=here / 'cert_c49_chain20_v12.json')
    args = parser.parse_args()
    raw = args.certificate.read_bytes()
    cert = json.loads(raw)
    base = cert['bases']['B']
    n, sets = base['n'], base['sets']
    validate(sets, n)
    if n != 20:
        raise ValueError('this task requires n=20')
    word = cert['word']
    if len(word) != 1 or word[0]['base'] != 'B':
        raise ValueError('expected a one-base chain')
    u, v = word[0]['u'], word[0]['v']
    if type(u) is not int or type(v) is not int or not (0 <= u < n and 0 <= v < n and u != v):
        raise ValueError('invalid ports')
    exponent = n - 2
    dist = sorted(Counter(s.bit_count() for s in sets).items())
    counts, m = matrix_from_sets(sets, u, v)
    print(f'Raw family: {len(sets)} sets; sizes {dist}; ports {(u, v)}', flush=True)
    cw, cpairs = containment(sets)
    print(f'Containment: witness={cw}; pairs={cpairs}', flush=True)
    sw, spairs, candidates = sunflower(sets, n)
    print(f'Sunflower: witness={sw}; pairs={spairs}; candidates={candidates}', flush=True)
    trace, det = polynomial_data(m)
    lower = Fraction('1.564773096082')
    relation, x, p = compare_rate(lower, m, exponent)
    csw, ccw, upper, ux, up = controls(sets, n, m, exponent)
    checks = [len(sets) == 6341, cw is None, sw is None, relation == -1,
              m == cert.get('M'), exponent == cert.get('S'),
              cert.get('claim') == '1.564773096082', (u, v) == (18, 19)]
    verdict = 'ACCEPT' if all(checks) else 'REJECT'
    rows = [
        ('Raw input', f'SHA-256 `{sha256(raw).hexdigest()}`'),
        ('Distinctness / range', f'{len(sets)} distinct integer masks; all in [0, 2^{n}); PASS'),
        ('Size distribution', '; '.join(f'{size}: {count}' for size, count in dist)),
        ('Sperner', f'{cpairs:,} unordered pairs, both containment directions; witness {cw}'),
        ('3-sunflower-free', f'{spairs:,} unordered pairs; {candidates:,} candidate third sets; witness {sw}'),
        ('Port counts (u,v)', '; '.join(f'{key}: {counts[key]}' for key in [(0,0),(0,1),(1,0),(1,1)])),
        ('Derived matrix', f'`{m}`; agrees with certificate: {m == cert.get("M")}'),
        ('Characteristic polynomial', f'p(t) = t^2 - {trace}t + {det}; discriminant {trace*trace-4*det}'),
        ('Exact lower bound', f'L = {fmt(lower)}; L^{exponent} > {trace}/2: {x > Fraction(trace,2)}; p(L^{exponent}) < 0: {p < 0}'),
        ('Planted sunflower control', f'Add masks 1,2,4 to actual family; REJECT, witness {csw}'),
        ('Planted containment control', f'Add mask 0 to actual family; REJECT, witness {ccw}'),
        ('Inflated bound control', f'U = {fmt(upper)}; U^{exponent} > {trace}/2: {ux > Fraction(trace,2)}; p(U^{exponent}) > 0: {up > 0}; REJECT'),
        ('Algorithm controls', 'All 256 families on {0,1,2} agree with direct pair/triple definitions; duplicate and out-of-range inputs rejected'),
        ('Lemma hypotheses / exponent', f'Disjoint relabelled copies; Sperner and 3-free verified; n_j = {exponent}j+2; rate exponent 1/{exponent}'),
    ]
    report = [f'VERDICT: {verdict}', '', '| Check | Observed result |', '|---|---|']
    report += [f'| {name} | {result} |' for name, result in rows]
    report += ['', '## Transfer matrix derived from the branches', '',
        'Write c_ij for the number of base members having u-membership i and v-membership j. '
        'Use a column vector x=(x_1,x_0)^T counting current members with the old output port present and absent. '
        'Branch I takes x_1 members with base u absent, so its contributions to the new present/absent states are c_01 x_1 and c_00 x_1. '
        'Branch II takes x_0 members with base u present, contributing c_11 x_0 and c_10 x_0. '
        'Thus x_new = [[c_01,c_11],[c_00,c_10]] x. '
        'These four classes give exactly the displayed certificate matrix in present-then-absent order. '
        'No certificate count, matrix entry, exponent, or claimed rate was used to construct it.', '',
        '## Exact root isolation', '',
        f'The larger eigenvalue is rho = ({trace} + sqrt({trace*trace-4*det}))/2. '
        'Both off-diagonal entries are positive, so the discriminant is positive and rho > trace/2. '
        'The polynomial is strictly increasing to the right of trace/2. '
        'The following reduced Fraction values, computed with integer arithmetic, isolate the root:', '',
        f'- L^{exponent} = `{fmt(x)}`.',
        f'- p(L^{exponent}) = `{fmt(p)}`.',
        f'- U^{exponent} = `{fmt(ux)}`.',
        f'- p(U^{exponent}) = `{fmt(up)}`.', '',
        f'Hence L^{exponent} < rho < U^{exponent}, and L < rho^(1/{exponent}) < U. '
        'The lower claim passes strictly and the inflated claim fails strictly. '
        'No floating-point calculations enter the checker or this evidence.', '',
        '## Exhaustiveness and the infinite chain', '',
        'For any fixed distinct A,B, equality of the three pairwise intersections forces C to contain every bit of A&B '
        'and exclude every bit of their symmetric difference. All other bits lie outside A|B and are free. '
        'The submask loop enumerates each such C exactly once, including the empty submask, excluding only C=A or C=B. '
        'Every unordered pair is tested, so every possible sunflower is covered (usually three times). '
        'Containment is separately tested in both directions for every unordered pair.', '',
        'The original [DEGKM paper, Lemma 4.4, printed p.129 / PDF page 12]'
        '(https://kostochk.web.illinois.edu/docs/old/jcta97degm.pdf) was checked directly. '
        'It applies to Sperner, 3-free families on disjoint ground sets, with one marked point in each. '
        'Its two branches match those in TASK.md and delete both marked points. '
        'Each step uses a freshly relabelled copy of B, so the ground sets are disjoint before gluing. '
        'The lemma preserves the two properties; v survives because u and v are distinct.', '',
        'For completeness, deletion is injective within each branch because the deleted port has fixed membership. '
        'Across branches an inclusion of output sets would force a strict inclusion in one of the original families, '
        'so the branches are disjoint and their union is Sperner. If a hypothetical sunflower uses one branch only, '
        'each factor projection has equal pairwise intersections. In a Sperner 3-free family such a triple, '
        'allowing repetitions, must have all projections identical: three distinct projections violate 3-freeness, '
        'and exactly two identical ones force containment of the repeated member in the third. '
        'Thus all three outputs would be identical. If two outputs use branch I and the third branch II, '
        'restore the base input port to the third projection. It occurs once and leaves all pairwise intersections equal. '
        'Three distinct base projections violate 3-freeness; a repeated projection forces strict containment. '
        'The other mixed case follows by interchanging the factors. This independently checks the lemma mechanism.', '',
        f'Begin with B on {n} points. Its state vector is '
        f'x_1=({counts[0,1]+counts[1,1]}, {counts[0,0]+counts[1,0]})^T, with both entries positive. '
        f'After j copies the vector is M^(j-1)x_1 and the ground size is '
        f'{n}j-2(j-1)={exponent}j+2. '
        'The positive matrix has a strictly positive eigenvector for rho. '
        'For example (M_01, rho-M_00)^T is positive, since rho exceeds both diagonal entries. '
        'Bound x_1 above and below by positive multiples of that vector and apply M repeatedly. '
        'The resulting family sizes are bounded above and below by positive constants times rho^(j-1). '
        f'Taking ({exponent}j+2)-th roots gives the limiting rate rho^(1/{exponent}). '
        'Since these are actual sunflower-free families, they imply the claimed lower bound for mu_3. '
        'For a full-limit definition of mu_3, let G(n) be the largest size of a Sperner 3-free family. '
        'Disjoint products preserve both properties by the projection argument above, so G is supermultiplicative '
        'and Fekete\'s lemma gives its exponential-rate limit. If F(n) is the largest unrestricted 3-free size, '
        'then G(n) <= F(n) <= (n+1)G(n): some cardinality layer of any family contains at least its size divided by n+1. '
        'Thus F and G have the same exponential-rate limit.', '',
        '## Concerns and reproducibility', '',
        'The unusual size-7 members are included in every check; no uniformity assumption was made. '
        'No unresolved concern affects the claim. This checks the supplied construction and rate, not novelty or optimality. '
        'Run `python3 check.py` beside the certificate to regenerate this report; Python standard library only. '
        'A failed family/property or certificate comparison produces REJECT and exit status 1; malformed input exits with an error.', '']
    (here / 'REPORT_openai_v12.md').write_text('\n'.join(report))
    print(f'VERDICT: {verdict}; M={m}; report={here / "REPORT_openai_v12.md"}', flush=True)
    return 0 if verdict == 'ACCEPT' else 1


if __name__ == '__main__':
    raise SystemExit(main())
