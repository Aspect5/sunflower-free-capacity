#!/usr/bin/env python3
"""INDEPENDENT exact checker for C_49 (3-sunflower-free capacity) certificates.

Entry 49a of the optimization-constants repository.

This file deliberately shares NO code with the searchers (search_*.py,
the search code): it imports nothing from them and copies no helper.  It uses the
standard library, mpmath (50-digit comparison) and numpy (vectorised exact
integer bit operations for the triple scan).  It re-reads the certificate from
disk and re-derives everything from the definition.

Definition enforced
-------------------
Three DISTINCT sets A,B,C form a 3-sunflower iff  A cap B = A cap C = B cap C.
F is sunflower-free (3-free) iff no three distinct members form a 3-sunflower.
Sets are integer bitmasks over the ground set [n] = {0,...,n-1}.
Equivalent form used by the scan: for distinct A,B put U = A|B, I = A&B; a
third member C (C != A, B) completes a sunflower iff (C & U) == I.

Format 1 (plain):  {"n": n, "sets": [masks]}
------------------------------------------
  * F 3-free and Sperner (no member contains another; every uniform family is
    Sperner)  ->  certified value = |F|^(1/n).
    Reason (DEGKM 1997 Lemma 3.1): the direct sum F^k on k disjoint copies of
    [n] is 3-free.  For three distinct members, look at block i: restrictions
    X,Y,Z have pairwise equal intersections.  If X = Y then X & Z = X, so
    X subset Z and Sperner forces X = Z.  So in every block the restrictions
    are all equal or pairwise distinct; the latter is a sunflower in F.  Hence
    all blocks agree and the members coincide.  f(kn) >= |F|^k.
  * F 3-free, NOT Sperner -> only its layers are used: certified value =
    (largest layer)^(1/n) (each layer is uniform, hence Sperner, and 3-free).
    Counterexample to using |F|: {{}, {1}} on n=1 squares to 2^[2], which
    contains the sunflower {}, {1}, {2}.

Format 2 (chain):  {"format": "c49-chain-v1",
                    "bases": {name: {"n": m, "sets": [masks]}},
                    "word":  [{"base": name, "u": int, "v": int}, ...]}
------------------------------------------------------------------------
  Every referenced base must be 3-free AND Sperner (checked explicitly).
  Gluing lemma (DEGKM 1997 Lemma 4.4; proof re-derived in the README): for
  Sperner 3-free A on X (a in X) and B on Y (b in Y), X, Y disjoint,
      G = {(C-a) u D : a in C in A, b notin D in B}
        u {C u (D-b) : a notin C in A, b in D in B}
  is Sperner and 3-free on (X u Y) - {a,b}; |G| = |A_a||B^b| + |A^a||B_b|.
  Chain: running family P with marked point z; a step (B,u,v) glues z to
  port u of a fresh copy of B and marks port v of that copy.  With
  p = (#members containing z, #members avoiding z) the new vector is
      p_in'  = p_in * #{D: u notin D, v in D} + p_out * #{D: u in D, v in D}
      p_out' = p_in * #{D: u notin D, v notin D} + p_out * #{D: u in D, v notin D}
  and the ground set grows by m_B - 2.  One period of the word multiplies p by
  P = M_k ... M_1.  The chain starts from the last word base with z = its v.
  The checker requires every entry of P and of the start vector to be > 0,
  so P^t p_0 >= c * rho^t * r (r the Perron vector, c > 0) and
      C_49 >= rho(P)^(1/S),  S = sum_i (m_i - 2).
  The acceptance test is EXACT: the checker picks a rational positive vector
  r ~ Perron vector and computes lam = min_i (P r)_i / r_i as a Fraction
  (Collatz-Wielandt: rho >= lam), then accepts iff Fraction(claim)^S <= lam.
  mpmath is used only to print rho^(1/S).

Usage
-----
  check_c49.py <certfile> --claim <float>   (exit 0 iff accepted)
  check_c49.py --selftest                    (planted bads + controls)

Baselines (absolute literals, report only): 1.551 DEGKM1997 published
(their proof gives 1.55118); 1.554 NS2017 arXiv, unpublished (Naslund);
upper bound 3/2^(2/3) = 1.8898815748.
"""

import sys
import os
import json
import tempfile
import itertools
from fractions import Fraction

import numpy as np
from mpmath import mp, mpf, power, sqrt as mpsqrt

mp.dps = 50


class Reject(Exception):
    pass


# ---------------------------------------------------------------- parsing
def _parse_family(obj, where):
    if not isinstance(obj, dict):
        raise Reject("%s is not an object" % where)
    if "n" not in obj:
        raise Reject("%s has no field 'n'" % where)
    n = obj["n"]
    if isinstance(n, bool) or not isinstance(n, int):
        raise Reject("%s: field 'n' is not an integer" % where)
    if n < 1 or n > 62:
        raise Reject("%s: field 'n' must be in [1,62], got %d" % (where, n))
    fam = None
    for key in ("sets", "subsets", "family", "masks", "members"):
        if key in obj:
            fam = obj[key]
            break
    if fam is None:
        raise Reject("%s has no family field (sets/subsets/family/masks/members)" % where)
    if not isinstance(fam, list) or len(fam) == 0:
        raise Reject("%s: family is not a non-empty list" % where)
    limit = 1 << n
    seen = {}
    for pos, x in enumerate(fam):
        if isinstance(x, bool) or not isinstance(x, int):
            raise Reject("%s: member %d is not an integer bitmask (%r)" % (where, pos, x))
        if x < 0:
            raise Reject("%s: member %d has negative bitmask %d" % (where, pos, x))
        if x >= limit:
            raise Reject("%s: member %d (bitmask %d) is not a subset of [%d]: "
                         "it uses a point outside the ground set" % (where, pos, x, n))
        if x in seen:
            raise Reject("%s: members %d and %d are duplicates (bitmask %d); "
                         "the family must consist of distinct sets" % (where, seen[x], pos, x))
        seen[x] = pos
    return n, list(fam)


def _popcount(x):
    return bin(x).count("1")


# ------------------------------------------------------------ exact scans
def _find_sunflower(masks):
    """Return indices (i,j,k) of a 3-sunflower or None.  Exact, near-quadratic.
    For each member A = masks[i]: every other member X has K = X & A and
    residual R = X & ~A.  {A, X, Y} is a sunflower iff X & A == Y & A (= K)
    and R_X & R_Y == 0 (then all three pairwise intersections equal K).  So we
    group the other members by K and look for a disjoint pair of residuals
    inside a group.  Valid for any family of distinct sets (no Sperner
    assumption).  Cross-checked against the cubic scan in the selftest."""
    a = np.array(masks, dtype=np.int64)
    m = len(a)
    for i in range(m):
        A = a[i]
        others = np.concatenate([np.arange(i), np.arange(i + 1, m)])
        K = a[others] & A
        R = a[others] & ~A
        order = np.argsort(K, kind="stable")
        Ks = K[order]
        starts = np.flatnonzero(np.r_[True, Ks[1:] != Ks[:-1]])
        ends = np.r_[starts[1:], len(Ks)]
        big = np.flatnonzero(ends - starts >= 2)
        for gi in big:
            idx = order[starts[gi]:ends[gi]]
            Rg = R[idx]
            for c0 in range(0, len(Rg), 1024):
                blk = Rg[c0:c0 + 1024]
                dis = (blk[:, None] & Rg[None, :]) == 0
                dis[np.arange(len(blk)), c0 + np.arange(len(blk))] = False
                if dis.any():
                    r, c = np.argwhere(dis)[0]
                    return (i, int(others[idx[c0 + r]]), int(others[idx[c]]))
    return None


def _find_sunflower_cubic(masks):
    """Return indices (i,j,k) of a 3-sunflower or None.  For every unordered
    pair i<j: U = A_i|A_j, I = A_i&A_j; any k != i,j with (A_k & U) == I is a
    sunflower (then A_k & A_i = I = A_k & A_j).  Exact int64 bit operations,
    vectorised over k and over chunks of j."""
    a = np.array(masks, dtype=np.int64)
    m = len(a)
    idx = np.arange(m)
    for i in range(m - 1):
        js = idx[i + 1:]
        for c0 in range(0, len(js), 256):
            jj = js[c0:c0 + 256]
            U = a[i] | a[jj]
            I = a[i] & a[jj]
            hit = (a[None, :] & U[:, None]) == I[:, None]
            hit[:, i] = False
            hit[np.arange(len(jj)), jj] = False
            if hit.any():
                r, k = np.argwhere(hit)[0]
                return (i, int(jj[r]), int(k))
    return None


def _find_containment(masks):
    """Return (i,j) with A_i a proper subset of A_j, or None."""
    a = np.array(masks, dtype=np.int64)
    m = len(a)
    for c0 in range(0, m, 256):
        blk = a[c0:c0 + 256]
        sub = (blk[:, None] & ~a[None, :]) == 0      # A_i subset A_j
        sub[np.arange(len(blk)), c0 + np.arange(len(blk))] = False
        if sub.any():
            r, j = np.argwhere(sub)[0]
            return (c0 + int(r), int(j))
    return None


def _verify_family(n, masks, where):
    tri = _find_sunflower(masks)
    if tri is not None:
        i, j, k = tri
        raise Reject("%s contains a 3-sunflower at indices (%d,%d,%d): masks %d,%d,%d "
                     "have equal pairwise intersections %d"
                     % (where, i, j, k, masks[i], masks[j], masks[k], masks[i] & masks[j]))
    return _find_containment(masks) is None


# ------------------------------------------------------------- format 1
def _verify_plain(obj, claim):
    n, masks = _parse_family(obj, "certificate")
    sperner = _verify_family(n, masks, "family")
    layers = {}
    for x in masks:
        w = _popcount(x)
        layers[w] = layers.get(w, 0) + 1
    if sperner:
        size, mode = len(masks), ("uniform" if len(layers) == 1 else "sperner(non-uniform)")
    else:
        r = max(sorted(layers), key=lambda t: layers[t])
        size, mode = layers[r], "non-sperner:largest-layer(r=%d)" % r
    value = power(mpf(size), mpf(1) / mpf(n))
    c = mpf(claim)
    if not (value >= c):
        raise Reject("certified value %s (%s, n=%d, counted=%d) < claim %s"
                     % (mp.nstr(value, 20), mode, n, size, mp.nstr(c, 20)))
    return ("OK plain n=%d |F|=%d mode=%s counted=%d value=%s >= claim=%s | beats_1.551=%s "
            "beats_1.554=%s" % (n, len(masks), mode, size, mp.nstr(value, 20), mp.nstr(c, 20),
                                "YES" if value > mpf("1.551") else "NO",
                                "YES" if value > mpf("1.554") else "NO"))


# ------------------------------------------------------------- format 2
def _step_matrix(masks, u, v):
    """[[#(u out, v in), #(u in, v in)], [#(u out, v out), #(u in, v out)]] :
    new (in, out) = M @ old (in, out) where old 'in' pairs with D avoiding u."""
    c = {(0, 0): 0, (0, 1): 0, (1, 0): 0, (1, 1): 0}
    for D in masks:
        c[((D >> u) & 1, (D >> v) & 1)] += 1
    return [[c[(0, 1)], c[(1, 1)]], [c[(0, 0)], c[(1, 0)]]]


def _matmul(A, B):
    return [[A[0][0] * B[0][0] + A[0][1] * B[1][0], A[0][0] * B[0][1] + A[0][1] * B[1][1]],
            [A[1][0] * B[0][0] + A[1][1] * B[1][0], A[1][0] * B[0][1] + A[1][1] * B[1][1]]]


def _verify_chain(obj, claim):
    bases = obj.get("bases")
    word = obj.get("word")
    if not isinstance(bases, dict) or not bases:
        raise Reject("chain certificate has no non-empty 'bases' object")
    if not isinstance(word, list) or not word:
        raise Reject("chain certificate has no non-empty 'word' list")
    parsed = {}
    for name in sorted(set(step.get("base") if isinstance(step, dict) else None for step in word)):
        if name not in bases:
            raise Reject("word references unknown base %r" % (name,))
        n, masks = _parse_family(bases[name], "base %r" % name)
        if not _verify_family(n, masks, "base %r" % name):
            raise Reject("base %r is not Sperner (a member contains another); the gluing "
                         "lemma needs Sperner bases" % name)
        parsed[name] = (n, masks)
    P = [[1, 0], [0, 1]]
    S = 0
    steps = []
    for pos, step in enumerate(word):
        u, v = step.get("u"), step.get("v")
        n, masks = parsed[step["base"]]
        for nm, x in (("u", u), ("v", v)):
            if isinstance(x, bool) or not isinstance(x, int) or not (0 <= x < n):
                raise Reject("word step %d: port %s=%r is not a point of base %r (n=%d)"
                             % (pos, nm, x, step["base"], n))
        if u == v:
            raise Reject("word step %d: ports u and v coincide (%d)" % (pos, u))
        M = _step_matrix(masks, u, v)
        P = _matmul(M, P)
        S += n - 2
        steps.append("%s(u=%d,v=%d,m=%d,M=%s)" % (step["base"], u, v, n, M))
    if min(P[0] + P[1]) <= 0:
        raise Reject("period matrix %s has a non-positive entry; Perron growth not certified" % (P,))
    n0, m0 = parsed[word[-1]["base"]]
    z = word[-1]["v"]
    p0 = (sum(1 for D in m0 if (D >> z) & 1), sum(1 for D in m0 if not (D >> z) & 1))
    if min(p0) <= 0:
        raise Reject("start vector %s has a zero entry" % (p0,))
    a, b, c, d = (mpf(P[0][0]), mpf(P[0][1]), mpf(P[1][0]), mpf(P[1][1]))
    rho = ((a + d) + mpsqrt((a - d) ** 2 + 4 * b * c)) / 2
    # rational Collatz-Wielandt lower bound
    r1 = Fraction(int(mp.floor(b * 10 ** 30)), 1)
    r2 = Fraction(int(mp.floor((rho - a) * 10 ** 30)), 1)   # eigvec (b, rho-a)
    if r1 <= 0 or r2 <= 0:
        raise Reject("could not form a positive test vector")
    lam = min((P[0][0] * r1 + P[0][1] * r2) / r1, (P[1][0] * r1 + P[1][1] * r2) / r2)
    cf = Fraction(claim)
    ok = cf ** S <= lam
    value = power(rho, mpf(1) / S)
    lam_root = power(mpf(lam.numerator) / mpf(lam.denominator), mpf(1) / S)
    if not ok:
        raise Reject("chain value rho^(1/S) = %s (rational CW bound %s, S=%d) < claim %s"
                     % (mp.nstr(value, 20), mp.nstr(lam_root, 20), S, repr(claim)))
    return ("OK chain S=%d steps=[%s] rho=%s value=%s cw_bound=%s >= claim=%s | beats_1.551=%s "
            "beats_1.554=%s" % (S, "; ".join(steps), mp.nstr(rho, 15), mp.nstr(value, 20),
                                mp.nstr(lam_root, 20), repr(claim),
                                "YES" if value > mpf("1.551") else "NO",
                                "YES" if value > mpf("1.554") else "NO"))


def verify(path, claim):
    if not os.path.isfile(path):
        raise Reject("certificate file not found: %s" % path)
    try:
        with open(path) as fh:
            obj = json.load(fh)
    except Exception as exc:
        raise Reject("certificate is not valid JSON: %s" % exc)
    if not isinstance(obj, dict):
        raise Reject("certificate JSON is not an object")
    try:
        float(claim)
    except Exception:
        raise Reject("claim %r is not a number" % (claim,))
    if obj.get("format") == "c49-chain-v1":
        return _verify_chain(obj, claim)
    if "format" in obj and obj["format"] not in ("plain", None):
        raise Reject("unknown format %r" % (obj["format"],))
    return _verify_plain(obj, claim)


# ---------------------------------------------------------------- selftest
def _write_tmp(obj):
    fd, path = tempfile.mkstemp(suffix=".json", prefix="c49selftest_")
    with os.fdopen(fd, "w") as fh:
        json.dump(obj, fh)
    return path


def _good_family():
    """8 transversals of the pairs {0,1},{2,3},{4,5}: 3-uniform, 3-free, n=6."""
    return {"n": 6, "sets": [(1 << a) | (1 << b) | (1 << c)
                             for a in (0, 1) for b in (2, 3) for c in (4, 5)]}


def _degkm_14():
    """DEGKM 1997 Sec. 4.2/4.3 families F (388) and L (352) on 14 points,
    rebuilt here from the paper's text, written independently of the searcher.
    W = {w1..w5, y}; H2 = {y,w_i,w_i+1}, {w_i,w_i-2,w_i+2} (mod 5)."""
    def brick_layers(off):
        pts = [off + i for i in range(6)]      # w1..w5 = off..off+4, y = off+5
        full = sum(1 << p for p in pts)
        h2 = []
        for i in range(5):
            h2.append((1 << pts[5]) | (1 << pts[i]) | (1 << pts[(i + 1) % 5]))
            h2.append((1 << pts[i]) | (1 << pts[(i + 3) % 5]) | (1 << pts[(i + 2) % 5]))
        return {0: [0], 1: [full & ~(1 << p) for p in pts], 2: h2, 3: [full & ~h for h in h2]}
    W1, W2 = brick_layers(0), brick_layers(6)
    X1, X2 = 1 << 12, 1 << 13

    def fam(quads):
        out = []
        for (i1, i2, j1, j2) in quads:
            for A in W1[i1]:
                for B in W2[i2]:
                    out.append(A | B | (X1 if j1 else 0) | (X2 if j2 else 0))
        return out
    F = fam([(1, 1, 0, 0), (2, 2, 1, 1), (1, 0, 1, 1), (0, 1, 1, 1),
             (1, 2, 1, 0), (3, 1, 1, 0), (1, 3, 0, 1), (2, 1, 0, 1)])
    L = fam([(1, 2, 0, 0), (2, 1, 0, 0), (2, 3, 1, 0), (3, 2, 0, 1),
             (1, 0, 1, 0), (0, 1, 0, 1), (3, 0, 1, 1), (0, 3, 1, 1)])
    return F, L


def _glue_explicit(A, na, a, B, nb, b):
    """Explicit Lemma-4.4 glue on ground [na-1 + nb-1] (a and b deleted)."""
    def drop(x, p):
        lo = x & ((1 << p) - 1)
        hi = x >> (p + 1)
        return lo | (hi << p)
    out = []
    for C in A:
        for D in B:
            if ((C >> a) & 1) + ((D >> b) & 1) == 1:
                out.append(drop(C, a) | (drop(D, b) << (na - 1)))
    return out


def _lemma_tests():
    """Brute-force the gluing lemma on small Sperner 3-free families and check
    the chain count recursion against an explicit two-step chain."""
    res = []
    fams = []
    g = _good_family()["sets"]
    fams.append((g, 6))
    k2 = [3, 5, 6]                                  # 2-subsets of [3]
    fams.append((k2, 3))
    c4 = [0b0111, 0b1011, 0b1101, 0b1110]          # 3-subsets of [4]
    fams.append((c4, 4))
    ok = True
    for (A, na), (B, nb) in itertools.product(fams, repeat=2):
        for a in range(na):
            for b in range(nb):
                G = _glue_explicit(A, na, a, B, nb, b)
                if len(set(G)) != len(G) or _find_sunflower(G) is not None \
                        or _find_containment(G) is not None:
                    ok = False
                Aa = sum(1 for C in A if (C >> a) & 1)
                Bb = sum(1 for D in B if (D >> b) & 1)
                if len(G) != Aa * (len(B) - Bb) + (len(A) - Aa) * Bb:
                    ok = False
    res.append(("gluing lemma brute force on %d family pairs x all ports" % (len(fams) ** 2), ok))
    import random as _r
    rng = _r.Random(12345)
    agree = True
    nsf = 0
    for trial in range(300):
        n = rng.randint(3, 7)
        fam = rng.sample(range(1 << n), rng.randint(3, min(14, 1 << n)))
        x1 = _find_sunflower(fam) is None
        x2 = _find_sunflower_cubic(fam) is None
        agree &= (x1 == x2)
        nsf += (not x2)
    res.append(("fast sunflower scan == cubic scan on 300 random families (%d with a sunflower)" % nsf,
                agree and 0 < nsf < 300))
    # explicit two-step chain on the 8-transversal family vs matrix recursion
    B, nb, u, v = g, 6, 0, 2
    P, npts, z = list(B), nb, v
    vec = [sum(1 for D in B if (D >> z) & 1), sum(1 for D in B if not (D >> z) & 1)]
    M = _step_matrix(B, u, v)
    ok2 = True
    for _ in range(2):
        P = _glue_explicit(P, npts, z, B, nb, u)
        z_new = (npts - 1) + (v if v < u else v - 1)
        npts = npts - 1 + nb - 1
        z = z_new
        vec = [M[0][0] * vec[0] + M[0][1] * vec[1], M[1][0] * vec[0] + M[1][1] * vec[1]]
        cnt = [sum(1 for D in P if (D >> z) & 1), sum(1 for D in P if not (D >> z) & 1)]
        if cnt != vec or _find_sunflower(P) is not None or _find_containment(P) is not None:
            ok2 = False
    res.append(("explicit 2-step chain matches M-recursion, is 3-free and Sperner", ok2))
    return res


def selftest():
    cases = []
    g = _good_family()
    F, L = _degkm_14()

    sun = {"n": 6, "sets": [0b000011, 0b000101, 0b001001, 0b010001]}
    cases.append(("plain: explicit 3-sunflower (common core {0})", sun, 1.0, "reject"))
    cases.append(("plain: duplicate member", {"n": 6, "sets": g["sets"] + [g["sets"][0]]}, 1.0, "reject"))
    true_v = power(mpf(8), mpf(1) / 6)
    cases.append(("plain: valid family, claim = true value + 1e-6", g, float(true_v + mpf("1e-6")), "reject"))
    out = {"n": 6, "sets": list(g["sets"])}
    out["sets"][0] = 1 << 6
    cases.append(("plain: member uses a point outside [n]", out, 1.0, "reject"))
    chain = {"n": 6, "sets": [0b000001, 0b000011, 0b000111, 0b001111]}
    cases.append(("plain: non-Sperner chain, claim uses |F|", chain,
                  float(power(mpf(4), mpf(1) / 6)), "reject"))
    cases.append(("plain: positive control (valid family, honest claim)", g,
                  float(mp.floor(true_v * 10 ** 12) / 10 ** 12), "accept"))
    cases.append(("plain: positive control (non-Sperner chain, largest-layer claim)", chain, 1.0, "accept"))
    fv = power(mpf(388), mpf(1) / 14)
    cases.append(("plain: DEGKM F (388 on 14, Sperner non-uniform) honest claim", {"n": 14, "sets": F},
                  float(mp.floor(fv * 10 ** 12) / 10 ** 12), "accept"))
    cases.append(("plain: DEGKM F claim = true value + 1e-6", {"n": 14, "sets": F},
                  float(fv + mpf("1e-6")), "reject"))
    # chain format: DEGKM L-F alternating chain (period on 24 points)
    ch = {"format": "c49-chain-v1", "bases": {"F": {"n": 14, "sets": F}, "L": {"n": 14, "sets": L}},
          "word": [{"base": "F", "u": 0, "v": 1}, {"base": "L", "u": 12, "v": 13}]}
    P = _matmul(_step_matrix(L, 12, 13), _step_matrix(F, 0, 1))
    a, b, c, d = map(mpf, (P[0][0], P[0][1], P[1][0], P[1][1]))
    chv = power(((a + d) + mpsqrt((a - d) ** 2 + 4 * b * c)) / 2, mpf(1) / 24)
    cases.append(("chain: DEGKM F/L word positive control (value %s)" % mp.nstr(chv, 12), ch,
                  float(mp.floor(chv * 10 ** 12) / 10 ** 12), "accept"))
    cases.append(("chain: DEGKM F/L word, claim = true value + 1e-6", ch, float(chv + mpf("1e-6")), "reject"))
    bad = json.loads(json.dumps(ch))
    bad["bases"]["F"]["sets"] = F[:-1] + [0b11]            # {0,1} is inside 5-sets -> not Sperner
    cases.append(("chain: non-Sperner base", bad, 1.0, "reject"))
    bad2 = json.loads(json.dumps(ch))
    bad2["bases"]["X"] = sun
    bad2["word"] = [{"base": "X", "u": 1, "v": 2}]
    cases.append(("chain: base containing a 3-sunflower", bad2, 1.0, "reject"))
    bad3 = json.loads(json.dumps(ch))
    bad3["word"][0]["u"] = 14
    cases.append(("chain: port outside the base ground set", bad3, 1.0, "reject"))
    bad4 = json.loads(json.dumps(ch))
    bad4["word"][0]["v"] = 0
    cases.append(("chain: u == v", bad4, 1.0, "reject"))
    bad5 = {"format": "c49-chain-v1", "bases": {"G": g}, "word": [{"base": "Y", "u": 0, "v": 1}]}
    cases.append(("chain: unknown base name", bad5, 1.0, "reject"))

    all_ok = True
    for name, obj, claim, want in cases:
        path = _write_tmp(obj)
        try:
            try:
                rep = verify(path, claim)
                got, detail = "accept", rep
            except Reject as exc:
                got, detail = "reject", str(exc)
        finally:
            os.unlink(path)
        ok = (got == want)
        all_ok = all_ok and ok
        print("[%s] want=%s got=%s | %s | %s" % ("PASS" if ok else "FAIL", want, got, name, detail[:160]))
    for name, ok in _lemma_tests():
        all_ok = all_ok and ok
        print("[%s] %s" % ("PASS" if ok else "FAIL", name))
    print("selftest %s (%d cases + lemma tests)" % ("PASS" if all_ok else "FAIL", len(cases)))
    return 0 if all_ok else 1


def main(argv):
    if "--selftest" in argv:
        return selftest()
    args = list(argv)
    if "--claim" not in args or len(args) < 3:
        print("usage: check_c49.py <certfile> --claim <float>  |  --selftest")
        return 2
    ci = args.index("--claim")
    try:
        claim = float(args[ci + 1])
    except (IndexError, ValueError):
        print("REJECT: --claim needs a float argument")
        return 1
    rest = [a for i, a in enumerate(args) if i not in (ci, ci + 1)]
    if len(rest) != 1:
        print("REJECT: expected exactly one certificate path, got %r" % (rest,))
        return 1
    try:
        print(verify(rest[0], claim))
        return 0
    except Reject as exc:
        print("REJECT: %s" % exc)
        return 1
    except Exception as exc:
        print("REJECT: unexpected error: %r" % (exc,))
        return 1


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
