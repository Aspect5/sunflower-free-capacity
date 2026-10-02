#!/usr/bin/env python3
"""Independent recheck of the C49 chain certificate (written independently of check_c49.py).

Checks, for the base family B (bitmasks on n points):
  S  Sperner: no A != B in the family with A subset of B (numpy, all ordered pairs).
  T  3-free: no 3 distinct sets with equal pairwise intersections.
     Method T1 (numpy): fix A; each other set S splits into a = S & A and r = S & ~A.
        B, C form a sunflower with A iff a_B == a_C and r_B & r_C == 0 (see README).
     Method T2 (C, sunflower_scan.c): plain triple loop, (C & (A|B)) == (A & B).
  R  Chain rate: M from port counts, exact rho via quadratic formula, rho^(1/(n-2)) to 30 digits.
Planted-bad inputs: a copy of B with one sunflower planted must fail T1 and T2;
a copy with a planted containment must fail S.

Usage: recheck_c49.py CERT [--plant]
"""
import json, sys, subprocess, os
from fractions import Fraction
import numpy as np
import mpmath

HERE = os.path.dirname(os.path.abspath(__file__))


def sperner_violations(sets):
    a = np.array(sets, dtype=np.int64)
    bad = 0
    for i in range(len(a)):
        # a[i] subset of a[j], j != i
        sub = (a[i] & a) == a[i]
        sub[i] = False
        bad += int(sub.sum())
    return bad


def sunflowers_t1(sets, n, stop_after=1):
    a = np.array(sets, dtype=np.int64)
    full = (1 << n) - 1
    found = []
    for i in range(len(a)):
        A = int(a[i])
        others = np.delete(a, i)
        inA = others & A
        out = others & (full ^ A)
        order = np.argsort(inA, kind="stable")
        inA_s, out_s = inA[order], out[order]
        bounds = np.flatnonzero(np.diff(inA_s)) + 1
        starts = np.concatenate(([0], bounds))
        ends = np.concatenate((bounds, [len(inA_s)]))
        for s, e in zip(starts, ends):
            if e - s < 2:
                continue
            g = out_s[s:e]
            disj = (g[:, None] & g[None, :]) == 0
            np.fill_diagonal(disj, False)  # r & r == 0 only if r == 0; one such set is not a pair
            if disj.any():
                x, y = np.argwhere(disj)[0]
                B = int(inA_s[s] | g[x]); C = int(inA_s[s] | g[y])
                found.append((A, B, C))
                if len(found) >= stop_after:
                    return found
    return found


def sunflowers_t2(sets):
    src = os.path.join(HERE, "sunflower_scan.c")
    exe = os.path.join("/tmp", "sunflower_scan_c49")
    subprocess.run(["gcc", "-O2", "-o", exe, src], check=True)
    inp = "\n".join(str(s) for s in [len(sets)] + list(sets)) + "\n"
    r = subprocess.run([exe], input=inp, capture_output=True, text=True, check=True)
    return r.stdout.strip()


def is_sunflower(A, B, C):
    return len({A, B, C}) == 3 and (A & B) == (A & C) == (B & C)


def chain_rate(sets, n, u, v):
    bu, bv = 1 << u, 1 << v
    c = {(vi, ui): 0 for vi in (1, 0) for ui in (1, 0)}
    for s in sets:
        c[(int(bool(s & bv)), int(bool(s & bu)))] += 1
    # (|P'_v|, |P'^v|) = M (|P_z|, |P^z|); G1 uses D with u out, G2 uses D with u in.
    M = [[c[(1, 0)], c[(1, 1)]], [c[(0, 0)], c[(0, 1)]]]
    tr = M[0][0] + M[1][1]
    det = M[0][0] * M[1][1] - M[0][1] * M[1][0]
    mpmath.mp.dps = 50
    rho = (tr + mpmath.sqrt(tr * tr - 4 * det)) / 2
    return M, tr, det, rho, rho ** (mpmath.mpf(1) / (n - 2))


def main():
    cert = json.load(open(sys.argv[1]))
    assert cert["format"] == "c49-chain-v1" and len(cert["word"]) == 1
    w = cert["word"][0]
    base = cert["bases"][w["base"]]
    n, sets = base["n"], [int(s) for s in base["sets"]]
    print("n", n, "|B|", len(sets), "distinct", len(set(sets)),
          "max mask ok", max(sets) < (1 << n), "sizes", sorted({bin(s).count('1') for s in sets}))
    sv = sperner_violations(sets)
    print("Sperner violations", sv)
    t1 = sunflowers_t1(sets, n)
    print("T1 sunflowers found", t1)
    t2 = sunflowers_t2(sets)
    print("T2", t2)
    bad = bool(sv) or bool(t1) or "triples=0" not in str(t2).split()   # fail closed (v1.1)
    M, tr, det, rho, rate = chain_rate(sets, n, w["u"], w["v"])
    print("M", M, "tr", tr, "det", det)
    print("rho", mpmath.nstr(rho, 30))
    print("rho^(1/%d)" % (n - 2), mpmath.nstr(rate, 30))
    claim = mpmath.mpf("1.5617847089505")
    print("rate - 1.5617847089505 =", mpmath.nstr(rate - claim, 5))
    if "--plant" in sys.argv:
        # planted sunflower: take A, B in B with A&B = K; add C = K | (fresh petal) where
        # the petal is a nonempty subset of points outside A|B, then C&A = C&B = K.
        # Drop any set that C contains or is contained in so the copy stays Sperner and the
        # only defect is the sunflower.
        full = (1 << n) - 1
        for i in range(len(sets)):
            for j in range(i + 1, len(sets)):
                A, B = sets[i], sets[j]
                rest = full ^ (A | B)
                if rest:
                    C = (A & B) | rest
                    break
            else:
                continue
            break
        assert is_sunflower(A, B, C)
        plant = [s for s in sets if not ((s & C) == s or (s & C) == C)] + [C]
        assert A in plant and B in plant
        print("PLANT sunflower C", C, "with", A, B, "copy size", len(plant),
              "Sperner viol", sperner_violations(plant))
        print("PLANT T1", sunflowers_t1(plant, n))
        print("PLANT T2", sunflowers_t2(plant))
        plant2 = sets + [sets[0] | (full ^ sets[0]) & -((full ^ sets[0])) ]
        print("PLANT containment Sperner viol", sperner_violations(plant2))
    if bad:
        print("RECHECK: FAIL"); sys.exit(1)
    print("RECHECK: OK")


if __name__ == "__main__":
    main()
