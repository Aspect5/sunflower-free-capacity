VERDICT: ACCEPT

| Check | Observed result |
|---|---|
| Raw input | SHA-256 `7472112ced9d98da8cf7b6a3d586b93ea66750e39989e7a6d6d5f52ffeb4e82e` |
| Distinctness / range | 6341 distinct integer masks; all in [0, 2^20); PASS |
| Size distribution | 7: 2; 10: 551; 11: 1108; 12: 3600; 13: 1080 |
| Sperner | 20,100,970 unordered pairs, both containment directions; witness None |
| 3-sunflower-free | 20,100,970 unordered pairs; 292,854,591 candidate third sets; witness None |
| Port counts (u,v) | (0, 0): 1080; (0, 1): 1993; (1, 0): 1908; (1, 1): 1360 |
| Derived matrix | `[[1993, 1360], [1080, 1908]]`; agrees with certificate: True |
| Characteristic polynomial | p(t) = t^2 - 3901t + 2333844; discriminant 5882425 |
| Exact lower bound | L = 782386548041/500000000000; L^18 > 3901/2: True; p(L^18) < 0: True |
| Planted sunflower control | Add masks 1,2,4 to actual family; REJECT, witness (1, 2, 1043712) |
| Planted containment control | Add mask 0 to actual family; REJECT, witness (0, 47071) |
| Inflated bound control | U = 1564773096083/1000000000000; U^18 > 3901/2: True; p(U^18) > 0: True; REJECT |
| Algorithm controls | All 256 families on {0,1,2} agree with direct pair/triple definitions; duplicate and out-of-range inputs rejected |
| Lemma hypotheses / exponent | Disjoint relabelled copies; Sperner and 3-free verified; n_j = 18j+2; rate exponent 1/18 |

## Transfer matrix derived from the branches

Write c_ij for the number of base members having u-membership i and v-membership j. Use a column vector x=(x_1,x_0)^T counting current members with the old output port present and absent. Branch I takes x_1 members with base u absent, so its contributions to the new present/absent states are c_01 x_1 and c_00 x_1. Branch II takes x_0 members with base u present, contributing c_11 x_0 and c_10 x_0. Thus x_new = [[c_01,c_11],[c_00,c_10]] x. These four classes give exactly the displayed certificate matrix in present-then-absent order. No certificate count, matrix entry, exponent, or claimed rate was used to construct it.

## Exact root isolation

The larger eigenvalue is rho = (3901 + sqrt(5882425))/2. Both off-diagonal entries are positive, so the discriminant is positive and rho > trace/2. The polynomial is strictly increasing to the right of trace/2. The following reduced Fraction values, computed with integer arithmetic, isolate the root:

- L^18 = `12066595278862709436583017244046660479630466854567658204891232507374477472057991968818915755902758073930105928164272736704832006434792217628082287134217265501731320011261463278412281789312317869684522201221337053521/3814697265625000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000`.
- p(L^18) = `-1067582697265777946033150318846161671600201778916217431613761419873606919060603623774196221603110588004684103746702484811988897240179131245905318585117090701826293993540937833634573275672828190316418834580529079633906314553409752622484902255607510410383400246355963538449404732279772517236043800446890520659739154228295495917659590313831201928604437501976924738466794241295953722901726124260011374910469296641981502559/14551915228366851806640625000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000`.
- U^18 = `3163185552818573065032271145489820164409961722920668719792717245674662621902283634829052865669183745976836000721838593912430464643323896680406568352622520105020756895715756080936536610573466305318973690441333474795146409/1000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000`.
- p(U^18) = `14888163112045126560316558866817835908802105808753590794744529435517323284945778527256187382016658262458384955069229421515766737708545499369724398688742703989642170249874735682939889181851836401707996234997041630278433396580013289632035128716096088465370618714931481381010938662778783527924539984917901845841307767755574308733500770587112538679657510001249215743307613211768593122033997997100353132414312716550677447543745595281/1000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000`.

Hence L^18 < rho < U^18, and L < rho^(1/18) < U. The lower claim passes strictly and the inflated claim fails strictly. No floating-point calculations enter the checker or this evidence.

## Exhaustiveness and the infinite chain

For any fixed distinct A,B, equality of the three pairwise intersections forces C to contain every bit of A&B and exclude every bit of their symmetric difference. All other bits lie outside A|B and are free. The submask loop enumerates each such C exactly once, including the empty submask, excluding only C=A or C=B. Every unordered pair is tested, so every possible sunflower is covered (usually three times). Containment is separately tested in both directions for every unordered pair.

The original [DEGKM paper, Lemma 4.4, printed p.129 / PDF page 12](https://kostochk.web.illinois.edu/docs/old/jcta97degm.pdf) was checked directly. It applies to Sperner, 3-free families on disjoint ground sets, with one marked point in each. Its two branches match those in TASK.md and delete both marked points. Each step uses a freshly relabelled copy of B, so the ground sets are disjoint before gluing. The lemma preserves the two properties; v survives because u and v are distinct.

For completeness, deletion is injective within each branch because the deleted port has fixed membership. Across branches an inclusion of output sets would force a strict inclusion in one of the original families, so the branches are disjoint and their union is Sperner. If a hypothetical sunflower uses one branch only, each factor projection has equal pairwise intersections. In a Sperner 3-free family such a triple, allowing repetitions, must have all projections identical: three distinct projections violate 3-freeness, and exactly two identical ones force containment of the repeated member in the third. Thus all three outputs would be identical. If two outputs use branch I and the third branch II, restore the base input port to the third projection. It occurs once and leaves all pairwise intersections equal. Three distinct base projections violate 3-freeness; a repeated projection forces strict containment. The other mixed case follows by interchanging the factors. This independently checks the lemma mechanism.

Begin with B on 20 points. Its state vector is x_1=(3353, 2988)^T, with both entries positive. After j copies the vector is M^(j-1)x_1 and the ground size is 20j-2(j-1)=18j+2. The positive matrix has a strictly positive eigenvector for rho. For example (M_01, rho-M_00)^T is positive, since rho exceeds both diagonal entries. Bound x_1 above and below by positive multiples of that vector and apply M repeatedly. The resulting family sizes are bounded above and below by positive constants times rho^(j-1). Taking (18j+2)-th roots gives the limiting rate rho^(1/18). Since these are actual sunflower-free families, they imply the claimed lower bound for mu_3. For a full-limit definition of mu_3, let G(n) be the largest size of a Sperner 3-free family. Disjoint products preserve both properties by the projection argument above, so G is supermultiplicative and Fekete's lemma gives its exponential-rate limit. If F(n) is the largest unrestricted 3-free size, then G(n) <= F(n) <= (n+1)G(n): some cardinality layer of any family contains at least its size divided by n+1. Thus F and G have the same exponential-rate limit.

## Concerns and reproducibility

The unusual size-7 members are included in every check; no uniformity assumption was made. No unresolved concern affects the claim. This checks the supplied construction and rate, not novelty or optimality. Run `python3 check.py` beside the certificate to regenerate this report; Python standard library only. A failed family/property or certificate comparison produces REJECT and exit status 1; malformed input exits with an error.
