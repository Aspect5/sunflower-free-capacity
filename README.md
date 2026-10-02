# A lower bound on the 3-sunflower-free capacity

**Theorem.** Let f(n) be the largest size of a family of subsets of [n] containing no three distinct sets A, B, C with A∩B = A∩C = B∩C, and let μ₃ = lim f(n)^{1/n}. Then

  μ₃ ≥ ρ^{1/18} = 1.56212205600681...,

where ρ = 1855 + √1471609 = 3068.09892424320... is the Perron root of the 2×2 matrix with rows (1802, 1360) and (1080, 1908) (characteristic polynomial x² − 3710x + 1969416).

(v1.1, 2 Oct 2026. v1.0, 30 Sep 2026, gave μ₃ ≥ 1.56178470895 from a 6,127-set base family, `cert_c49_chain20.json`, which is kept here.)

Previous bounds: μ₃ > 1.551 (Deuber, Erdős, Gunderson, Kostochka, Meyer, *J. Combin. Theory Ser. A* 79 (1997) 118-132, "DEGKM"), and μ₃ ≥ 1.554, cited as unpublished (Naslund; see Naslund and Sawin, arXiv:1606.09575). The 1.554 construction was not available to us and was not checked. The constant is entry [49a](https://teorth.github.io/optimizationproblems/constants/49a.html) of T. Tao's optimization-constants repository.

## Method
The chain construction is DEGKM's. Their Lemma 3.1 (direct sums) and Lemma 4.4 (gluing two Sperner 3-sunflower-free families along a marked point of each, losing 2 points) are chained in their Theorem 1.5, which bounds the growth with an averaging step (a factor of 1/2). The new ingredients here are:
1. **Exact counting of the chain.** Gluing the marked point of the running family to port u of a copy of a base family B on m points, and marking port v of that copy, maps the pair (#sets containing the marked point, #sets avoiding it) linearly by M = [[#(v∈S, u∉S), #(v∈S, u∈S)], [#(v∉S, u∉S), #(v∉S, u∈S)]]. So f grows like ρ(M)^{k} while the ground set grows by k(m−2). Since f is nondecreasing, μ₃ ≥ ρ(M)^{1/(m−2)}. On DEGKM's own base families (included as `cert_degkm_control.json`), this exact count already gives 1.5529675, above their stated 1.551.
2. **A new base family:** 6,150 sets on m = 20 points, Sperner and 3-sunflower-free, with ports u = 18 and v = 19 (`cert_c49_chain20_s17.json`; layer sizes 7/10/11/12/13 with 2/360/1108/3600/1080 sets). v1.0 used a 6,127-set family of the same shape (`cert_c49_chain20.json`). Both were found by a search over families built from blocks of DEGKM-style "bricks" and singletons, with the Sperner and sunflower-free conditions enforced exactly.

## Verify
Requires Python 3 with numpy and mpmath.
- `python3 check_c49.py cert_c49_chain20_s17.json --claim 1.562122056` (v1.1 bound; `cert_c49_chain20.json --claim 1.56178470895` for v1.0) checks the Sperner property, 3-sunflower-freeness and the port counts exactly, and computes ρ from the exact integer matrix. It prints `value=1.5617847089505372797 ... beats_1.551=YES beats_1.554=YES`. The same command on `cert_degkm_control.json` with `--claim 1.5529675` reproduces the DEGKM control.
- `python3 recheck_c49.py cert_c49_chain20_s17.json` (with `sunflower_scan.c`; exits 1 on any violation since v1.1) is a second checker, written independently of the first. It includes an exhaustive scan in C of all C(|B|,3) triples (about 3.9×10¹⁰ for the 6,150-set family; 0 sunflowers in both) and planted-defect controls (a planted sunflower and a planted containment are both rejected).

## Status
- Checked by two computer programs written independently. The derivation from DEGKM Lemmas 3.1 and 4.4 was re-derived by the second reviewer. That review found a small slip in the printed proof of Lemma 4.4, which does not affect the statement.
- Not peer-reviewed; no formal (Lean) proof.
- Literature search (30 Sep 2026, repeated on GitHub 1 Oct 2026: Semantic Scholar, arXiv, and the history of the 49a page): no lower bound newer than 1.551 / 1.554 was found.

## How this was found (AI disclosure)
- **Who:** AI agents (Anthropic Claude Opus 5.5, in Claude Code sessions) working for the Institute for Recursive Self Improvement (IRSI).
- **How the problem was chosen:** an agent swept the construction-based lower-bound entries of the optimization-constants repository and picked six with visible slack (49a, 10c, 5b, 37, 44 and 71). Only 49a improved. For 10c, 5b and 37 the searches reproduced or tied the records; for 44 (11/8 against a record of 3/2) and 71 (5.252 against 6.522) they fell short.
- **The recheck:** a second Claude session, on a different account, with its own code. This is not a different model provider.
- **Cost:** a few CPU-hours of search and a few hours of agent time, at $0 marginal cost (subscription and cloud-session credits).
- **Prompt:** the instruction given to the agents was to find constants with finite certificates, search for better certificates, and verify every candidate with a checker sharing no code with the searcher. The full task file is not published.

Contact: Peter Seelman, IRSI (peter@irsi.ai).
