# A lower bound on the 3-sunflower-free capacity

**Theorem.** Let f(n) be the largest size of a family of subsets of [n] containing no three distinct sets A, B, C with A∩B = A∩C = B∩C, and let μ₃ = lim f(n)^{1/n}. Then

  μ₃ ≥ ρ^{1/18} = 1.56477309608283...,

where ρ = (3901 + √5882425)/2 = 3163.18555281243... is the Perron root of the 2×2 matrix with rows (1993, 1360) and (1080, 1908) (characteristic polynomial x² − 3901x + 2333844).

(v1.2, 17 Oct 2026, from a 6,341-set base family, `cert_c49_chain20_v12.json`. Earlier versions are kept here: v1.1, 2 Oct 2026, gave μ₃ ≥ 1.56212205601 from a 6,150-set family, `cert_c49_chain20_s17.json`; v1.0, 30 Sep 2026, gave μ₃ ≥ 1.56178470895 from a 6,127-set family, `cert_c49_chain20.json`.)

Previous bounds: μ₃ > 1.551 (Deuber, Erdős, Gunderson, Kostochka, Meyer, *J. Combin. Theory Ser. A* 79 (1997) 118-132, "DEGKM"), and μ₃ ≥ 1.554, cited as unpublished (Naslund; see Naslund and Sawin, arXiv:1606.09575). The 1.554 construction was not available to us and was not checked. The constant is entry [49a](https://teorth.github.io/optimizationproblems/constants/49a.html) of T. Tao's optimization-constants repository.

## Method
The chain construction is DEGKM's. Their Lemma 3.1 (direct sums) and Lemma 4.4 (gluing two Sperner 3-sunflower-free families along a marked point of each, losing 2 points) are chained in their Theorem 1.5, which bounds the growth with an averaging step (a factor of 1/2). The new ingredients here are:
1. **Exact counting of the chain.** Gluing the marked point of the running family to port u of a copy of a base family B on m points, and marking port v of that copy, maps the pair (#sets containing the marked point, #sets avoiding it) linearly by M = [[#(v∈S, u∉S), #(v∈S, u∈S)], [#(v∉S, u∉S), #(v∉S, u∈S)]]. So f grows like ρ(M)^{k} while the ground set grows by k(m−2). Since f is nondecreasing, μ₃ ≥ ρ(M)^{1/(m−2)}. On DEGKM's own base families (included as `cert_degkm_control.json`), this exact count already gives 1.5529675, above their stated 1.551.
2. **A new base family:** 6,341 sets on m = 20 points, Sperner and 3-sunflower-free, with ports u = 18 and v = 19 (`cert_c49_chain20_v12.json`; layer sizes 7/10/11/12/13 with 2/551/1108/3600/1080 sets). Compared with the v1.1 family, nearly all of the gain is in the class with v∈S, u∉S (1802 → 1993 sets). The family was found by an exact search over block structures ("types") followed by a local search over individual sets, with exact incremental counts of containments and sunflowers.

## Verify
- `python3 check_openai_v12.py` (Python standard library only; about 15 seconds) checks the v1.2 certificate from scratch: all 20,100,970 pairs for containment, an exhaustive sunflower search over all 292,854,591 candidate third sets, the port counts and matrix, and the bound by exact rational root isolation (1.564773096082 < ρ^{1/18} < 1.564773096083). It runs planted-sunflower, planted-containment and inflated-bound controls, and writes `REPORT_openai_v12.md`, which also re-derives the gluing lemma and the limit argument. It exits 0 on ACCEPT and 1 on REJECT.
- `python3 check_c49.py cert_c49_chain20_v12.json --claim 1.564773096` (requires numpy and mpmath) checks the Sperner property, 3-sunflower-freeness and the port counts exactly, and computes ρ from the exact integer matrix. It prints `value=1.5647730960828312992 ... beats_1.551=YES beats_1.554=YES`. The same command with `cert_c49_chain20_s17.json --claim 1.562122056` or `cert_c49_chain20.json --claim 1.56178470895` checks the earlier versions, and with `cert_degkm_control.json --claim 1.5529675` reproduces the DEGKM control.
- `python3 recheck_c49.py <certificate>` (with `sunflower_scan.c`; exits 1 on any violation) is a third checker, written independently of `check_c49.py`. It includes an exhaustive scan in C of all C(|B|,3) triples and planted-defect controls.

## Status
- v1.2 is checked by programs written independently by models from two different providers: `check_openai_v12.py` was written from the definitions and the DEGKM paper by OpenAI's GPT-6.1 Sol (via Codex), without access to the other checkers; `check_c49.py` and `recheck_c49.py` were written by separate Claude sessions. All accept.
- The derivation from DEGKM Lemmas 3.1 and 4.4 has been re-derived by two reviewers. One found a small slip in the printed proof of Lemma 4.4, which does not affect the statement.
- Not peer-reviewed; no formal (Lean) proof.
- Literature search (30 Sep and 1 Oct 2026: Semantic Scholar, arXiv, GitHub and the history of the 49a page; repeated by an outside AI referee in October 2026): no lower bound newer than 1.551 / 1.554 was found.

## How this was found (AI disclosure)
- **Who:** AI agents (Anthropic Claude Opus 5.5, in Claude Code sessions) working for the Institute for Recursive Self Improvement (IRSI). The v1.2 independent checker was written by OpenAI GPT-6.1 Sol.
- **How the problem was chosen:** an agent swept the construction-based lower-bound entries of the optimization-constants repository and picked six with visible slack (49a, 10c, 5b, 37, 44 and 71). Only 49a improved. For 10c, 5b and 37 the searches reproduced or tied the records; for 44 (11/8 against a record of 3/2) and 71 (5.252 against 6.522) they fell short.
- **Cost:** about 33 CPU core-hours of search for v1.2 and a few hours of agent time, at $0 marginal cost (subscription and cloud-session credits).
- **Prompt:** the instruction given to the agents was to find constants with finite certificates, search for better certificates, and verify every candidate with a checker sharing no code with the searcher. The full task file is not published.

Contact: Peter Seelman, IRSI (peter@irsi.ai).
