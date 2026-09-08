# Translation round: the anti-period reduced to a Schur-coefficient palindrome

Scratch, local-only. Gates: `_rb_trans.py` (target + Lever B d(q) catalogue),
`_rb_trans2.py` (36-weight alphabet + e_n build + content), `_rb_trans3.py`
(alternant-pinning attempt), `_rb_trans4.py` (symmetric-function gates A1-A3 + full-37).
Run from repo root, `PYTHONIOENCODING=utf-8`. Exact over ℤ.

## Target (reconfirmed, 0 fails all k)
`F_k(m+44) = −F_k(m)` for both partners in-wall `|·| ≤ 32−2k` (⇔ R_k vanishes on the window
⇔ the two-row law). In n-units (`n = 18−m/2`, `c_k(n) = (−1)ⁿ F_k(36−2n)`): the reflection is the
**palindrome `c_k(n) = c_k(14−n)` about n = 7**, and with the proven reciprocity (anti-palindrome
`c_k(36−n) = −c_k(n)` about n = 18) it is the **anti-period `c_k(n+22) = −c_k(n)`**. 18−7 = 11, 2·11 = 22.

## Lever B — the q-period route: DEAD, now catalogued across all (k,m)
`d_{k,m}(q) = g_{k,m}(q) − g_{k,44−m}(q) = q⁵·(q²−1)·h_{k,m}(q²)`, h a polynomial in q² with
**non-cyclotomic roots** — no roots of unity, so no period-22 specialisation closes it.
Near-wall examples (degree-2 in q²): k=0 → (4q²−1), k=2 → (11q²−3), k=4 → (9q²−5)
[roots q² = 1/4, 3/11, 5/9]. The factor (q²−1) makes `d(1)=0` automatic, so "theorem ⇔ (q²−1)|d(q)"
is a true restatement but **circular as a proof** (d(1)=0 IS the theorem). Confirms + extends the
prior round's negative.

## Lever A — the slice reformulation: the anti-period is a symmetric-function statement
The t₁^m-slice of X factorises (REDUCTION ROUND) as `X_m = Δ̂'·M²·(−1)ⁿ·e_n(w)`, and in-window
`F_k(m) = content₅(X_m; tail_k)`, `tail_k = (2k+10,8,6,4,2)`. The **36-weight alphabet is explicit**:
```
w = { 2·e_u : u = 1..5 }                            (E_V, 5 axis weights)
  ∪ { −2·L' : L' ∈ {±1}⁵ \ {(1,1,1,1,1)} }          (E_U, 31 hypercube vertices minus the apex)
```
i.e. the full ±hypercube of vertices **minus its apex −2·(1⁵)** (the deleted top weight), plus the
five axis weights. `e_n(w) = [xⁿ]∏(1 + x·w_i)`.

**GATED exact over ℤ (`_rb_trans4.py`):**
- **A1** `e_n(w)` is **S₅-symmetric** (0/312, n=0..12) — because both E_V and E_U are S₅-orbits.
  Consequently the *naive* full-W(C₅) content of e_n is 0; the content is carried entirely by the
  slice's Vandermonde alternant Δ̂'·M², making `L(f) := content₅(Δ̂'·M²·f; tail)` a **Weyl / Schur
  read-off** with `L(e_n) = c_k(n)`.
- **A2** deletion recurrence **`e_n(w₃₆) = e_n(w₃₇) − apex·e_{n−1}(w₃₆)`**, apex = −2·(1⁵) (0/12) —
  the deleted-sheet mechanism as a clean symmetric-function recurrence.
- **A3** reciprocity **`e_{36−n}(w) = prod(w)·e_n(1/w)`**, prod(w) = (4,4,4,4,4) (0/13) — the PROVEN
  reciprocity, now confirmed as a Laurent-polynomial identity of the alphabet.
- full-37 reciprocity `e_{37−n}(w₃₇) = prod₃₇·e_n(1/w₃₇)`, prod₃₇ = (2,2,2,2,2) = M² (0/14).

**Reformulation (the round's product):** the affine anti-period is EQUIVALENT to the Schur-coefficient
palindrome `[s_{λ*(k)}] e_n(w) = [s_{λ*(k)}] e_{14−n}(w)` (⇔ `c_k(n) = c_k(14−n)`), where λ*(k) is the
partition selected by tail_k under the Weyl read-off. The problem is thereby moved out of
deleted-sheet sine-products into the **Schur expansion of e_n over an explicit S₅-symmetric,
hypercube-minus-apex alphabet**.

## Residual, sharpened
Prove `c_k(n) = c_k(14−n)` (the palindrome about 7). The natural mechanism is A2: if the **undeleted
full-37 system** carries the palindrome, the deletion recurrence transports it and its failure term is
exactly R_k (the defect). Open bookkeeping: pin (δ, μ0) of the alternant Δ̂'·M² to make λ*(k) and
`[s_{λ*}]` fully explicit (an overcounting factor — S₅-stabilizer — was seen in the naive
`_rb_trans3.py` implementation; not needed for the reformulation's validity, needed to run the Schur
route to the end).

## Net
Lever B closed (no q-period). Lever A converts the sole open residual (the affine coroot translation)
into a concrete symmetric-function palindrome over an explicit alphabet, with reciprocity + the
deletion mechanism gated exactly. The theorem stands as a FACT (0 fails all k); as a PROOF, the finite
reflection (negation) and reciprocity are done, and the remaining reflection/translation is now the
Schur palindrome `c_k(n)=c_k(14−n)`.
