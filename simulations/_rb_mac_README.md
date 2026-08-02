# Route 1 (the twisted-affine genre lead): A_12^(2) worked out concretely

Scratch, local-only (`_rb_mac*.py`). Run from repo root, `PYTHONIOENCODING=utf-8`. The two-row
reflection law `n_(j,k)=n_(10-j,k)` (⇔ `mu1 -> 22-mu1`, `m=2mu1 -> 44-m`) has the pinned symmetry
signature: infinite dihedral `<s: mu1->-mu1 (char -1), s0: mu1->22-mu1 (char +1)>`, translation
`22Z` carrying `-1`, wall `mu1=11`. This round tests the resumption point's genre lead: **is X (or the
two-row slice) a finite specialization of a twisted-affine A_12^(2) denominator/numerator identity?**

Verdict up front: **PARTIAL.** The Lie-type / character signature is the twisted A_12^(2) genre
exactly (so(13) horizontal, spin weights, mixed `(-,+)` characteristic). But X is provably **not** the
finite specialization of the A_12^(2) Weyl-Kac / Macdonald denominator, and the operative wall/period
arithmetic (11/22/44) is the Langlands-partner C_6 = Sp(12) level-4 scale, **not** A_12^(2)'s native
scale 13. Concrete refutation reasons + numeric checks below.

## (1) The twisted system, derived and pinned: A_2n^(2), n=6 (= A_12^(2))

Source for the concrete denominator: Rosengren-Schlosser, *Elliptic determinant evaluations and the
Macdonald identities for affine root systems*, arXiv:math/0505213. **Dictionary fact (load-bearing):**
Macdonald's *reduced* affine root system **BC_n IS the twisted Kac-Moody A_2n^(2)** (Macdonald A<->A^(1),
B<->B^(1), B^v<->A_{2n-1}^(2), C<->C^(1), C^v<->D_{n+1}^(2), **BC<->A_2n^(2)**, D<->D^(1)). So the paper's
`W_{BC_6}` is exactly the A_12^(2) Weyl-Kac denominator.

| datum | value for A_12^(2) (n=6) | how known |
|---|---|---|
| horizontal finite subalgebra | **so(13) = B_6** (grade-0 fixed algebra of the order-2 automorphism of sl(13)=A_12) | Kac; the 64 sheets = spin weights of so(13) (F133 proof §5) |
| finite root system | **non-reduced BC_6** = {±e_i (short, B), ±2e_i (long, C), ±e_i±e_j} | contains BOTH so(13) short roots and sp(12) long roots |
| finite Weyl group | W(B_6)=W(C_6), hyperoctahedral, order 2^6·6! | matches the W(C_6) antisymmetry of X |
| dual Coxeter number h^v | **2n+1 = 13** | quasi-period exponent + Macdonald series scale (both = 2n+1), verified `_rb_mac3.py` |
| affine real roots | 3 lengths: `±e_i+kd` (from θ(x_i)), `±2e_i` (from θ(px_i^2;p^2)), `±e_i±e_j` (from θ(x_ix_j^±)) | read off the product form below |

**The A_12^(2) = BC_6 Weyl-Kac denominator (verbatim, math/0505213):**
```
   W_{BC_n}(x) = prod_i theta(x_i) theta(p x_i^2; p^2) * prod_{i<j} x_i^{-1} theta(x_i x_j) theta(x_i/x_j)
   theta(x;p) = prod_{k>=0} (1 - p^k x)(1 - p^{k+1}/x)          (so theta(x;0) = 1-x)
   BC_n theta quasi-period:  f(px) = p^{-n} x^{-(2n+1)} f(x),  f(1/x) = -x^{-1} f(x)   [exponent 2n+1]
```
**The Macdonald identity (series == product, verbatim):**
```
   (p;p)_inf^n W_{BC_n}(x) = sum_{m in Z^n} sum_{sigma in S_n} sgn(sigma)
        prod_i x_i^{(2n+1) m_i} p^{(2n+1) C(m_i,2)+n m_i} ((x_i p^{m_i})^{sigma(i)-n} - (x_i p^{m_i})^{n+1-sigma(i)})
```
The `(2n+1)=13` in both the quasi-period and the series is the object's native level scale.

## (2) The candidate, and its refutation WITH REASONS

**Candidate.** X = Dhat * Prod_31 is the finite (classical, p->0) specialization of `W_{BC_6}`, so the
two-row law is the finite shadow of the A_12^(2) affine reflection.

**Refuted, three independent structural reasons:**

- **R1 - factor count / homogeneity (structural, exact).** X = `prod_{u<v} sin((x_u-x_v)/2)` (15
  differences) * `prod_31 sin(L.x/2)` (31 spin sheets) = **46 half-angle sine factors**. The classical
  limit of `W_{BC_6}` is the ordinary **B_6 Weyl denominator** `prod_u sin(x_u/2) * prod_{u<v}
  sin((x_u+x_v)/2) sin((x_u-x_v)/2)` = 6 single + 15 sum + 15 diff = **36 factors**. 46 != 36: different
  homogeneity degree, so X cannot be proportional to it. (`_rb_mac2.py` T-A)
- **R2 - the spinor product is absent from the denominator (the deepest reason).** The A_12^(2)
  denominator is built from `theta(x_i)`, `theta(px_i^2)` and PAIRWISE `theta(x_i x_j^±)` only - rank-<=2
  factors. It contains **no spinor factor** (no product over the 2^{n-1} sign-sheets `L.x/2`, which are
  genuinely six-fold). So the "B_6 spin weight system appears in the finite specialization, with a 1^6
  deletion" hypothesis is **false for the A_12^(2) denominator**: the spin weights never appear, hence
  there is no spinor factor from which to delete the top sheet. The spinor product is a different
  (dual-pair / super-denominator) object, not a Kac-Moody denominator. Numerator route also fails: a
  Weyl-Kac numerator A_{Lambda+rho} is a single alternant (determinant), not a 46-fold product.
- **R3 - level/period mismatch (arithmetic).** Measured: wall `mu1=11`, period `mu1->mu1+22`, `m->m+44`.
  - Sp(12)=C_6 dual reading: wall = `l + h^v(C_6) = 4 + 7 = 11`  => level 4, period 2*11=22. **MATCHES.**
  - A_12^(2) native reading: `h^v = 2n+1 = 13`; `11 = l+13 => l = -2` (negative); 13 does not divide
    22 or 44, and the Macdonald series scale is 13, not 11. So the reflection is **not** at an A_12^(2)
    integral level; the operative scale is the Langlands-partner C_6=Sp(12) side (where the n_lambda
    characters already live). A_12^(2) supplies the Weyl group and the non-reduced B/C bridge, not the
    period.

**What genuinely MATCHES the twisted genre (the partial-credit half):**
- The finite Weyl group W(C_6), the horizontal so(13), and the spin-weight sheets are exactly the
  A_12^(2) data (F133 proof §5 already names so(13)/spin and its Langlands dual Sp(12)).
- The mixed `(-,+)` character is the native twisted/non-reduced signature: the affine wall belongs to a
  root of different length (the long C-root 2e_1) than the finite walls (short), so X's extra affine
  symmetry `s0` comes with the OPPOSITE sign to the finite `s` - precisely `(char s, char s0)=(-1,+1)`.
  An untwisted numerator would be `(-1,-1)`; a plain theta `(+1,+1)`. Only the twisted/non-reduced side
  produces `(-1,+1)`, and A_12^(2) is the unique affine system bridging so(13) (weights) and sp(12)
  (characters) via the non-reduced BC_6. **Genre correct; specific identity not.**

## (3) Numerical verifications (exact / machine precision, from below)

- **D_inf structure on live F_k(m)** (`_rb_mac1.py`, exact over Z via dict `Xc`):
  `s`  (m->-m,  char -1): `F_k(-m) = -F_k(m)`      **0 mismatches** (all k, m in [-60,60]).
  `s0` (m->44-m, char +1): `F_k(44-m) = F_k(m)`     **0 of 36 in-wall pairs**.
  translation (m->m+44, char -1): `F_k(m+44)=-F_k(m)` **0 of 36 in-wall pairs**.
  Confirms the pinned `(-,+)` signature directly on the data. (Also: `F_k(22)` nonzero only for odd k
  (=-2,-3): the affine center is even, consistent with char(s0)=+1.)
- **A_12^(2) = BC_6 denominator is the genuine Weyl-Kac denominator** (`_rb_mac3.py`): Macdonald identity
  series == product at nome p=0.06, worst rel dev **3.6e-14** over 3 points (S_n sum collapsed to a 6x6
  determinant by multilinearity; sum over Z^6 truncated at |m|<=8).
- **Classical limit = B_6 denominator** (`_rb_mac2.py` T-C): `W_{BC_6}|_{p->0} / D_{B_6}^{det}` is
  **constant** across 8 random points (spread **4.2e-13**), i.e. the p->0 limit IS the B_6 Weyl
  denominator (determinant form 1.2b).
- **X is NOT that denominator** (`_rb_mac2.py` T-B): `X / D_{B_6}^{sine}` is **non-constant** across 8
  points (spread **1.4e-1**). Together with R1's factor count: X != finite spec of `W_{BC_6}`.
- **l=2 break guard.** The two-row law holds `l in {0,1,3}` and breaks at three-row `l=2` (schur round,
  0/36 two-row; 8/16 at l=2). A correct identity MUST break at l=2. The C_6 level-4 affine-reflection
  reading (rho_4-tail resonance at l=2) is consistent with this; the A_12^(2) denominator candidate was
  refuted before this guard even applies (R1-R3), so no false positive was accepted.

## (4) VERDICT

**PARTIAL.**
- MATCH: the genre. Horizontal so(13)=B_6, spin-weight sheets, finite Weyl group W(C_6), and the mixed
  `(-,+)` characteristic are exactly the twisted A_12^(2) signature (the unique affine system bridging
  so(13) weights and sp(12)=C_6 characters through the non-reduced BC_6). The `(-,+)` odd-even
  characteristic is genuinely the twisted/non-reduced fingerprint, as the lead predicted.
- NO-MATCH: the specific identity. X is provably not a finite specialization of the A_12^(2) Weyl-Kac /
  Macdonald denominator (R1 factor count 46 vs 36; R2 the denominator has no spinor factor at all, so the
  spin system and the 1^6 deletion have no home there; a numerator is a single alternant not a product).
  The operative wall/period (11/22/44) is the Langlands-partner C_6=Sp(12) **level 4** (4+7=11), not the
  A_12^(2) native scale 13 (which reads level -2). X is the deleted-top-weight **spinor discriminant** -
  a dual-pair / super-denominator object (King spin<->symplectic kin, F133 proof §5), not a Kac-Moody
  affine denominator. The twisted-affine lead correctly identifies the symmetry group and character but
  does not supply the closed identity for the 143 n_lambda.

Files: `_rb_mac1.py` (live D_inf harvest), `_rb_mac2.py` (factor-count + classical-limit + X-vs-B_6),
`_rb_mac3.py` (A_12^(2) Macdonald series==product, 3.6e-14).
