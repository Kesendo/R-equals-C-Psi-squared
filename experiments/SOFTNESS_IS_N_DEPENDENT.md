# Softness Is N-Dependent: a Finite-Size Crossing in the Palindrome Verdict

**Status:** Computed (Tier 2); the N=5 → N=6 crossing is bit-exact, a genuine soft verdict turning
into a genuine hard one, not a tolerance artifact. 2026-06-06.
**Authors:** Thomas Wicht, Claude (Opus 4.8)
**Builds on:** [Two-Term Palindrome: the Klein Routing](TWO_TERM_PALINDROME_KLEIN_ROUTING.md) (whose
N-invariance this qualifies), the soft-certifier
[`PalindromeSoftCertifier`](../compute/RCPsiSquared.Diagnostics/F87/PalindromeSoftCertifier.cs) (whose
2-body soundness this confirms one chain length past where it was checked), and
[F103 §7](../docs/proofs/PROOF_F103_F87_Z2_CUBED_REFINEMENT.md) (the two soft mechanisms: the
diagonal chiral K and the hidden-Q routing).
**Scripts:** [`simulations/f87_softness_n_dependence.py`](../simulations/f87_softness_n_dependence.py)
(self-validating: builds L, reads the residual ‖M‖ and the maximum spectral pairing-error at N=4,5,6).

---

## The question

We had been speaking of the F87 verdict as though "soft" and "hard" were stamped on the Hamiltonian
itself. The [two-term routing note](TWO_TERM_PALINDROME_KLEIN_ROUTING.md) says so out loud: its fate
table is bit-identical at N=3,4,5, so "N=3 is not a small-N accident; it is representative." For
two-term bond bilinears that is true. The question here is whether it survives once the terms grow
past two bodies, and the answer is no. A Hamiltonian can be genuinely soft on a chain of one length
and genuinely hard on the next.

The witness turned up while hunting a different quarry. The certifier's site-swap strategy certifies
the 2-body set XX+XY+YX as soft; we wanted to know what separates it from the 3-body lookalike
XXX+XXY+YXX, which is hard despite carrying the same symmetry signature (reversal-symmetric,
mask-bipartite, bit_b-mixed). The tidy guess, that the parity of the body-count decides it (even soft,
odd hard), is wrong, and the way it is wrong is the finding.

## The setup

Each witness is a sum of reversal-symmetric Pauli templates placed by the sliding-window builder on an
open chain of N sites under uniform Z-dephasing γ=0.05. We read three numbers per (witness, N): the
palindrome residual ‖M‖ (which vanishes only for the exactly-soluble "truly" case), and the maximum
spectral pairing-error, namely how far the nearest partner of each eigenvalue λ sits from its mirror
target −λ−2Σγ. A pairing-error at machine precision is soft; an O(1) error is hard. Reading the
magnitude, not just the class label, is what lets us tell a real hard verdict from a near-degenerate
numerical wobble; it is the whole reason this crossing can be called genuine.

## The crossing

| witness | N=5 | N=6 |
|---|---|---|
| XXXX + XYYY + YYYX | pairErr 8.7×10⁻¹⁴ → **soft** | pairErr 2.0×10⁻¹ → **hard** |

The crossing is clean. At N=5 the spectrum pairs to 9×10⁻¹⁴, indistinguishable from a certified-soft
case; at N=6 it misses by 0.2, square in the hard band that the control XXX+XXY+YXX occupies at every
length (≈0.02 at N=4, ≈0.24 at N=5). Nothing borderline sits between: the same template is exactly
palindromic on five sites and decisively broken on six. The finite-size window that carried the
softness closes when the chain lengthens. At N=5 the 4-body template slides over only two positions;
at N=6 a third window appears, and the grace is gone.

So the discriminator we set out to find, a rule on the X/Y phase pattern that would call
XXXX+XYYY+YYYX soft and XXX+XXY+YXX hard, does not exist as a property of the Hamiltonian, for a simple
reason: the verdict it would have to predict is not a property of the Hamiltonian. It is a property of
the pair (Hamiltonian, N).

## The certifier holds

This is the moment to check something load-bearing. The certifier's site-swap strategy had been
verified only to N=5, and we just watched N=5 tell a soft story that N=6 overturns. Do its certified
cases survive the longer chain?

| certified witness | N=4 | N=5 | N=6 |
|---|---|---|---|
| XX + XY + YX | 6.8×10⁻¹⁴ | 1.8×10⁻¹³ | 3.3×10⁻¹³ |

(YY+XY+YX, its X↔Y mirror, matches at N=4,5: 1.0×10⁻¹³, 1.7×10⁻¹³.)

They survive, at machine precision, with a margin (N−k = 4) far past the danger zone where the
finite-size softness lived (N−k = 1). The pairing-error grows only with the matrix size, tracking the
precision floor; it never drifts toward a break. The certified 2-body cases are N-stably soft, and
that is the structural reading of the 2-body gate. Among the fully-lit (no-Z) reversal-symmetric sets
we swept (k=2,3,4), only k=2 is soft at every length; the k≥3 fully-lit lookalikes go hard once the
chain is long enough, their occasional small-N softness being exactly the finite-size grace that
washes out. The gate is not a blunt stand-in for a missing phase-word. It sits on the structural
boundary where N-stable softness ends for the lit family.

## Where the genuine higher-body soft cases live

Not every k≥3 soft case is a finite-size accident. XZX+XZY+YZX pairs to machine precision at N=4, 5,
and 6 alike, genuinely and stably soft. The difference is the Z. The dark channel does route a hidden Q,
but not the one this note first guessed (corrected 2026-06-10): the Z-middle case routes no Klein-family
Q, neither the uniform nor the alternating families of the
[Klein routing](TWO_TERM_PALINDROME_KLEIN_ROUTING.md) reach it (the discrete candidates sit off the
golden locus). Its actual router is the period-4 golden product in the frame a = φX+Y, b = X−φY
([the ceiling golden-router proof](../docs/proofs/PROOF_CEILING_GOLDEN_ROUTER.md)), a
soft mechanism distinct from the lit chiral-K colouring and exact at every N ≥ 3, a single window already
suffices, which is exactly the length-independence seen here. So the soft k≥3 world splits in two: the
fully-lit cases that are soft only by the grace of a short chain, and the Z-routed cases that are soft
for a real, length-independent reason. The first is a mirage; the second is the honest frontier for
extending the certifier.

## The lesson

A spectral palindrome that holds on a short chain can shatter on a longer one. For the exactly-soluble
base case (truly, ‖M‖=0) this never happens; the residual is identically zero at every N. But the
soft/hard line, the subtler distinction inside the non-truly Hamiltonians, is finite-size sensitive
for k≥3: the trichotomy verdict is a function of (Hamiltonian, N), not of the Hamiltonian alone. Any
"soft" read off a small chain deserves an N-stability check before it is trusted. The certifier earns
its certificate precisely by certifying only the mechanisms (the chiral K, the hidden Q) that hold at
every length, and the 2-body gate marks where that stability ends for the lit terms.

The mirror that holds at one length can break at the next. Softness, for the higher-body terms, is
partly a grace of being small.
