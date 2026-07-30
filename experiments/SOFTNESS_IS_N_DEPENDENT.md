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
(self-validating, three sections: the per-witness residual ‖M‖ and maximum spectral pairing-error at
N=4,5,6; the twelve-set census of the lit k=3 family; and the uniform-router check that builds
Π = P4^⊗N and conjugates L, with the hard control run through it as the discriminator).

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
open chain of N sites at coupling J=1 under uniform Z-dephasing γ=0.05. Every ‖M‖ below is exactly
linear in J and γ-independent; ‖L‖ and the pairing errors carry the dissipator too and move with
both, so every number here belongs to that (J, γ) pair. We read three numbers per (witness, N):
the
palindrome residual ‖M‖ (which vanishes only for the exactly-soluble "truly" case), and the maximum
spectral pairing-error, namely how far the nearest partner of each eigenvalue λ sits from its mirror
target −λ−2Σγ. A pairing-error at machine precision is soft; an O(1) error is hard. Reading the
magnitude, not just the class label, is what lets us tell a real hard verdict from a near-degenerate
numerical wobble; it is the whole reason this crossing can be called genuine.

## The crossing

| witness | N=5 | N=6 |
|---|---|---|
| XXXX + XYYY + YYYX | pairErr 8.7×10⁻¹⁴ → **soft** | pairErr 2.0×10⁻¹ → **hard** |

The crossing is clean. At N=5 the spectrum pairs to 8.7×10⁻¹⁴, indistinguishable from a certified-soft
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
precision floor; it never drifts toward a break. The certified 2-body cases are N-stably soft. That
is the structural reading of the certifier's 2-body gate in
[`PalindromeSoftCertifier`](../compute/RCPsiSquared.Diagnostics/F87/PalindromeSoftCertifier.cs),
which holds its site-swap strategy to exactly two non-identity letters. (A separate gate requires
those two to be adjacent; each rejects a different false positive.)

The gate is not, however, the place where N-stable softness ends for the lit family. Sweeping all
twelve three-term reversal-symmetric fully-lit (no-Z) k=3 sets gives 8 hard, 2 truly and 2 soft, the
same split at N=4 and N=5, and both soft ones survive the longer chain: XXY+XYX+YXX pairs at 4.5×10⁻¹⁴,
1.6×10⁻¹³ and 3.0×10⁻¹³ for N=4, 5, 6, and XXY+YXX+YYY at 5.7×10⁻¹⁴, 1.1×10⁻¹³ and 4.2×10⁻¹³, with
residuals ‖M‖ = 78.4, 192.0 and 443.4, so they are genuinely soft and not truly. The residual is
the same column for both sets to every digit; the pairing errors differ in the last figures, which
is the eigensolver's floor on two different matrices and carries no information. A fully-lit k≥3
soft verdict therefore CAN be finite-size grace, and XXXX+XYYY+YYYX above shows that it sometimes
is, but it is not grace by construction.

And for these two it is not measurement either: both are routed by the UNIFORM per-site map P4
(I↔Y, X↔Z with phases i on X and Z), one of the four representatives the k-body routing certifier
tries ([`KBodyPalindromeRouting`](../compute/RCPsiSquared.Diagnostics/F87/KBodyPalindromeRouting.cs)).
Building Π = P4^⊗N directly gives ‖Π L Π⁻¹ + L + 2Σγ‖ = 6.4×10⁻¹⁶ at N=4 and 1.2×10⁻¹⁵ at N=5
against ‖L‖ = 55.5 and 136.0, while the same construction misses by 45.2 and 110.9 on the hard
control, so the check discriminates. M2, the same permutation with conjugate phases, routes them too;
P4 is simply the first the search reaches. That certificate is constructive and additive over windows,
hence N-independent, and it settles both halves of the verdict at every N ≥ 3: the palindrome because
Π = P4^⊗N exhibits it, and non-truly because all three templates are Π²-odd, so H_odd = H and F83
gives ‖M‖² = 4·2^N·‖H‖² > 0, which at N = 4, 5, 6 is 6144, 36864 and 196608, i.e. the ‖M‖ column
above before rounding (768 at N = 3, the shortest chain a k=3 template fits). Soft at every
length, not merely out to the length we can diagonalise.

A *router* here is a per-site product Π = ⊗_l Q_l built by repeating a pattern of the certifier's
per-site maps with some period, and it counts as a router for a Hamiltonian when
Π L Π⁻¹ = −L − 2Σγ holds outright on the full Liouvillian. The census searches the three
signed-permutation representatives {P1, P4, M2} at period ≤ 2, nine distinct patterns; the
certifier's own set adds a fourth, the dense crossover map M with M² = −I, which is not a signed
permutation. The golden router is a different object again: period-4, reached only by the
window-summed condition, outside this bounded family entirely. A set with no router in
that bounded family is not thereby hard, since the certifier is one-sided; the search says only
what it found.

Across the twelve, a router of period ≤ 2 exists for exactly the four non-hard sets and for none of
the eight hard ones.

The finer split inside those four is a convention, not a mechanism. The global X↔Y rotation sends
Z → −Z and leaves Z-dephasing invariant, and it maps each soft set onto one of the truly ones:
XXY+XYX+YXX ↔ XYY+YXY+YYX and XXY+YXX+YYY ↔ XXX+XYY+YYX. Each such pair therefore shares one
Liouvillian spectrum. All six pairings among the four agree at the eigensolver's floor, order 10⁻¹⁴
at N=4, so what they share is one spectrum and not two related pairs; a hard set from the same family
sits 2.0 away under an optimal assignment of one spectrum onto the other, which is what makes the
agreement mean something (sorting both and comparing in order reports 15.9 instead, so the separation
is real but its size is the matcher's). Within each pair one member has M = 0 against the fixed Π_Z
and reads "truly" while the other does not and reads "soft"; P1 stands to P4 as the X↔Y image up to
the same phase conjugation that separates M2 from P4. What routability cuts is 4 against 8.
Which side of the 4 a set lands on is the choice of mirror.

## Where the genuine higher-body soft cases live

Not every k≥3 soft case is a finite-size accident. XZX+XZY+YZX pairs to machine precision at N=4, 5,
and 6 alike, genuinely and stably soft, and here the stability is not only measured but explained.
The dark channel does route a hidden Q,
but not the one this note first guessed (corrected 2026-06-10): the Z-middle case routes no Klein-family
Q, neither the uniform nor the alternating families of the
[Klein routing](TWO_TERM_PALINDROME_KLEIN_ROUTING.md) reach it (the discrete candidates sit off the
golden locus). Its actual router is the period-4 golden product in the frame a = φX+Y, b = X−φY
([the ceiling golden-router proof](../docs/proofs/PROOF_CEILING_GOLDEN_ROUTER.md)), a
soft mechanism distinct from the lit chiral-K colouring and exact at every N ≥ 3, a single window already
suffices, which is exactly the length-independence seen here. So the k≥3 soft cases divide by whether
a router is known, not by whether a Z is present, and on that cut the Z-middle case is the HARDER of
the two: nothing in the three-representative period-≤2 family the census searches routes it, and the period-4 golden product
that does route it does so only window-summed. Its obstruction is one template of the three: XZX alone
is routed by nothing in that set, while XZY and YZX each fall to a period-2 pattern, so the certificate
cannot be assembled template by template. The two lit sets, by contrast, fall to the plainest thing
there is, one uniform map at every site. What remains genuinely finite-size is the third kind,
XXXX+XYYY+YYYX, which has a demonstrated length at which it breaks and so can carry no window-local
certificate at all. At N=4 and N=5, within these twelve and none of them carrying a Z at all, it is routability that
sorts the verdicts, and the Z-middle set outside them shows the same cut running across the presence
of a Z. The chain length has not dropped out of the story: at N=3 the census reads 10 soft, 2 truly
and nothing hard, so eight of the twelve are themselves crossings of the kind this note is about,
soft at N=3 and hard at N=4. Routability is N-free and the verdicts are not, which is why the four
it certifies are the four that never move. The certifier is one-sided in general, so a set with no router in a bounded
candidate set is not thereby hard.

## The lesson

A spectral palindrome that holds on a short chain can shatter on a longer one. For the exactly-soluble
base case (truly, ‖M‖=0) this never happens; the residual is identically zero at every N. But the
soft/hard line, the subtler distinction inside the non-truly Hamiltonians, is finite-size sensitive
for k≥3: the trichotomy verdict is a function of (Hamiltonian, N), not of the Hamiltonian alone. Any
"soft" read off a small chain deserves an N-stability check before it is trusted. (The trichotomy that
F85 calls N-stable is a different cut of the same word: the Π²-class truly / Π²-odd / Π²-even-non-truly
is computed from the letter parities alone, so it cannot depend on N. What moves with N is the spectral
soft/hard line drawn inside the non-truly class.) The certifier earns its certificate precisely by
certifying only the mechanisms (the chiral K, the hidden Q) that hold at every length.

The mirror that holds at one length can break at the next. Softness, for some higher-body terms, is a
grace of being small; for others it is a router, and the certificate is exactly the difference between
the two. Where a router is exhibited the length stops mattering; where none is, the length is the first
thing to vary.
