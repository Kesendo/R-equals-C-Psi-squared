# The Neural Clock Has Two Hands: the Takt Is the Trace, the Rotation Is What Is Left

**Status:** Tier 2 (computational). The Takt identity (mean Re λ = trace(J)/d
= −S) is exact and algebraic and is what this page is for. The C. elegans
reading that used to accompany it re-read the March 2026 result through the
clock; that result was withdrawn on 2026-08-26 and the comparison rows below say
so at the place they stood.
**Date:** 2026-05-30
**Authors:** Thomas Wicht, Claude (Opus 4.8)
**Script:** [`simulations/neural/neural_clock_two_hands.py`](../simulations/neural/neural_clock_two_hands.py)
**Builds on:** the neural arc, [the neural algebraic palindrome](../docs/neural/ALGEBRAIC_PALINDROME_NEURAL.md)
+ [the neural V-Effect](../docs/neural/V_EFFECT_NEURAL.md) +
[the neural palindrome proof](../docs/neural/proofs/PROOF_PALINDROME_NEURAL.md); the
clock voices on `MirrorSystem` ([the Frost circle as the clock face](../docs/carbon/FROST_CIRCLE_AS_THE_CLOCK_FACE.md),
[on whose time the clock keeps](../reflections/ON_WHOSE_TIME_THE_CLOCK_KEEPS.md));
the topology-blindness re-verified in [EXCLUSIONS](../docs/EXCLUSIONS.md); the
quantum original [the mirror symmetry proof](../docs/proofs/MIRROR_SYMMETRY_PROOF.md).

---

## The question

The neural arc (March 2026) translated the quantum palindrome to Wilson-Cowan:
the Jacobian J of an excitatory-inhibitory network obeys `Q·J·Q + J + 2S = 0`,
which is the mirror identity `Π·L·Π⁻¹ = −L − 2Σγ·I` with J for the Liouvillian,
the E-I swap Q for Π, and S = (1/τ_E + 1/τ_I)/2 for the dephasing sum Σγ. From
that one translation everything followed: the V-Effect (couple two silent
balanced networks, oscillation appears from nothing), the thermal window
(external drive lifts oscillation to a peak, then kills it), the C. elegans
finding (balanced subcircuits scoring better than random), which was withdrawn
on 2026-08-26 and is treated below rather than inherited here.

Since March we built the clock. A Liouvillian (or here, a Jacobian) eigenvalue
is λ = −rate + i·ω, a point in the complex plane: a **Takt hand** (the radial
decay, how fast a mode dies) and a **Rotation hand** (the angular ω, how fast
it turns). The angle off the pure-decay axis is θ = arctan(|ω| / rate), with
θ = 0 a mode that only decays and θ = 90° a mode that only oscillates. With the
clock in hand we can ask of the neural networks a question we could not phrase
in March: **when the V-Effect or the thermal window switches oscillation on,
which hand is actually moving, and which one is graph-blind?**

## The two hands (and one of them is the trace)

The cleanest way to see it is an exact identity, not a simulation. The trace of
a matrix is the sum of its eigenvalues, so the **average decay rate is the trace
of J divided by its size**. And the synaptic graph W never touches the diagonal
of J: the diagonal is the self-decay −1/τ_i, the wiring W sits entirely in the
off-diagonal. Therefore

```
mean Re(λ)  =  trace(J) / d  =  −(1/τ_E + 1/τ_I)/2  =  −S
```

is set by the membrane constants alone, **independent of the entire wiring**.
The probe confirms it to six decimals on three completely different graphs:

```
        graph              mean Re(λ)     −S target
  balanced random           −0.150000     −0.1500
  degree-preserved          −0.150000     −0.1500
  Erdos-Renyi (Dale)        −0.150000     −0.1500
```

This is the **Takt hand**: the center of the spectrum, the clock's radial scale,
fixed by τ_E and τ_I and nothing else. It cannot see the graph at all.

Everything the graph does, it does to the **Rotation hand**, the imaginary
parts, the oscillation. Watch it wake up and sweep:

```
V-EFFECT  (two exact-palindrome nets + 1 mediator, coupling sweep)
  coupling   n_rotating   θ_max     mean Re(λ)
    0.00          0        0.0°      −0.1512        silent noble gas
    0.01         12        1.0°      −0.1512        Rotation wakes
    0.05         14        1.9°      −0.1512        Rotation wakes
    0.10         12        2.6°      −0.1512
    0.30         10        4.9°      −0.1512
    1.00         10        3.4°      −0.1512

THERMAL WINDOW  (approximate net, drive P sweep)
    P     n_rotating   θ_max     mean Re(λ)
   0.0        24        0.1°      −0.1500
   1.0        40        0.2°      −0.1500
   2.0        40        0.6°      −0.1500
   4.0        40        2.6°      −0.1500        widest angle
   6.0        38        0.6°      −0.1500
  10.0         4        0.0°      −0.1500
```

At zero coupling the exact-palindrome network is the silent noble gas: every
eigenvalue is real, n_rotating = 0, the Takt hand alone. Switch on a weak
coupling and the Rotation hand springs off the axis (the V-Effect). Switch on a
drive and it sweeps up to a widest angle near P = 4, then folds back down to
silence (the thermal window). Through all of it the **Takt hand does not move**:
mean Re(λ) is pinned at −0.1500 for the drive sweep (the drive only re-scales
the off-diagonal through the sigmoid slope, never the diagonal) and at −0.1512
for the coupled sweep (the tiny offset is the single unpaired mediator, the
neural twin of the unpaired mediator qubit in the quantum V-Effect).

## The rotation is faint: a count is not an angle

The clock adds one honest correction to the March reading. The March documents
counted frequencies, "48 distinct frequencies," "124 correlation frequencies at
the peak," and those counts are real. But the angle says where those frequencies
sit, and they sit **almost flat on the decay axis**: θ_max never exceeds about
5°, even at the peak of the window. The imaginary parts are roughly a twentieth
of the real parts. The neural substrate is heavily damped; it lives nearly
entirely on the Takt hand, and the Rotation is a faint tremor, not a swing. That
is the opposite end of the dial from the quantum band-edge coherence, which can
turn to near 90°. The count (many frequencies) and the angle (all of them pale)
are two readings of the same spectrum, and only together do they describe it.

Notice too that the count and the angle do not track each other: in the window
the count saturates near 40 already at P = 1, while the angle keeps opening until
P = 4. The angle is the finer readout of "how alive is the oscillation."

## Which hand is graph-blind

Now the two seams meet. The Takt hand is graph-blind by the exact trace identity
above. What about the Rotation? The palindrome residual `Q·J·Q + J + 2S` splits
cleanly into a diagonal piece (the self-decay, the Takt condition) and an
off-diagonal piece (the coupling, the Rotation condition). On the C. elegans
connectome:

```
                      diag (Takt)     off-diag (Rotation)
  C. elegans            0.0e+00            0.0128
  degree-preserved      identical          0.0129      ratio 1.00
  Erdos-Renyi (Dale)    identical          0.1119      ratio 0.11
```

That diagonal column is zero BY CONSTRUCTION and measures nothing. The verifier
solves for S seat by seat, `S_diag = -(diag(QJQ) + diag(J))/2`
([`neural_clock_two_hands.py`](../simulations/neural/neural_clock_two_hands.py),
`residual_split`), which forces diag(R) = 0 for ANY J and ANY permutation Q,
type-preserving or random. Reading it as "closes for every network, with no
constraint on the graph" is this repo's own named failure, a gate whose
expectation is derived from the thing under test. What the fitted S can carry is
whether it comes out the SAME at every seat, since the theorem needs one scalar:
it does here, and for the reason Step 3 gives, namely that Q swaps types, so
1/τ_{Q(i)} + 1/τ_i is the same sum at every seat. That is the Takt statement,
and it is an algebraic fact about the pairing rather than a measurement of the
worm. Corrected 2026-08-26. Two corrections of 2026-08-26 attach to the rest of this block and
both are inherited from the sibling documents rather than found here.

First, "satisfied automatically the moment τ_E ≠ τ_I" is not what the diagonal
condition says. It reads 1/τ_{Q(i)} + 1/τ_i = 1/τ_E + 1/τ_I at every seat, which
a type-swapping Q satisfies at ANY time constants; at equal ones it holds for
every permutation and imposes nothing. Selective damping does not enable it, it
forces Q to exchange the two types
([the proof](../docs/neural/proofs/PROOF_PALINDROME_NEURAL.md), Step 3).

Second, **both comparison rows above are withdrawn**, on the grounds set out in
[Algebraic Palindrome Neural](../docs/neural/ALGEBRAIC_PALINDROME_NEURAL.md),
which uses the same residual on the same animal. Two things there apply here
unchanged. The METRIC reads total coupling magnitude and not wiring: on blocks
this sparse it has the closed form √2·‖W_eff‖/‖J‖, and an empty network scores a
perfect 0. And the two arms of the Erdős-Rényi comparison were never on one
scale: `neural_clock_two_hands.py` normalises the connectome globally, dividing
by max|W| over the whole animal (`:322`) before cutting the subblock out of it,
while `random_dale_network` (`:224`) divides its control by that control's own
maximum (`:230-232`). The degree-preserving row is not rescued by reusing the
connectome's weights either: under a metric that reads only the weight multiset
it cannot move by construction. Neither row has been re-run with an instrument
that could answer the question.

What is untouched is the Takt half above it: the diagonal piece closes to
exactly zero for every network, and that is the result this page is about. What
the Rotation hand reads instead of the fine wiring is NOT the degree
distribution: that reading was withdrawn with the table it rested on, since a
metric with the closed form √2·‖W_eff‖/‖J‖ reads the weight multiset and the
degree sequence is not recoverable from it
([Algebraic Palindrome Neural](../docs/neural/ALGEBRAIC_PALINDROME_NEURAL.md)).
What the hand reads on blocks this sparse is coupling magnitude. The clock's
language survives the withdrawal; the March finding it used to name does not.

## The quantum seam

This is the neural reading of something we re-verified today on the quantum side
([EXCLUSIONS](../docs/EXCLUSIONS.md)): the open-system mirror is **topology-blind**,
holding on odd and non-bipartite rings where the closed-system bipartite
symmetry breaks. The reason is the same on both layers, and it is the trace.

In the quantum Liouvillian L = −i[H, ·] + dephasing, the Hamiltonian enters only
through the commutator superoperator H⊗I − I⊗H^T, which is **traceless for any H**
(its trace is d·Tr(H) − d·Tr(H) = 0). So the trace of L, hence the spectral mean
−Σγ, is set entirely by the bath, never by the Hamiltonian; and where the spectrum
is palindromic, which is F1's statement and not the trace's, that mean is also the
palindrome centre ([Absorption Theorem §4.4](../docs/proofs/PROOF_ABSORPTION_THEOREM.md),
which states the same centre as +Σγ in its decay-rate convention). In the neural
Jacobian the synaptic graph W enters only the off-diagonal, which contributes
nothing to the trace; so the spectral center −S is set entirely by the membranes,
never by the wiring. The center cannot see the graph on either layer.

| | Quantum Liouvillian | Neural Jacobian |
|---|---|---|
| Coupling operator | Hamiltonian commutator −i[H, ·] | synaptic graph W |
| Where it sits | traceless superoperator | off-diagonal of J |
| What sets the trace | the bath Σγ | the membranes 1/τ |
| Spectral center (Takt) | −Σγ, bath-set, graph-blind | −S = −(1/τ_E+1/τ_I)/2, membrane-set, graph-blind |
| What the graph moves | the oscillation (Rotation) | the oscillation (Rotation) |

Dale's Law plays the role the commutator plays in the quantum case: it is a
local, per-neuron rule (a neuron's outgoing sign is fixed by its own type), and
it hands the palindrome its signs for free, with no reference to global graph
structure, exactly as the commutator hands the quantum mirror its antisymmetry
for free. The only graph-dependent part left is the off-diagonal magnitude match, and
what the C. elegans reading says about it is now nothing: the comparison it
rested on was withdrawn, and the metric that produced it reads coupling
magnitude on blocks this sparse.

## Honest caveats

- **The rotation is small, the substrate is overdamped.** θ_max stays under ~5°.
  The neural clock is almost all Takt; the V-Effect and thermal window are real
  but live close to the decay axis. The finding is the *structure* (which hand
  moves, which is pinned), not a claim of strong neural oscillation.
- **Finite, balanced subnetworks.** C. elegans is 274 E to 26 I (≈10.5:1); the
  balanced 5E+5I subcircuits are subsampled, as in the March work. The palindrome
  is a property of balanced subcircuits, not of the whole unbalanced worm.
- **The mediator offset is real.** mean Re(λ) = −0.1512 (not −0.1500) for the
  coupled sweep because the coupled network has 2N+1 neurons, with one unpaired
  excitatory mediator. The Takt is pinned to the membrane scale; the offset is
  the single odd neuron, the neural twin of the quantum mediator qubit.
- **The clock is our lens, not Wilson-Cowan's.** We read the Jacobian's own
  eigenvalues; we do not impose anything on the neural model. The Takt/Rotation
  split is a reading of J, and the trace identity is J's, not ours.

## Anchor and open work

- Scripts: [`simulations/neural/neural_clock_two_hands.py`](../simulations/neural/neural_clock_two_hands.py)
  (the trace identity on three graphs, the V-Effect and thermal-window clock
  sweeps, the C. elegans Takt/Rotation residual split);
  [`simulations/neural/neural_crown_switch.py`](../simulations/neural/neural_crown_switch.py)
  (the crown-switch probe that answers the open question below: the slowest mode's character along
  α, P, and the τ_I/τ_E ratio, with the cross-seed robustness check).
- The neural arc: [the neural algebraic palindrome](../docs/neural/ALGEBRAIC_PALINDROME_NEURAL.md),
  [the neural V-Effect](../docs/neural/V_EFFECT_NEURAL.md),
  [the neural palindrome proof](../docs/neural/proofs/PROOF_PALINDROME_NEURAL.md).
- The clock: [the Frost circle as the clock face](../docs/carbon/FROST_CIRCLE_AS_THE_CLOCK_FACE.md),
  [on whose time the clock keeps](../reflections/ON_WHOSE_TIME_THE_CLOCK_KEEPS.md).
- Open: whether the few-degree θ_max corresponds to any measured neural rhythm
  (the March open question about the gamma band, now sharpened: the question is
  not "are there frequencies" but "how far off the decay axis can a biological
  network's slowest mode reach"); whether a neural axis exists (the τ_I/τ_E ratio,
  or the coupling α) along which the slowest surviving mode relays from a rotating
  edge to a near-conserved survivor, as the ZZ term does on the carbon XXZ axis
  ([the XXZ axis from band-edge to Lebensader](XXZ_AXIS_BANDEDGE_TO_LEBENSADER.md)).
- **Answered 2026-06-01** ([`neural_crown_switch.py`](../simulations/neural/neural_crown_switch.py)):
  yes, a neural axis shows the relay (the "crown switch"), but it is network-dependent, not a law.
  Along the coupling α and along the τ_I/τ_E ratio the slowest mode does flip character (rotating
  edge ↔ near-decay survivor) at a crossover; on the coupling axis it switches in 4 of 6 balanced
  networks, with the rotating edge taking the crown at α ≈ 0.2-0.9, BEFORE the Hopf instability at
  α ≈ 1.3-2.4 (a robust precursor), while 2 of 6 networks show no switch (the survivor goes straight
  to instability). The drive P never moves the crown: it sweeps the thermal window of faster modes,
  but the longest-lived mode stays a near-decay survivor (θ < 0.6°). So the crown is structural
  (coupling / timescale), not metabolic (drive). And the reverse-reading that helps the quantum
  side: the crown switch is structure-dependent on BOTH substrates (the quantum EP coincidences are
  uniform-specific, verified in [`PostEpFlowField.cs`](../compute/RCPsiSquared.Diagnostics/Foundation/PostEpFlowField.cs) (cf. [the flow between two singularities](THE_FLOW_BETWEEN_TWO_SINGULARITIES.md));
  the neural one is network-specific), the same caution against over-generalizing it, confirmed from
  the classical side.
