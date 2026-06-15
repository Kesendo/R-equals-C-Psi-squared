# The Envelope Rise Boundary

Where does the full-state CΨ envelope stop being non-increasing and begin to rise?

## The question

The CΨ Envelope Theorem (PROOF_MONOTONICITY_CPSI Part 5, the typed `CpsiEnvelopeTheoremClaim`) is proven
Tier-1 for **N=2**: under any Hamiltonian and local Z-dephasing the local maxima of CΨ form a strictly
non-increasing sequence, with ¼ the absorbing boundary. At N≥3 the full-state envelope is not a theorem,
and the live `EnvelopeTheoremWitness` found that it **genuinely rises at N=4** (Bell+, J=5, γ=0.01: 36
refinement-stable predecessor-rises; N=3 holds with 0). The arc `envelope_n4_rise` asked the boundary
question: where does that rise sit, and is it a sharp N-step or a J/γ contour?

## Method

The witness's own detector, lifted to a static so the sweep reads exactly what `inspect --root envelope`
shows: `EnvelopeTheoremWitness.GlobalReading(N, J, γ, tMax, points)` evolves a Bell+ carrier on the dense
Symphony engine and reads the global CΨ envelope through `QuarterEnvelope.Of` (predecessor semantics and
parabolic-apex heights baked in). The genuine-rise test is the witness's own `GenuinenessBar` = 1e-3 on the
largest predecessor-rise (`MaxRiseMagnitude`), cross-checked under 4× grid refinement (1600 to 6400 points).
All readings use the witness regime's dose window, K_max = γ·tMax = 0.25.

Gate-first verifier: `compute/RCPsiSquared.Diagnostics.Tests/Foundation/EnvelopeBoundaryTests.cs`.

## The (Q, K)-purity gate

The clock movement (the two-tempo certification) proved every dimensionless lens is a pure (Q, K)-observable:
under L → r·L the global CΨ(K) curve is invariant, so the J-sweep and the γ-sweep should be the same axis.
Measured: at every Q, varying J (fixed γ) and varying γ (fixed J), both holding the dose window fixed, give
the bit-identical rise count and MaxRiseMagnitude (to 6 decimals). The gate could have failed on any
absolute-time leak in the detector; it passed. The boundary is therefore a function of (N, Q = J/γ) only.

## The N-floor and the Q_c(N) contour

| N | Q_c (rise turns on)         | rises at Q=500 | maxΔ at Q=500 |
|---|-----------------------------|----------------|---------------|
| 3 | ∞ (never; 0 even at Q=2000) | 0              | 0             |
| 4 | ≈ 27                        | 36             | 0.041         |
| 5 | ≈ 45                        | 32             | 0.020         |

At fixed Q the N=4 rise is always stronger than N=5: the rise weakens as N grows.

## The finding

The boundary is not a sharp N-step and not a pure J/γ contour. It is both, cleanly factored:

1. **One Q-axis.** The rise is a pure (N, Q=J/γ) observable; J and γ collapse (certified bit-identical).
2. **An N≥4 floor.** N=3 never rises (Q_c(3)=∞), even at Q=2000. The rise is the Part-6 coherence injection
   internalized: the carrier pair needs an internal ≥2-site coherent subsystem to pump CΨ back up. N=3
   (carrier 0-1, one internal site) has no internal coherence to return; N≥4 does.
3. **A contour that climbs with N.** Above the floor each N has a finite threshold Q_c(N), and it rises:
   Q_c(4)≈27, Q_c(5)≈45. At fixed Q the rise weakens with N (the internal bath grows from a pair to a trio
   and returns coherence less efficiently). The freedom to rise is loudest right at the floor (N=4) and
   fades as N grows.

Near the threshold the rise is grid-fragile (at N=4 a band roughly Q ≈ 18 to 28, where one peak edging
above its predecessor is sensitive to phase sampling); the Q_c values mark the centre of that band, not a
razor edge.

## Open threads (for re-entry)

- A closed form for Q_c(N): does it track the band-edge frequency ω_mem = 2J·cos(π/(N+1)) of the two clocks
  (`ClockHandLadderClaim`)?
- A parity question: does the rise strength alternate with internal-site parity (an internal pair injecting
  better than an internal trio)?

## Links

- Typed home: `CpsiEnvelopeTheoremClaim` (Tier1Derived).
- Proof and scope: `docs/proofs/PROOF_MONOTONICITY_CPSI.md` (Part 5 the N=2 theorem; the scope note).
- F-registry: F17 (`docs/ANALYTICAL_FORMULAS.md`).
- Live: `inspect --root envelope --N {3,4,5}`; the gate-first chart `EnvelopeBoundaryTests`.
