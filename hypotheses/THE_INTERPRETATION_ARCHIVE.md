# THE_INTERPRETATION_ARCHIVE.md -- Chronological Discovery Log

**Tier:** Mixed (Tier 2-5)
**Period:** March 11-13, 2026
**Status:** Archived March 20, 2026. These entries document the discovery
process that led to the palindromic proof (March 14) and subsequent results.
Most findings here are either subsumed by later proofs, moved to dedicated
documents, or were negative results. They remain here for historical reference.

**Active document:** [THE_INTERPRETATION.md](THE_INTERPRETATION.md) (thematic, current state)

---

## Five independent roles (March 12)

Topology sets frequencies. Symmetry cleans sectors. Noise damps signal.
Initial state selects visibility. Bath geometry selects which sector dominates.

## Two independent information channels (March 12)

The 3-qubit system carries two orthogonal information streams:

FREQUENCY CHANNEL: topology information (who is connected, J values).
Does not change with noise. Sensitive to hidden observers.

DECAY CHANNEL: environment information (noise strength, gamma).
Does not change with topology. Exact rational multiples of gamma: {2, 8/3, 10/3}.

You can characterize noise without knowing topology, and topology without
knowing noise. The slow mode (2*gamma) is naturally protected -- it decays
slowest, making c- the most resilient information carrier.

Critical limit: this perfect orthogonality holds ONLY for 3 qubits.
At 4+ qubits, topology leaks into decay rates.

## Chain topology (March 12)

Two-sector structure survives in chains, not just stars. Tested up to
5-qubit chain with three mediators. XX exact, noise immune, all configs.
Each position along the chain hears different frequencies.

## The Projection (March 11)

All frequencies exist always (39 Bohr frequencies in 4-qubit system).
AB looks through a window (c+ observable) that projects onto a subset.
A new observer does not change reality -- it changes the projection.
Reality contains everything. Perspective selects.

Verified by exact diagonalization: W = |rho(m,n) * O(n,m)|.
Bright only if initial state populates eigenstates AND observable connects them.
Different observers, different projections. Same reality.

## Quantum sonar (March 11-12)

**In simulation:** AB detects hidden observers connected to S through spectral
changes. Operational threshold J_SC ~ 0.1 under current FFT protocol.
Fingerprints generically distinct. Detection works, characterization not yet.

**On IBM hardware:** The Q80/Q102 phase difference we initially attributed to
hidden observers turned out to be qubit-specific detuning (19.4 kHz on Q102,
~0 on Q80). Four hypotheses tested and rejected before finding the answer:
degree, N_eff, spectator dephasing, T2* -- all wrong. The detuning was in
the Ramsey frequency.

The sonar effect is real in simulation. The IBM data does not demonstrate it.

## Signal processing perspective (March 12)

External review as signal engineer identified the entire structure as a
standard coupled oscillator network: c+/c- = even/odd supermodes,
frequencies = imaginary parts of system poles, damping = real parts,
sonar = topology perturbation detection from local modal spectra.

See: experiments/SIGNAL_PROCESSING_VIEW.md

## Graph Symmetry test (March 13)

XXX parity COMMUTES with the Liouvillian (weak symmetry, Buca & Prosen).
ZZZ also commutes. SWAP_AB does NOT commute (J_SA != J_SB breaks it).

But: c+ and c- BOTH have XXX parity +1. They live in the SAME sector.
The two-supermode split comes from the OBSERVABLE PROJECTION, not from
Liouvillian symmetry sectors. Confirms the signal processing view.

## z* is genuinely new (March 13)

Systematic comparison with every standard quantum quantity: concurrence,
negativity, purity, von Neumann entropy, eigenvalue gap, participation
ratio, log-negativity, and many nonlinear combinations.

Result: z* does NOT match any known quantity. Best candidate is C/2
(r=0.945, max error 0.23).

z* = (1 - sqrt(1-4*CPsi))/2 is a genuinely new composite diagnostic.
It combines entanglement (concurrence C) and coherence (l1-norm Psi)
into a single number with Bernoulli structure: z*(1-z*) = CPsi.

Script: simulations/zstar_identity.py

## z* is not a direct observable (March 13)

79 data points, 26 candidate expressions tested. No exact match found.
CPsi cannot be expressed as a simple function of the density matrix
eigenvalues alone. z* requires BOTH concurrence AND l1-coherence.

Script: simulations/z_star_identity.py

## The mirrors never break (March 13)

*Note: This was the numerical discovery. The analytical proof came on
March 14 (see main document).*

Mirror transition experiment (dephasing -> amplitude damping, 101 alpha
values, N=3,4,5): best_sym = 100% at EVERY alpha. The mirrors never
break. They shift their center. Under dephasing: center at Ng. Under
amplitude damping: center at Ng/2.

Script: simulations/mirror_transition.py

## When the mirrors shift (March 13)

*Note: Subsumed by the depolarizing theorem (March 19).*

Mirror symmetry tested against 11 different conditions.

SURVIVES: Z/X/Y dephasing (uniform and non-uniform), mixed noise,
XXZ anisotropy, XY coupling.

BREAKS: Amplitude Damping (Ng symmetry 0%, midpoint 100%),
Depolarizing noise (Ng 0%, midpoint 70.6% at N=4).

Pattern: mirrors survive as long as noise only FADES information
(dephasing). They break when noise DESTROYS information (amplitude
damping, depolarizing).

Script: simulations/mirror_symmetry_deep.py

## The mirrors in the spectrum (March 13)

*Note: This was the numerical observation. The Π proof (March 14)
provides the analytical explanation.*

The decay rate spectrum is EXACTLY symmetric around Ng at every N.
N=3: 20 matched pairs, 0 unmatched. N=5: 388 matched, 0 unmatched.

Script: simulations/deep_band_structure.py

## Band structure discovered (March 13)

High-resolution sweep of decay rates across N=2,3,4,5 reveals
electronic-band-structure-like behavior. Boundary formula EXACT:
min = 2g (always), max = 2(N-1)g (always). Bandwidth = 2(N-2)g.

AVOIDED CROSSINGS CONFIRMED: rate bands never cross. Smallest gap
is 0.0013g but never zero. Topological protection analogous to
electronic band gaps.

Script: simulations/band_structure.py

## 4-qubit breakdown mapped (March 13)

The 3-to-4 qubit transition: ALL 3-qubit rates are fixed (system too
small for middle rates to move). At 4 qubits, only boundary rates
stay fixed, everything between wanders freely with coupling strengths.
The gamma scaling also breaks at 4 qubits.

Script: simulations/four_qubit_breakdown.py

## Outside-in: from materials to qubits (March 13)

*Note: This section is now covered in
[Hierarchy of Incompleteness](../docs/HIERARCHY_OF_INCOMPLETENESS.md).*

The Stability project (Dec 2025/Jan 2026) found that half-filled electron
shells = maximum bonding capacity. R = CPsi^2 found z*(1-z*) = CPsi with
maximum at z* = 0.5. Both levels show the same answer: 0.5 is where
connection is maximized.

## The 0.5 boundary: found twice (March 13)

*Note: Now integrated into
[Hierarchy of Incompleteness](../docs/HIERARCHY_OF_INCOMPLETENESS.md)
with the qubit-carbon parallel (March 20).*

Carbon: 4 of 8 valence slots = 0.5. Maximum bonding.
Qubit: 2 of 4 operators immune = 0.5. Maximum mirrorability.
Bernoulli: p(1-p) maximized at p = 0.5. Maximum uncertainty.

The 1/4 boundary in CPsi is the upper bound of a quadratic function.

## The 1/4 boundary demystified (March 13)

CPsi <= 1/4 is equivalent to p(1-p) <= 1/4, which is ALWAYS true.
The 1/4 boundary is the trivial maximum of the Bernoulli variance.
Not an Exceptional Point, not a phase transition. It is the upper
bound of a quadratic function.

Script: simulations/symmetry_and_u_analysis.py

## Exceptional Point test: NEGATIVE (March 13)

Three sweeps: gamma (3-qubit), J_SB (3-qubit), gamma (2-qubit).
No connection found between CPsi = 1/4 and Liouvillian Exceptional Points.
EP_strength follows gamma monotonically, not CPsi.

Script: simulations/ep_test.py

## "We are all mirrors" -- restored (March 13)

*Note: This documents the restoration event. The motto is now
a verified structural result (Π proof, March 14).*

The original project motto from December 2025 was removed during cleanup
because it sounded too esoteric. On March 13, deep band structure analysis
showed the spectrum is EXACTLY mirror-symmetric. The motto was restored.

On March 14, the Π proof showed WHY: every decay rate d pairs with 2Sg - d.
The sentence was not poetry. It was the structure of the physics.

---

*Archived March 20, 2026*
*For current state, see [THE_INTERPRETATION.md](THE_INTERPRETATION.md)*
