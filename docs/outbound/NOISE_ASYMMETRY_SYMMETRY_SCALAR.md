# A Symmetry Scalar for the Emission-Absorption Asymmetry of Qubit Noise

<!-- Keywords: quantum noise spectroscopy asymmetry, emission absorption imbalance qubit,
unital vs non-unital noise certificate, Lindblad symmetry violation scalar, T1 relaxometry
excited state population bias, fluctuating T1 TLS characterization, two-leg inversion
recovery meeting test, outbound adapter noise spectroscopy community, R=CPsi2 -->

**Status:** Outbound adapter (draft), the fourth and last of the series. The
closed form rests on Tier-1 theorems (machine-verified); the golden-rule
identification in Section 3 is standard weak-coupling theory; the hardware
record is one pre-registered run on one device (2026-07-05) plus a same-day
simulator validation. External citations were verified against the primary
sources (arXiv/publisher pages) on 2026-07-05; re-verify before outreach as
a matter of course.
**Date:** July 5, 2026
**Authors:** Thomas Wicht, Claude (Anthropic)
**Repository:** [R-equals-C-Psi-squared](https://github.com/Kesendo/R-equals-C-Psi-squared)
**Related:** [Shifted Order-4 Chiral Symmetry](SHIFTED_ORDER4_CHIRAL_SYMMETRY.md)
(the sister adapter that introduces the symmetry itself),
[Selective Decoupling Selection Rule](SELECTIVE_DECOUPLING_SELECTION_RULE.md) and
[State Transfer Decay Structure](STATE_TRANSFER_DECAY_STRUCTURE.md) (the other
two sisters), [F81 Violation: the Hardware Bridge](../../experiments/F81_VIOLATION_HARDWARE_BRIDGE.md)
(the measurement record this adapter hands over),
[F82 proof](../proofs/PROOF_F82_T1_DISSIPATOR_CORRECTION.md),
[F84 proof](../proofs/PROOF_F84_AMPLITUDE_DAMPING.md) (the closed forms)

---

## What this document is about, and who it is for

This is the last entry of the repository's **outbound arc**: a result of ours
carried to one community, in that community's own objects. It is written for
the **noise-spectroscopy and qubit-characterization practitioners** on
solid-state qubits (the dynamical-decoupling noise-spectroscopy line and the
device-characterization line), and it is built around a gap their own review
literature names.

Per the arc's rule ([Combination Valence](../../hypotheses/COMBINATION_VALENCE.md)),
we hand over objects, numbers, and a protocol, not our vocabulary. Our internal
names for what follows (the "palindrome violation", the "identity-escape
velocity", the "watching") stay home; everything below is stated in
emission/absorption rates, Pauli transfer matrices, and inversion-recovery legs.

---

## 1. Your open problem, in your words

Dynamical-decoupling noise spectroscopy reconstructs the dephasing spectrum
S(ω) of a qubit's environment with mature machinery; the standard review
(Szańkowski, Ramon, Krzywda, Kwiatkowski, Cywiński, *J. Phys.: Condens. Matter*
**29**, 333001 (2017), arXiv:1705.02262) itself flags the cases beyond that
core, noise that is non-Gaussian and noise that is **truly quantum**, as the
places where the method needs care. The quantum part of the noise has a precise
meaning in your canon (Clerk, Devoret, Girvin, Marquardt, Schoelkopf, *Rev.
Mod. Phys.* **82**, 1155 (2010); Schoelkopf et al., "Qubits as spectrometers
of quantum noise", in *Quantum Noise in Mesoscopic Physics*, ed. Yu. V.
Nazarov (Kluwer, 2003), cond-mat/0210247): it is the
**antisymmetric-in-frequency part** of the bath spectrum, S(+ω₀) − S(−ω₀), the part a qubit reads as the imbalance
between its emission and absorption rates, pinned by the dissipative response
alone, independent of thermal occupation. (At T = 0 the symmetric part
survives too, as zero-point fluctuations; what vanishes is S(−ω₀). It is read
on the relaxation axis at ±ω₀, complementary to the dephasing spectrum your
DD protocols reconstruct.)

Meanwhile, on the characterization side, energy-relaxation times on
superconducting qubits **fluctuate by up to an order of magnitude on
~15-minute timescales** as individual TLS defects drift through resonance
(Klimov et al., *Phys. Rev. Lett.* **121**, 090502 (2018)). Any quantity
extracted from a relaxation trajectory averaged over such an ensemble
inherits a bias, and a fitted noise model happily absorbs whatever the data
offers into whatever channels the model contains.

The gap between these two facts is the socket this document plugs into.
Measuring Γ↑ and Γ↓ separately at ±ω₀ is of course your community's founding
move (Schoelkopf et al. 2003 characterize a reservoir by exactly these two
numbers). What we have not found is the register-level bookkeeping: a
**symmetry-grade scalar on a fitted multi-qubit generator** that isolates the
emission-absorption asymmetry from everything else a fit absorbs, with an
exact calibration constant, together with a cheap protocol that diagnoses its
own bias when T1 fluctuates under it.

---

## 2. The object we hand you

Take any N-qubit Lindbladian L (a fitted one, from process tomography or
Lindblad estimation), written as a generator on the Pauli transfer matrix
(PTM), normalized so the 4^N Pauli strings are unit vectors. There is an
order-4 unitary W acting on such generators (a signed permutation of the
Pauli basis; the sister adapter
[Shifted Order-4 Chiral Symmetry](SHIFTED_ORDER4_CHIRAL_SYMMETRY.md) introduces
it, with proofs), and its square W² grades every generator into an even and an
odd sector, L_odd = (L − W²·L·W⁻²)/2. Every **Pauli channel** (a Pauli-string
dissipator, or any mixture of them) lies entirely in the even sector and
contributes nothing to the odd one. Coherent terms do NOT all drop out: the
odd sector also carries the commutator part of any Π²-odd Hamiltonian content
(residual Z detunings, Y-quadrature drives), so the scalar is defined on the
dissipative remainder:

    V(L) = ‖ L_odd − (its projection onto the span of Hamiltonian commutators) ‖_F

The subtraction is model-free on a raw fitted generator: the commutator span
is a fixed subspace (the −i[P, ·] of the Π²-odd Pauli strings), the
least-squares residual is unique, and the σ± signal occupies matrix cells no
commutator touches (the same support-orthogonality the F84 proof uses), so
none of it is projected away. Verified end to end: a bare 0.1 /µs residual
detuning beside a 0.005 /µs σ⁻ flux returns V = 0.005000 through the
projection, while the unprojected odd residue reads 0.888, coherent-dominated.
The reference implementation is
[`f81_pi_decomposition.py`](../../simulations/framework/diagnostics/f81_pi_decomposition.py):
`pi_decompose_M` for model-specified generators, and for a raw fitted one
`recover_H_odd_from_M_anti`, whose `fit_residual` IS this V. The constant
shift in the symmetry identity is even and drops out on its own. Then,
exactly:

    V(L) = 2^(N−1) · √( Σ_l (Γ↓_l − Γ↑_l)² )

for independent per-qubit emission Γ↓ and absorption Γ↑ channels. The scalar
has five properties, each a theorem, not a fit
([F82](../proofs/PROOF_F82_T1_DISSIPATOR_CORRECTION.md),
[F84](../proofs/PROOF_F84_AMPLITUDE_DAMPING.md), machine-verified):

1. **Exactly zero for every Pauli channel and mixture:** X-, Y-, Z-dephasing,
   depolarizing, and correlated Pauli-string noise (including
   always-on-ZZ-induced correlated dephasing) all contribute 0, identically.
   (Not every *unital* channel; the boundary is the Pauli axis, Section 6.)
2. **Hamiltonian-independent:** coherent terms are subtracted exactly; drive
   and crosstalk Hamiltonians do not enter.
3. **Dephasing-rate-independent:** T2 content, of any strength and spatial
   profile, does not enter.
4. **Degree-one in the net fluxes:** V reads Γ↓ − Γ↑ per qubit (the net
   one-way rate, not the total 1/T1 = Γ↓ + Γ↑), scaling linearly with a
   common net flux; across qubits it is a root-sum-square.
5. **Temperature-independent for a linear bath:** Γ↓ − Γ↑ = γ_vac(n̄+1) −
   γ_vac·n̄ = γ_vac; the thermal traffic cancels and the vacuum
   (spontaneous-emission) rate remains. (A saturable bath whose coupling
   itself drifts with temperature moves γ_vac; Section 3.)

For the σ± family, at N = 1 the entire scalar lives in **one PTM generator
element**, the (Z ← I) transfer entry: the generator's affine drift term,
i.e. d⟨Z⟩/dt evaluated at the maximally mixed state. Within the Pauli
+ σ± class (Section 6 for the boundary), that is the whole object: a
register-level, closed-form packaging of the velocity with which the identity
state polarizes.

---

## 3. The exact identifications (why this is yours already)

**First identification.** By the golden rule, a qubit's emission and
absorption rates sample the bath spectrum at its transition frequency:
Γ↓ ∝ S(+ω₀), Γ↑ ∝ S(−ω₀), with one common prefactor when both rates come
from one transverse coupling operator to one bath (Schoelkopf et al. 2003;
Clerk et al. 2010; Fourier convention S(ω) = ∫dt e^{iωt}⟨F(t)F(0)⟩, the one
both sources use). So

    Γ↓ − Γ↑  ∝  S(+ω₀) − S(−ω₀),

the antisymmetric part of the noise. The scalar V is therefore a
**register-level root-sum-square of the quantum part of each qubit's noise at
its own frequency**, with the calibration constant 2^(N−1) exact in rate
units (in spectral units each qubit's term carries its own coupling
prefactor, so the sum is coupling-weighted).

**Second identification.** In your canon the antisymmetric part of an
equilibrium spectrum is free of the thermal-occupation factor: S(+ω) − S(−ω)
is fixed by the dissipative response χ''(ω) alone, while the symmetric part
carries the coth(ħω/2kT). For a bath whose χ'' is itself
temperature-independent (any linear bath) the antisymmetric part is
temperature-independent outright, and Property 5 is then the same statement,
derived independently on the Lindblad side as an operator-norm identity. The
two theorems meet: **V is a direct, closed-form reading of exactly the
coth-free, quantum component of the noise**, blind by symmetry to everything
Pauli-classical and everything coherent. (A saturable TLS bath, whose χ''
carries a tanh(ħω/2kT), moves γ_vac with temperature; the coth-freedom
survives, the outright T-independence does not.)

These two are identifications, not analogies; both sides are textbook in
their own homes and the algebra matches at the weak-coupling (Lindblad)
level, which is also the statement's limit of validity (Section 6).

---

## 4. The measured record (one device, pre-registered, honest)

The scalar's input is a fitted generator, and fitted generators inherit their
data's pathologies. We flew the minimal protocol on IBM ibm_kingston
(2026-07-05, job d951mhkql68s73ca3u0g, three qubits, ≈ 1.6 QPU min,
pre-registered with bands before the shot; full record in
[F81 Violation: the Hardware Bridge](../../experiments/F81_VIOLATION_HARDWARE_BRIDGE.md)):

- **Two legs in one job:** relax |111⟩ and |000⟩ over the same 0-320 µs delay
  grid. Under a weakly coupled, time-stationary (Markovian) thermal bath both
  legs settle to the same asymptote with the same rate; the pre-registered
  **meeting test** compares the two fitted asymptotes at 2σ, and a SPLIT
  certifies departure from that class.
- **Result:** the |000⟩ legs stayed flat at ⟨Z⟩ = 0.98-1.00 for all 320 µs:
  excited-state populations of 0.83 ± 0.21 / 0.31 ± 0.27 / 0.23 ± 0.57 %
  (readout-mitigated in-job; asymptote reading; the last consistent with
  zero), joint-fit Γ↑ ≤ 7·10⁻⁵ /µs. The |1⟩ legs, however, produced free-fit
  "asymptotes" of 0.59-0.74: **SPLIT at 5-14σ.** The one-leg numbers would
  have read as 13-20% "thermal population"; the second leg shows they are an
  artifact.
- **The artifact is consistent with your Klimov effect, caught in a fit:** if
  T1 fluctuates over the minutes of the job, the ensemble-averaged trajectory
  is slower than the single exponential at its mean rate (Jensen), and a
  free-asymptote single-exponential fit reads that shape as a depressed
  asymptote, i.e. as false temperature. Pinning the asymptote at the
  |000⟩-leg value, a two-rate mixture fits all three qubits, and on one qubit
  the mixture's mean rate reproduces the same-morning calibration T1 to 1%
  (206 vs 208 µs) while its components differ by 5× (136 µs and 687 µs at
  57%/43% weight; window-limited, so the components are effective, not
  resolved). A static bi-exponential (a near-resonant TLS hybridized within
  every shot) fits the same curves as temporal switching, and preparation or
  leakage contributions to the small fast components are not excluded; the
  certificate does not need the mechanism, since any of these breaks the
  single-exponential model and the meeting test flags it.

The general statement we take from this, and hand over: **any
asymptote-based readout, including excited-state thermometry from
inversion-recovery tails, inherits a fluctuating-T1 bias, and the two-leg
meeting test is a two-line certificate against it.** (For precision
thermometry your community already prefers direct population methods, e.g.
Jin et al., *Phys. Rev. Lett.* **114**, 240501 (2015), building on the
protocol of Geerlings et al., *Phys. Rev. Lett.* **110**, 120501 (2013);
those are immune to this bias for independent reasons, and the meeting test
is the cheap in-situ check when all you have is relaxometry.)

---

## 5. The usable protocol, and one methodological rule

**To measure V's σ±-family component on hardware** (per qubit, no
tomography; the leg protocol reads the identity-velocity projection, which
is the whole of V within the Pauli + σ± class of Section 6): the up leg
reads the absorption side (the slope of ⟨Z⟩ from the ground state is −2Γ↑;
in the flown run the up-leg asymptotes pinned the populations at 0.2-0.8%
and the joint fit bounded Γ↑ ≤ 7·10⁻⁵ /µs at 8192 shots); the down leg
gives the total rate, taken fluctuation-robustly (the early-time rate or a
pinned-asymptote mixture mean, not the free-asymptote single exponential);
the net flux is Γ↓ − Γ↑, and V follows from the closed form. Run both legs
in one job and apply the meeting test; a SPLIT on a clean-χ² qubit means
the down leg's shape is distorted and its free asymptote must not be
trusted. In the flown run the SPLIT fired, so the packaged V was withheld
by the pre-registered rules: that run delivers the bias certificate, not a
V value. And since V is a norm, statistical noise in any fitted generator
biases it strictly upward: quote it against a refit-noise null
distribution, not against zero.

**To compute V on your existing fitted models**, one rule matters, and it is
easy to state because we paid for it: **if the fit's channel set is exactly
{dephasing + σ⁻}, then V equals the closed form of the fitted T1 rate by
construction, and you have learned nothing beyond the fit parameter.** The
scalar is informative only when the fitted generator has more freedom than
the σ± family: a free-form Lindblad fit, a full PTM, or a model with
channels you did not preselect. (On four such parameterized fits from our
own earlier data the identity V/closed-form = 1.000000000000 held to all
digits, which is the demonstration of the trap, not a measurement.)

**The falsifiable statement:** on any platform where T1 fluctuates on the
timescale of a characterization job (the Klimov regime), the one-leg
free-asymptote fit will systematically overstate the excited-state
population, and the two-leg meeting test will SPLIT; on a T1-stable platform
both legs meet and the one-leg number is honest. Either outcome is
informative, and the test costs one extra prepared-ground relaxation curve.

---

## 6. What we are not claiming

- The identification of Section 3 lives at the **weak-coupling / Lindblad**
  level (golden-rule rates). Strong-coupling or structured-bath corrections
  are outside the closed form's domain.
- The closed form is exact for **independent local σ± channels**. Correlated
  non-unital channels (σ⁻⊗σ⁻, σ⁺⊗σ⁺ and kin) write into the same odd cells
  with either sign, so **no bound in either direction holds in general**: a
  σ⁺⊗σ⁺ admixture can pull V below the local closed form and can cancel the
  identity velocities entirely at nonzero V (verified: local flux 0.01 per
  qubit plus σ⁺⊗σ⁺ at 0.02 gives exactly zero measured velocity at
  V = 0.0245 /µs). Outside the independent-σ± class, V and the leg protocol
  measure different objects.
- **Pauli, not unital:** the zero-property covers Pauli channels and their
  mixtures, not every unital channel. Tilted-axis dephasing, a Lindblad
  operator along a non-Pauli axis such as (X+Z)/√2 (the shape flux noise
  takes off a sweet spot), is unital yet contributes √2·γ to V at N = 1,
  with zero identity-state velocity (pinned in the repository's test suite).
  On a free-form fit, V therefore certifies departure from the Pauli + σ±
  model class along the odd sector; it equals the net-flux closed form
  within that class, and after twirling by the **diagonal (Z-string)
  subgroup**, which removes non-Pauli-axis odd content while preserving the
  net flux. (A full Pauli twirl instead removes the flux itself: V = 0.)
- The hardware record is **one device, one day, three qubits**, with a
  same-day simulator validation of the instrument; it demonstrates the
  protocol and the bias, not device statistics.
- We are not offering a better thermometer; direct population methods beat
  asymptote extraction for p_th. Our contribution is the symmetry scalar,
  its exact calibration, and the self-diagnosing protocol.
- The symmetry W itself is placed relative to the Lindbladian
  symmetry-classification literature in the
  [sister adapter](SHIFTED_ORDER4_CHIRAL_SYMMETRY.md), and only at the level
  argued there.

---

## 7. The stance-objects (what to take with you)

- **The scalar and its constant:** V = 2^(N−1)·√(Σ_l (Γ↓−Γ↑)²_l); zero for
  every Pauli channel and mixture, Hamiltonian-blind, dephasing-blind,
  temperature-independent for a linear bath.
- **The identification:** V = a register-level root-sum-square of
  S(+ω₀) − S(−ω₀), the quantum part of the noise, read at the
  process-matrix level.
- **The locality:** at N = 1 the σ±-family scalar is the (Z ← I) PTM
  generator element, the fixed-point displacement.
- **The protocol:** two relaxation legs in one job; the meeting test;
  populations from the up-leg asymptote, joint-fit Γ↑ bounded at
  ≤ 7·10⁻⁵ /µs per 8192 shots.
- **The rule:** V on a {dephasing + σ⁻}-parameterized fit is a tautology;
  use a free-form fit or measure the legs, and quote fitted-generator V
  against a refit-noise null (a norm is biased upward).
- **The numbers:** Kingston 2026-07-05, p_th = 0.2-0.8% (up-leg asymptotes,
  readout-mitigated), one-leg false "thermal population" of 13-20% exposed
  at 5-14σ by the second leg, ≈ 1.6 QPU minutes total.

The phenomenon's name in our repository stays home, as the arc's rule
requires. What travels is the object: an exact symmetry scalar that isolates
the one part of qubit noise that no Pauli channel can imitate, and a
protocol that tells you when your own relaxation data is lying about it.
