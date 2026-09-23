# A bond and a seat under fixed illumination

**Status:** local `(1,1)`-block producer at the repository's `γ₀=0.05`,
`J_hop=0.075` carrier-book point, 2026-09-23; finite-sample numerical finding,
no new F number or hardware confirmation.
**Run:** `python simulations/handshake_bond_seat_fixed_gamma.py`.
**Retained output:** [γ₀=0.05 producer run](../simulations/results/handshake_bond_seat_repo_gamma005_run.txt).
**Controls:** `python -m pytest -q simulations/test_handshake_bond_seat_fixed_gamma.py`.
**Coherent comparison:** [the bond/seat γ=0 calculation](THE_BOND_AND_THE_SEAT.md).
**Finite-shot continuation:** [one-million-shot known-defect classification](THE_BOND_AND_THE_SEAT_WITH_FINITE_SHOTS.md).

The [F124 proof](../docs/proofs/PROOF_HANDSHAKE_TRANSITION_INVARIANT.md)
owns the bond-to-mode matrix, [F157](../docs/ANALYTICAL_FORMULAS.md) the
blind-seat node count, and [Handshake Geometry](../hypotheses/HANDSHAKE_GEOMETRY.md)
uses the same bond operator. The [Q anchor map](../docs/Q_REGIME_ANCHORS.md)
supplies `γ₀=0.05`, `J=0.075`, while the [glossary's coupling table](../docs/GLOSSARY.md)
fixes which `J` convention that means. The [open arc](../compute/RCPsiSquared.Core/OpenArcs/OpenArcsRegistry.cs)
already separates dissipative response from Hamiltonian mixing; the
[Confirmations registry](../simulations/framework/confirmations.py) does not
make this bond/seat numerical run a hardware confirmation. The
[caught-errors ledger](../docs/CAUGHT_ERRORS.md) records the factor-of-two
trap between the two XY coupling books.

## System and question

Take one excitation on a seven-site open XY chain, prepared in the
highest-energy standing wave `ψ₁` of the unperturbed hopping Hamiltonian
`H₀ = J_hop A_path = (J_hop/2) Σ_b(X_b X_{b+1}+Y_b Y_{b+1})`.
A small additive change `δJ_b` acts on one bond `b` through
`V_b = |b⟩⟨b+1| + |b+1⟩⟨b|`. Read the population at site `j` over time.
The question is which independent small bond changes that local time trace
can distinguish when dephasing remains present.

Here `γ₀=0.05` is **fixed** and `J_hop=0.075`, hence
`Q=J_hop/γ₀=1.5`. This is the repository's convenient numerical anchor,
not a measured value of a universal constant. The producer's `J_hop` is the
single-excitation hopping and the coefficient of `(XX+YY)/2`. A script using
`J_Pauli·(XX+YY)` would need `J_Pauli=0.0375` for the same Hamiltonian; using
`0.075` in that second book would instead give hopping `0.15` and carrier
`Q=3`. The primary profile puts the same `γ₀` on every site.
Changing only units would map this run to `γ₀'=1`, `J_hop'=1.5`,
`t'=0.05t`: the response to an **additive** hopping change then scales by
`20`. The scale-control test checks that relation on both dephasing profiles;
`γ₀=1, J=1` would be `Q=1`, not this point.
As a separate comparison to F157's one-watched-site model, the producer also
puts that same `γ₀` on the centre site alone. That profile is a different
choice of dephasing channels, not a different value of `γ₀`. The total rate
`Σγ` changes from `0.35` (all seven sites) to `0.05` (centre only), so this
comparison does not isolate placement at a matched total rate.

The density matrix stays in the single-excitation `(1,1)` block: its cells
are `|a⟩⟨b|`, so the block has `N²` complex coordinates. Its equation is

```text
ρ̇ = -i[H₀ + δJ_b V_b, ρ] + D(ρ),
D(ρ)_{ab} = -2(γ_a + γ_b) ρ_{ab}  for a ≠ b,   D(ρ)_{aa} = 0.
```

For uniform illumination `γ_a=γ₀`, every off-diagonal pays `−4γ₀`.
For the centre-only comparison `γ_3=γ₀` and the other site rates vanish.
This is the convention used by [Cone](../compute/MirrorWorld/Cone.cs) and the
[single-excitation Haken-Strobl block](../simulations/coherence_horizon_se_block.py).

The producer computes `∂ρ_jj(t)/∂δJ_b` by a Fréchet derivative of this
`N²×N²` propagator. An independent matrix-shaped tangent ODE integrates
`ρ̇=Lρ` and `σ̇=Lσ−i[V_b,ρ]` with `σ(0)=0`. It shares neither the Kronecker
generator construction nor the matrix-exponential derivative.

## Result at the declared time samples

At `γ₀=0.05`, `J_hop=0.075` the ten physical times are
`t=[2,4,6,8,10,12,15,18,24,30]`, spanning `γ₀t=0.1..1.5`.
Each readout matrix has ten time rows and six bond-direction columns.
Subdividing every interval into four gives 37 time samples. The numerical
ranks at sites `j=0..6` are `[6,6,6,3,6,6,6]` for both uniform and
centre-only dephasing, on **both** grids. The rank tolerance is
`64 ε_machine·max(n_times,6)·σ_max`; the smallest retained singular value
across these matrices is `5.143×10⁷` times its respective tolerance. On the
ten-time grid the weakest noncentral singular value is only `7.309×10⁻⁶`
of its matrix's leading singular value. Thus all six directions clear numerical roundoff,
while their recovery under finite measurement noise is not established.
The producer's PASS also requires the smallest retained-to-leading ratio on
either grid to exceed `√ε_machine = 1.490×10⁻⁸`; the observed minimum is
`7.309×10⁻⁶`. This is a numerical-separation check, not an instrument-noise
budget.

An independent matrix-shaped tangent ODE agrees with four selected Fréchet
slopes at `t=10,20`; its largest residual is `0.046` times the machine-error
scale `ε_machine N (|t|(1+|J_hop t|)+|s_Fréchet|+|s_ODE|)`, against budget
`16`. Halving the rate **only** in the Fréchet generator makes that check fail.
The separate explicit-site-`Z` Lindblad ODE in the test file also reproduces
all `10×6×7×2 = 840` sampled slopes within its integration-error budget,
without calling the producer's Kronecker generator or Fréchet derivative.

The centre cannot distinguish the two bonds within any reflection pair:
reflection fixes the centre site, prepared `ψ₁`, and both dephasing profiles.
It therefore caps the centre's continuous-trace rank at three. The tested
sample reaches three. Across the ten times and both profiles, the largest
centre mirror residual is `1.08×10⁻¹⁵`. At `t=10`, site `j=2` separates
the two end bonds by `0.184782` in slope under uniform illumination.
The PASS check requires that separation to exceed the same roundoff-scale
budget used for the centre reflection residual; deleting it while preserving
the total bond slope makes the check fail.

A uniform change of all bonds no longer has zero response. At that same
sample the sum of the six slopes at site `j=2` is `−0.192335351` for uniform
illumination and `−0.520299408` for the centre-only channel. With `γ₀`
fixed, changing every hopping changes `J/γ₀`, so the population trajectory
changes even under a uniform bond perturbation. The earlier
[`γ₀=1` comparison](../simulations/results/handshake_bond_seat_fixed_gamma_run.txt)
is still runnable with `--unitized-control`; its `Q=1,2,10` values are
different points from this `Q=1.5` run.

## What happens to the F124–F157 connection

[F124](../docs/proofs/PROOF_HANDSHAKE_TRANSITION_INVARIANT.md) gives the full
Hamiltonian bond-to-mode transition matrix. [F157](../docs/ANALYTICAL_FORMULAS.md)
counts sine-mode nodes at a watched site. Together they give the exact
`γ=0` continuous-trace rank
`N−gcd(j+1,N+1)−1`, which at `N=7` is `[5,4,5,2,5,4,5]`.
That remains a valid **mathematical calibration** and is produced by the
[coherent-limit script](../simulations/handshake_bond_seat_readout.py).
An F157 blind mode can have zero *effective* decay at a watched site while
`γ₀` stays positive. That is different from deleting the dissipator from
every density-matrix cell, which is what the `γ=0` calibration does.

At fixed positive `γ₀` the density generator acts on mode *pairs* and mixes
the Hamiltonian channels. The coherent node count therefore does not count
the bond-to-population response channels in this regime. In particular,
the centre's sampled rank is three rather than two, and the noncentral sites
reach six in the declared fixtures. This finding does not alter F157's
blind-state theorem on its own object.

The repository's [primordial γ₀ hypothesis](../hypotheses/PRIMORDIAL_GAMMA_CONSTANT.md)
motivates holding `γ₀` fixed; this producer does not test whether it is a
universal constant. Nor does it compute the dissipative PTF `α` profile or
certify recovery from a finite, noisy instrument. A structural next step is
to derive a fixed-`γ₀` observability criterion for the continuous trace and
identify where the numerical rank can change with `Q` or the channel profile.
The linked finite-shot continuation tests a narrower thirteen-hypothesis,
known-magnitude classifier under ideal Z-position counts at this same
`γ₀/J_hop` point; it does not turn this continuous-trace rank into an
unrestricted or hardware bond-recovery result.
