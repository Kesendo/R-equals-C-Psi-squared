# A bond and a seat under fixed illumination

**Status:** local `(1,1)`-block producer at fixed positive `γ₀`, 2026-09-23;
finite-sample numerical finding, no new F number or hardware confirmation.
**Run:** `python simulations/handshake_bond_seat_fixed_gamma.py`.
**Retained output:** [fixed-γ producer run](../simulations/results/handshake_bond_seat_fixed_gamma_run.txt).
**Controls:** `python -m pytest -q simulations/test_handshake_bond_seat_fixed_gamma.py`.
**Coherent comparison:** [the bond/seat γ=0 calculation](THE_BOND_AND_THE_SEAT.md).

## System and question

Take one excitation on a seven-site open XY chain, prepared in the
highest-energy standing wave `ψ₁` of the unperturbed hopping Hamiltonian
`H₀ = J A_path`. A small additive change `δJ_b` acts on one bond `b` through
`V_b = |b⟩⟨b+1| + |b+1⟩⟨b|`. Read the population at site `j` over time.
The question is which independent small bond changes that local time trace
can distinguish when dephasing remains present.

Here `γ₀ > 0` is **fixed**. In units where `γ₀=1`, the producer compares
`J=1,2,10`, hence `Q=J/γ₀=1,2,10`. This changes the hopping, not the
dephasing constant. The primary profile puts the same `γ₀` on every site.
As a separate comparison to F157's one-watched-site model, the producer also
puts that same `γ₀` on the centre site alone. That profile is a different
choice of dephasing channels, not a different value of `γ₀`.

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

Ten times are sampled. At `Q=10` their `γ₀t` values are
`[0.02,0.04,0.07,0.10,0.15,0.23,0.35,0.50,0.70,1.0]`; for `Q=1,2`
the observation times are rescaled so the same `Jt` values are compared while
`γ₀` stays fixed; each grid is printed in the retained output. For each
readout site the resulting matrix has ten rows (times)
and six columns (bond changes). Its numerical ranks are:

| `Q` | Uniform illumination, sites `j=0..6` | Centre-only channel, sites `j=0..6` |
|---:|:---|:---|
| 1 | `[6,6,6,3,6,6,6]` | `[6,6,6,3,6,6,6]` |
| 2 | `[6,6,6,3,6,6,6]` | `[6,6,6,3,6,6,6]` |
| 10 | `[6,6,6,3,6,6,6]` | `[6,6,6,3,6,6,6]` |

The smallest retained singular value in these 42 matrices is approximately
`3.14×10⁹` times the declared roundoff rank tolerance
`64 ε_machine·max(10,6)·σ_max` of each response matrix. Rank here counts
independent small-bond combinations visible to the sampled population trace;
it does not imply each individual bond can be identified. This is a robust
rank reading for these samples, not a theorem for every `Q`, time set, or
measurement noise level. The independent tangent ODE agrees with four
selected Fréchet slopes; its largest residual is `0.650` times the
machine-error scale
`ε_machine N (|t|(1+|Jt|)+|s_Fréchet|+|s_ODE|)`, against budget `16`.
Halving the rate **only** in the
Fréchet generator makes that gate fail.

The centre cannot distinguish the two bonds within any reflection pair:
reflection fixes the centre site, prepared `ψ₁`, and both dephasing profiles.
It therefore caps the centre's continuous-trace rank at three. The tested
sample reaches three. Across all ten `Q=10` times and both profiles, the
largest centre mirror residual is `6.94×10⁻¹⁷`. At `Q=10`, `γ₀t=0.1`, site
`j=2` separates the two end bonds by `0.005450` in slope under uniform
illumination.

A uniform change of all bonds no longer has zero response. At that same
sample the sum of the six slopes at site `j=2` is `−0.000951315` for uniform
illumination and `−0.000614574` for the centre-only channel. With `γ₀`
fixed, changing every hopping changes `J/γ₀`, so the population trajectory
changes even under a uniform bond perturbation.

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
