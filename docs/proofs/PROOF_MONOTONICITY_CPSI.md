<!-- QUARTER-CURRENT -->
# Named CΨ Decays, and Where the Old Monotonicity Story Broke

Current reading: Parts 1–3 prove the named Bell+ channel formulas; exact local
Markov counterexamples rule out the universal pointwise and absorbing-boundary
statements.  The autonomous N=2 successive-local-maxima question remains open.

**Status:** Tier 1 for the named Bell+/channel formulas in Parts 1–3 and the instantaneous Pauli invariance
in Part 7. The former universal pointwise, absorbing-boundary, and local-control package is false. The
autonomous N=2 successive-local-maxima claim is open: the old argument did not prove it, and the exact
counterexamples below do not settle it.
**Date:** 2026-03-22; current-truth repair 2026-09-14 (history remains in git)
**Authors:** Thomas Wicht, Claude (Anthropic), Codex
**Reference formulas:** [F25](../ANALYTICAL_FORMULAS.md) (Bell+ Z closed form),
[F26](../ANALYTICAL_FORMULAS.md) (Bell+ Pauli closed forms), and
[F27](../ANALYTICAL_FORMULAS.md) (named-channel K values). The old universal F28 absorber reading is
withdrawn; [Subsystem Crossing](PROOF_SUBSYSTEM_CROSSING.md) retains the conditional convergence result.

---

## What survived the sharper look

The old story imagined a ball that could only roll downhill. That image is right for several named
trajectories, but CΨ is not a Lyapunov function for arbitrary local Markovian evolution. A local Hamiltonian
can turn population into computational-basis coherence, a Pauli pulse can change the following derivative
while leaving the instantaneous value fixed, and a local Lindblad semigroup can cross upward through 1/4.

The useful structure remains. Parts 1–3 preserve exact Bell+ formulas for named Z, Pauli, and amplitude-
damping channels. Part 4 preserves the familiar |01⟩ strong-coupling approximation as a named trajectory,
not a universal template. Part 5 records the failed proof step and exact counterexamples. Part 6 states the
conditional result that really is available: convergence to a fixed state below 1/4 implies eventual stay
below. Part 7 keeps instantaneous N-qubit Pauli invariance, without promoting it to trajectory invariance.

This correction makes the landscape more interesting, not less: named rivers still run downhill, while
the general dynamics can bend, turn, and cross. The finite rise atlas asks where those turns live.

---

## Exact positive statements

For the named Bell+ preparations and channel models below,

`CΨ(t) = Tr(ρ(t)²) · L₁(ρ(t))/(d-1)`

has the displayed closed form and negative derivative for t>0. These are trajectory theorems with named
initial states and generators. They are not a theorem for every state, every Hamiltonian, or every local
Markovian channel.

---

## Part 1: Pure Z-Dephasing

### Setup

Bell+ = (|00⟩ + |11⟩)/√2, Lindblad: L_k = √γ σ_z^(k) for k = 0, 1.

Under Z-dephasing, diagonals are preserved, off-diagonals decay:

```
ρ(t) = [[1/2,  0,  0,  f/2],
         [0,    0,  0,  0  ],
         [0,    0,  0,  0  ],
         [f/2,  0,  0,  1/2]]
```

where **f = e^{-4γt}** (each Z operator contributes 2γ to off-diagonal decay).

### CΨ in closed form

- **Purity:** C = Tr(ρ²) = 2·(1/2)² + 2·(f/2)² = (1 + f²)/2
- **L₁ coherence** (sum of absolute values of all off-diagonal elements)**:** |ρ₀₃| + |ρ₃₀| = f
- **ψ_norm:** Ψ = f/3 (d = 4, so d-1 = 3)
- **CΨ = C·Ψ = f(1 + f²)/6**

### Derivative

df/dt = -4γf, so by chain rule:

```
dCΨ/df = d/df [f(1+f²)/6] = (1 + 3f²)/6

dCΨ/dt = (dCΨ/df)(df/dt) = [(1 + 3f²)/6] · (-4γf)

        = -2γf(1 + 3f²)/3
```

### Sign

For f > 0 (all finite t) and γ > 0:
- f > 0 ✓
- (1 + 3f²) > 0 ✓ (always)
- γ > 0 ✓

**Therefore dCΨ/dt < 0 strictly for all t > 0. QED (Z-dephasing).**

### Crossing point

CΨ = 1/4 when f(1 + f²) = 3/2. Newton's method gives f* ≈ 0.8612.

t_cross = -ln(f*)/(4γ) = 0.1495/(4γ) → **K = γ·t_cross = 0.0374**

Numerical verification: K_Z = 0.0374 ± 0.0000. ✓

---

## Part 2: General Pauli Channels

### Setup

Local noise with rates (γ_x, γ_y, γ_z) on each qubit. Lindblad operators:
L_k^(i) = √γ_k · σ_k^(i) for k ∈ {x,y,z}, i ∈ {0,1}.

Bell+ stays Bell-diagonal. In the correlation representation:

```
ρ(t) = (I⊗I + c₁ σ_x⊗σ_x + c₂ σ_y⊗σ_y + c₃ σ_z⊗σ_z) / 4
```

where for Bell+ initial state:
- c₁(t) = e^{-αt}, with α = 4(γ_y + γ_z)
- c₂(t) = -e^{-βt}, with β = 4(γ_x + γ_z)
- c₃(t) = e^{-δt}, with δ = 4(γ_x + γ_y)

### CΨ in closed form

**Purity:** C = (1 + c₁² + c₂² + c₃²)/4 = (1 + e^{-2αt} + e^{-2βt} + e^{-2δt})/4

**L₁ coherence:** In computational basis, the off-diagonals are:
- |ρ₀₃| = |ρ₃₀| = |c₁ - c₂|/4 = (e^{-αt} + e^{-βt})/4
- |ρ₁₂| = |ρ₂₁| = |c₁ + c₂|/4 = |e^{-αt} - e^{-βt}|/4

L₁ = (|c₁-c₂| + |c₁+c₂|)/2 = max(e^{-αt}, e^{-βt})

(Using the identity (a+b+|a-b|)/2 = max(a,b) for a,b > 0.)

**ψ_norm:** Ψ = max(e^{-αt}, e^{-βt}) / 3

### Without loss of generality: α ≤ β

Then e^{-αt} ≥ e^{-βt} for all t ≥ 0, so L₁ = e^{-αt}.

Define u = e^{-αt}, v = e^{-βt}, w = e^{-δt}:

```
CΨ = u(1 + u² + v² + w²) / 12
```

### Derivative

```
dCΨ/dt = [du/dt · (1+u²+v²+w²) + u · (2u·du/dt + 2v·dv/dt + 2w·dw/dt)] / 12

       = [-αu(1+u²+v²+w²) + u(-2αu² - 2βv² - 2δw²)] / 12

       = -u/12 · [α(1+u²+v²+w²) + 2αu² + 2βv² + 2δw²]

       = -u/12 · [α + 3αu² + (α+2β)v² + (α+2δ)w²]
```

### Sign

Every coefficient in the bracket is ≥ 0:
- α ≥ 0
- 3α ≥ 0
- α + 2β ≥ 0
- α + 2δ ≥ 0

And every variable u², v², w² > 0 for finite t. The bracket is zero
only if α = β = δ = 0 (no noise). For any nonzero noise:

**dCΨ/dt < 0 strictly for all t > 0. QED (General Pauli).**

### K values for special cases

| Channel | α | β | δ | K = γ_eff · t_cross |
|---------|---|---|---|---------------------|
| Pure Z (γ) | 4γ | 4γ | 0 | 0.0374 |
| Pure X (γ) | 0 | 4γ | 4γ | 0.0866 |
| Pure Y (γ), re-sorted | 0 | 4γ | 4γ | 0.0866 |
| Depolarizing (γ/3 each) | 8γ/3 | 8γ/3 | 8γ/3 | 0.0440 |

**K_X = K_Y = ln(2)/8 = 0.0866433…; K_Z = 0.0374 is the odd one out.** The
discriminating fact is whether the *l₁-coherence norm* L₁ = Σ_{i≠j}|ρ_ij| =
max(|c_x|, |c_y|), the max(u, v) above, is decoherence-free or decays. It is the
norm and not any single element: under pure X the |00⟩⟨11| element runs ½ → ¼ while
the norm stays at 1, what leaves that element reappearing in |01⟩⟨10|. The XX
correlation is pinned under pure X; the YY correlation is pinned under pure Y; in **both**
cases L₁ stays nonzero, so L₁ = max(u,v) = 1 and CΨ = (1+v²)/6, crossing ¼ at
v² = ½ ⟹ K = ln(2)/8 = 0.0866433… Under pure Z **both** XX and YY decay, so L₁ → 0,
CΨ = u(1+u²)/6, and K_Z = 0.0374. (Note: F26 with γ_y only gives the physical rates
(α,β,δ) = (4γ, 0, 4γ); since β = 0 < α this **violates the WLOG α ≤ β**, so it must
be re-sorted to α = 0, giving L₁ = e^{−αt} = 1, exactly pure X's form. Dropping the
re-sort silently yields the pure-Z form and a wrong K_Y = K_Z; F27 in the registry
carries the same trap warning.)
All four are exact; what differs is their degree, and that is what decides how they
are written. K_X = K_Y = ln(2)/8 is a rational multiple of ln 2, because the pinned
coherence makes the crossing condition v² = ½ elementary, and it is carried as the
expression, never as a decimal. K_Z and K_depol are logs of cubic irrationals:
with u = e^{−4K} for pure Z and u = e^{−8K/3} for depolarizing, the crossings are
u³ + u − 3/2 = 0 and u³ + u/3 − 1 = 0, each with negative discriminant, hence one
real root reachable by real radicals (no casus irreducibilis), so

    K_Z     = −¼·ln(∛(3/4 + √(0.599537…)) + ∛(3/4 − √(0.599537…))) = 0.03735013…
    K_depol = −⅜·ln u*,  u* the real root of u³ + u/3 − 1        = 0.04395476…

The decimals 0.0374 and 0.0440 are four-place readings of those closed forms, not
a confession that none exists.

---

## Part 3: Amplitude Damping

### Setup

L_k = √γ |0⟩⟨1|^(k) for k = 0, 1. Non-unital: fixed point is |00⟩.

With q = e^{-γt}, p = 1-q:

```
ρ(t) = [[(1+p²)/2,  0,     0,     q/2  ],
         [0,         pq/2,  0,     0    ],
         [0,         0,     pq/2,  0    ],
         [q/2,       0,     0,     q²/2 ]]
```

### CΨ in closed form

**Purity:** C = a² + 2b² + d² + 2(q/2)²
where a = (1+p²)/2 = (2-2q+q²)/2, b = pq/2 = (1-q)q/2, d = q²/2.

C = (2-2q+q²)²/4 + (1-q)²q²/2 + q⁴/4 + q²/2

which collapses to the exact closed form (verified numerically to
machine precision at 7 sample points):

**C = (q² − q + 1)²**

**L₁ coherence:** Only ρ₀₃ and ρ₃₀ are nonzero off-diagonal:
L₁ = 2 · |q/2| = q

**ψ_norm:** Ψ = q/3

**CΨ = C(q) · q/3**

### Key observation

Both C and Ψ are functions of q = e^{-γt} only. Since dq/dt = -γq:

```
dCΨ/dt = (dCΨ/dq)(dq/dt) = (dCΨ/dq)(-γq)
```

We need dCΨ/dq > 0 (CΨ increases with q, i.e., decreases as q decays).

With C = (q² − q + 1)², CΨ = q(q² − q + 1)²/3, and the derivative
factors exactly:

```
3 · dCΨ/dq = (q² − q + 1)(5q² − 3q + 1)
```

Both quadratic factors have negative discriminant (1 − 4 < 0 and
9 − 20 < 0), so each is strictly positive for all real q, and
dCΨ/dq > 0 everywhere. (C(q) itself is NOT monotone in q; it dips to
9/16 at q = 1/2 and recovers. Only the product CΨ = q·C/3 is monotone,
which is exactly what the theorem needs.)

**Therefore dCΨ/dt = (positive)(−γq) < 0 for all t > 0, γ > 0. QED.**

The crossing q(q² − q + 1)² = 3/4 gives q* = 0.90219, so
K_AD = −ln(q*) = 0.1029.

### Numerical verification

K_AD = 0.1029 ± 0.0000 (CV = 0.0%). Heisenberg coupling J has zero
effect (Bell+ is eigenstate of H_Heisenberg).

---

## Summary

| Named Bell+ channel family | Result on that trajectory | K value | Basis |
|----------------------------|---------------------------|---------|-------|
| Pure Z-dephasing | strictly decreasing | 0.0373501… | exact (Part 1) |
| Pure X-noise | strictly decreasing | ln(2)/8 = 0.0866433… | exact (Part 2) |
| Pure Y-noise | strictly decreasing | ln(2)/8 = 0.0866433… | exact (Part 2) |
| Depolarizing | strictly decreasing | 0.0440… | exact (Part 2) |
| General Pauli rates `(γ_x,γ_y,γ_z)` | strictly decreasing in the stated ordering | varies | exact (Part 2) |
| Amplitude damping | strictly decreasing | 0.1029… | exact (Part 3) |
| Combined AD + Z | decreasing in 124 named configurations | varies | finite numerical catalogue |

Nothing in this table promotes the named Bell+ result to arbitrary initial states or generators, and it does
not make 1/4 an absorbing set for all Markovian dynamics.

---

## Earlier finite extensions (March 22, 2026)

**Script:** [monotonicity_remaining.py](../../simulations/monotonicity_remaining.py)

### General initial states (Test A)

Nineteen states were sampled (4 Bell, 5 product, 10 Haar-random). In those finite runs, every state that
started above 1/4 crossed below, and no increasing sequence of sampled peaks was resolved. That is a
catalogue observation, not a universal quantifier or a proof of the successive-local-maxima claim.

### Collective noise (Test B)

The named local and collective Z/X runs gave the same CΨ trajectories on Bell+. Anti-correlated Z noise
`Z₁-Z₂` leaves Bell+ in a decoherence-free subspace. These finite and symmetry-specific facts do not
classify all collective noise.

### N > 2 subsystems (Test C)

For the sampled GHZ (N=3,4,5) and W (N=3,4) preparations, the selected two-qubit reductions started below
1/4, stayed below on the sampled window, and approached zero. These named reduced-state runs do not inherit
a universal N=2 envelope theorem and do not establish an all-subsystem law.

## Part 4: Named strong-coupling |01⟩ approximation

### Setup

This section follows the named |01⟩ preparation under Heisenberg J + Z-dephasing γ. The state stays in the
{|01⟩, |10⟩} subspace. Define a = ρ_{01,01} (population), v = Im(ρ_{01,10})
(the only nonzero off-diagonal component, since Re = 0 by symmetry).

### Equations of motion

```
da/dt = -4Jv
dv/dt = -4γv - 2J(1 - 2a)
```

### Solution (damped oscillation)

With x = a - 1/2, the characteristic equation is λ² + 4γλ + 16J² = 0:

λ = -2γ ± 2i√(4J² - γ²) ≡ -2γ ± iω

For J >> γ (typical regime): ω ≈ 4J.

```
a(t) = 1/2 + (1/2)e^{-2γt}[cos(ωt) + (2γ/ω)sin(ωt)]             (exact)
v(t) = [J/√(4J²-γ²)]e^{-2γt}sin(ωt) ≡ V₀e^{-2γt}sin(ωt)        (exact)
```

The coefficient `2γ/ω` is fixed by the initial data `a(0)=1`, `v(0)=0`: it makes
`a'(0)=0` and satisfies both equations above. In the named `J≫γ` approximation the sine correction is
small, leaving `a(t)≈1/2+(1/2)e^{-2γt}cos(ωt)` and `ω≈4J`.

### CΨ for |01⟩

In the full 4×4 basis, only ρ_{01,10} and ρ_{10,01} are nonzero off-diagonal:

- **Purity:** C = 2a² - 2a + 1 + 2v² = 1/2 + 2(x² + v²)
- **L₁ coherence:** L₁ = 2|v|
- **Ψ:** Ψ = 2|v|/3
- **CΨ = [1/2 + 2(x² + v²)] × 2|v|/3**

### Approximate peak curve in the regime J ≫ γ

At the peaks of |sin(ωt)| (where |v| is maximal and cos(ωt) ≈ 0):

```
x² + v² ≈ e^{-4γt} [(1/4)cos²(ωt) + V₀²sin²(ωt)]
```

Since V₀ ≈ 1/2 for J >> γ: x² + v² ≈ (1/4)e^{-4γt}

**At the named approximate peaks:**

```
CΨ_max(t) ≈ [1/2 + (1/2)e^{-4γt}] × (2V₀/3)e^{-2γt}
```

### Derivative of the approximate peak curve

```
dCΨ_max/dt = (V₀/3) e^{-2γt} [-2γ(1 + e^{-4γt}) - 4γe^{-4γt}]
           = (V₀/3) e^{-2γt} [-2γ - 6γe^{-4γt}]
           < 0   for all γ > 0, t ≥ 0.
```

Thus this strong-coupling approximation decreases for the named |01⟩ trajectory. It is useful intuition,
not a proof for every J/γ, every initial state, or the autonomous N=2 peak sequence.

---

## Part 5: The historical envelope argument and its missing step

The historical claim said that, for every two-qubit state under arbitrary H and local Z-dephasing, the
successive local maxima of CΨ are non-increasing. Its spectral argument does not establish that claim.

A stable linear generator can expand a trajectory into decaying modes, but this gives only a decaying
upper bound on distance from a stationary subspace. It does not make a nonlinear, basis-dependent
functional monotone. In particular:

1. In the interaction picture, noncommuting H rotates the dephasing operators. Computational-basis
   off-diagonal entries need not obey independent scalar decays.
2. A decreasing bound `CΨ(t) ≤ B(t)` does not order the values of CΨ, or its local maxima, at two times.
3. CΨ is not a coordinatewise increasing function of absolute Liouvillian-mode amplitudes; modes can
   interfere and can create computational-basis coherence from populations.
4. Consecutive maxima need not have a common period or “similar oscillatory phase.”

So the former Step 5 assumed the ordering it needed to prove. The autonomous N=2 successive-local-maxima
claim remains an interesting question, but it is presently unproved.

### Exact pointwise and control counterexamples

Use `d=4`, `CΨ=P·L₁/3`, and first-site Z dephasing
`D(ρ)=(Z₁ρZ₁-ρ)/4`. Let

`ρ₀=((I+X/2+Z/2)/2) ⊗ |0⟩⟨0|`, `H=Y⊗I`.

This is a positive, trace-one density matrix, `CΨ(ρ₀)=1/8`, and no previously zero off-diagonal entry is
used in the one-sided derivative. Exact evaluation gives

`(P'(0), L₁'(0), CΨ'(0)) = (-1/8, +3/4, +1/6)`.

Thus **CΨ'(0) = +1/6** under a time-independent local Markovian generator. The H→0 mutation gives

`(P'(0), L₁'(0), CΨ'(0)) = (-1/8, -1/4, -1/12)`,

so **CΨ'(0) = -1/12** and the same exact gate distinguishes coherent creation from dephasing decay.

Two controls sharpen the basis dependence:

- A local Hadamard sends `|00⟩` from `CΨ=0` to `CΨ=1/3`.
- An active first-site Z pulse preserves the instantaneous `CΨ=1/8`, yet under the same laboratory H it
  changes the derivative from `+1/6` to `-1/3`. In a passive frame, transforming H to `Z₁HZ₁` as well,
  the derivative remains `+1/6`.

These examples disprove universal pointwise monotonicity and the shortcut from instantaneous Pauli
invariance to trajectory invariance. A positive initial derivative is not by itself a pair of successive
local maxima, so it **does not settle the successive-local-maxima question**.

### What the finite atlas says

`EnvelopeTheoremWitness` now exposes a finite rise atlas. Named Bell+ runs at N=4 and N=5 resolve
predecessor rises above a reporting bar, while one named N=3 run resolves none. The retained same-(N,Q,K)
pair agrees to six decimals. Those are finite numerical observations: they neither decide the autonomous
N=2 peak question nor yield an all-Q/all-N classification, an absence claim, or a mechanism. See
[The Finite Envelope-Rise Atlas](../../experiments/ENVELOPE_RISE_BOUNDARY.md).

---

## Part 6: Crossing is not a Markovian/non-Markovian watershed

The old page assigned upward crossing to bath memory. The fixed local Lindblad semigroup below crosses
upward without memory, so Markovianity is not the dividing line. Earlier structured-bath runs remain
interesting named trajectories, not evidence for an exclusive mechanism.

### Earlier structured-bath sightings

**Script:** [non_markovian_revival.py](../../simulations/non_markovian_revival.py)

A named structured-bath model (2 system qubits + 1 bath qubit in |+⟩) produced sampled subsystem-CΨ
revivals above 1/4 after an earlier downward crossing:

| J_SB | γ_B | Max Revival | Crossings ↑ | Sustained |
|------|-----|-------------|-------------|-----------|
| 5.0 | 0.01 | 0.3001 | 97 | 0.5 |
| 5.0 | 0.05 | 0.3009 | 17 | 0.3 |
| 5.0 | 0.50 | **0.3035** | 3 | 0.1 |
| 2.0 | 0.01 | 0.2731 | 37 | 1.4 |
| 0.5 | 0.01 | 0.2566 | 11 | 5.0 |

Largest value in this finite table: **CΨ = 0.3035**. It is one model catalogue, not a general statement
about non-Markovian dynamics or eventual behavior.

### An exact upward crossing under a fixed local semigroup

Take two local jumps with H=0,

`L₁=|+⟩⟨-|⊗I`, `L₂=I⊗|0⟩⟨1|`.

Their Liouvillian has one-dimensional kernel. Its unique stationary density is
`ρ*=|+0⟩⟨+0|`, and each dissipator separately annihilates that state. Hence `CΨ(ρ*)=1/3`, already above
1/4. Replacing the coherent-axis first jump by computational-axis damping `|0⟩⟨1|⊗I` changes the unique
target to `|00⟩⟨00|`, with CΨ=0; this mutation fails the same stationary-target gate.

There is also an explicit path, not merely an endpoint argument. Start at `|00⟩⟨00|` and write `q=e^{-t}`.
The exact trajectory obeys

`CΨ(t)=((1-q)(q²-q+2))/6`,

and

`dCΨ/dt=q(3q²-4q+3)/6 > 0` for `0<q≤1`.

At `q=1/8`,

`CΨ=847/3072=1/4+79/3072`.

So one fixed, local, time-independent Markovian semigroup carries a state from below 1/4 to above it. This
refutes the universal absorbing-set claim. The curve has no finite interior local maxima, so it does not
settle the successive-local-maxima question either.

### The conditional statement that survives

Let `ρ(t)` be a continuous trajectory with `ρ(t) → ρ*`. Because CΨ is continuous, if
`CΨ(ρ*) < 1/4` and the trajectory starts above 1/4, then it eventually stays below 1/4. It need not cross
only once and need not be monotone on the way. Named basis-aligned T1/T2/depolarizing models may use this
conclusion only under their stated convergence assumptions. See
[Subsystem Crossing](PROOF_SUBSYSTEM_CROSSING.md), which also records a separable primitive-CPTP fixed-state
counterexample with CΨ≈0.2935.

The repaired picture has two independent questions:

- What does a named generator do along a named trajectory? Parts 1–3 answer several such cases exactly.
- Where does a trajectory converge? A target below 1/4 gives eventual stay-below; a target above 1/4 can
  instead attract upward.

Neither question turns 1/4 into a universal one-way membrane. That loss of a slogan opens a better search:
classify generators, stationary states, and finite peak sequences separately.

---

## Part 7: Instantaneous N-qubit Pauli invariance (March 26, 2026)

### Theorem

CΨ is exactly invariant under the full N-qubit Pauli group at the instant of conjugation. For any density
matrix rho and any N-qubit Pauli operator U (tensor product of {I, X, Y, Z}):

    CΨ(U rho U+) = CΨ(rho)     (exact, not approximate)

### Proof sketch

CΨ = Tr(rho^2) x L1(rho) / (d-1). Purity Tr(rho^2) is invariant under
all unitaries (standard result). L1 coherence (sum of absolute values of
off-diagonal elements) is invariant under the Pauli group because Pauli
operators permute computational basis states with phase factors, and the
absolute value absorbs the phases.

### Numerical verification

Tested on partially decohered 2-qubit states with all 16 Pauli group
elements (I,X,Y,Z)^2: delta CΨ = 0.00e+00 in every case. Non-Pauli
unitaries (Ry, Rx, Hadamard, CNOT, random U(4)) change CΨ by up to
-0.24 (H x I on a typical state).

### What instantaneous invariance does—and does not—say

An active Pauli pulse leaves CΨ unchanged at that instant. It may still change the subsequent laboratory-
frame trajectory because the state has changed while H and the dissipators have not. The exact Z-pulse
example in Part 5 flips `CΨ'(0)` from `+1/6` to `-1/3` under the same H. Only a passive change of frame,
which conjugates the generator too, preserves the derivative.

Consequently, this algebraic identity does not establish trajectory invariance for dynamical decoupling,
does not forbid future upward crossing, and does not make local non-Pauli unitaries harmless. The local
Hadamard example `|00⟩→|+0⟩` changes CΨ from 0 to 1/3. Earlier periodic-pulse simulations remain named
protocol runs; they cannot be universalized from Pauli invariance alone.

Subsystem oscillations in coupled systems remain legitimate measured phenomena. They simply are not the
only route to an upward CΨ crossing, because Part 6 supplies a fixed local Markovian route as well.

---

## References

### Sibling proofs in the 1/4-boundary roadmap

- [Uniqueness Proof](UNIQUENESS_PROOF.md): Within the assumed normalized recurrence/power family, purity motivates α=2 but does not derive the recurrence; physical selection remains open. Layer 1 pins the discriminant-zero coordinate of the chosen α=2 normal form.
- [Subsystem Crossing](PROOF_SUBSYSTEM_CROSSING.md): the conditional convergence implication and the primitive-CPTP counterexample
- [Proof Roadmap Quarter Boundary](PROOF_ROADMAP_QUARTER_BOUNDARY.md): the seven-layer master roadmap; this proof is Layer 5

### F-formula registry

- [F25, F26, F27 in ANALYTICAL_FORMULAS](../ANALYTICAL_FORMULAS.md): Bell+ closed forms (Z and general Pauli) and K-values per channel, derived here
- F28 (historical fixed-point absorber): withdrawn as a universal claim; the conditional replacement lives in [Subsystem Crossing](PROOF_SUBSYSTEM_CROSSING.md)

### Scripts

- [generalized_pauli_channels.py](../../simulations/generalized_pauli_channels.py): 124/124 channel configurations verified
- [amplitude_damping_test.py](../../simulations/amplitude_damping_test.py): non-unital channel
- [non_markovian_revival.py](../../simulations/non_markovian_revival.py): transient revivals (Part 6)
- [monotonicity_remaining.py](../../simulations/monotonicity_remaining.py): Test A/B/C extensions (general states, collective noise, N>2 subsystems)

### Related experiments and hypotheses

- [Information Geometry](../../experiments/INFORMATION_GEOMETRY.md): Bures-geodesic picture of dCΨ/dt for Bell+ (states the Hamiltonian moves are off the geodesic)
- [Temporal Sacrifice](../../experiments/TEMPORAL_SACRIFICE.md): a named CΨ-heartbeat experiment; instantaneous Pauli invariance alone does not explain its whole trajectory
- [V-Effect Palindrome](../../experiments/V_EFFECT_PALINDROME.md): finite frequency-bin and F87 routing censuses under explicitly different generators; no mode ancestry or coupling-only complexity mechanism is inferred
- [Resonance Not Channel](../../hypotheses/RESONANCE_NOT_CHANNEL.md): the resonator framework around the named coupled-system runs
