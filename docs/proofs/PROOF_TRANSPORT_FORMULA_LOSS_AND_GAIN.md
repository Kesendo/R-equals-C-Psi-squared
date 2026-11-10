# PROOF: The Transport Formula Under Loss and Gain, and the Seventh Centred Trace at Second Order in the Rates

**Status:** Tier 1 derived by hand, for hopping bonds with ZZ on any graph under loss and gain per site, without fields or dephasing. Theorem A (the transport formula) holds for every Hamiltonian that conserves the excitation number. Theorem B gives the part of degree two in the rates of the seventh centred trace in closed form; that the trace has no parts beyond degree two and four rests on F137's three maps, cited and not shown again; its part of degree four is not derived. Not typed. Gate: [`simulations/seventh_moment_transport.py`](../../simulations/seventh_moment_transport.py), integers, every comparison exact.
**Date:** 2026-11-10
**Registry:** [F137](../ANALYTICAL_FORMULAS.md#f137), whose fifth centred trace is the lowest order of Theorem A.
**Authors:** Thomas Wicht, Claude
**Builds on:** F137 in [the formula registry](../ANALYTICAL_FORMULAS.md#f137) (the setting, the three maps that sort the traces by their degrees, the fifth centred trace); [`simulations/fifth_moment_flux_area.py`](../../simulations/fifth_moment_flux_area.py) (the generator the gate imports).

Under loss and gain an odd centred trace Tr((L + c/2)^m) that is not zero rules out the pairing of the spectrum about −c/2. F137 gives the fifth in closed form: flux through a triangle of hopping bonds times the area of the triangle its three rate points span. This proof states the formula behind that result on its own, reads it as a picture, and takes it one order further: the part of the seventh trace that is of degree two in the rates, on any graph, with ZZ.

<!-- STAGE 0 SWEEP RECORD (for us, not for the reader: hidden on GitHub, read by every session, agent and hook)

Named in repo primitives: the odd centred traces Tr((L + c/2)^m) of a generator with a number-conserving H under loss and gain per site (F137's setting without fields or dephasing), at second order in the rates.

Swept for Theorem A: the registry holds it, as step (3) of F137's derivation of the fifth centred trace, with the same f and the same reading (an excitation against a hole); FluxAreaFifthMomentClaim carries it in its doc comment and in its node on the derivation. This file gives it a proof of its own and adds nothing to its statement. Elsewhere docs/proofs uses Duhamel's formula in J_b (PROOF_F90_F86C2_BRIDGE, PROOF_F86B_OBSTRUCTION, a mention in PROOF_F89_PATH_D_CLOSED_FORM) and, in a rate, in PROOF_DEPHASING_FRONT_RENEWAL's refill term under Z dephasing; none of them for loss and gain.

Swept for Theorem B, the seventh trace under loss and gain in closed form: the registry (F137: the recorded 430080 of the temperature gradient, the sevenths 1708448/9 with fields and −263277/400 with X dephasing, both outside this class and single values; its pointer to WHICH_BONDS; F117 and F120, odd power sums under Z dephasing; no closed form at the seventh order with heating), docs/proofs (MIRROR_SYMMETRY_PROOF's thermal bullets, third moment only; PROOF_MOMENT_TOWER_PUMP_CHANNEL, Z dephasing; PROOF_F87_WINDOWED_MONOMIAL_CONVERSE §3, the bra walk and the ket walk round an odd cycle under Z dephasing; PROOF_WEIGHT1_DEGENERACY Theorem E, a triangle threaded by flux at the Zeno end; PROOF_ANTILINEAR_TRIANGLE §7, the weighted transpose at one ratio and transfers both ways), experiments (WHICH_BONDS_PAIR_UNDER_COOLING: the seventh in closed form on the open chain under cooling alone, of degree two in the rates and carrying the anisotropy |w|², which is zero for bonds that conserve the excitation number; its rings and its triangle, where under cooling alone the parts of degree two or less vanish; THERMAL_BREAKING, the break on the triangle off one line, no seventh; THE_CUBE_SITE_BY_SITE, a fifth with dephasing), the typed layer in both halves (CooledBondPowerSumsClaim and its witness, cooling alone; FluxAreaFifthMomentClaim and its witness, the fifth; ExactThermalGenerator.OddPowerSums, the generic route; no claim for a seventh under loss and gain), MirrorWorld (Sextant.OddTraceObstruction and Sextant.FluxAreaFifth; a test records that the temperature gradient is silent at the fifth and obstructs at the seventh; no closed form beyond the fifth), the gates fifth_moment_flux_area.py (the 430080) and cube_moves_gate.py, the OpenArcs registry (nothing on odd centred traces under loss and gain), the glossary (nothing), CAUGHT_ERRORS (nothing on fifth or seventh traces), fw.Confirmations (nothing), framework/diagnostics/odd_traces.py (the traces modulo primes, no closed form). No store holds the part of degree two of the seventh trace, its weights, the ring of five or the ring of four.

Checked for adjacency: Theorem B at cooling alone against WHICH_BONDS (γ↑ = 0 puts every rate point on one axis and the part of degree two is zero on any graph, as its count says for the parts of degree two or less; its seventh on the chain lives on |w|², outside this class; note that its s and a are parts of the coupling matrix, not rates); the part of degree four against the same document's triangle (under cooling alone and one class bond repeated it is zero there; with gain it is not, the 430080); the two ways round a loop against PROOF_F87_WINDOWED_MONOMIAL_CONVERSE §3 (the same picture under Z dephasing; F137's record already names the pair); the amplitudes f against F152's (0,1) block, the hop amplitudes as a one-excitation generator (f at lowest order is one entry of its exponential; not joined further); the weights' dependence on the origin of the rate plane against F137's open question on flux reversal (silent on a line where the three weights are equal; not decided). Theorem B is not typed.
-->

## 1. Setting and statements

N sites, H = Σ over bonds of J(X_aX_b + Y_aY_b) + D_z(X_aY_b − Y_aX_b) + Δ_ab·Z_aZ_b on any graph, a bond being free to carry ZZ without hopping. Write t_ab = J + iD_z for the amplitude of the hop a → b, t_ba its conjugate, and t_ab = 0 where no bond hops. With σ⁺ the raising operator and Z = −1 on the excited state, σ⁺σ⁻ = (1 − Z)/2, this is H = Σ_{a≠b} τ_ab·σ⁺_aσ⁻_b + Σ Δ_ab·Z_aZ_b with τ_ab = 2·t_ba, the coupling of the hop b → a. (τ with two indices is always this coupling; τ alone is the time.) Loss σ⁻_l at rate γ↓_l and gain σ⁺_l at rate γ↑_l on every site, no fields, no dephasing, c = Σ_l (γ↓_l + γ↑_l), and as in F137

  s_l = (γ↓_l + γ↑_l)/2,  δ_l = (γ↓_l − γ↑_l)/2,  L + c/2 = −i·ad_H + Σ_l (s_l·𝓗_l + δ_l·𝓐_l),

with 𝓗(ρ) = σ⁺ρσ⁻ + σ⁻ρσ⁺ and 𝓐(ρ) = σ⁻ρσ⁺ − σ⁺ρσ⁻ + ½{Z, ρ} on one site.

**Theorem A (transport formula).** Let H conserve the excitation number (any such H, fields and further terms included). For every pair of sites l, k the coefficient of s_l·δ_k in Tr e^{τ(L + c/2)} is

  τ·∫₀^τ (|f⁺_lk(x)|² − |f⁻_lk(x)|²) dx,   f^±_lk(x) = Tr(σ^±_l·e^{−iHx}·σ^∓_k·e^{−iH(τ−x)}),

the traces over the Hilbert space. Its coefficient of τ^m/m! is the coefficient of s_l·δ_k in Tr((L + c/2)^m). It is antisymmetric in (l, k) and zero at l = k.

**Theorem B (the seventh centred trace at second order in the rates).** In the setting above, Tr((L + c/2)⁷) has parts of degree two and four in the rates and no others. The part of degree two is

  −168·4^N · Σ_{l<k} β_lk · (γ↓_l·γ↑_k − γ↓_k·γ↑_l),

that is, 336·4^N·β_lk is the coefficient of s_l·δ_k, with β_lk = −β_kl and

  β_lk = Σ_j Im(t_lk·t_kj·t_jl) · (W_lkj + Q_lkj)
       − 15·Σ_{i,j} |t_li|²·Im(t_ik·t_kj·t_ji) − 15·Σ_{i,j} |t_kj|²·Im(t_lj·t_ji·t_il)
       + 15·Σ_{a,b,d} Im(t_lk·t_kd·t_db·t_ba·t_al) − 15·Σ_{m,i,j} Im(t_lm·t_mk·t_kj·t_ji·t_il)
       + 5·Σ_{i,j} Im(t_lk·t_kj·t_ji·t_il)·(Δ_lk − Δ_li − Δ_ij − Δ_jk)
       + 5·Σ_{m,j} Im(t_lm·t_mk·t_kj·t_jl)·(Δ_lj + Δ_jk),

every inner sum running over the ordered tuples of sites that are distinct and not l or k, and

  W_lkj = 30·S − 10·(d_l + d_k) − 5·d_j − 2·|t_lk|²,   S = Σ_bonds |t|²,  d_v = Σ_w |t_vw|² over all sites w,

  Q_lkj = 15·Σ Δ² (ZZ bonds that do not touch l, k, j) + 5·Σ Δ² (one end at l or k, the other off the three)
        + 10·Σ Δ² (one end at j, the other off the three)
        + 5·(Δ_lj² + Δ_jk²) + (10⁄3)·Δ_lk² + 2·Δ_lj·Δ_jk − (8⁄3)·Δ_lk·(Δ_lj + Δ_jk).

Read by lines. **The triangle** (l, k, j) carries the imaginary part of its loop times a weight that counts every bond of the graph by how it touches the triangle. For a hopping bond: 30·|t|² if it does not touch, 25 hanging at the third corner j, 20 hanging at l or k, 15 for the triangle's edges through j, 8 for the edge l–k; in one rule, 30 for every bond, minus 10 for each touch of l or k, minus 5 for a touch of j, and the edge l–k two less. A ZZ bond counts through the quadratic form Q (a bond may carry both): off the triangle's edges 15·Δ², minus 10 for a touch of l or k, minus 5 for a touch of j, and on the edges the last four terms, which follow no such rule. **The copy** (second line): where l is joined to a corner i of a triangle (i, k, j), the pair (l, k) reads that triangle's loop Im(t_ik·t_kj·t_ji), l in i's place, with the weight −15·|t_li|², the sign opposite to that of a triangle through l itself; the same for k, and so l pairs with both other corners. **The pentagon** (third line): a ring of five through l and k counts +15 where they are neighbours and −15 where they are two steps apart, which on one ring v₀ → v₁ → … → v₄ → v₀ adds up to 2520·4^N·Im(loop) times the pentagram of the five rate points minus their pentagon: Σ_n (γ↓_n·γ↑_{n+2} − γ↓_{n+2}·γ↑_n) minus Σ_n (γ↓_n·γ↑_{n+1} − γ↓_{n+1}·γ↑_n), the sites counted the way the loop runs, each sum twice a signed area. **The square** (last two lines) enters only with ZZ: a ring of four with flux, times a linear form in the Δ on its edges, the edge l–k against the other three where l and k are neighbours, one way from l to k against the other where they are opposite (the ordered pairs (m, j) and (j, m) together give the loop times Δ_lj + Δ_jk − Δ_lm − Δ_mk); ZZ on the diagonals does not enter.

On one triangle alone the weights W + Q of the three pairs are equal when the three |t| are equal and the three edges carry one Δ; then the part of degree two is, like the fifth trace, flux times the area of the rate triangle and does not change when the rate points are moved by one vector. With unequal weights it changes under such a move: it sees the origin of the rate plane. The temperature gradient that F137 records on the unit triangle at flux π/2 (γ = 2, n̄ = (1, 2, 3), Tr((L + c/2)⁷) = 430080) is such a case of equal weights on collinear points: its part of degree two is zero and the 430080 is all of degree four, which this proof does not derive. The two settings at N = 4 that F137 records from gate M11 as breaking above the fifth trace are decided here: on the unit triangle at flux π/2 and ratio ¼ with a fourth site at ratio 1⁄9 joined by a hopping bond, and on the ring of four at Δ = 7⁄10, flux π/2 and four ratios, the seventh trace is not zero and all of degree two, −44800 and −20380416⁄25, the ring through its ZZ alone (the ring of four of the fifth trace's gate, another setting, has −68239360 and zero without its ZZ).

## 2. Proof of Theorem A

Write G = −i·ad_H and expand Tr e^{τ(G + ε·𝓗_l + η·𝓐_k)} to first order in ε and in η. By Duhamel's formula and the cyclicity of the trace the coefficient of εη is

  τ·∫₀^τ Tr(𝓗_l·e^{xG}·𝓐_k·e^{(τ−x)G}) dx,

the trace now over the operators on the Hilbert space. Each of the four factors is a sum of maps ρ ↦ AρB, and the trace of ρ ↦ AρB over the operators is Tr A·Tr B; e^{xG} is ρ ↦ e^{−iHx}·ρ·e^{iHx}. So every term is a product of a trace over the ket side and one over the bra side. The term ½{Z, ρ} of 𝓐_k puts nothing that changes the excitation number beside the σ^± that 𝓗_l puts on each side, and the trace of a single σ^±_l between factors that conserve the number is zero. Of the remaining four terms the two that raise at l and at k on the same side vanish for the same reason. What is left is σ⁺_lρσ⁻_l with σ⁻_kρσ⁺_k, and σ⁻_lρσ⁺_l with −σ⁺_kρσ⁻_k:

  Tr(σ⁺_l e^{−iHx} σ⁻_k e^{−iH(τ−x)})·Tr(σ⁻_l e^{iH(τ−x)} σ⁺_k e^{iHx}) − (the same with σ⁺ and σ⁻ exchanged),

and the second trace of each product is the complex conjugate of the first. That is |f⁺_lk|² − |f⁻_lk|². The coefficient of τ^m/m! of Tr e^{τ(L + c/2)} is Tr((L + c/2)^m), and of it the terms with one 𝓗_l and one 𝓐_k. By cyclicity f⁺_kl(x) = f⁻_lk(τ − x) and f⁻_kl(x) = f⁺_lk(τ − x), so the integral is antisymmetric in (l, k) and zero at l = k. ∎

**The picture.** f⁺_lk is the amplitude that an excitation put in at l is taken out at k, f⁻_lk the same for a hole; a rate pair (l, k) reads the one in its square against the other. With one side as the ket and its conjugate as the bra: the excitation has to be carried from l to k on the ket side and on the bra side, and where the two take different ways round a loop of the graph their interference carries the flux.

## 3. The time circle

From here on the setting of §1. The spin flip X⊗N conjugates H to its complex conjugate H̄ and exchanges σ⁺ with σ⁻, so with y = τ − x

  f⁺_lk = Σ_n (−i)ⁿ·g_n(x, y),   f⁻_lk = Σ_n (−i)ⁿ·ḡ_n(x, y),   g_n = Σ_{p+q=n} A_pq·x^p·y^q/(p!·q!),   A_pq = Tr(σ⁺_l H^p σ⁻_k H^q),

ḡ_n being g_n with every coupling conjugated. Every g_n is a sum of products of couplings times real polynomials in x and y, so ḡ_n is its complex conjugate, and the terms (n, n′) and (n′, n) of |f⁺|² − |f⁻|² together are 2i·((−i)ⁿ·i^n′ − (−i)^n′·iⁿ)·Im(g_n·ḡ_n′); g₀ = 0 for l ≠ k. At order five, the order of the seventh trace:

  |f⁺_lk|² − |f⁻_lk|² at order five = 4·[Im(g₁·ḡ₄) − Im(g₂·ḡ₃)].   (1)

To evaluate g_n, place the operators on a circle of circumference τ = x + y in the order in which they act, the word read from the right: σ⁺_l at 0, then the arc of length y with q terms of H, σ⁻_k at y, then the arc of length x with p terms, back to 0; the polynomial x^p·y^q/(p!·q!) is the volume of the ordered points on the two arcs. The trace is a product over the sites. A site that carries no σ^± gives a factor 2. A site that carries some is excited from each of its σ⁺ forward to its next σ⁻ and gives 1 if they alternate round the circle and 0 otherwise. A term Δ_ab·Z_aZ_b is diagonal and multiplies by ±Δ_ab according to the two sites' states where it stands.

**A simple path.** For h hops along l = v₀ → v₁ → … → v_h = k through sites all different, with the hop v_{i−1} → v_i at the point p_i, the site v_i is excited on the arc from p_i forward to p_{i+1} (p₀ = 0, p_{h+1} = y). Every position of the hops is allowed, so with κ = 2^{N−2}

  g_h ∋ κ·2^{1−h}·τ_{v_h v_{h−1}}⋯τ_{v₁v₀}·τ^h,

independent of x: g₁ = κ·τ_kl·τ and the two-step part of g₂ is (κ/2)·Σ_j τ_kj·τ_jl·τ².

**The hard core.** Where a site is visited more than once it carries four or more σ^±, and only the positions in which they alternate are allowed. All the numbers of Theorem B at Δ = 0 come from these fractions (§5).

**ZZ.** With z terms Z_aZ_b among the n, the points of the hops fixed, the terms ZZ integrate to ⟨Φ^z⟩/z!, where Φ = Σ_bonds Δ_ab·∮ z_a(u)·z_b(u) du is the ZZ energy integrated round the circle, z_a = −1 where site a is excited and +1 where it is not, and ⟨·⟩ averages over the two states of every site that carries no σ^±, for which ⟨z_b⟩ = 0 and ⟨z_b·z_b′⟩ = δ_bb′. So a word of h hops along a simple path with z terms ZZ contributes

  κ·2^{1−h}·(the couplings of the path)·τ^{h+z}·F_hz(y/τ),   F_h0 = 1,   (2)

F_hz the integral of ⟨Φ^z⟩/z! over the positions of the hops at τ = 1, each hop's point uniform on the circle. Two formulas carry every average below: ζ_a = ∮ z_a du is 1 minus twice the length of the arc on which a is excited, and ∮ z_a·z_b du is 1 minus twice the two lengths plus four times the length on which both are excited.

## 4. The fifth trace

At order three, |f⁺|² − |f⁻|² = −4·Im(g₁·ḡ₂). The part of g₂ with one hop and one ZZ is a real multiple of τ_kl and drops out against g₁ = κ·τ_kl·τ; the two-step part gives −2κ²·τ³·Σ_j Im(τ_kl·conj(τ_kj·τ_jl)) = −4^N·τ³·Σ_j Im(t_lk·t_kj·t_jl), the direct hop on the one side against the detour over a third site on the other. This is step (3) of F137's derivation, and the rest is there.

## 5. The seventh trace without ZZ

Five hops in (1). Their product has to close a loop with an imaginary part: the triangle (l, k, j) with one bond walked there and back besides, a triangle elsewhere reached over one bond, or a ring of five. Write I_j = Im(t_lk·t_kj·t_jl); then Im(τ_kl·conj(τ_kj·τ_jl)) = 8·I_j, and |τ_e|² = 4·|t_e|².

**The triangle with a bond e walked there and back.** In g₃ the word is the direct hop with the bond, in g₄ the detour over j with the bond. A bond that shares a site with the path leaves one site more without a σ^±, which doubles the trace, and on the shared site four σ^± have to alternate, which only a fraction of the orders do. With e the bond and κ·τ^n the common factor, the coefficient of the couplings is

| position of e | in g₃ beside the direct hop | in g₄ beside the detour over j |
|---|---|---|
| touches none of the path's sites | 1⁄4 | 1⁄8 |
| one end on the path, the other off it | 1⁄6 (2 of 6 orders, at l or k) | 1⁄12 (8 of 24, at l, k or j) |
| the edge l–j or j–k | (one end at l or k: 1⁄6) | 1⁄24 (the hop on that edge stands twice; 2 words) |
| the edge l–k | (x³ + y³)/6 | φ₁/2 |

with φ₁ = (x⁴ + y⁴)/6 + (x³y + xy³)/2 at τ = 1 (of the words with all four terms on one arc 4 alternate, with three on one arc and one on the other 3, with two and two none; a word with p and q terms on the arcs weighs x^p·y^q/(p!·q!), and three sites carry σ^±, so the factor is κ/2). On the edge l–k only the orders with all three terms on one arc alternate in g₃, as the two-site problem shows directly: there f⁺ = −i·(τ_kl/|τ_kl|)·(sin |τ_kl|x + sin |τ_kl|y). One more word enters g₄ on the edge l–k: the direct hop twice and the other two edges of the triangle walked back, τ_kl·τ_kl·τ_jk·τ_lj, with the coefficient φ₂/2, φ₂ = (x⁴ + y⁴)/12 + (x³y + xy³)/6 + x²y²/2 (2, 1 and 2 words), and against g₁ it carries −8·I_j·|τ_kl|².

Put into (1), with the two-step part of g₂ against g₃ and g₁ against g₄, the coefficient of 8·I_j·|τ_e|²·κ² in ∫₀¹(Im g₁ḡ₄ − Im g₂ḡ₃) dy is: for a bond off the triangle 1⁄8 + 1⁄8 = 1⁄4; hanging at l or k 1⁄12 + 1⁄12 = 1⁄6; hanging at j 1⁄12 + 1⁄8 = 5⁄24 (beside the direct hop a bond at j touches nothing); an edge through j 1⁄24 + 1⁄12 = 1⁄8; the edge l–k ½·∫(φ₁ − φ₂) + ½·∫(x³ + y³)/6 = 1⁄40 + 1⁄24 = 1⁄15. The coefficient of s_l·δ_k in the seventh trace is 7!·4 times the integral, and 7!·4·κ²·8·4 = 40320·4^N = 336·4^N·120, so the five fractions are the weights 30, 20, 25, 15 and 8 of |t_e|². The direct hop with a triangle other than (l, k, j) walked once round comes with both directions of the loop at the same polynomial, the two being exchanged by renaming the triangle's sites off the path, and so has no imaginary part against g₁; round (l, k, j) itself the two directions are the two words of the edge l–k above, with different polynomials.

**The copy and the pentagon.** In −Im(g₂·ḡ₃) the two-step path l → m → k stands against the three-step path l → i → j → k, each simple, with coefficients κ/2 and κ/4. At m = i the product is |τ_il|² times the loop of the triangle (i, k, j), which is the copy with l in i's place, its sign the minus of (1), and at m = j the same for k; at m different from both it is the ring of five with l and k two steps apart. The ring of five through the bond l–k is g₁ against the four-step path in g₄, coefficient κ/8, with the plus of (1). Each comes to 7!·4·κ²·(1/8)·32 = 5040·4^N = 336·4^N·15 times the imaginary part of the product of the five amplitudes, with the signs of Theorem B.

## 6. The seventh trace with ZZ

With ZZ the five terms of H in (1) can be three hops and two ZZ, on a triangle, or four hops and one ZZ, on a ring of four; five hops are §5 unchanged, and one hop with ZZ has no imaginary part. Write P = xy at τ = 1, so ∫₀¹P dy = 1⁄6 and ∫₀¹P² dy = 1⁄30. The integrals (2), each over the positions of the hops on the circle:

- **One hop, one ZZ.** The two sites are in opposite states on the arc that holds the hop's point and in equal states on the other: F₁₁ = −Δ_lk·(1 − 4P).
- **One hop, two ZZ.** ⟨Φ²⟩ = Δ_lk²·(∮ z_l z_k)² + Σ_b (Δ_lb·ζ_l + Δ_kb·ζ_k)² + Σ Δ² over the bonds off l and k, b running over the other sites. Over the hop's point ζ_l² and ζ_k² average to 1⁄3 and ζ_l·ζ_k to 2P − 1⁄3, whose integral over y is zero:
  F₁₂ = ½·[Δ_lk²·(1 − 4P) + Σ_b ((Δ_lb² + Δ_kb²)/3 + 2·Δ_lb·Δ_kb·(2P − 1⁄3)) + Σ_off Δ²].
- **Two hops over j, one ZZ.** For two sites that follow each other on the path, ∮ z z averages to −1⁄3 whatever y is: for the first two sites of a path, excited from 0 to the hop point a and from a to the next point b, it is 1 − 2b where a < b and 2b − 1 where b < a, and each half integrates to −1⁄6; for l and k it averages to 2P − 1⁄3:
  F₂₁ = −(Δ_lj + Δ_jk)/3 + Δ_lk·(2P − 1⁄3).
- **Two hops over j, two ZZ.** Here ζ_l, ζ_j and ζ_k are uncorrelated, and the three integrals Z_e = ∮ z z of the triangle's edges have the second moments V below, each integrated piece by piece over the six orders of the two hop points and the point y:
  F₂₂ = ½·[Σ_{e,e′} V_ee′·Δ_e·Δ_e′ + Σ_b (Δ_lb² + Δ_jb² + Δ_kb²)/3 + Σ_off Δ²],
  V_lj,lj = V_jk,jk = 1⁄3,  V_lk,lk = 1⁄3 − 4P/3 + 8P²/3,  V_lj,jk = 4P²,  V_lj,lk = V_jk,lk = −2P/3 + 4P²/3,
  b now off l, j and k. As a check, Σ_e Z_e depends only on the number of excitations on the two arcs, and is 3x − y, −1 or 3y − x on the six orders; its second moment 1 − 4P + 16P² is the sum of the entries of V.
- **Three hops l → i → j → k, one ZZ.** Neighbours on the path −1⁄3; the two pairs one step further apart have independent arcs and average to zero; l and k as before:
  F₃₁ = −(Δ_li + Δ_ij + Δ_jk)/3 + Δ_lk·(2P − 1⁄3).

**The triangle.** In (1), g₁ against the detour with two ZZ, the detour against the direct hop with two ZZ, and the direct hop with one ZZ against the detour with one:

  ∫₀¹ (Im g₁ḡ₄ − Im g₂ḡ₃) dy ∋ 4κ²·Σ_j I_j·∫₀¹ (F₂₂ + F₁₂ − F₁₁·F₂₁) dy,

and 7!·4·4κ² = 336·4^N·15, so Q_lkj is 15 times the integral. Term by term: a bond off the three sites ½ + ½, so 15; one end at l or k 1⁄6 + 1⁄6, so 5; one end at j 1⁄6 + ½ (for the direct hop j is a spectator and the bond is off the path), so 10; two bonds from one outside site to two of the three have no mixed term after the integral over y; Δ_lj² and Δ_jk² 1⁄6 + 1⁄6, so 5; Δ_lk² 1⁄10 + 1⁄6 − 2⁄45 = 2⁄9, so 10⁄3; Δ_lj·Δ_jk 2⁄15, so 2; Δ_lk·Δ_lj and Δ_lk·Δ_jk −1⁄15 − 1⁄9 = −8⁄45, so −8⁄3.

**The square.** Where l and k are neighbours on a ring l → k → j → i → l: g₁ against the three-step path with one ZZ, and the direct hop with one ZZ against the three-step path, 4κ²·Im(t_lk·t_kj·t_ji·t_il)·∫(F₃₁ − F₁₁) dy, the integral being (Δ_lk − Δ_li − Δ_ij − Δ_jk)/3. Where they are opposite, the two-step path over m against the two-step path over j with one ZZ, −4κ²·Im(t_lm·t_mk·t_kj·t_jl)·F₂₁ summed over the ordered pairs (m, j); the loop changes sign when m and j are exchanged and the term Δ_lk of F₂₁ does not, so it drops out, and (Δ_lj + Δ_jk)/3 remains. With 7!·4·4κ²/3 = 336·4^N·5 these are the last two lines of Theorem B.

**A second way to one number.** At one Δ on the three edges of a triangle alone, N = 3, H_ZZ is a function of the excitation number, 3Δ for 0 or 3 excitations and −Δ for 1 or 2, and commutes with the hopping. Sorting f⁺ by the number of excitations on its two arcs, f⁺ = e^{iΔτ}·[f° + (e^{−4iΔx} − 1)·U(y) + (e^{−4iΔy} − 1)·U(x)] with f° the amplitude at Δ = 0 and U(u) the entry (k, l) of e^{−iuτ̂}, τ̂ the matrix of the couplings τ_ab (no excitation on the arc x beside one on the arc y, or two beside three). The order five of |f⁺|² − |f⁻|² at Δ² is 256·Δ²·xy·I_j at τ = 1, which gives Q = 10·Δ², and the form above gives 5 + 5 + 10⁄3 + 2 − 16⁄3 = 10.

**The degrees.** By F137's three maps every trace of a power is even in the rates, and its degree in H plus its degree in δ is even. At the seventh power the part of degree zero is Tr((−i·ad_H)⁷) = 0, the spectrum of ad_H being symmetric about zero, and the part of degree six has one term of H, which is zero by the argument F137 gives for the first order in H at the fifth power, that argument not counting the rates. Both the maps and that argument are F137's and are not shown again here. A term of degree two has five terms of H, so its degree in δ is odd: it has one s and one δ, and Theorem A reads it. That leaves degree four, three terms of H with four rates. ∎

## 7. The gate

[`simulations/seventh_moment_transport.py`](../../simulations/seventh_moment_transport.py) works with M = 2L + c and the rates s = γ↓ + γ↑, a = γ↓ − γ↑ of the fifth trace's gate, where the coefficient of s_l·a_k in Tr(M⁷) is 10752·4^N·β_lk (10752 = 336·32, the powers of two of M⁷ against the halved rates). Traces modulo ten primes joined by the Chinese remainder theorem, the degrees in the rates separated by scaling them, every comparison exact.

- **(T)** Theorem A order by order: the Taylor coefficients of the integral, a sum of products of the traces A_pq in exact rationals, against the coefficient of s_l·a_k read off the generator, at m = 5 and m = 7, every ordered pair, on six settings at N = 3 and 4: four with triangles, a ring of four with ZZ, and one without a loop, which is silent.
- **(D)** no part of degree zero or six.
- **(L)** Theorem B against the generator on a triangle with ZZ, a triangle with a hopping pendant, a triangle with ZZ only on bonds to a fourth site, a ring of four with flux and ZZ on edges and diagonals, K₄, two graphs at N = 5 and seven random graphs.
- **(V)** the gate holds Theorem B as 21 integer coefficients, those of Q and of the rings of four at three times their value and each term of the two linear forms on its own; each of the 21, moved by one, is noticed by a row.
- **(R)** the temperature gradient of F137, whose 430080 is all of degree four; a ring of five against 2520·4^N·Im(loop)·(pentagram − pentagon); the triangle with a fourth site and the two rings of four of §1's last paragraph, where the seventh trace is not zero, has no part of degree four and is given by Theorem B, one ring's being zero without its ZZ; and the triangle of §8 with three unequal weights.
- **(C)** controls: with a longitudinal field the seventh trace is no longer even in the rates; an XX term moves the part of degree two off Theorem B.

The gate checks Theorem A and Theorem B. The polynomials of §5 and §6 are the derivation and are not asserted one by one.

## 8. Open

- **The part of degree four** of the seventh trace (three terms of H, four rates). It decides the temperature gradient above and is zero on the three settings at N = 4 of (R); Theorem A does not reach it, a fourth Duhamel order would.
- **F137's open question**, whether the spectrum can be even under flux reversal at a flux other than 0 and π with ratios that differ between the sites. Three rate points on a line that misses the origin are such ratios; there the part of degree two is zero where the three weights are equal and is not identically zero where they differ, though it has zeros there too: with hoppings 1, 1 + i and 2i the weights are 98, 91 and 77, and at γ↓ = (3, 5, 2), γ↑ = γ↓ − 2 the part of degree two is zero and Tr((L + c/2)⁷) = 322560 is all of degree four. The question is not decided here.
- **Not typed.** Theorem B has no claim and no live witness in the typed layer; the fifth trace has (`FluxAreaFifthMomentClaim`).
