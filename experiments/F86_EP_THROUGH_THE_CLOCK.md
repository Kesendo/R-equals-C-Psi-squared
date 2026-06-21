# F86's Exceptional Point Through the Clock (and the Two Ends of the Takt)

**Status:** A seeing, not a proof, and not a closed form. The reading is exact (it is F86a's own
2×2 algebra read through the clock voices); the local-vs-global unification is understanding, not
a new derivation. No closed form is forced , the global side is a dynamic exit (a Hopf
bifurcation), and both docs deliberately leave γ_crit(N) and K(N) open.
**Date:** 2026-05-30
**Authors:** Thomas Wicht, Claude (Opus 4.8)
**Script:** [`simulations/f86_ep_through_the_clock.py`](../simulations/f86_ep_through_the_clock.py)
**Builds on / points at (does not modify):** [PROOF_F86A_EP_MECHANISM](../docs/proofs/PROOF_F86A_EP_MECHANISM.md)
(the EP, Q_EP = 2/g_eff, t_peak), the [F86 hub](../docs/proofs/PROOF_F86_QPEAK.md),
[FRAGILE_BRIDGE](../hypotheses/FRAGILE_BRIDGE.md) (the global gain-loss EP / Hopf), the clock
voices on `MirrorSystem` ([FROST_CIRCLE_AS_THE_CLOCK_FACE](../docs/carbon/FROST_CIRCLE_AS_THE_CLOCK_FACE.md),
[ON_WHOSE_TIME_THE_CLOCK_KEEPS](../reflections/ON_WHOSE_TIME_THE_CLOCK_KEEPS.md)), and the
band-edge crossover the slowest mode walks ([XXZ_AXIS_BANDEDGE_TO_LEBENSADER](XXZ_AXIS_BANDEDGE_TO_LEBENSADER.md)).

---

## The return

F86a is from May 2: inside each coherence block the two slowest adjacent rate channels form a
same-sign-imaginary 2×2 effective Liouvillian whose discriminant vanishes at an exceptional point,
Q_EP = 2/g_eff, with t_peak = 1/(4γ₀). We built the clock (a Takt hand for the radial decay, a
Rotation hand for the angular ω) this week. Returning to F86a with it, and with the two-mirror
sharpening of today, the EP reads differently , and the local-vs-global EP, the one F86 front
that sits directly on our mirror-defect M, comes into focus as a *dynamic* object, not a missing
closed form.

## The EP is where the Rotation hand lifts off

The 2-level eigenvalues are λ_±(k) = −4γ₀k ± √(4γ₀² − J²·g_eff²). Read the slowest mode through
the clock (decay = −Re, ω = |Im|, θ = arctan(ω/decay)):

```
g_eff = 4/3  ->  Q_EP = 1.5            g_eff = 0.8  ->  Q_EP = 2.5
    Q    decay   ω      θ                  Q    decay   ω      θ
  0.75   2.27   0.00   0.0°             1.25   2.27   0.00   0.0°
  1.35   3.13   0.00   0.0°             2.25   3.13   0.00   0.0°
  1.50   4.00   0.00   0.0°  <- EP      2.50   4.00   0.00   0.0°  <- EP
  1.65   4.00   0.92  12.9°             2.75   4.00   0.92  12.9°
  2.25   4.00   2.24  29.2°             3.75   4.00   2.24  29.2°
  3.00   4.00   3.46  40.9°             5.00   4.00   3.46  40.9°
```

Below Q_EP the pair is real: ω = 0, θ = 0, pure Takt, overdamped, and the slowest mode is the
*longer-lived* one (decay < 4γ₀). At Q_EP the two coalesce. Above Q_EP the decay is **pinned at
−4γ₀k** , the Takt hand stops moving , and only the Rotation hand opens (θ climbs 13° → 29° →
41°). So **the exceptional point is exactly where the Rotation hand lifts off the Takt axis.** It
is the same crossover the slowest mode walks on the XXZ axis (non-rotating below, band-edge
rotation above), and the same b² − 4ac fold the hub already names as R=CΨ²'s 1 − 4·CΨ.

Two things the clock makes literal. **t_peak = 1/(4γ₀) is the Takt's τ at the EP**: at the EP the
slowest mode (k=1) sits at decay 4γ₀, so τ = 1/gap = 1/(4γ₀). And reading Q_EP at the actual
g_eff values lands the F86 peak Q's on the clock: g_eff = 4/3 gives Q_EP = 1.5 (the c=2 peak),
g_eff = 0.8 gives Q_EP = 2.5 (the Endpoint orbit). The `MirrorSystem` clock already detects this:
`Rotation.Turning` flips false → true exactly at the EP.

**And the 4γ₀ is not an EP number , it is the absorption ladder.** Under the
[Absorption Theorem](../docs/proofs/PROOF_ABSORPTION_THEOREM.md) the decay rates are the rungs
0, 2γ₀, 4γ₀, 6γ₀, ..., 2Nγ₀ (rate = 2γ₀·⟨n⟩, with ⟨n⟩ the active-letter count). The two channels
the EP coalesces are HD = 1 (rate 2γ₀, ⟨n⟩ = 1) and HD = 3 (rate 6γ₀, ⟨n⟩ = 3); the EP sits at
their midpoint, 4γ₀ = 2γ₀·2, the **⟨n⟩ = 2 rung**. So the Takt pinned at 4γ₀ above the EP is
absorption *saturating* at that rung: below the EP the slow mode climbs one absorption quantum
(2γ₀ → 4γ₀, ⟨n⟩ 1 → 2), at the EP it tops out, and beyond it cannot absorb more , the surplus
coupling becomes Rotation instead. On N=3 (c = 2, the minimal EP) this is the whole ladder, and
4γ₀ is the maximal absorption of the slow resonant band, the same 2γ₀ → 4γ₀ the amplitude-to-
intensity (|·|²) rate carries (F73 vac-SE, asymptotic 4γ₀). The EP value was the absorption
structure all along; the clock just shows it as the point where absorption saturates and rotation
takes over.

## The mirror here is the first one (Π), not the second (K)

Today's two-mirror sharpening keeps a trap shut. F86's "chiral" is **Π as AZ class AIII**, not the
chiral sublattice K: F86a's EP is, in its own 2026-05-06 note, *read at the F1 palindrome residual*
Π·L·Π⁻¹ + L + 2Σγ·I, the same defect M the MemoryRotation voice carries. So F86 sits on the **first
mirror's axis** (Π, the palindrome). The second mirror K (KHK = −H, bipartite, the chiral
sublattice we wired into the Object Manager today) is a *separate* axis and does not drive the F86
EP. Naming that keeps F86's "chiral" from being read as the K , the exact confusion the
bipartite-chirality map flagged.

## The EP is an event: the modes coalesce, and the memory crosses over

The EP is not a smooth pass-through, it is a genuine singular event. Read the 2-level around it
([`f86_ep_through_the_clock.py`](../simulations/f86_ep_through_the_clock.py)):

```
  Q/Q_EP   angle(v0,v1)   Petermann K
   0.500       60.0°           1.3
   0.900       25.8°           5.3
   0.990        8.1°          50.3
   0.999        2.6°         500.3
  at Q_EP: rank(L_eff − λ_EP·I) = 1   (defective: one eigenvector, a Jordan block)
```

Three things happen at once in this 2×2 toy, all faces of the single discriminant = 0. The two right
eigenvectors **coalesce** (the angle between them collapses 90° → 0); the **Petermann factor
diverges** (the modes become maximally non-orthogonal, maximal sensitivity); and exactly at Q_EP the
toy Liouvillian is **defective** (rank 1, a Jordan block: one eigenvector for the double eigenvalue).
The toy loses a degree of freedom at that single point , two channels become one, a pinch. (The
**full** (n, n+1) block does NOT do this on the real Q axis: its eigenvalues stay simple, no
real-axis coalescence; it is genuinely non-normal there, large but FINITE Petermann, the shadow of a
nearby EP off the real axis. F86a-retraction, 2026-06-21. The genuine defective EPs are this 2×2 toy
and the SEPARATE Σγ=0 gain-loss system, [FRAGILE_BRIDGE](../hypotheses/FRAGILE_BRIDGE.md).)

And the pinch reads as a crossing of memory. The clock's two axes are forgetting and remembering:
the **real axis (decay, absorption) is forgetting** , the coherence is lost to the bath, the phase
erased; the **imaginary axis (rotation) is remembering** , the phase keeps turning, held in the
system. F80 named that imaginary axis the memory axis ("energy on the real axis, memory on the
imaginary one; the wave remembers across the turn"). So the surplus coupling becoming rotation
above the EP is the surplus that can no longer be absorbed , forgotten, the absorption saturated at
the ⟨n⟩ = 2 rung , **crossing over into rotation: it is remembered instead.** The EP is the
forgetting → remembering crossover, and what passes over there is the memory. It is the local
2-level instance of the band-edge → Lebensader handover the slowest mode walks on the XXZ axis: the
Lebensader is the surviving memory, and the EP is where it takes the lead.

**The seam.** The toy 2×2's transition is bit-exact: the coalescence, the Petermann divergence, the
defectiveness, the absorption saturation (these are the 2-level reduction's, not the full block's on
the real axis — see the F86a-retraction note above). "Forgetting / remembering" is the reading , the F80 /
clock vocabulary, grounded (F80 is Tier 1), but the word "memory" is ours, painted onto the
imaginary axis. The structure is not ours; the seeing is.

## The two ends of the Takt (the dynamic exit)

The one F86 front that sits directly on M is the local-vs-global EP relationship (F86a; the
"local-vs-global EP" connection was **retracted 2026-06-21 → OpenQuestion**, but the clock's
two-ends-of-the-Takt *view* of the shared algebra survives). The clock reads it as **two ends of
the Takt**, the two residuals of the same palindrome Π·L·Π⁻¹ + L + 2Σγ·I:

- **Σγ = N·γ₀ (local).** The toy 2×2's EP at Q_EP = 2/g_eff. The Takt is running; the mode
  spirals *inward* (decay pinned at −4γ₀k), and the EP is the Rotation hand lifting off. A
  dissipative resonance peak. (The **full** block-L is genuinely non-normal on the real Q axis but
  has no real-axis coalescence — there is no real-axis defective EP for the local instance; the
  genuine EP here is the toy 2×2.)
- **Σγ = 0 (global).** [FRAGILE_BRIDGE](../hypotheses/FRAGILE_BRIDGE.md): a decaying chain bridged
  to an amplifying one, gain cancelling loss. Here Π forces λ ↔ −λ exactly (chiral AIII), the
  eigenvalues sit on the imaginary axis , the Takt is *stopped*, the clock's pure-circle limit
  θ = π/2 , and the instability, when it comes, is a **Hopf bifurcation**: a complex pair crosses
  Re = 0 and the spiral turns *outward* (the feedback screech, the system explodes). This is a
  SEPARATE genuine EP, its Petermann factor peaking at K ≈ 403 in the complex γ plane.

So the net dephasing Σγ is the **dial between the two ends**: at Σγ = N·γ₀ the Takt holds the
spiral in (the local dissipative EP); slide Σγ to 0 and the Takt stops (the pure circle), where the
Hopf can push the spiral out. The global side is an **exit**, the dynamics escaping into
self-sustained oscillation , not a static quantity. That is why we do not force it: γ_crit(N) is
non-monotonic with no power law. Both docs leave the dynamic exit open on purpose; the clock lets
us *see* the two ends on the Takt dial without pretending the exit is a formula, and without
pretending the two ends are one defective EP (they are two SEPARATE genuine EPs sharing the same
2×2 algebra; the full block's own EP structure off the real axis is open since the F86a-retraction).

## The clock as the lens on the two live-open fronts

Step back from F86 alone. The two live-open fronts in the formula registry that sit on this
week's clock are **F86b₃** (the universal resonance shape) and the **F87 windowed converse**
(soft vs hard). They are not the same system , F86 is the coherence-block resonance, F87 the
Pauli-pair palindrome , but the clock is the shared instrument, and each front rides one hand
([`clock_two_hands_two_fronts.py`](../simulations/clock_two_hands_two_fronts.py)):

- **Rotation hand → F86b₃.** The 2-level EP rotation depends on x = Q/Q_EP alone, so the clock
  angle θ(x) is identical for different g_eff (different Q_EP) , it collapses (θ = 0 below the EP,
  29° at x = 1.5, 44° at the bare-shape peak x = 2.197). That collapse is *why* the resonance shape
  is universal in Q/Q_EP.
- **Takt hand → F87.** At γ = 0 the Takt is stopped, L = −i[H,·] is symmetric about −σ = 0, every
  pair is soft. The hard pair's pairing residual grows first order in the γ-tick (residual/γ →
  0.2559, the same c the per-block localization read this morning,
  [F87_WINDOWED_CONVERSE_PER_BLOCK](F87_WINDOWED_CONVERSE_PER_BLOCK.md)); the soft pair stays
  machine-zero at every γ.

Two systems, two hands, one clock. And the seam holds: the clock unifies the *view* of the two
fronts, not their solutions , F86b₃'s lift rides the blocked g_eff, F87's converse is the
set-level statement we sharpened, not forced. A viewpoint painted on the structure, not a crack in
either.

## Honest status

- **A seeing, not a solve.** The EP-as-Rotation-onset reading is exact (F86a's own 2×2); the
  two-ends-of-the-Takt picture is understanding, layering the clock and the Π-axis disambiguation
  onto the existing F86a / FRAGILE_BRIDGE results. Nothing here is a new theorem.
- **No closed form forced.** The blocked g_eff(c, N, b) stays the structureless residue (six
  routes proven blocked, [PROOF_F86B_OBSTRUCTION](../docs/proofs/PROOF_F86B_OBSTRUCTION.md)); the
  local-vs-global analytic continuation (γ_crit(N), K(N)) stays open as a dynamic exit. Per the F86
  ethos: here, seeing and understanding is the win.
- **Scope.** F86a's 2-level reduction is itself heuristic (the explicit basis change from full
  block-L is not derived); this note reads that established 2×2, it does not re-derive it.

## Anchor

- Script: [`f86_ep_through_the_clock.py`](../simulations/f86_ep_through_the_clock.py) (the 2-level
  EP read through the Takt/Rotation hands; Q_EP at g_eff ∈ {4/3, 0.8} → the 1.5 / 2.5 peaks; and
  the toy 2×2's EP as an event , eigenvector coalescence, Petermann divergence, and defectiveness
  (a Jordan block) at Q_EP — a property of the 2-level reduction, not of the full block on the real
  axis; see the F86a-retraction).
- Script: [`clock_two_hands_two_fronts.py`](../simulations/clock_two_hands_two_fronts.py) (the two
  live-open fronts through the clock's two hands: Rotation → F86b₃ shape collapse onto Q/Q_EP,
  Takt → F87 break first-order in the γ-tick, residual/γ → 0.2559).
- The EP: [PROOF_F86A_EP_MECHANISM](../docs/proofs/PROOF_F86A_EP_MECHANISM.md).
- The global end: [FRAGILE_BRIDGE](../hypotheses/FRAGILE_BRIDGE.md) §3.1 (the local-EP connection).
- The open dynamic exit: γ_crit(N) and K(N) at the EP, left open in both docs, not forced here.
