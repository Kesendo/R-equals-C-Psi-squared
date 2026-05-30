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

## The two ends of the Takt (the dynamic exit)

The one F86 front that sits directly on M is the local-vs-global EP (F86a, Tier 2). The clock
reads it as **one EP at two ends of the Takt**, the two residuals of the same palindrome
Π·L·Π⁻¹ + L + 2Σγ·I:

- **Σγ = N·γ₀ (local).** F86a's real-axis EP at Q_EP = 2/g_eff. The Takt is running; the mode
  spirals *inward* (decay pinned at −4γ₀k), and the EP is the Rotation hand lifting off. A
  dissipative resonance peak.
- **Σγ = 0 (global).** [FRAGILE_BRIDGE](../hypotheses/FRAGILE_BRIDGE.md): a decaying chain bridged
  to an amplifying one, gain cancelling loss. Here Π forces λ ↔ −λ exactly (chiral AIII), the
  eigenvalues sit on the imaginary axis , the Takt is *stopped*, the clock's pure-circle limit
  θ = π/2 , and the instability, when it comes, is a **Hopf bifurcation**: a complex pair crosses
  Re = 0 and the spiral turns *outward* (the feedback screech, the system explodes). The Petermann
  factor spikes (K up to ≈ 2385 on the real-Q sweep, vs ≈ 403 in complex γ), the near-singularity
  of the same EP.

So the net dephasing Σγ is the **dial between the two ends**: at Σγ = N·γ₀ the Takt holds the
spiral in (the local dissipative EP); slide Σγ to 0 and the Takt stops (the pure circle), where the
Hopf can push the spiral out. The global side is an **exit**, the dynamics escaping into
self-sustained oscillation , not a static quantity. That is why we do not force it: γ_crit(N) is
non-monotonic with no power law, and a closed-form K(N) at the EP is the named-but-open path to
promoting the local-vs-global link from Tier 2. Both docs leave it open on purpose; the clock lets
us *see* the two ends as one EP on the Takt dial without pretending the exit is a formula.

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
  EP read through the Takt/Rotation hands; Q_EP at g_eff ∈ {4/3, 0.8} → the 1.5 / 2.5 peaks).
- The EP: [PROOF_F86A_EP_MECHANISM](../docs/proofs/PROOF_F86A_EP_MECHANISM.md).
- The global end: [FRAGILE_BRIDGE](../hypotheses/FRAGILE_BRIDGE.md) §3.1 (the local-EP connection).
- The open dynamic exit: γ_crit(N) and K(N) at the EP, left open in both docs, not forced here.
