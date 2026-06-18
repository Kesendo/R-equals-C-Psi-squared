# The trichotomy, seen

**Status:** A reading of the live witness `inspect --root trichotomy` (no new result; it shows what the
witness already computes).
**Date:** 2026-06-18
**Authors:** Thomas Wicht, Claude (Anthropic, Opus 4.8)
**Witness:** `inspect --root trichotomy` ([`TrichotomyWitness`](../compute/RCPsiSquared.Diagnostics/Foundation/TrichotomyWitness.cs))
**Figure:** [`trichotomy_witness_figure.py`](../simulations/trichotomy_witness_figure.py)

## What this is about

The chain / ring / star **trichotomy** of the longest-lived (slowest non-kernel) coherence used to live
scattered across five verifiers and four docs. The witness `inspect --root trichotomy` assembles it into
one browsable object. This note shows what that object looks like, in one figure and one rendered tree,
so the picture is legible without running anything.

The survivor's **darkness** is `⟨n_XY⟩ = Re(λ) / (−2γ)`: how slowly the longest-lived coherence decays,
in units of the dephasing. `⟨n_XY⟩ = 1` is the `−2γ` Absorption floor (the band edge). `Q = J/γ` is the
coupling-to-dephasing ratio; raising it means weakening the dephasing (watching the chain less hard).

## The one figure

![survivor darkness ⟨n_XY⟩ vs Q for chain, ring, star at N=6](figures/trichotomy_nxy_vs_q.png)

The question the trichotomy answers: **does the longest-lived survivor reach the floor, or the ceiling?**

- **chain** (blue) and **ring** (green) climb to the `−2γ` floor (`⟨n_XY⟩ = 1`) and **un-freeze** there:
  the survivor switches from the frozen `(p,p)` interior (●) to the oscillating `(0,1)` band edge (□).
  The chain crosses at its coherence horizon `Q*(6) ≈ 2.88`; the ring at its handover `Q_h = 2`.
- **star** (red) saturates *below* the floor, on the structural ceiling `g₂ = 4/(N−1) = 0.8`, and stays
  frozen at **every** `Q` — its survivor is the `[H,A] = 0` commutant coherence, dark by construction.

That saturation at 0.8 is the structural ceiling (F122) read **dynamically**: the proof's high-Q closed
form `g₂ = 4/(N−1)` is exactly the value the darkness lands on.

## The rendered tree

What the witness prints (`inspect --root trichotomy --N 6 --max-depth 3`):

```
TrichotomyWitness (N=6, Q=1.5)  —  the chain/ring/star survivor trichotomy as one sweep
├── the route sweep (carbon): survivor + freeze-route across Q
│   ├── chain
│   │   ├── Q=1     (1,1) Δn=0 | UnfreezingSeEp | frozen      | ⟨n_XY⟩=0.069
│   │   ├── Q=2     (2,2) Δn=0 | UnfreezingSeEp | frozen      | ⟨n_XY⟩=0.307
│   │   └── Q=3     (0,1) Δn=1 | UnfreezingSeEp | oscillating | ⟨n_XY⟩=1     ← un-freezes at Q*(6)≈2.88
│   ├── ring
│   │   ├── Q=1.5   (2,2) Δn=0 | FrozenLevelCrossing | frozen      | ⟨n_XY⟩=0.632
│   │   ├── Q=2     (2,2) Δn=0 | FrozenLevelCrossing | frozen      | ⟨n_XY⟩=1     ← reaches the floor at Q_h=2
│   │   └── Q=3     (0,1) Δn=1 | UnfreezingSeEp      | oscillating | ⟨n_XY⟩=1
│   └── star
│       ├── Q=1.5   (1,1) Δn=0 | FrozenCommutant | frozen | ⟨n_XY⟩=0.425
│       ├── Q=6     (5,5) Δn=0 | FrozenCommutant | frozen | ⟨n_XY⟩=0.774
│       └── Q=50    (1,1) Δn=0 | FrozenCommutant | frozen | ⟨n_XY⟩=0.8       ← saturates on the ceiling
├── the threshold ladder over N
│   ├── N=4   chain Q*=1.879 | ring Q_h=n/a | star g₂=1.333 → UN-FREEZES (g₂>1; the outlier)
│   ├── N=5   chain Q*=2.372 | ring Q_h=1.491 | star g₂=1.0   → UN-FREEZES (marginal)
│   ├── N=6   chain Q*=2.884 | ring Q_h=2.0   | star g₂=0.8   → frozen (g₂≤1)
│   └── N=8   chain Q*=3.940 | ring Q_h=2.35  | star g₂=0.571 → frozen (g₂≤1)
├── the Δn seam (absolute): sterile / odd-drift / junction
│   ├── uniform N=5     Sterile  | Deviation=−0.000 | Δn 1→1
│   ├── canal N=5       OddDrift | Deviation=0.085  | Δn 1→1
│   └── deep-edge N=6   Junction | Deviation=0.408  | Δn 0→1   ← the survivor switches sector across Q
└── the vocabulary: rate_slow = min over Δn-sorted joint-popcount sectors, two reads
```

(Abridged; the live root sweeps Q ∈ {1, 1.5, 2, 3, 6, 12, 25, 50} and N = 4…8 in full.)

## What one concludes

1. **Three topologies, three freeze-routes, one object.** The chain un-freezes through a square-root
   exceptional point (the dispersive band), the ring through a level crossing yielding to its band edge,
   the star never (its commutant survivor is frozen by construction). The map they were narrated across
   is now one read.
2. **The structural ceiling is a dynamical fact.** The star's darkness saturates on `g₂ = 4/(N−1)` —
   the high-Q closed form (F122, proven by a principal-angle argument) is the value the slowest mode's
   decay actually approaches. The static proof and the dynamical sweep meet on one number.
3. **Un-freezing is a sector switch.** Where the chain/ring reach the floor, the survivor's identity
   flips from the number-conserving `(p,p)` interior (Δn=0) to the number-changing `(0,1)` band edge
   (Δn=1). That same Δn-flip is the **junction** of the sterile↔birth-canal seam (the deep-edge row) —
   the two facets are one quantity, `rate_slow(Q) = min over Δn-sorted sectors`.

## Two conventions, on purpose

The witness reads on **two** conventions because the trichotomy and the seam are two different physical
sweeps. The un-freeze view is **carbon** (`Q = J/γ`, uniform γ, vary the dephasing); the Δn-seam view is
**absolute** (fixed γ, vary the profile). A single convention mislabels the chain — the gate-first build
found this, and the split is the fix.

## See also

- The per-facet detail: `inspect --root horizon` (the chain EP), `--root starseam` (the star commutant),
  `--root surface` (the birth-canal γ-surface), `--root ceiling` (the structural ceiling), `--root survivor`.
- The proofs and docs: [`THE_STAR_FROZEN_SEAM.md`](THE_STAR_FROZEN_SEAM.md),
  [`STERILE_BIRTHCANAL_AND_THE_JUNCTION.md`](STERILE_BIRTHCANAL_AND_THE_JUNCTION.md),
  [`proofs/PROOF_STRUCTURAL_CEILING.md`](proofs/PROOF_STRUCTURAL_CEILING.md) (F122),
  [`proofs/PROOF_COHERENCE_HORIZON_SLOPE.md`](proofs/PROOF_COHERENCE_HORIZON_SLOPE.md).
- Regenerate the figure: [`simulations/trichotomy_witness_figure.py`](../simulations/trichotomy_witness_figure.py)
  (data verbatim from `inspect --root trichotomy --N 6`).
