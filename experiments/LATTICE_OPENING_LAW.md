# The lattice opening law: the heavier sock minus the living spook

**Date:** 2026-07-16 (evening, found playing)
**Status:** verified from below. Gate: [`simulations/lattice_opening_law.py`](../simulations/lattice_opening_law.py) (N = 2, 3, all PASS, exit 0). C# pin: `LatticeTests.The_Opening_Law_Holds_On_The_Cat_Pair` in `compute/MirrorWorld.Tests`.
**Where it came from:** playing the bridged lattice (`compute/MirrorWorld/Lattice.cs`, the Klein V₄ of watchings) with the SpookyAction seed, in an interactive toy built the same evening. Tom read the curves and said the closed form approaches the lattice curve; chasing that sentence produced the general law. The first headwork guess ("the Bell seed closes the lattice for all t") was WRONG; the measured number 0.432 = (1 − e⁻²)/2 exposed the closed form.

## Setup

The bridged lattice runs four worlds side by side: e = ρ(t) under the normal rule
(rate −2Σ_l γ_l per disagreeing site), and the one-sided reading L = X^N ρ(t) under
the turned rule (rate −2Σ_l γ_l per AGREEING site). The bridge L(t)[i,j] = e(t)[~i,j]
is exact, entry for entry, at every tick; this holds for ANY site-resolved watching
profile γ_l, one-sided included (gate T0, machine zero; the wrong rule breaks at 0.5).

Seed the **cat pair** ψ(θ) = cos θ·|0…0⟩ + sin θ·|1…1⟩. This is SpookyAction's
sock-drawer skeleton: two definite records on the diagonal plus the spook, the k = N
coherence |0…0⟩⟨1…1| paying the full rate −2Γ, Γ = Σ_l γ_l
(see [SPOOKY_ACTION_TRANSLATED](../docs/quantum/SPOOKY_ACTION_TRANSLATED.md)).

Define the **lattice opening** as the entry-wise distance between the two worlds,

    opening(t) = max_ij | e(t)[i,j] − L(t)[i,j] |.

## The law

    opening(t) = max(cos²θ, sin²θ) − cos θ · sin θ · e^(−2Γt)

**In words: opening = the heavier sock's weight minus the LIVING spook.** Two
contributions, cleanly separated:

- the seed's **chirality floor** max(cos²θ, sin²θ): timeless, set at preparation,
  untouched by the watching;
- the **spook term** cos θ·sin θ·e^(−2Γt): dies under the watching at the full
  k = N rate and closes exactly the gap it owns, never more.

Corners: at θ = 45° the floor is ½ and the opening is purely the dead fraction of
the spook, (1 − e^(−2Γt))/2: at t = 0 all four lattice vertices are ONE matrix
(ψ(45°) is an X^N eigenstate), and the lattice stays TWO worlds forever (e = LR and
L = R exactly, gate T3): the Bell lattice is split only by this spook-death meter.
At θ = 0 there is no spook and the lattice stands open at 1, constant.

## Why it is exact, and blind to J

The cat sector is **H-dead**: an excitation-conserving hop needs an excitation next
to a hole, and |0…0⟩ has no excitation while |1…1⟩ has no hole, so the XY handshake
annihilates both ends at every N; the ZZ bond gives both ends the SAME diagonal
energy, so the spook collects no phase either. The e-world trajectory is therefore
pure dephasing: populations frozen at cos²θ / sin²θ, spook = cos θ·sin θ·e^(−2Γt).
The bridge then gives the opening entry-wise, and with c² + s² = 1 and c·s·k ≤ ½ the
maximum entry is max(c², s²) − c·s·k. Gate T2 pins the consequence from below: the
opening trajectories at J = 0.4 and J = 1.7, and with zz = 0.8, are bit-identical
(deviation exactly 0.0).

Because the derivation only uses "no hop can touch the cat ends" and "both ends
carry the same diagonal energy", the law is
J-free, topology-free, and holds for every watching profile γ_l through the single
number Γ = Σ_l γ_l, one more sighting of the Absorption Theorem's site-sum.

## What this is, next to what we had

The committed SpookyAction results (pages blind, spook paying −2Γ, one-sided
watching leaving the unwatched page at I/2) are all inputs here, not re-derived.
New is the LATTICE reading: the distance between a world and its one-sided mirror
reading is a meter of how much spook is still alive, with the
chirality of the seed as its floor. The turn trades spook and sock drawer (the
spook cell sits on L's and R's DIAGONAL), so the same law read from the turned side
says: the turned worlds' STRUCTURE decays toward the untouched chirality gap.

Not minted as an F-number; candidate if it earns siblings. The site-resolved turned
rule (gate T0) is the one genuinely new lattice ingredient beyond the committed
uniform-γ `Lattice.cs`; it stays gate-level here (the C# `Lattice` object keeps
uniform watching until a second consumer wants the profile).

## Files

- Gate: [`simulations/lattice_opening_law.py`](../simulations/lattice_opening_law.py)
  (T0 site-resolved turned rule + wrong-rule discriminator, T1 the law across θ for
  one- and two-sided watching, T2 J/ZZ blindness, T3 the Bell corner; N = 2, 3).
- C# pin: `The_Opening_Law_Holds_On_The_Cat_Pair` in
  `compute/MirrorWorld.Tests/LatticeTests.cs` (uniform γ, N = 3, Γ = Nγ).
- The lattice object: `compute/MirrorWorld/Lattice.cs` (the Klein V₄ of watchings,
  every edge exact; commit eca401b).
- The meaning neighbours: [SPOOKY_ACTION_TRANSLATED](../docs/quantum/SPOOKY_ACTION_TRANSLATED.md)
  (the spook, the pages, no signalling), `docs/proofs/` F1/Absorption for the −2Γ rate.
- The thread the review left behind (a longitudinal field breaks the bridge but not
  the law) is chased in [LATTICE_H_THREAD](LATTICE_H_THREAD.md): X^N as a third F131
  mirror, the mixed two-field pencil, and the doubly-mirrored zeros.
