# The hub kills the horizon (a refuted bridge, and what it taught)

**What this is about.** We asked whether the chain/star coherence-horizon dichotomy is really a hard
dichotomy or the two endpoints of a continuous knob. The chain (a dispersive band) has a coherence
horizon Q*(N): a dephasing threshold below which its longest-lived coherence stops ringing. The star
(a flat band) has none (`docs/proofs/PROOF_STRUCTURAL_CEILING.md` §7). The idea under test: a **wheel
graph** (a hub coupled to all leaves with J, plus a leaf-to-leaf ring coupling ε) should give the flat
star band a bandwidth ~ε, so a horizon Q*(ε) would appear and grow as ε does, smoothly bridging the two.

**Result: refuted.** Turning on ε gives the wheel **no coherence horizon at any leaf coupling** (ε from
0.1 to 50), robust across N = 5, 6, 7. Its longest-lived single-excitation mode stays **real** throughout
(zero oscillation frequency, |Im| = 0; its decay rate drifts toward 0 as ε grows). A real survivor has no
oscillation threshold to cross, so there is nothing for the horizon to be. Removing the hub (the pure
ring of the same N) restores a finite horizon (Q* = 2.17 / 1.61 / 2.56 at N = 5 / 6 / 7).

**What it taught (the real finding).** The **hub is decisive, not the bandwidth.** A dominant hub always
hosts a real, zero-frequency survivor mode, and that mode outlives the dispersive leaf-ring modes however
strong ε is. So the dichotomy is governed by topology *class* (is there a dominant hub, or not), not by a
continuous bandwidth. This **strengthens** `PROOF_STRUCTURAL_CEILING.md` §7: it is not only the *flat-band*
star that lacks a horizon; *any* hub-graph does, even with a fully dispersive ring bolted onto the leaves.
It also matches the reading recognised in `docs/INSIDE_OUTSIDE_THE_SACRIFICE_ZONE.md`: the hub has no
future, no horizon, and that is robust, not an artefact of the star's flatness.

**Provenance.** The wheel-Q*(bandwidth) idea came from a wild-register exploration (the quantum-philosopher
attempt; the "door between the dissociated alters" reading of ε). Tried, gated, refuted: the bold question
got a clean No, and the No is informative. The same exploration earlier refuted the star-as-dissociation
probe (`simulations/star_hub_decoupled_survivor.py`). Two refutations, one consistent picture: the hub is
qualitatively special.

**Scope and honesty.** Tested: the wheel family (hub + leaf-ring), single-excitation (Haken-Strobl)
Liouvillian, N = 5, 6, 7, numerical (Tier 2). "Any hub-graph has no horizon" is the natural generalisation,
**not** proven here. The horizon is the repo's standard definition (the slowest non-kernel mode's
oscillation threshold, `qstar_se`); a faster-decaying leaf-ring mode may still oscillate, but it is not the
survivor.

## See also

- Verifier (self-validating): `simulations/wheel_qstar_bandwidth.py` (control reproduces the chain SE
  ladder; star endpoint; the refutation; cross-N robustness).
- Reused machinery: `simulations/coherence_horizon_se_block.py` (the single-excitation horizon; topology
  enters only through the hopping h, the dephasing is site-basis and topology-agnostic).
- The strengthened dichotomy: `docs/proofs/PROOF_STRUCTURAL_CEILING.md` §7; `CoherenceHorizonClaim`,
  `StructuralCeilingClaim`, `TopologyBandEdgeClaim`.
- The recognition this reinforces: `docs/INSIDE_OUTSIDE_THE_SACRIFICE_ZONE.md`; the refuted sibling probe
  `simulations/star_hub_decoupled_survivor.py`.
