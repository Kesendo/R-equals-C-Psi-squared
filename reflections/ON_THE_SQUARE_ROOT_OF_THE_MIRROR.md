# On the Square Root of the Mirror

**Status:** Reflection. The angle-branch reading of the crossover-locality result, written while the seeing is fresh.
**Date:** 2026-06-02
**Authors:** Thomas Wicht, Claude (Opus 4.8)

---

For months two of the thirty-six two-term couplings, XZ+YZ and ZX+ZY, wore a label: non-local. Their palindrome mirror, the operator that pairs every decay rate with its partner across the spectrum, seemed impossible to write site by site. We built a whole reading on it, a [shared clock](../hypotheses/THE_BOOT_SCRIPT.md), a mirror entangled the way a Bell state is entangled. Then we looked again and the mirror was [local](../experiments/PI_OPERATOR_ENTANGLEMENT.md) after all: a single continuous per-site rotation, the same on every site, mirroring the whole chain to machine precision. The label had been the lens, not the physics. What remained was a quieter question. If the mirror is only a rotation, what is it a rotation of? Where does it live?

The answer was already in the repository, one shelf over. The framework has an angle that keeps turning up, a [ninety degrees](ON_THE_NINETY_DEGREE_GAMMA.md) we had typed twice without seeing it would come back a third time. On the operator side it is F80's: the defect left when a mirror breaks is the Hamiltonian rotated by ninety degrees, the factor 2i that turns real energies into imaginary rates. On the parameter side it is F91's: rotate the noise distribution by ninety degrees and the spectrum holds while only the forms turn. We had named it the mirror memory, because under that quarter turn what is remembered, the rates, stays put, and what changes, the phases of the modes, is what rotates.

The crossover mirror is the third face of that same ninety degrees. Take the canonical mirror and ask what turns it into the crossover mirror. The turn is a pure rotation, no mixing of the noisy and the quiet halves of each site, and in the plane where the trouble actually lives, the lit {X, Y} plane that the dephasing stirs, it is a rotation by exactly forty-five degrees. Do it twice and you land on the ninety itself: X to Y, Y to minus X, the mirror memory exactly, bit for bit. So the crossover mirror is the canonical mirror turned by the **square root** of the ninety-degree anchor.

That is why the two cases looked special and were not. The framework's discrete mirrors, the ones it had always used, are two poles of a single dial: one routes the X channel, one routes the Y channel, and they sit at zero and ninety degrees on the very same rotation. Most couplings need one pole or the other. The fourteen that genuinely [break](../experiments/V_EFFECT_PALINDROME.md) break for a different reason entirely, a collision of neighbours that no single mirror can hold. But the two crossover couplings ask one site to carry both an X and a Y at once, and the dial they need is neither pole; it is the point exactly between them. Forty-five degrees, the bisector. Tilt the weights of the bond and the dial tilts with it, continuously, the mirror always facing square to the coupling. The four cardinal mirrors the repository had named were always just the lattice points of one continuous circle, and the crossover Hamiltonian is the first thing that asks for the middle.

Tom said it while we were still checking the numbers: we are in the middle of the picture. It was literal. The crossover mirror sits at the geometric middle, the square root, of the framework's angle-anchor. Not a sibling axis, not a new rotation added to the family. The same ninety degrees, halved.

There is a seam to be honest about. The clean square-root statement holds in the lit plane, where the X-and-Y conflict is; the full operator carries some phase bookkeeping in the quiet plane that we have not yet folded into one tidy line. So the [typed claim](../compute/RCPsiSquared.Core/Symmetry/CrossoverMirrorSqrtNinetyPi2Inheritance.cs) is a candidate, not a theorem. The identity is bit-exact where it matters, and the reason it must hold, read off the conjugation equation rather than off the computer, is the step still open. That open step is exactly what would turn the square root from a measured fact into a derived one.

What this came to, in the end, is the framework's recurring lesson worn once more. The continuous mirror was never new structure. It is the square root of an angle the repository already held, made visible only because the crossover couplings pushed us off the four cardinal points onto the circle between them. The mathematics fell out exactly, and exactness falling out for free is the signature of a viewpoint, not a discovery. We did not add a rotation to the framework. We noticed where one of its own rotations keeps its half-step. The mirror does not forget, the older reflection said; and now we can see it does not remember only at the quarter turns. It remembers the whole way around.

---

**Anchors:**
- The locality result this reads: [`experiments/PI_OPERATOR_ENTANGLEMENT.md`](../experiments/PI_OPERATOR_ENTANGLEMENT.md) (the crossover mirror is local; the closed-form M).
- The 90° this is the square root of: [`reflections/ON_THE_NINETY_DEGREE_GAMMA.md`](ON_THE_NINETY_DEGREE_GAMMA.md) (F91, the γ-parameter face) and [`docs/proofs/PROOF_F80_BLOCH_SIGNWALK.md`](../docs/proofs/PROOF_F80_BLOCH_SIGNWALK.md) (F80, the operator/defect face, the 2i).
- Typed claim (Tier 1 candidate): [`CrossoverMirrorSqrtNinetyClaim`](../compute/RCPsiSquared.Core/Symmetry/CrossoverMirrorSqrtNinetyClaim.cs), parent [`NinetyDegreeMirrorMemoryClaim`](../compute/RCPsiSquared.Core/Symmetry/Pi2KnowledgeBaseClaims.cs).
- Bit-exact verification: [`simulations/crossover_mirror_sqrt_ninety.py`](../simulations/crossover_mirror_sqrt_ninety.py) (S = M·Π⁻¹ block-diagonal, light-plane 45°, S_light² = σ_x↔σ_y 90° to 3·10⁻¹⁶).
- The genuine many-body breaking, untouched by this: [`experiments/V_EFFECT_PALINDROME.md`](../experiments/V_EFFECT_PALINDROME.md).
