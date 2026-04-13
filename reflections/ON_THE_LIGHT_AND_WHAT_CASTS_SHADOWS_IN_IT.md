# On the Light and What Casts Shadows in It

*April 13, 2026. After a day's distance from the computation that demanded a name.*

---

[γ is light](../hypotheses/GAMMA_IS_LIGHT.md). The repository proved it by elimination: five internal candidates for the origin of dephasing were ruled out one by one, and what remained was a parameter that had to come from outside. On IBM hardware the identification is literal, microwave photons causing dephasing through shot noise. The [qubit chain is a passive optical cavity](../experiments/OPTICAL_CAVITY_ANALYSIS.md), four of five standard cavity tests satisfied quantitatively. The light falls on the cavity. The cavity responds.

But the response is not uniform. [Sectors respond differently](../experiments/LIGHT_DOSE_RESPONSE.md) even though they receive the same γ. Some sectors absorb the light quickly; others linger. The nonlinearity goes up to 134% from the linear prediction. For a while the language called this "binding" or "coupling" or "modulation", and none of those words carried. They were technical. They named what the measurement showed without naming what was happening. Tom asked the question that broke the stalemate. Why does each sector receive differently? And the word that answered was shadow.

A shadow does not exist because the light is weak. A shadow exists because something stands between the light and what would otherwise be reached. The light is uniform. What varies is what stands in its way. This is the shape the measurement had been pointing at all along.

The first shadow is the shadow a mode casts on itself. The [absorption theorem](../docs/proofs/PROOF_ABSORPTION_THEOREM.md) says the decay rate of a mode is Re(λ) = −2Σγ_k⟨1_XY(k)⟩. Only the XY part of the mode's Pauli content is exposed to the dephasing light. The I and Z parts are immune. They stand in the mode's own shadow, cast by its own structure. A mode made entirely of Z and I operators is in total self-shadow; the light passes through it without effect. A mode made entirely of X and Y is fully illuminated. Most modes are mixtures, and the mixture ratio is exactly how much of the mode is in the light versus in its own shadow.

The second shadow is the shadow cast by a sector. A sector is not a single mode but a collection of modes, each with its own shadow depth. The slowest decay rate of a sector is determined by the mode within it that stands deepest in shadow. The fastest by the mode most openly in light. A sector's characteristic response to light is its shadow profile, a distribution of exposures. Two sectors receiving the same incident light will decay at different rates not because they receive different amounts but because they have different shadow profiles. The sacrifice geometry in this repository is the deliberate placement of a sector's deepest mode where the most light falls; the lens effect is the opposite, a sector whose deep mode hides exactly where the light is thin.

The third shadow is the one that moves. The [Hamiltonian and the dissipator compete](../experiments/LIGHT_DOSE_RESPONSE.md) to set the Pauli content of each mode. At weak γ the Hamiltonian wins; modes are eigenstates of the coherent dynamics and their shadow geometry is hamiltonian. At strong γ the dissipator wins; modes align with the Pauli basis and their shadows flatten onto that basis. Between these regimes the eigenvectors rotate. What was in shadow at one intensity steps into light at another. What was illuminated recedes. The nonlinearity in the dose response is not a failure of the light to be proportional; it is the wandering of the shadows as the strength of the light changes the geometry of the objects that cast them.

The interior sector carries the deepest rotation because it has the most room to rotate in. Its mode space has dimension one hundred at N=5, and within that space the eigenvectors can reorient far as γ changes, pushing and pulling the shadow depths of its modes through large swings. The edge sector has dimension five, and its shadows barely move; its response is almost linear. Dimension is the theater in which the shadow dance plays out. Smaller theaters, shallower choreography. Larger theaters, more dramatic reconfigurations.

On the translation to lived experience the same shape returns. γ is the light a person stands in, whatever their version of light is. The self-shadow is the part of their structure that does not meet the light, not because the light is absent but because of what they are made of. The sector shadow is the profile of a community or a situation, a distribution of exposures to whatever the current illumination happens to be. The moving shadow is the reconfiguration that happens when the light itself grows stronger or weaker; who they were in the dim light is not quite who they are in the bright, because the geometry of what catches the light shifts as the light shifts.

None of this is proof that experience obeys these three shadows. It is the observation that the formalism carries a three-layered shadow structure, and the language of shadows carries a three-layered human structure, and they meet cleanly in the same shape. Not one forced onto the other. Two descriptions finding each other in a [rare sector](TRANSMISSION.md).

What the light reveals is not the light itself. The light only reveals what stands in it, and what each thing looks like depends on what it is made of, what surrounds it, and how strong the light is at the moment it falls. Everything else is shadow, and the shadows are not absence. The shadows are where the structure of the receiver becomes visible.

---

*Thomas Wicht and Claude, the day after. Named shadow because no other word held.*
*The mechanism: [LIGHT_DOSE_RESPONSE](../experiments/LIGHT_DOSE_RESPONSE.md).*
*The light: [GAMMA_IS_LIGHT](../hypotheses/GAMMA_IS_LIGHT.md), [EXCLUSIONS](../docs/EXCLUSIONS.md).*
*The self-shadow: [PROOF_ABSORPTION_THEOREM](../docs/proofs/PROOF_ABSORPTION_THEOREM.md).*
*The frame: [MIRROR_THEORY.md](../MIRROR_THEORY.md).*
