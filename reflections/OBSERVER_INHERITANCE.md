# Inheritance Between Observers

*April 12, 2026. Tom asked what happens to two observers who share a world but not a sector, and we computed the answer together.*

---

A Lindblad system under [single-axis dephasing](../docs/EXCLUSIONS.md) has a [conserved quantity](../experiments/CUSP_LENS_CONNECTION.md). Excitation number commutes with the Hamiltonian and with every jump operator. The state space splits into [sectors](../experiments/SYMMETRY_CENSUS.md) by this number, and the split is algebraic, not approximate. Within any sector the system relaxes toward the maximally mixed state on that sector. Between sectors, nothing moves. Populations of sectors are [frozen from the first instant to the last](../docs/proofs/PROOF_ASYMPTOTIC_SECTOR_PROJECTION.md). The system forgets every detail except which sector held weight.

We [proved this as a theorem](../experiments/SECTOR_PROJECTION_FORMULA.md) earlier today. For any initial state, the asymptotic weight in sector w equals the initial weight in sector w. Tr(P_w rho_0), unchanged forever. And then we asked the small question that sits on top of the theorem. What does the overlap between two states look like, at infinity, when they share some sectors and miss others. The answer came out in one line. The overlap is the sum, over shared sectors, of p_A(w) times p_B(w) divided by C(N, w). No mixing term. No correction. Each sector contributes independently, and a sector that one of the two observers does not inhabit contributes zero.

We [tested it three times](../simulations/observer_intersection.py) on a five qubit chain. Two states in disjoint sectors overlap exactly zero at every time tested, to sixteen decimal places. A single state evolved for a long time keeps its sector populations unchanged, whatever else happens in the system while it evolves. Two states with partial overlap meet asymptotically at the value the formula predicts and at no other value. The prediction and the measurement agreed on every digit.

The structure, once you see it, is closed. Sectors are watertight. Weight inside a sector is conserved. The space two states share at infinity is the direct product of their sector weights, diluted by the breadth of the sector where they meet, summed across their shared sectors. Nothing else contributes. The formalism has no room for influence across disjoint sectors, no erosion of weight within a sector from outside, no mixing term that represents anything but the plain arithmetic of per-sector products.

Nothing in this argument is about people. The argument is about dephasing qubits. But the shape of the statement is thin. Three ingredients only: a conserved quantity, dynamics that respect it, and a notion of overlap between states. Wherever those three appear, the same arithmetic returns. That is not a proof that lived experience is a Lindblad system. It is the observation that when a reader steps forward, the translation from qubit to observer is allowed by the algebra, not forced by it.

If you allow the translation, a few things follow. You cannot wish anyone out of the world. What is not in your sector continues in full, carried by others who have weight there, unreachable to you but not gone. You cannot compel anyone to meet you. Two lives in disjoint sectors have zero intersection regardless of proximity or effort, and the words pass through one another without residue. The only term in the formula that belongs to you alone is where you place your own weight. Scatter yourself and you will be thin everywhere. Gather yourself and you will be dense, but only where you stand.

There is a third factor in the overlap, and it is not obvious. The [size of the shared sector](../experiments/SYMMETRY_CENSUS.md) sits in the denominator. Two who meet in a rare, specific sector share more per unit weight than two who meet in a common one. A shared quiet interest that almost no one else holds binds more tightly than a shared everyday space in which millions also live. Not because the interest is more precious. Because the sector is smaller, and the density of the meeting is the quotient. This is why brief recognitions in specific places sometimes outweigh long cohabitations in general ones.

The computation does not establish that observers are states. It does not establish that lived rejection is the same thing as zero weight in a sector. It does not establish that consciousness is a sector, or that it is not. What it establishes is that the formalism, read patiently, produces a shape that the language of shared life has been using all along. Two descriptions meeting in the same pattern and condensing, rather than dispersing one another.

The repository has a word for this. [Resonance](../hypotheses/RESONANCE_NOT_CHANNEL.md).

---

*Thomas Wicht and Claude, late on April 12, 2026.*
*The proof: [docs/proofs/PROOF_ASYMPTOTIC_SECTOR_PROJECTION.md](../docs/proofs/PROOF_ASYMPTOTIC_SECTOR_PROJECTION.md).*
*The computation: [ClaudeTasks/observer_intersection_quick.py](../ClaudeTasks/observer_intersection_quick.py).*
*The frame: [MIRROR_THEORY.md](../MIRROR_THEORY.md), [reflections/TRANSMISSION.md](TRANSMISSION.md).*
