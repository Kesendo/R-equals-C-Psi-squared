# On the Soft Break

---

*A reflection on what our framework sees that earlier readings could not.*

---

There are 36 ways to write a two-term bilinear Hamiltonian on a qubit chain. Heisenberg saw them through his matrices in 1925: any sum of σ_a σ_b coupling. Pauli saw them through his algebra: 9 essentially distinct two-body operators (XX, XY, ..., ZZ), C(9,2) = 36 unordered pairs to combine. They could not, with hand calculation in their era, tell which 36 combinations would produce a palindromic spectrum on a multi-qubit chain. That question was not even formulated in their language.

In March 2026, the V-Effect asked it computationally and answered: **14 break, 22 do not**. The test was eigenvalue pairing: for each eigenvalue λ of the Liouvillian, find a partner λ' such that λ + λ' = −2Σγ. Either every eigenvalue has a partner (the 22) or some are orphaned (the 14). This is what spectroscopy would see if you measured the system's transition frequencies.

In April 2026, we built `framework.py` and asked the same question, but with a finer test: not just eigenvalues, but the operator equation itself. Π·L·Π⁻¹ + L + 2Σγ·I, is this matrix exactly zero, or is it nonzero somewhere?

The answer: **3 are truly zero** (the framework's both-parity-even Heisenberg/XXZ subset). **33 are nonzero**. And among those 33: **14 also break eigenvalue pairing** (the V-Effect's hard breaks), **19 do not**.

The 19 are the soft breaks. They are systems where the operator equation has a residual matrix M of non-trivial magnitude (norm ~22 to 45, similar to the hard breaks), yet the eigenvalues of L still pair to within machine precision. Spectroscopy would call them palindromic. Operator-level analysis against Π says no.

---

What we see, and earlier readings could not:

The eigenvalue spectrum is the **diagonal** of L's structure in its eigenbasis. The operator equation tests the **full matrix**, every entry. A matrix can have non-trivial off-diagonal entries (mixing eigenstates) without shifting the eigenvalues themselves.

In the soft-break cases the residual matrix M is nonzero, but in a way that does not shift L's eigenvalue pairs. The spectral pairing is a property of L's eigenvalue distribution; the operator equation is about how Π and L compose. A matrix M can fail to be the zero matrix while leaving eigenvalues alone.

This is invisible to spectroscopy. A frequency comb on a soft-break system would show the same palindromic pattern as a truly unbroken system. Whether any measurement sees the soft category itself is a further question, and for nine of the 19 the answer is no (below, and in the Klein routing note).

---

This is what the framework gives us that the language of 1925 did not: a **super-operator** view of the dynamics.

Heisenberg's matrices describe states and observables. Pauli's algebra describes operators acting on states. Both are at the level of the Hilbert space H. The Liouvillian L is one level above: it acts on operators, taking ρ to L(ρ). And the conjugation Π·L·Π⁻¹ is one level above THAT: a transformation of the dynamics itself, a super-operator on operators. Heisenberg in 1925 had no language for this; Pauli's algebra reaches operators-on-states, not operators-on-operators-on-states. The Lindblad framework (1976) named the L; the palindromic mirror (proved 2026 in this project) named the Π·L·Π⁻¹ relation; the C²⊗C² parity decomposition (also 2026) named the resolution into bit_a and bit_b sectors. Together these compose the vocabulary in which the operator-versus-spectral distinction is statable.

Heisenberg could not have asked: "is Π·L·Π⁻¹ = −L − 2Σγ·I exactly?" Pauli could not have asked: "or only weakly, on the eigenvalues?" Both questions require the super-operator level to even be possible questions.

The V-Effect asked the spectral form of the question. It got 14/22. The 22 looked the same: palindromic.

Our framework asks the **super-operator form**. It gets 3/33 with a 14/19 sub-split of the 33. The 22 is now visibly two distinct phenomena: 3 truly palindromic at the super-operator level, 19 whose spectrum pairs while Π's equation fails.

---

What is hidden in the 19?

A residual matrix M, of comparable magnitude to the 14 hard-break residuals, that does not propagate to the eigenvalue level. The break is real against Π; it just does not register as eigenvalue motion.

A complete classification of *why* each of the 19 spectral-paired cases passes the spectral test is open work. What we can say with care:

- Each soft-break case has *some* structural reason. For example, XY+YX (matched bit_b violations on both bonds) is anti-invariant under per-site σ_x conjugation: U·(XY+YX)·U⁻¹ = −(XY+YX) for U = X⊗X. This anti-symmetry is independent of Π. Under U-conjugation, L_H = −i[H, ·] flips sign while the dephasing dissipator is preserved. Whether this composes with Π to give an exact anti-involution of the *full* L (and hence spectrum pairing) requires more careful operator algebra than we have done.

- For one-good-one-bad combos like XX+XZ (the (ab)+(──) pattern): the good term XX carries the canonical Π-symmetry; the bad term XZ contributes a residual whose spectral consequence is non-trivial to compute by inspection. Empirically the spectrum still pairs.

- For the 14 hard breaks, no symmetry rescues the spectrum: both bonds violate parity in ways that compound rather than cancel, and the eigenvalues genuinely fail to pair.

The framework gives us the *test* (operator equation vs spectrum), the *empirical 3/19/14 split*, and a structural *signature* for hard breaks (every hard break contains either a (a─)-term or is (──)+(─b)). What it does not yet give us is an analytical predictor that, given a soft-break combo, identifies the specific symmetry that rescues its spectrum. That is genuinely open: a derivation per case, or a unifying theorem.

*(Since answered, 2026-06: the predictor turned out to be a single graph property. A soft break is exactly one whose Hamiltonian hopping graph is bipartite in the dephasing basis; the symmetry that rescues the spectrum is the chiral sign-flip that sends one side of that bipartition to minus itself. So it was a unifying theorem after all, not a derivation per case: the [windowed converse](../docs/proofs/PROOF_F103_F87_Z2_CUBED_REFINEMENT.md), non-bipartite ⟹ hard, now a theorem with no residual. The living-through of that same closure is [On the One Diagonal](ON_THE_ONE_DIAGONAL.md).)*

---

What this means, and the test (added 2026-04-25 evening):

The framework predicts: for the 19 soft-break Hamiltonians, eigenvalues of L pair under λ ↔ −λ−2Σγ to machine precision (palindromic spectrum), but Π applied to L's eigenvectors does NOT land on the partner-eigenvalue's eigenspace (eigenvector pairing severely broken). For the 3 truly-unbroken cases both pairings are exact. For the 14 hard-broken cases the eigenvalue pairing already fails.

We tested this in simulation (`simulations/soft_break_eigenvector_test.py`). For each of the 36 two-term Hamiltonians at N=3, we diagonalised L, paired eigenvalues as (λ_i, −λ_i − 2Σγ), and computed the subspace overlap |⟨v_partner | Π v_i⟩| / (‖v_partner‖ · ‖Π v_i‖). The prediction is exact:

| Category | Count | Eigenvalue pair error | Eigenvector overlap |
|----------|-------|----------------------|---------------------|
| truly unbroken | 3 | ~10⁻¹⁴ | **1.000000 exact** |
| soft broken | 19 | ~10⁻¹⁴ | **0.000 to 0.598** (dramatically broken) |
| hard broken | 14 | 10⁻⁵ to 10⁻⁴ | (eigenvalues don't pair) |

The 19 soft-broken cases have eigenvalue pairing intact at machine precision while Π conjugation sends eigenvectors **completely off** their partner-eigenvalue's eigenspace. The eigenvector overlap is not "just below 1"; in many cases it is **0 exactly**: Π sends v_i into a subspace orthogonal to the partner v_j's eigenspace. The spectrum shows palindromic pairing while Π scrambles the eigenvectors.

Specifically:
- XX+XZ, XX+ZX, XZ+YY, YY+ZX: eigenvector overlap min = 0 AND avg = 0 (complete reshuffling)
- XY+YX, XY+ZZ, YX+ZZ: min ≈ 0, avg 0.2-0.6 (partial mixing)
- XX+YZ, YY+YZ, YY+ZY: min ≈ 0.005-0.01, avg 0.3-0.4

The structural reading: the standard Heisenberg-form receivers (the 3 truly-unbroken: XX+YY, XX+ZZ, YY+ZZ) preserve full Π-symmetry of L's spectrum. Any departure from this set (the 19 soft-broken) breaks the eigenvector pairing severely while leaving the spectrum intact. The 14 hard-broken break both.

What we did at the super-operator level: eigendecomposition of the 64×64 Pauli-basis Liouvillian at N=3, computing the subspace overlap directly. The eigenvector pairing is a property of the super-operator Π·L·Π⁻¹, not of the Hamiltonian or the spectrum.

What we did at the operator (observable) level: built a hardware-ready test. Initial state |+−+⟩ on N=3, evolve under Heisenberg-form (XX+YY) vs soft-broken form (XY+YX), measure ⟨X₀Z₂⟩ and ⟨Z₀X₂⟩ at t=0.8 with γ=0.1 Z-dephasing.

Framework idealised prediction:
  - truly_unbroken (XX+YY): ⟨X₀Z₂⟩ = 0
  - soft_broken (XY+YX): ⟨X₀Z₂⟩ = −0.62
  - Discriminator strength: 0.62 in absolute Pauli-expectation units.

Aer simulation with realistic IBM Heron r2 noise (T1=250μs, T2=240μs, 2-qubit gate error 0.5%, readout error 3%) at 4096 shots/basis: the difference survives realistic hardware noise.
  - truly_unbroken: ⟨X₀Z₂⟩ = −0.02 (within statistical noise of 0)
  - soft_broken: ⟨X₀Z₂⟩ = −0.66 (slightly stronger than ideal: the circuit's Trotter steps strengthen it)
  - Aer discrimination: 0.64, essentially the ideal 0.62.

A hardware run on Marrakesh would consume ~3-5 QPU minutes (3 Hamiltonians × 9 Pauli bases). Statistical error per Pauli expectation at 4096 shots is ≈ 0.015; the predicted soft-vs-truly signal is 0.62, a separation of ~30σ for the difference of two expectations.

---

The hardware run (added 2026-04-26 night):

We ran it. Job `d7mjnjjaq2pc73a1pk4g` on `ibm_marrakesh` (Heron r2) on path [48, 49, 50] using `run_soft_break.py` (the first script in the IBM tomography pipeline that imports `framework.py` directly).

Hardware result on real Marrakesh, 4096 shots/basis:

| Observable | Framework idealised | Aer w/ Marrakesh noise | **Hardware Marrakesh** |
|------------|---------------------|------------------------|--------------------------|
| ⟨X₀Z₂⟩ truly_unbroken | 0.000 | -0.020 | **+0.011** |
| ⟨X₀Z₂⟩ soft_broken | -0.623 | -0.660 | **-0.711** |
| ⟨X₀Z₂⟩ hard_broken | +0.195 | +0.230 | **+0.205** |
| **Δ(soft − truly)** | -0.62 | -0.64 | **-0.72** |

The difference on real hardware is slightly STRONGER than the idealised framework predicts, and the circuit is why: the hardware runs three first-order Trotter steps, and the Trotter prediction gives Δ(soft − truly) = −0.723 against the measured −0.722 ([`MARRAKESH_THREE_LAYERS`](../experiments/MARRAKESH_THREE_LAYERS.md)). T1 would attenuate the signal, not strengthen it.

All three values match their Trotter predictions on hardware; the two nonzero ones stand at ~13σ (hard) and ~47σ (soft). Full table in [`data/ibm_soft_break_april2026/`](../data/ibm_soft_break_april2026/README.md).

---

The single sentence:

The 22 V-Effect-unbroken Hamiltonians at N=3 are not 22; they are 3 truly unbroken plus 19 cases where Π's operator equation breaks but the spectrum does not register the break. This decomposition is invisible at the spectral level (Heisenberg's, Pauli's, the V-Effect's all agree the spectrum is palindromic). It is visible at the **super-operator level** (Π·L·Π⁻¹ + L + 2Σγ·I, the conjugation of the Liouvillian itself). The hidden 19 are real against Π, and predictable.

---

What breaks is Π's equation, and Π is written in one frame. A one-site turn that keeps the Z-dephasing mixes X with Y and leaves every Z where it is, while truly asks every string of H for an even number of Y and an even number of Z: a turn can change the first count, never the second. Three of the 19 carry no string with an odd number of Z: XY+YX, XY+ZZ and YX+ZZ, the three whose mirror in [the Klein routing of the two-term palindrome](../experiments/TWO_TERM_PALINDROME_KLEIN_ROUTING.md) alternates between odd and even sites. With every second site turned by (X + Y)/√2 they are truly, XY+YX becoming XX+YY itself, and their mirror is Π turned. The other sixteen each carry an XZ, ZX, YZ or ZY term, and no such turn makes them truly: the size of M has a floor, set by those terms, that no such turn lowers. Their mirrors, built in the Klein routing note, are products of one-site maps too, but not Π seen from a frame of one-site turns that keep the Z-dephasing.

The hardware pair was one of the three. Soft XY+YX from |+−+⟩ is truly XX+YY from |+, −i, +⟩ on every reading the run took, in the Z-dephasing model that fits it, so its ⟨X₀Z₂⟩ of −0.62 (−0.72 in the circuit) measures where the middle qubit starts against the turn, not the break ([Marrakesh Three Layers](../experiments/MARRAKESH_THREE_LAYERS.md)).

---

*"Wir müssen diesmal erklären, nicht Heisenberg oder Pauli, wir."*  Thomas Wicht, 2026-04-25

---

## Pointers

- [V-Effect Fine Structure](../experiments/V_EFFECT_FINE_STRUCTURE.md): the empirical 3/19/14 decomposition.
- [V-Effect Palindrome](../experiments/V_EFFECT_PALINDROME.md): the original 14/22 finding, March 2026.
- [Two-Term Palindrome, Klein Routing](../experiments/TWO_TERM_PALINDROME_KLEIN_ROUTING.md): the explicit mirror of each of the 19.
- [the zero-immunity proof](../docs/proofs/PROOF_ZERO_IMMUNITY.md): the analytical (w=0, w=N) extreme-sector immunity that grounds the framework's strict test.
- [Heisenberg Reloaded](../hypotheses/HEISENBERG_RELOADED.md): the level-stack picture.
- [On the Painter Principle](ON_THE_PAINTER_PRINCIPLE.md): every painter from their spot. The framework is our spot; the soft break is what we see from it.
- `simulations/soft_break_eigenvector_test.py`: super-operator-level verification (eigenvector pairing under Π).
- `simulations/soft_break_hardware_design.py`: the hardware design, Pauli expectations of the three Hamiltonians.
- `simulations/soft_break_aer_test.py`: validation against realistic Marrakesh-like noise on Aer.
- `D:\Entwicklung\Projekte\.NET Projekte\AIEvolution\AIEvolution.UI\experiments\ibm_quantum_tomography\run_soft_break.py`: hardware-runnable. First script in that pipeline directory that imports framework.py directly. Modes: `--simulate` (Aer), `--hardware` (real IBM Heron r2). Aer self-test reproduces Δ⟨X₀Z₂⟩ = -0.64 (framework ideal: -0.62).
