# IBM Marrakesh Soft-Break Hardware Test, April 2026

Hardware run of one two-term Pauli-pair Hamiltonian from each class of
the F87 trichotomy (truly, soft, hard), against the framework's
predictions for every two-qubit Pauli expectation on the two end
qubits. The soft Hamiltonian is the truly one with its middle qubit
turned (see the Reading), so the run tests the predictions, not the
soft class.

Run on `ibm_marrakesh` (Heron r2) on path [48, 49, 50] using the framework-
based pipeline at `D:\...\ibm_quantum_tomography\run_soft_break.py`. This
is the FIRST script in that pipeline directory that imports `framework.py`
directly to construct Hamiltonians and observables.

## Files

| File | Description |
|------|-------------|
| `soft_break_ibm_marrakesh_20260426_001101.json` | Raw counts and reconstructed 2-qubit Pauli expectations for 3 Hamiltonians × 9 tomography bases on (q0=48, q2=50). |

## Experiment summary

- **Backend:** ibm_marrakesh (Heron r2). Same backend as 2026-04-25 K-partnership run.
- **Date:** 2026-04-26, 00:11 UTC.
- **Job ID:** `d7mjnjjaq2pc73a1pk4g`.
- **Path:** [48, 49, 50] (3 contiguous Marrakesh qubits, same chain prefix as the K-partnership 2026-04-25).
- **Initial state:** |+−+⟩ X-Néel (Hadamard on q0, q2; Hadamard+Z on q1).
- **Three Hamiltonians:** `truly_unbroken` = XX+YY (XY-model, Heisenberg's both-parity-even subset), `soft_broken` = XY+YX (matched bit_b violations), `hard_broken` = XX+XY (mixed parity violations).
- **Evolution:** uniform J = 1, t = 0.8, n_trotter = 3 (first-order Trotter via Qiskit `PauliEvolutionGate`).
- **Tomography:** 9 Pauli bases on (q0, q2), 4096 shots/basis. Total 27 circuits.
- **QPU cost:** ~3 minutes wall-clock (well within budget).

## Results: hardware vs prediction

| Observable | Continuous Lindblad (γ=0.1 Z-only) | Trotter n=3 (γ=0.1 Z-only) | Aer w/ Marrakesh-like noise | **Hardware Marrakesh** |
|------------|------------------------------------|----------------------------|------------------------------|-------------------------|
| ⟨X₀Z₂⟩ truly_unbroken | +0.000 | +0.000 | -0.020 | **+0.011** |
| ⟨X₀Z₂⟩ soft_broken | -0.623 | -0.723 | -0.660 | **-0.711** |
| ⟨X₀Z₂⟩ hard_broken | +0.195 | +0.213 | +0.230 | **+0.205** |
| **Δ(soft − truly)** | **-0.62** | **-0.72** | **-0.64** | **-0.72** |

The hardware Δ(soft − truly) matches the **Trotter n=3** prediction to
within 0.0014, not the continuous-Lindblad idealization. T1 thermal
relaxation does not produce the hardening: it monotonically *attenuates*
the soft signal (γ_T1=0.5 gives Δ = -0.44 at γ_Z=0.1, further from
hardware). Quantification in
[`simulations/marrakesh_t1_amplification_test.py`](../../simulations/marrakesh_t1_amplification_test.py).

The actual mechanism is **Trotter discretization at δt = 0.267**: with
‖H‖ ~ J·N ~ 3, the small-step condition ‖H·δt‖ ≪ 1 fails. The Trotter
column is the circuit the hardware ran: the runner builds one
`SparsePauliOp` with the bond loop outside the term loop, and Qiskit's
default Lie-Trotter synthesis keeps that order (the framework's
`predict_signature_table` is gated against it in
[`test_predict_signature_table.py`](../../simulations/framework/tests/workflows/test_predict_signature_table.py)).
The two orders differ by one swap, the second term on bond (0,1)
against the first on bond (1,2). For soft (XY+YX) those meet on the
shared site as X and X and commute, so the order is invisible: Trotter moves ⟨X₀Z₂⟩ from −0.62 to −0.72, onto the hardware.
For hard (XX+XY) they do not, and the order is what the bias depends on:
bonds outside gives +0.213 against the hardware's +0.205, so Trotter
alone matches ⟨X₀Z₂⟩ in all three categories.

Joint optimization of (γ_Z, γ_T1) over all 45 hardware
observable-Hamiltonian pairs converges to γ_T1 ≈ 0 with γ_Z = 0.105
(total RMS 0.039); the data does not support adding T1.

On the Y-containing observables (⟨Y₀Z₂⟩, ⟨I₀Y₂⟩, ⟨X₀Y₂⟩) continuous Lindblad,
Trotter and hardware agree in sign where the models are nonzero (truly's
⟨Y₀Z₂⟩, hard's ⟨I₀Y₂⟩ and ⟨X₀Y₂⟩); Trotter changes the magnitudes. The
other entries are exact zeros of both models, and soft's ⟨Y₀Z₂⟩ = +0.098
is the largest hardware reading among them: XY+YX keeps a real state real,
so every odd-Y reading of soft is zero for any real channel, up to a small
floor from noise inside the compiled frames and the readout's leak. The
same reading, on path [4, 5, 6] with the gates of the 2026-04-30 run, came
out nonzero again on 2026-10-29
([the reality-window flight](../../experiments/REALITY_WINDOW_HARDWARE_PREDICTION.md)),
and the same reading stood at the same sign on Kingston [13, 14, 15]
on 2026-11-16 (raw +0.067), on a day whose readout broke the guards on
two of its qubits, the odd-Y readings there not seen to turn with
the frame of any of the three path qubits, at that flight's sensitivity
and within the limits stated in
[the driver's frame](../../experiments/DRIVE_FRAME_HARDWARE_PREDICTION.md) §7.

## Other discriminating Pauli expectations

The three Hamiltonians separate on hardware in at least one of the 16
two-qubit Pauli expectations:

| Pauli (P_0, P_2) | truly_unbroken | soft_broken | hard_broken |
|------------------|-----------------|-------------|-------------|
| ⟨X₀Z₂⟩ | +0.011 | **-0.711** | +0.205 |
| ⟨Z₀X₂⟩ | -0.002 | **-0.479** | -0.042 |
| ⟨Y₀Z₂⟩ | +0.583 | +0.098 | +0.022 |
| ⟨Z₀Y₂⟩ | +0.187 | +0.017 | -0.004 |
| ⟨X₀X₂⟩ | -0.010 | -0.018 | **+0.212** |
| ⟨Y₀Y₂⟩ | +0.472 | +0.436 | **-0.011** |
| ⟨Z₀Z₂⟩ | +0.163 | +0.204 | +0.032 |

`truly` and `soft` differ most in ⟨X₀Z₂⟩, ⟨Y₀Z₂⟩ and ⟨Z₀X₂⟩. `hard` separates from both
in ⟨X₀X₂⟩ and ⟨Y₀Y₂⟩. Statistical error per Pauli expectation at 4096 shots
is ≈ 0.015. On ⟨X₀Z₂⟩, hard's +0.205 stands at 13σ and soft's −0.711 at
47σ.

## Reading

The run measured what the framework's Z-dephasing model, run in the
circuit's Trotter order, predicts for the three Hamiltonians. It does not
read the soft class. XY+YX is XX+YY with the middle qubit turned:
V = (X + Y)/√2 on q1 exchanges X and Y and sends Z to −Z, so it carries
XX+YY onto XY+YX bond by bond, keeps the Z-dephasing, and leaves q0, q2
and every reading on them alone. Under Z-dephasing, the model the joint
fit above selects (γ_T1 ≈ 0), soft from |+−+⟩ is therefore truly from
|+, −i, +⟩ (V|−⟩ is |−i⟩ up to a phase), on all 16 two-qubit Paulis and
on ⟨n⟩ (the mean excitation number on q0 and q2), continuous and in the
circuit's Trotter order, to 3·10⁻¹⁶: the soft − truly difference is one
dynamics at two middle-qubit states. The residual M = Π·L·Π⁻¹ + L + 2Σγ·I
is nonzero for XY+YX and zero for XX+YY because Π is written in one
frame; Π turned by V is an exact mirror of XY+YX
([MARRAKESH_THREE_LAYERS](../../experiments/MARRAKESH_THREE_LAYERS.md),
Layer 3).

With amplitude damping the identification with XX+YY fails, because the X
in V turns decay on q1 into pumping. The quarter turn S on q1 keeps σ⁻ and
carries XY+YX onto YY − XX, again built from truly strings: both arms are
then truly chains, XX+YY from |+−+⟩ and YY − XX from |+, −i, +⟩, and no
longer one dynamics.

This run does NOT establish:

- That the predictions hold for the other soft or hard combinations; one
  Hamiltonian per class was flown.
- A complete error budget. With Trotter in the circuit's order (which
  fully explains Δ_hw vs Δ_continuous) and γ_Z = 0.105, residual
  per-observable RMS 0.039 remains; what carries it (readout, gate errors,
  anything beyond Z-dephasing) is not separated here.

## Reproducing the run

```bash
cd "D:\Entwicklung\Projekte\.NET Projekte\AIEvolution\AIEvolution.UI\experiments\ibm_quantum_tomography"
python run_soft_break.py --hardware --backend ibm_marrakesh --path 48,49,50 --shots 4096
```

The run-time on Marrakesh queue + execution: ~3 minutes total.

## Pointers

- [V_EFFECT_FINE_STRUCTURE](../../experiments/V_EFFECT_FINE_STRUCTURE.md): the 3/19/14 split from which the three Hamiltonians are drawn.
- [ON_THE_SOFT_BREAK](../../reflections/ON_THE_SOFT_BREAK.md): the super-operator framing.
- [PROOF_ZERO_IMMUNITY](../../docs/proofs/PROOF_ZERO_IMMUNITY.md): analytical anchor for the (w=0, w=N) extreme-sector immunity that grounds the framework's strict test.
- [framework/](../../simulations/framework/): the lean cockpit package; primitives used to build the Hamiltonians and predict observables.
- [soft_break_eigenvector_test.py](../../simulations/soft_break_eigenvector_test.py): super-operator-level verification (eigenvector pairing).
- [soft_break_aer_test.py](../../simulations/soft_break_aer_test.py): Aer with Marrakesh-like noise.
- [marrakesh_t1_amplification_test.py](../../simulations/marrakesh_t1_amplification_test.py): T1 attenuates rather than amplifies the soft signal; Trotter n=3 fully accounts for the hardening.
- [f80_ibm_soft_break_check.py](../../simulations/f80_ibm_soft_break_check.py): F80's M-spectrum for the three Hamiltonians.
- `D:\Entwicklung\Projekte\.NET Projekte\AIEvolution\AIEvolution.UI\experiments\ibm_quantum_tomography\run_soft_break.py`: the hardware-runnable pipeline.
