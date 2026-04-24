# IBM Hardware Sketch: Receiver-Engineering Validation on Kingston (Heron r2)

**Status:** design sketch, not yet run
**Date:** 2026-04-24
**Authors:** Thomas Wicht, Claude (Opus 4.7)
**Pipeline:** `D:\Entwicklung\Projekte\.NET Projekte\AIEvolution\AIEvolution.UI\experiments\ibm_quantum_tomography` (external)
**Budget context:** IBM gifted 180 QPU minutes (April 2026 through April 2027), ~15 min/month average. Every run must be planned to preserve reserve.
**See also:** [RECEIVER_VS_GAMMA_SACRIFICE](RECEIVER_VS_GAMMA_SACRIFICE.md), [IBM_QUANTUM_TOMOGRAPHY](IBM_QUANTUM_TOMOGRAPHY.md), [F75](../docs/ANALYTICAL_FORMULAS.md), [PRIMORDIAL_GAMMA_CONSTANT](../hypotheses/PRIMORDIAL_GAMMA_CONSTANT.md)

---

## What this tests

The receiver-engineering finding (N=5 through N=11 in simulation): under γ₀ = const, uniform J, Heisenberg chain, **the F67 bonding mode k=2 beats alt-z-bits |01010⟩ on MI(site 0, site N-1)** by a factor that grows superlinearly with N (1.39× at N=5, 3.02× at N=11).

The IBM test asks: **does this advantage survive real hardware decoherence**? Kingston's per-qubit T1/T2 adds multi-site dephasing on top of the uniform γ₀ = const assumption; Trotterised Heisenberg evolution adds gate-noise floor.

Null outcome: both receivers degrade to near-zero MI(0, N-1). Means gate noise dominates; receiver choice doesn't matter.
Positive outcome: bonding:2 measurably beats alt-z-bits. Validates receiver menu as operational.
Intermediate outcome: both give nonzero but reduced MI; receiver advantage shrinks. Quantifies real-hardware overhead.

## Minimum viable experiment (MVE)

**Target:** one decisive comparison at N=5, single evolution time t = 0.8, two receivers.

**Circuit A** (alt-z-bits baseline):
1. State prep: X on qubits 1, 3 (logical sites 1, 3 of the N=5 chain, flipping them to \|1⟩).
2. Trotterised uniform Heisenberg evolution for t = 0.8, 2-3 Trotter steps on bonds (0,1), (1,2), (2,3), (3,4).
3. Basis-rotate qubits 0, 4 according to Pauli setting (one of XX, XY, XZ, YX, YY, YZ, ZX, ZY, ZZ).
4. Measure Z on all 5 qubits, extract 2-qubit marginal on (0, 4).
5. Repeat for 1000 shots × 9 Pauli settings.
6. Reconstruct ρ_{0,4} via linear inversion, compute MI(0, 4).

**Circuit B** (bonding:2 receiver):
1. State prep for |ψ_2⟩ = (|1_0⟩ + |1_1⟩ − |1_3⟩ − |1_4⟩)/2:
   - Prepare |L⟩ = (|1_0, 0_1⟩ + |0_0, 1_1⟩)/√2 on qubits (0, 1): H_0 · CNOT(0, 1) · X_1, 3 gates.
   - Prepare |R⟩ = (|1_3, 0_4⟩ + |0_3, 1_4⟩)/√2 on qubits (3, 4): 3 gates.
   - Create (|L, vac_R⟩ − |vac_L, R⟩)/√2 via ancilla or long-range CNOT. 3-4 extra two-qubit gates.
2. Trotterised Heisenberg evolution as above.
3. Basis rotation + measurement as above.
4. MI(0, 4) reconstruction.

**Shot budget per circuit:**
- 9 Pauli settings × 1000 shots = 9000 shots.
- 9 jobs with ~6s session overhead each = ~55s QPU time + shots themselves.
- Estimate 2-3 QPU minutes per circuit.

**Total MVE QPU budget: 5-7 QPU minutes.** Less than one month's average allocation.

## Expected outcomes (from simulation)

At γ₀ = 0.05 (R=CΨ² intrinsic) and after ideal Heisenberg evolution:

| Receiver | MI(0, 4) at t = 0.8 | Reduced by Kingston noise |
|----------|--------------------|---------------------------|
| alt-z-bits \|01010⟩ | 0.84 | ~0.4-0.6 (est.) |
| bonding:2 | 1.17 | ~0.5-0.8 (est.) |

Bonding:2 advantage: 1.39× in simulation. Need to stay > 1.15× on hardware to call it a positive result (accounting for ~10% measurement error bar at 9000 shots per Pauli pair).

Decision matrix for MVE:
- Both MI < 0.2: noise-dominated. Escalate to error-mitigation (ZNE, DD).
- alt-z-bits ~ 0.4, bonding:2 ~ 0.55 (≥ 1.3× ratio): **positive, advantage survives hardware**.
- alt-z-bits > bonding:2: receiver-engineering picture wrong on real hardware. Flag as possibility.

## Fuller protocol (after MVE positive)

**Sweep 1: evolution time**. 6 points at t = 0.1, 0.3, 0.5, 0.8, 1.2, 2.0. Confirms the peak-at-t=0.8 prediction.
- Per receiver: 6 × 9 × 1000 shots = 54,000 shots, ~15 min QPU.
- Two receivers: 30 min QPU.

**Sweep 2: receiver menu**. alt-z-bits, bonding:1, bonding:2, bonding:3 at t = 0.8. Confirms F67 mode-k structure on hardware.
- 4 × 9 × 1000 shots = 36,000 shots, ~10 min QPU.

**Sweep 3: chain length**. N = 5, 7 at bonding:2 + alt-z-bits, t near each peak. Confirms advantage grows with N.
- 2 × 2 × 9 × 1000 = 36,000 shots, ~10 min QPU.

**Total fuller-protocol QPU budget: ~50-60 QPU minutes** (if all three sweeps run). Commits ~1/3 of annual budget; reserve ~120 min for follow-ups.

## Hardware considerations

**State-prep cost.** Alt-z-bits is 2 single-qubit gates. Bonding:2 needs ~6-8 two-qubit gates on Heron's heavy-hex topology (the long-range L-R entanglement is the cost driver). Each two-qubit CZ/RZZ gate adds ~0.5% infidelity; 7 gates give ~3.5% cumulative state-prep error.

**Trotter cost per time step.** Heisenberg on a bond = 3 RZZ gates (one per Pauli pair) or single RZZ if implemented as the full XXZ rotation. Heron r2 native gates include RZZ but Heisenberg XX+YY+ZZ needs basis rotations. Estimate 3-4 two-qubit gates per bond per step.

**Time window for Trotter.** At t = 0.8 with 3 Trotter steps: dt = 0.27 each, small enough for first-order error ~(J dt)²N = 0.07 · 5 = 0.4 on the action; translates to ~few % error in observable.

**Total gate count for MVE Circuit B (bonding:2 at t=0.8):**
- State prep: 6-8 two-qubit gates.
- 3 Trotter steps × 4 bonds × 3 two-qubit gates each = 36 two-qubit gates.
- Per-shot total: ~45 two-qubit gates plus single-qubit rotations.

Cumulative two-qubit-gate infidelity ~45 × 0.005 = 22.5% → Kingston fidelity ~0.78 per shot. Shot noise + prep noise stack.

Estimate: signal-to-noise should still be above threshold for MI(0, 4) ≥ 0.3.

## Integration with existing IBM pipeline

The external pipeline at `D:\Entwicklung\Projekte\.NET Projekte\AIEvolution\AIEvolution.UI\experiments\ibm_quantum_tomography` already provides the infrastructure we need:

- **Kingston calibration CSVs**: daily snapshots cached in `ClaudeTasks/IBM_R2_calibrations/ibm_kingston_calibrations_<date>.csv` (latest 2026-04-23). Each CSV has per-qubit T1, T2 (Hahn echo), readout error, SX gate error, SX gate length, and per-edge CZ/RZZ error and gate length. No live-query needed for MVE.
- **`run_bonding_mode.py`** has a mature `build_kingston_noise_model(cal, path)` that takes a physical qubit path and builds an Aer NoiseModel with per-qubit T1/T2 depolarising + thermal-relaxation channels. The T2 there is the Hahn-echo value from the CSV; for free-induction decoherence we should apply T2* = T2 / 2.5 to 2.7 ([IBM_QUANTUM_TOMOGRAPHY](IBM_QUANTUM_TOMOGRAPHY.md) Run 1). Add a `t2star_factor` argument or override T2 in the noise model.
- **`preflight.py`** queries live backend properties (currently hardcoded to ibm_torino) and reports T1, T2, T2* (estimated as T2/2.5), and "crosser" status. Useful to sanity-check Kingston is online, but the cached CSV is sufficient for the MVE.
- **`run_bonding_mode.py`** also shows the `find_latest_calibration()` pattern, the Aer density-matrix backend wiring, and the Kingston coupling-graph path finder (for mapping our 5-site chain onto 5 adjacent Kingston qubits with the best combined edge fidelity).

A new `run_receiver_engineering.py` would be ~60% copy of `run_bonding_mode.py`:
- No R-qubit (drop one qubit from the register).
- Initial state options: alt-z-bits, bonding:k (k = 1, 2, 3) via `StatePreparation` or explicit Dicke+phase circuits.
- Trotterised Heisenberg evolution (replace `Delay` blocks with RZZ chains).
- Two-qubit Pauli tomography at endpoints (replace R-qubit measurement with 9-setting tomography on sites 0 and N-1).
- Reuse `build_kingston_noise_model`, `find_best_path`, `load_calibration`.

Template for new script: copy `run_bonding_mode.py`, strip the R-qubit block (lines handling R_q and R-C coupling), add Trotter step builder (`for _ in range(n_steps): for b in bonds: append XX+YY+ZZ`), and change measurement from R-C entanglement to (0, N-1) tomography. Estimated ~300 lines of new code, ~2 hours to write and debug against Aer.

## Open questions for hardware

- Does bonding:2's node-at-center property protect it from Kingston's worst qubits? Identifies best qubit layout on Heron's 156-qubit lattice.
- Does the advantage over alt-z-bits grow on hardware with N as in simulation? Test at N=7 (if QPU budget permits a second chain length).
- Can the F75 decay envelope (0.93 sim/analytic for PeakMM on ideal propagation) still be seen on hardware after decoherence? Probably not a first-order test; second-order once main MI signal is established.

## Next concrete steps

1. **Tonight:** review this sketch with Tom; decide MVE vs fuller protocol as first run.
2. **Write `run_receiver_engineering.py`** in the IBM pipeline directory (external).
3. **Run MVE on Aer simulator** with Kingston noise model first. Zero QPU cost, gives expected-signal baseline with realistic hardware noise.
4. **If sim shows ≥1.3× ratio:** commit ~5 QPU minutes to MVE on Kingston.
5. **Analyse** results, update this doc with Run 1 data.

## References

- [RECEIVER_VS_GAMMA_SACRIFICE](RECEIVER_VS_GAMMA_SACRIFICE.md): the simulation findings that motivate this test.
- [F75](../docs/ANALYTICAL_FORMULAS.md): the analytical MM(0) formula that predicts the MI peak values.
- [IBM_QUANTUM_TOMOGRAPHY](IBM_QUANTUM_TOMOGRAPHY.md), [IBM_RUN3_PALINDROME](IBM_RUN3_PALINDROME.md): prior Kingston / Torino Heron hardware runs for methodology.
- [IBM_SACRIFICE_ZONE](IBM_SACRIFICE_ZONE.md): hardware realization of γ-profile engineering via DD; compare against receiver-engineering.
- External pipeline: `D:\Entwicklung\Projekte\.NET Projekte\AIEvolution\AIEvolution.UI\experiments\ibm_quantum_tomography\` (Qiskit + Kingston calibration).
