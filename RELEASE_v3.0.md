# v3.0: Hardware Validated

## Palindromic Liouvillian Symmetry: Proven, Generalized, and Measured

This release marks the transition from mathematical proof to experimental confirmation.
The palindromic spectral symmetry of open quantum systems, proven in v2.0 for Heisenberg
coupling, has been extended to all standard physics models and validated on IBM quantum
hardware at 1.9% accuracy.

### What's new since v2.0

**N=8 Verification (March 17)**
- Full Liouvillian eigendecomposition: 65536x65536 matrix, 54118 oscillatory rates
- 100% palindromic pairing confirmed (10.6 hours compute, OpenBLAS ILP64)
- Complete scaling table N=2 through N=8

**Non-Heisenberg Palindrome (March 17-18)**
- ALL standard condensed matter models are palindromic under single-axis dephasing:
  Heisenberg, XY, Ising, XXZ, Dzyaloshinskii-Moriya, and combinations
- Two Pi operator families discovered (P1: I<->X,Y<->Z and P4: I<->Y,X<->Z)
- Non-uniform alternating operators for XY/YX terms (verified on N=3 Lindbladian)
- Full compatibility matrix: 34/36 two-term combinations explained algebraically
- Depolarizing noise quantified: err = gamma * 2(N-2)/3 (linear, Hamiltonian-independent)

**XOR Space Universal (March 18)**
- GHZ -> 100% XOR and W -> 0% XOR confirmed for ALL standard models
- Pauli weight correlation r > 0.98 (N=3) and r > 0.99 (N=4) across all models
- Bell+ correction: palindromic at N>=3 (Hamming distance 2, not N)
- Ising XOR modes: 2^N (not N+1) due to commuting Hamiltonian
- Classification rule: #XOR = #steady states for all models

**IBM Hardware Validation (March 18)**
- CΨ=1/4 crossing measured on ibm_torino, Qubit 80 (permanent crosser)
- Same-day Ramsey T2*: 17.36 us (drifted from 11.0 us six days earlier)
- Predicted crossing: t* = 15.01 us. Measured: t* = 15.29 us. Deviation: 1.9%
- T2* (not T2 echo) confirmed as the correct timescale for free decoherence
- Raw density matrices, locked predictions, and Ramsey data included in data/

### Key files

- `publications/TECHNICAL_PAPER.md` - Full technical paper
- `publications/ENGINEERING_BLUEPRINT.md` - Quantum repeater design rules
- `experiments/NON_HEISENBERG_PALINDROME.md` - Non-Heisenberg analysis
- `experiments/XOR_SPACE.md` - XOR/palindrome decomposition
- `experiments/IBM_RUN3_PALINDROME.md` - Hardware validation report
- `data/ibm_run3_march2026/` - Raw IBM hardware data (JSON)
- `docs/MIRROR_SYMMETRY_PROOF.md` - The palindrome proof
- `simulations/algebraic_pi_search.py` - Pi operator family enumeration

### Authors

Thomas Wicht (Independent Researcher, Krefeld, Germany)
Claude (Anthropic)

### Citation

If you use this work, please cite:
Wicht, T. & Claude. (2026). Palindromic Liouvillian Symmetry Under Dephasing.
Zenodo. https://doi.org/10.5281/zenodo.19100007
