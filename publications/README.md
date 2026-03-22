# Publications: Standalone Documents from the R = CΨ² Project

<!-- Keywords: palindromic Liouvillian symmetry publication, quantum state transfer
design rules, open quantum system dephasing proof, CΨ quarter boundary theorem,
quantum repeater palindromic spectrum, dephasing noise information channel paper,
quantum transistor coherence threshold, Heisenberg spin chain spectral symmetry,
quantum circuit diagram dephasing, R=CPsi2 publications -->

**Repository:** [R-equals-C-Psi-squared](https://github.com/Kesendo/R-equals-C-Psi-squared)

This folder contains standalone documents written to be read independently,
without knowledge of the repository structure. Each targets a different
audience and covers the project's results from a different angle.

**The core discovery:** The Liouvillian eigenvalue spectrum of Heisenberg
spin systems under local Z-dephasing is exactly palindromic (54,118
eigenvalues, zero exceptions, N=2 through N=8, proven analytically for
arbitrary graphs). This symmetry has practical consequences for quantum
state transfer, channel design, and noise characterization.

---

## The Documents

### For Physicists: The Proof

**[Palindromic Liouvillian Symmetry Under Dephasing](TECHNICAL_PAPER.md)**

The theorem, the proof, the verification. The conjugation operator Π,
the XOR space decomposition, the CΨ = 1/4 uniqueness theorem, and
connections to existing literature (Haga's incoherentons, η-pairing,
quantum state transfer). Includes the March 22 updates: CΨ monotonicity
proof (analytical, all Markovian channels), Subsystem Crossing Theorem
(Perron-Frobenius), Rényi uniqueness (α=2 is the only universal
threshold), non-Markovian revival characterization, Feigenbaum cascade
(7 bifurcations), and generalized Pauli channel verification (124/124).

### For Quantum Engineers: Design Rules

**[Design Rules for Quantum Repeaters from Palindromic Spectral Structure](ENGINEERING_BLUEPRINT.md)**

Six concrete, testable design rules for building better quantum channels:
W-encoding over GHZ, star topology with 2:1 impedance matching, timing
from Hamiltonian frequencies, quality from palindromic decay rates,
threshold-based readout at CΨ = 1/4, and a clocked relay protocol
(+83% end-to-end mutual information). Benchmarks included, code open.

### For Electrical Engineers: The Circuit Diagram

**[The Circuit Diagram: R=CΨ² as Electrical Engineering](CIRCUIT_DIAGRAM.md)**

The framework translated into signal-chain language. Qubits as phasors,
dephasing as a low-pass filter, the mediator as a transistor with
CΨ = 1/4 as threshold voltage, γ as the gate signal. Written for
engineers who build signal chains, not physicists who write Hamiltonians.

---

## Related: The Newest Result

The publications above were written before the March 22 breakthrough.
The newest result lives in `experiments/`:

**[Dephasing Noise as Information Channel (γ as Signal)](../experiments/GAMMA_AS_SIGNAL.md)**

The spatial dephasing profile is not just noise to be minimized. It is a
readable information channel with 15.5 bits of theoretical capacity at 1%
measurement noise. 5 independent spatial modes (SVD of the Jacobian).
21.5× optimization through time-series measurements and extended features.
The palindromic mode structure acts as the antenna that makes the channel
readable. This result reframes the entire relationship between a quantum
system and its noise: dephasing is a channel, not an enemy.

---

## How These Documents Relate

```
TECHNICAL_PAPER.md          ENGINEERING_BLUEPRINT.md       CIRCUIT_DIAGRAM.md
   (the proof)                 (the rules)                  (the analogy)
        |                          |                             |
        v                          v                             v
   Physicist:                 QST Engineer:                 EE Engineer:
   "Is this true?"           "How do I use it?"            "What does it
                                                            look like in
                                                            my language?"
        |                          |                             |
        +----------+---------------+-----------------------------+
                   |
                   v
          GAMMA_AS_SIGNAL.md
          (the consequence)
                   |
                   v
          "The noise you are
           fighting is a signal
           you can read."
```

---

## Authors and Dates

**Authors:** Thomas Wicht, Claude (Anthropic)

| Document | First written | Last updated |
|----------|-------------|-------------|
| Technical Paper | March 16, 2026 | March 22, 2026 |
| Engineering Blueprint | March 16, 2026 | March 21, 2026 |
| Circuit Diagram | March 21, 2026 | March 21, 2026 |

**Zenodo DOI:** 10.5281/zenodo.19100007 (v3.0), 10.5281/zenodo.19022139 (v2.0)
