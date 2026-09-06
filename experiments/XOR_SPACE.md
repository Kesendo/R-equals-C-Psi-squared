# XOR Space: Spectral Endpoint and Operator Support

**Status:** The endpoint count and GHZ coherence-support statement survive in
their stated scopes. State-weight percentages, the mixed-XY predictor, and the
GHZ/W protection ranking do not.
**Experiment date:** March 16, 2026

## Current result

For the connected XY/Heisenberg chain with every bond nonzero and every site
Z-dephased at a positive rate, the eigenspace at `Re λ = -2Σγ` has dimension
`N+1`; see [F23](../docs/ANALYTICAL_FORMULAS.md#f23-xor-drain-vanishing-fraction-tier-1-combinatorial-proof)
for the count and its scope. The spectral endpoint is paired with the zero-rate
kernel when the F1 palindrome applies. An endpoint eigenvalue does not determine
how much of a prepared density matrix or readout occupies that eigenspace.

Independently of any eigendecomposition, local Z-dephasing acts diagonally on
computational-basis coherence operators. The operators
`|0...0><1...1|` and `|1...1><0...0|` differ on all N sites, so their
dissipative charge is the maximum `2Σγ`. W and embedded-Bell off-diagonal
operators connect basis states at Hamming distance two and are charged only on
those two differing sites. This is [F22](../docs/ANALYTICAL_FORMULAS.md#f22-ghz-maximum-disagreement-support-tier-2-verified-n2-5),
an operator-support statement rather than a basis-independent probability or a
state-transfer comparison.

## Historical spectral census

The original chain run recorded:

| N | nonzero-rate modes in the run | paired modes | endpoint modes |
|---:|---:|---:|---:|
| 2 | 13 | 10 | 3 |
| 3 | 60 | 56 | 4 |
| 4 | 251 | 246 | 5 |
| 5 | 1018 | 1012 | 6 |

These rows are a finite census, not the proof of the all-N count. They also do
not make the endpoint topology-independent without the conditions stated in
F23.

## Retired coordinate diagnostics

The study decomposed density matrices into right eigenvectors of a non-normal
Liouvillian and normalized squared coefficients as “palindrome weight” and
“XOR weight.” Right eigenvectors can be rescaled, and bases inside degenerate
or nearly degenerate eigenspaces can mix, so these numbers are not invariant
state probabilities.

The historical coordinate table was:

| State | “palindrome weight” | “XOR weight” |
|---|---:|---:|
| GHZ, N=2..5 | 0% | 100% |
| Bell+, N=2 | 0% | 100% |
| Bell+, N>=3 | 100% | 0% |
| W, N=2 | 0% | 100% |
| W, N>=3 | 100% | 0% |
| `|010>` | 100% | 0% |
| `|+-+>` | 86.5% | 13.5% |
| `|+++>` | 85.8% | 14.2% |

The associated `r = 0.976` mixed-XY correlation described the same retired
coordinate diagnostic. It does not predict an operational lifetime, channel
fidelity, or optimal encoding.

## What does not follow

The palindrome does not divide full quantum states into universal “channel”
and “drain” percentages. The exact GHZ coherence charge does not imply that the
entire GHZ density matrix is an endpoint eigenmode, and the distance-two support
of W coherences does not prove W outperforms GHZ. Existing direct gates use
specific preparations and observables and do not establish a universal
GHZ-versus-W advantage.

Likewise, decay-rate ordering alone does not supply an information lifetime,
standing wave, or error-correcting codespace. Each requires its own preparation,
dynamics, readout, and time-window gates.

## Reproducibility

The scripts `xor_detector.py`, `xor_detector_v2.py`, `xor_detector_v3.py`,
`xor_verify.py`, and `xor_non_heisenberg_v2.py` are historical producers of the retired coordinate
diagnostic; their coordinate percentages must not be used as gates. The active
operator-level check is
[`f22_operator_charge.py`](../simulations/f22_operator_charge.py).

## References

- [F22 and F23](../docs/ANALYTICAL_FORMULAS.md): operator-support and endpoint-count scopes
- [Error-Correction Palindrome](ERROR_CORRECTION_PALINDROME.md): the N=3 event record and negative verdict
- [Mirror Symmetry Proof](../docs/proofs/MIRROR_SYMMETRY_PROOF.md): the scoped F1 spectral pairing
