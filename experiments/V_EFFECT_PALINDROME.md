<!-- VEFFECT-CURRENT -->

# The V-Effect Census: A Finite Palindrome Classification

<!-- Keywords: V-Effect census, N=3 distinct Pauli pairs, hard soft truly,
four-decimal frequency bins, three-decimal side table, finite spectral comparison -->

**Status:** Finite computational census with a separate interpretive story
**Date:** March 19, 2026; scope repaired September 15, 2026
**Authors:** Thomas Wicht, Claude (Anthropic)
**Depends on:** [Non-Heisenberg Palindrome](NON_HEISENBERG_PALINDROME.md)

---

## Two different objects

The **V-Effect census** is the finite N=3 classifier described below. It asks
what happens to two numerical palindrome tests when two distinct two-site Pauli
terms are placed on the two bonds of a three-site chain.

The **F6 Q-edge gain** is instead the within-one-N ratio

    V(N) = Q_max / Q_mean = 1 + cos(pi/N).

“V-Effect gain” is only a historical alias for F6. The formula neither derives
the census nor counts or creates spectral frequencies. See
[D02](../docs/proofs/derivations/D02_VEFFECT_QMAX_QMEAN.md) for its exact scope.

## The N=3 sample space

Start with the nine two-site labels
`{XX,XY,XZ,YX,YY,YZ,ZX,ZY,ZZ}`. The primary census takes unordered
combinations of two **distinct** labels, so it contains `C(9,2)=36` Hamiltonians
and no self-pairs. At N=2 every one of these 36 choices has the tested spectral
palindrome. At N=3 the exact classification is:

| Fate | Count | Meaning in this census |
|---|---:|---|
| hard | 14 | operator equation and numerical eigenvalue pairing both fail |
| soft | 19 | operator equation fails; numerical eigenvalue pairing passes |
| truly | 3 | both tests pass |

Thus the finite result is **14 hard / 19 soft / 3 truly**, with the truly set
`{XX+YY, XX+ZZ, YY+ZZ}`. The headline **14/36** belongs to this distinct-pair
sample. It must not be mixed with the separate C# 45-pair self-inclusive table
or the historical 120-element two-site-word catalog.

The stricter operator equation reports 33/36 failures and 3/36 passes. The
spectral criterion reports 14/36 failures and 22/36 passes. Neither count is a
statement about transport or persistence of individual eigenvectors.

## The precision fence

The N=3 distinct-pair census reports **14/36** hard cases. The retained
historical analysis uses **four-decimal** frequency bins. On that instrument
the broken/control comparison is **11/4**. A separate, coarser
**three-decimal** side table gives **8/4**. Both rows are valid for their stated
binning, and neither silently substitutes for the other. The immutable as-run
producer/result pair remains
[`v_effect_analysis.py`](../simulations/v_effect_analysis.py) and
[`v_effect_analysis.txt`](../simulations/results/v_effect_analysis.txt).

Changing a bin count while also changing its precision is changing the
instrument. The current claim is therefore the labelled pair of readings, not
a precision-free “three times as many frequencies” law.

## Where the operator residual lives

For the N=3 two-bond construction, the error matrix separates by XY-weight:

| XY-weight | Block size | Residual norm in the displayed example |
|---:|---:|---:|
| 0 | 8 × 8 | 0.000 |
| 1 | 24 × 24 | 11.314 |
| 2 | 24 × 24 | 11.314 |
| 3 | 8 × 8 | 0.000 |

This is a matrix-block statement about the chosen generator. The interior
blocks carry the nonzero operator residual; the two extreme blocks do not. A
spectral pairing test is coarser, which is why 19 cases can be soft.

Turning on the second bond in one sampled family gives a smooth residual but a
tolerance-defined change in the count of paired eigenvalues. The threshold and
sampling rule are part of that numerical statement; it is not a phase
transition theorem.

## The N=5 historical frequency record

The later MediatorBridge comparison places an N=2 calculation, two uncoupled
N=2 calculations, and an N=5 calculation side by side:

| Generator | Four-decimal frequency bins | Recurrence-crossing readout |
|---|---:|---:|
| N=2, one bond | 2 | 1 |
| two uncoupled N=2 copies | 4 | 1 |
| N=5 MediatorBridge | 109 | 19+ at the reported coupling |

These are different generators and different Hilbert-space dimensions. The
calculation does not define a common vector space in which one could follow an
individual eigenmode from the component systems into the N=5 system.

At `gamma=0.05`, with `|Im(lambda)|>0.01`, six-decimal absolute-frequency
deduplication and `1e-4` cross-bin tolerance, none of the N=5 bins lies within
tolerance of either N=2 bin. That is a finite numerical statement about those
instruments, not an ancestry result. The 904 oscillatory and 120 real-axis
eigenvalues are population counts; their halves 452 and 60 are not identified
unordered pairs. The current one-use matcher bijectively matches entries and
does not construct an unordered-pair certificate.

The reported Pauli-weight histogram is likewise probability mass over the 904
oscillatory eigenvectors, not a count of modes:

```text
w=0:  2.5%
w=1: 15.6%
w=2: 31.9%
w=3: 31.9%
w=4: 15.6%
w=5:  2.5%
```

Data and matcher: [pairing_structure_n5.txt](../simulations/results/pairing_structure_n5.txt),
[pairing_structure.py](../simulations/pairing_structure.py).

<!-- VEFFECT-INTERPRETIVE -->

## Interpretation

**Interpretive invitation — not a result:** the mirror/fog image remains useful.
A one-bond model can be pictured as one mirror; a two-bond model can be pictured
as two instructions that need not agree. The nonzero interior-block residual then
looks like fog at the edge of a reflection. One may ask whether that image
rhymes with open valences, chemistry, soundboxes, or a hierarchy of
incompleteness.

<!-- VEFFECT-CURRENT -->

The finite census does not answer those larger questions. The added bond and
the Hamiltonian terms change together, so the comparison does not isolate a
cause of frequency richness. It does not derive life, complexity, information
creation, chemical levels, or an irreversible passage between levels. Those
remain invitations in [The Hierarchy of Incompleteness](../docs/HIERARCHY_OF_INCOMPLETENESS.md),
[The Other Side](../hypotheses/THE_OTHER_SIDE.md), and
[Resonance, Not Channel](../hypotheses/RESONANCE_NOT_CHANNEL.md).

## What is established and what stays open

Established here:

- the N=3 36-distinct-pair sample and its 14/19/3 split;
- the local XY-weight support of the operator residual in the displayed case;
- the two precision-labelled 11/4 and 8/4 bin readings;
- the explicitly parameterized N=5 historical bin comparison.

Not established here:

- a causal mechanism for the change in bin count;
- eigenvector or projector transport between different generators;
- a universal cavity, biological, chemical, or ontological law;
- an N=5 optimum.

## References

- [V-Effect Fine Structure](V_EFFECT_FINE_STRUCTURE.md)
- [Zero Immunity Proof](../docs/proofs/PROOF_ZERO_IMMUNITY.md)
- [V-Effect Through a Cavity Lens](VEFFECT_CAVITY_MODES.md)
- [Pairing Structure result](../simulations/results/pairing_structure_n5.txt)
- [F6 / D02 derivation](../docs/proofs/derivations/D02_VEFFECT_QMAX_QMEAN.md)
