# D09: Palindrome-Paired Band SFF

**Derives:** equal frequency SFFs for palindrome-paired decay-rate bands.
**From:** F1's full-spectrum multiset pairing.
**Status:** Exact corollary for mirrored bands; the numerical producer uses average-light bins, not invariant XY-weight eigenvalue sectors.

## Statement and scope

Let λ=-d+iω run through the full Liouvillian spectrum with algebraic multiplicity.
Under uniform per-site dephasing γ, F1 pairs

    λ ↦ −λ−2Nγ, hence d ↦ 2Nγ−d and ω ↦ −ω.

For a decay-rate interval I, define its reflected interval I'={2Nγ−d : d∈I}.
The bands B_I and B_I' contain equal numbers M_I of eigenvalues. Their raw
frequency SFFs, unnormalized or normalized by M_I², agree for every real t:

    K_I(t) = |Σ_{λ∈B_I} exp(i Im(λ)t)|² / M_I² = K_I'(t).

The normalized statement requires nonempty bands. Boundary membership must be
reflected consistently, including endpoint inclusions and numerical tolerances.

## Proof

The spectral palindrome supplies a multiplicity-preserving bijection from B_I
to B_I'. It negates every frequency. Thus the partner trace amplitude is the
complex conjugate of the original trace amplitude; their modulus squares agree.
Equal cardinalities give the same normalization. No eigenvector weight-sector
assumption is used.

## What the producer bins

[spectral_form_factor.py](../../../simulations/spectral_form_factor.py) selects
eigenvalues by |d−2wγ|<γ. The reflected bin is centred at 2(N−w)γ.
These are palindrome-paired decay-rate bands. At real Hamiltonian parameters,
the Absorption Theorem identifies d=2γ⟨n_XY⟩, so they are average-light bins.

A Hamiltonian can mix Pauli strings of different XY weight. Its Liouvillian
eigenvectors therefore need not live in fixed XY-weight eigenspaces. The pure
Pauli weight complement under Π is not a license to assign every mixed
eigenmode a fixed integer weight.

## Exact zero-frequency endpoint bands

For a connected uniform chain at strictly positive uniform dephasing, the N+1 stationary modes and their palindrome
partners at λ=−2Nγ, every frequency is zero. Their trace amplitude is N+1 for
all t. The unnormalized SFF is the constant (N+1)²; the normalized SFF is 1
for all t, not an impulse at t=0. This statement concerns the exact endpoint
eigenspaces, not the entire I/Z or XY-weight-N Pauli subspaces, nor a finite-width
bin containing additional eigenvalues. At γ=0 the two endpoints coincide with
λ=0 and merge into the generally larger Hamiltonian-commutator kernel, so the
N+1 endpoint multiplicity statement does not apply there.

## Executable checks

The producer reports mirrored decay-rate bins at N=3−5. The frequency-kernel
tests in [test_sff_window_summary.py](../../../simulations/tests/test_sff_window_summary.py)
check constant zero-frequency sums against an oscillatory two-frequency control.
A pair contributes 2cos(ωt) to the trace amplitude; squaring a sum of pairs
produces doubled and cross frequencies in the SFF.
