# D6: Spectral Gap and Mixing Time

**What this derivation is about:** The spectral gap (the decay rate of the slowest non-stationary mode, which sets the timescale for the system to forget its initial state) is exactly 2γ, independent of chain length, coupling strength, or topology. This means the mixing time (how long until the system reaches its equilibrium) scales only as N·ln(4)/(2γ), growing linearly with system size.

**Source formulas:** 1 (palindrome), 3 (decay rate bounds)
**Tier:** 1 (from palindrome proof)
**Status:** VERIFIED (N=2-5, deviation < 1e-14)

## Derivation

From F3: the minimum nonzero decay rate is 2*gamma (w=1 modes).

    Spectral gap = min{|Re(lambda)| : lambda != 0} = 2*gamma

This is the rate of the slowest decaying mode. It determines:
- The slowest decoherence timescale: tau_slow = 1/(2*gamma)
- The mixing time (convergence to steady state):

    t_mix <= (1/gap) * ln(dim) = ln(4^N) / (2*gamma) = N*ln(4) / (2*gamma)

The gap is independent of the Hamiltonian H and the topology.
It depends only on the per-qubit dephasing rate gamma.

## Numerical verification

| N | Gap (numerical) | Gap (formula) | Deviation |
|---|----------------|--------------|-----------|
| 2 | 0.1000000000   | 0.1000000000 | < 4e-16   |
| 3 | 0.1000000000   | 0.1000000000 | < 2e-15   |
| 4 | 0.1000000000   | 0.1000000000 | < 6e-15   |
| 5 | 0.1000000000   | 0.1000000000 | < 1e-14   |

(gamma = 0.05 for all tests)

Script: [`simulations/verify_derivations.py`](../../../simulations/verify_derivations.py)

## Replaces

Numerical search for the smallest nonzero eigenvalue. The gap
is always 2*gamma, so mixing time estimates need only N and gamma.
