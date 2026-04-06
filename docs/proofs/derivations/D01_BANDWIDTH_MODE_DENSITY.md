# D1: Bandwidth and Mode Density

**What this derivation is about:** The oscillatory modes of a spin chain span a finite frequency range (the bandwidth). This derivation shows that the bandwidth approaches 8J (eight times the coupling strength) as the chain grows, and that modes pile up at the band edges (van Hove singularities), not in the middle. Both results follow algebraically from the cosine dispersion relation.

**Source formulas:** 2 (w=1 dispersion)
**Tier:** 1 (algebraic consequence of F2)
**Status:** VERIFIED (N=2-5, max frequency error < 5e-9)

## Derivation

From F2: omega_k = 4J * (1 - cos(pi*k/N)), k = 1, ..., N-1.

    omega_1   = 4J * (1 - cos(pi/N))
    omega_{N-1} = 4J * (1 + cos(pi/N))

    BW = omega_{N-1} - omega_1
       = 4J * [(1 + cos(pi/N)) - (1 - cos(pi/N))]
       = 8J * cos(pi/N)

As N -> infinity: cos(pi/N) -> 1, so BW -> 8J.

Mode density (continuous limit): k/N -> q, omega(q) = 4J(1-cos(pi*q)).

    d(omega)/dq = 4J*pi*sin(pi*q)
    rho(omega) = |dq/d(omega)| = 1 / (pi * sqrt(omega * (8J - omega)))

This is the 1D tight-binding density of states (the standard nearest-neighbour hopping model in solid-state physics) with van Hove singularities (points where the density of modes diverges because the dispersion curve is flat) at the band edges omega = 0 and omega = 8J.

## Numerical verification

| N | BW (numerical) | BW (formula) | Error |
|---|----------------|--------------|-------|
| 2 | 0 (1 mode)     | 0            | 0     |
| 3 | 4.0000         | 4.0000       | <1e-15|
| 4 | 5.6569         | 5.6569       | <5e-9 |
| 5 | 6.4721         | 6.4721       | <3e-9 |

Script: [`simulations/verify_derivations.py`](../../../simulations/verify_derivations.py)

## Replaces

Numerical mode density estimation for the w=1 sector.
