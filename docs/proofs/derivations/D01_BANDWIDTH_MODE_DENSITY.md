# D1: Bandwidth and Mode Density

**Source formulas:** 2 (w=1 dispersion)
**Tier:** 1 (algebraic consequence of formula 2)
**Status:** VERIFIED (N=2-5, max frequency error < 5e-9)

## Derivation

From formula 2: omega_k = 4J * (1 - cos(pi*k/N)), k = 1, ..., N-1.

    omega_1   = 4J * (1 - cos(pi/N))
    omega_{N-1} = 4J * (1 + cos(pi/N))

    BW = omega_{N-1} - omega_1
       = 4J * [(1 + cos(pi/N)) - (1 - cos(pi/N))]
       = 8J * cos(pi/N)

As N -> infinity: cos(pi/N) -> 1, so BW -> 8J.

Mode density (continuous limit): k/N -> q, omega(q) = 4J(1-cos(pi*q)).

    d(omega)/dq = 4J*pi*sin(pi*q)
    rho(omega) = |dq/d(omega)| = 1 / (pi * sqrt(omega * (8J - omega)))

This is the 1D tight-binding density of states with van Hove
singularities at the band edges omega = 0 and omega = 8J.

## Numerical verification

| N | BW (numerical) | BW (formula) | Error |
|---|----------------|--------------|-------|
| 2 | 0 (1 mode)     | 0            | 0     |
| 3 | 4.0000         | 4.0000       | <1e-15|
| 4 | 5.6569         | 5.6569       | <5e-9 |
| 5 | 6.4721         | 6.4721       | <3e-9 |

Script: simulations/verify_derivations.py

## Replaces

Numerical mode density estimation for the w=1 sector.
