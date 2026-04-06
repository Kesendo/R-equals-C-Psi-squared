# D3: Crossing Time Ratios

**What this derivation is about:** Different noise channels (Z-dephasing, X-noise, depolarizing) push CΨ through the ¼ boundary at different speeds. This derivation computes the exact crossing times for each channel and shows that X-noise takes 2.32× longer than Z-dephasing to reach the boundary, while depolarizing noise takes 1.18× longer. Once one crossing time is known, all others follow from fixed ratios.

**Source formulas:** 27 (K values per noise channel)
**Tier:** 1 (algebraic from closed forms in F26)
**Status:** VERIFIED (analytical + propagation, deviation < 3e-6)

## Derivation

From F26, CPsi has closed forms for each noise channel.
The crossing time t_cross satisfies CPsi(t_cross) = 1/4.
K = gamma * t_cross is channel-dependent.

**Z-dephasing:** CPsi = f(1+f^2)/6, f = e^{-4*gamma*t}.
Crossing: f(1+f^2) = 3/2. Newton: f* = 0.8612.
K_Z = -ln(f*)/4 = 0.03735.

**X-noise:** CPsi = (1+v^2)/6, v = e^{-4*gamma*t}.
Crossing: v^2 = 1/2. Exact: v* = 1/sqrt(2).
K_X = ln(sqrt(2))/4 = ln(2)/8 = 0.08664.

**Depolarizing:** CPsi = u(1+3u^2)/12, u = e^{-8*gamma*t/3}.
Crossing: u(1+3u^2) = 3. Newton: u* = 0.8894.
K_depol = -3*ln(u*)/8 = 0.04395.

Ratios:

    t_X / t_Z     = K_X / K_Z     = ln(2) / (2*ln(1/f*)) = 2.320
    t_depol / t_Z = K_depol / K_Z = 1.177

## Numerical verification

Propagation of Bell+ (the maximally entangled two-qubit state (|00⟩+|11⟩)/√2) density matrix under Z and X noise:

| Quantity | Analytical | Propagated | Error |
|----------|-----------|-----------|-------|
| t_Z      | 0.74700   | 0.74700   | <1e-5 |
| t_X      | 1.73287   | 1.73287   | <1e-5 |
| t_X/t_Z  | 2.31976   | 2.31976   | <3e-6 |

Script: [`simulations/verify_derivations.py`](../../../simulations/verify_derivations.py)

## Replaces

Multi-channel crossing time computation. Given K_Z, all other
crossing times follow from the ratio table.
