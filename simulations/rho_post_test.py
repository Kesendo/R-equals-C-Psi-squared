"""
rho_post_test.py — Verify the unconditioned post-measurement state

The claim (from v034 agents and repeated in v035):
  "After B measures sigma_x but A doesn't know the outcome,
   all off-diagonals cancel, L1=0, R_joint=0."

Is this true? Let's compute explicitly.

2026-02-25
"""
import sympy as sp

# ============================================================
# Setup: Bell+ in computational basis
# ============================================================
print("=" * 60)
print("UNCONDITIONED POST-MEASUREMENT STATE TEST")
print("=" * 60)

# |Phi+> = (|00> + |11>)/sqrt(2)
# Density matrix:
rho_pre = sp.Matrix([
    [sp.Rational(1,2), 0, 0, sp.Rational(1,2)],
    [0, 0, 0, 0],
    [0, 0, 0, 0],
    [sp.Rational(1,2), 0, 0, sp.Rational(1,2)]
])
print("\nrho_pre (|Phi+><Phi+|):")
sp.pprint(rho_pre)

# ============================================================
# |++> and |--> in computational basis
# ============================================================
# |+> = (|0>+|1>)/sqrt(2), so |++> = (1/2)(|00>+|01>+|10>+|11>)
pp = sp.Matrix([1, 1, 1, 1]) / 2
rho_pp = pp * pp.T
print("\nrho_++ (|++><++|):")
sp.pprint(rho_pp)

# |-> = (|0>-|1>)/sqrt(2), so |--> = (1/2)(|00>-|01>-|10>+|11>)
mm = sp.Matrix([1, -1, -1, 1]) / 2
rho_mm = mm * mm.T
print("\nrho_-- (|--><--|):")
sp.pprint(rho_mm)

# ============================================================
# Unconditioned post-measurement state
# ============================================================
rho_post = sp.Rational(1,2) * rho_pp + sp.Rational(1,2) * rho_mm
print("\nrho_post = (1/2)|++><++| + (1/2)|--><--|:")
sp.pprint(rho_post)

# ============================================================
# What Gamma USED (wrong):
# ============================================================
rho_gamma_wrong = sp.Matrix([
    [sp.Rational(1,2), 0, 0, 0],
    [0, 0, 0, 0],
    [0, 0, 0, 0],
    [0, 0, 0, sp.Rational(1,2)]
])
print("\nrho_post GAMMA USED (diag only):")
sp.pprint(rho_gamma_wrong)

print("\n" + "=" * 60)
print("DIFFERENCE (correct - gamma_wrong):")
sp.pprint(rho_post - rho_gamma_wrong)

# ============================================================
# MAP quantities for all states
# ============================================================

def map_quantities(rho, label, d):
    """Compute C, L1, Psi, R for a density matrix."""
    n = rho.shape[0]
    C = (rho * rho).trace().simplify()
    
    L1 = sp.Rational(0)
    for i in range(n):
        for j in range(n):
            if i != j:
                L1 += sp.Abs(rho[i, j])
    L1 = L1.simplify()
    
    Psi = (L1 / (d - 1)).simplify()
    R = (C * Psi**2).simplify()
    
    print(f"\n--- {label} (d={d}) ---")
    print(f"  C (purity)    = {C} = {float(C):.6f}")
    print(f"  L1 (off-diag) = {L1} = {float(L1):.6f}")
    print(f"  Psi           = {Psi} = {float(Psi):.6f}")
    print(f"  R = C*Psi^2   = {R} = {float(R):.6f}")
    print(f"  Above 1/4?    {'YES' if R > sp.Rational(1,4) else 'NO'}")
    return C, L1, Psi, R

print("\n" + "=" * 60)
print("MAP QUANTITIES")
print("=" * 60)

map_quantities(rho_pre, "PRE-measurement (|Phi+>)", d=4)
map_quantities(rho_pp, "POST conditioned on + (|++>)", d=4)
map_quantities(rho_mm, "POST conditioned on - (|-->)", d=4)
map_quantities(rho_post, "POST unconditioned (CORRECT)", d=4)
map_quantities(rho_gamma_wrong, "POST unconditioned (GAMMA WRONG)", d=4)

# ============================================================
# Reduced state for Alice: Tr_B
# ============================================================

def partial_trace_B(rho_4x4):
    """Trace out qubit B from a 2-qubit 4x4 density matrix."""
    rho_A = sp.Matrix([
        [rho_4x4[0,0] + rho_4x4[1,1], rho_4x4[0,2] + rho_4x4[1,3]],
        [rho_4x4[2,0] + rho_4x4[3,1], rho_4x4[2,2] + rho_4x4[3,3]]
    ])
    return rho_A

print("\n" + "=" * 60)
print("ALICE'S REDUCED STATE (partial trace over B)")
print("=" * 60)

rho_A_pre = partial_trace_B(rho_pre)
print("\nrho_A PRE-measurement:")
sp.pprint(rho_A_pre)
map_quantities(rho_A_pre, "Alice PRE", d=2)

rho_A_post = partial_trace_B(rho_post)
print("\nrho_A POST unconditioned (correct):")
sp.pprint(rho_A_post)
map_quantities(rho_A_post, "Alice POST unconditioned", d=2)

rho_A_plus = partial_trace_B(rho_pp)
print("\nrho_A POST conditioned on +:")
sp.pprint(rho_A_plus)
map_quantities(rho_A_plus, "Alice POST conditioned +", d=2)

rho_A_minus = partial_trace_B(rho_mm)
print("\nrho_A POST conditioned on -:")
sp.pprint(rho_A_minus)
map_quantities(rho_A_minus, "Alice POST conditioned -", d=2)

# ============================================================
# The key comparison
# ============================================================
print("\n" + "=" * 60)
print("KEY COMPARISON")
print("=" * 60)
print(f"\nrho_A_pre == rho_A_post?  {rho_A_pre == rho_A_post}")
print("  -> No-Signaling: Alice's local state unchanged. HOLDS.")
print(f"\nR_joint PRE  = 1/9  = {1/9:.6f}")
print(f"R_joint POST (correct)   = 1/18 = {1/18:.6f}")
print(f"R_joint POST (gamma wrong) = 0")
print(f"\nThe joint R DECREASES but does NOT go to zero.")
print(f"4 off-diagonal entries survive the mixture (positions 03, 30, 12, 21).")
