"""
Derive the dwell-time prefactor from Pauli sector weights.

Question: Can the Bell+ dwell prefactor 1.080088 be derived from the
Pauli sector weight distribution, without the explicit CPsi(t) formula?

Key result: For Bell+ (only k=0 and k=2 sectors),
    dCPsi/dt = -2*gamma*Psi(1 + 6W_2)
    prefactor = 2(1 + 2W_2) / (1 + 6W_2)
where W_2 is the k=2 sector weight at the crossing.
"""

import numpy as np
from scipy.optimize import brentq
import os

OUT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")
os.makedirs(OUT_DIR, exist_ok=True)

lines = []  # collect output for text file


def pr(s=""):
    print(s)
    lines.append(s)


# ============================================================
# Step 1: Bell+ Pauli decomposition under Z-dephasing
# ============================================================

pr("=" * 72)
pr("STEP 1: Bell+ Pauli decomposition under Z-dephasing (N=2)")
pr("=" * 72)
pr()
pr("rho(0) = (1/4)(II + XX - YY + ZZ)")
pr("Under Z-dephasing with rate gamma per qubit:")
pr("  a_P(t) = a_P(0) * exp(-2*gamma*w(P)*t)")
pr("  where w(P) = XY-weight (number of X or Y factors)")
pr()
pr("Bell+ Pauli coefficients:")
pr("  a_II = 1     (w=0, stationary)")
pr("  a_ZZ = 1     (w=0, stationary)")
pr("  a_XX = f     (w=2, f = exp(-4*gamma*t))")
pr("  a_YY = -f    (w=2)")
pr("  All others = 0")
pr()

# Verify numerically
gamma = 1.0
t_test = 0.5
f_test = np.exp(-4 * gamma * t_test)

# Build rho(t) explicitly
II = np.eye(4)
XX = np.array([[0, 0, 0, 1], [0, 0, 1, 0], [0, 1, 0, 0], [1, 0, 0, 0]])
YY = np.array([[0, 0, 0, -1], [0, 0, 1, 0], [0, 1, 0, 0], [-1, 0, 0, 0]])
ZZ = np.diag([1, -1, -1, 1])

rho_bell = (II + f_test * XX - f_test * YY + ZZ) / 4

# Check purity
C_num = np.trace(rho_bell @ rho_bell)
C_formula = (1 + f_test**2) / 2
pr(f"Numerical check at t={t_test}, f={f_test:.6f}:")
pr(f"  Purity (matrix):  {C_num:.10f}")
pr(f"  Purity (formula): {C_formula:.10f}")
pr(f"  Match: {abs(C_num - C_formula) < 1e-12}")

# L1 coherence
off_diag = np.abs(rho_bell) - np.abs(np.diag(np.diag(rho_bell)))
L1_num = np.sum(off_diag)
Psi_num = L1_num / 3
Psi_formula = f_test / 3
pr(f"  Psi (matrix):  {Psi_num:.10f}")
pr(f"  Psi (formula): {Psi_formula:.10f}")
pr(f"  Match: {abs(Psi_num - Psi_formula) < 1e-12}")


# ============================================================
# Step 2: CPsi(t) from sector weights
# ============================================================

pr()
pr("=" * 72)
pr("STEP 2: CPsi from sector weights")
pr("=" * 72)
pr()
pr("Sector weights (W_k = (1/d) * sum |a_P|^2 for w(P)=k):")
pr("  W_0 = (1/4)(1 + 1) = 1/2     (constant: {II, ZZ})")
pr("  W_1 = 0                        (no k=1 content)")
pr("  W_2 = (1/4)(f^2 + f^2) = f^2/2  (decaying: {XX, YY})")
pr()
pr("Purity: C = W_0 + W_2 = 1/2 + f^2/2 = (1+f^2)/2")
pr("Coherence: Psi = f/3 = sqrt(2*W_2)/3")
pr("CPsi = C * Psi = (1/2 + W_2) * sqrt(2*W_2)/3")
pr()

# Show the relationship W_2(Psi)
pr("State-specific relation: W_2 = f^2/2 = (3*Psi)^2/2 = 9*Psi^2/2")


# ============================================================
# Step 3: dCPsi/dt at the cusp from weights
# ============================================================

pr()
pr("=" * 72)
pr("STEP 3: dCPsi/dt at the cusp from sector weights")
pr("=" * 72)
pr()
pr("General: dW_k/dt = -4*gamma*k*W_k (exponential decay per sector)")
pr("  dW_0/dt = 0 (stationary)")
pr("  dW_2/dt = -8*gamma*W_2")
pr()
pr("dC/dt = dW_0/dt + dW_2/dt = -8*gamma*W_2")
pr("dPsi/dt = -4*gamma*f/3 = -4*gamma*Psi (since Psi = f/3)")
pr()
pr("dCPsi/dt = dC/dt * Psi + C * dPsi/dt")
pr("       = -8*gamma*W_2 * Psi + (1/2 + W_2) * (-4*gamma*Psi)")
pr("       = -4*gamma*Psi * [2*W_2 + 1/2 + W_2]")
pr("       = -4*gamma*Psi * (1/2 + 3*W_2)")
pr("       = -2*gamma*Psi * (1 + 6*W_2)")
pr()
pr("THIS IS THE WEIGHT-BASED FORMULA:")
pr("    dCPsi/dt = -2*gamma * Psi * (1 + 6*W_2)")
pr()

# Verify numerically
f_cross = brentq(lambda f: f * (1 + f**2) - 1.5, 0.5, 1.0)
t_cross = -np.log(f_cross) / (4 * gamma)
W2_cross = f_cross**2 / 2
Psi_cross = f_cross / 3
C_cross = 0.5 + W2_cross

dcpsi_formula25 = -2 * gamma * f_cross * (1 + 3 * f_cross**2) / 3
dcpsi_weights = -2 * gamma * Psi_cross * (1 + 6 * W2_cross)

pr(f"At the crossing (CPsi = 1/4):")
pr(f"  f_cross = {f_cross:.10f}")
pr(f"  W_2     = {W2_cross:.10f}")
pr(f"  C       = {C_cross:.10f}")
pr(f"  Psi     = {Psi_cross:.10f}")
pr(f"  CPsi      = {C_cross * Psi_cross:.10f} (should be 0.25)")
pr()
pr(f"  dCPsi/dt (Formula 25):       {dcpsi_formula25:.10f}")
pr(f"  dCPsi/dt (weight formula):   {dcpsi_weights:.10f}")
pr(f"  Match: {abs(dcpsi_formula25 - dcpsi_weights) < 1e-14}")
pr()

# The prefactor
pr("At the crossing, Psi = 1/(4C) = 1/(2 + 4*W_2). Substituting:")
pr()
pr("  |dCPsi/dt| = 2*gamma * (1+6*W_2) / (2+4*W_2)")
pr("  K_dwell = 2*delta / |dCPsi/dt| = delta * (2+4*W_2) / (gamma*(1+6*W_2))")
pr()
pr("  PREFACTOR = (2 + 4*W_2) / (1 + 6*W_2)")
pr()

prefactor_weight = (2 + 4 * W2_cross) / (1 + 6 * W2_cross)
prefactor_direct = 2 / abs(dcpsi_formula25)  # = 2/1.851701

pr(f"  Prefactor (weight formula): {prefactor_weight:.10f}")
pr(f"  Prefactor (direct):         {prefactor_direct:.10f}")
pr(f"  Match: {abs(prefactor_weight - prefactor_direct) < 1e-10}")
pr()
pr(f"  The prefactor 1.080088 IS a pure function of W_2 at the crossing.")


# ============================================================
# Step 4: |+>^2 comparison
# ============================================================

pr()
pr("=" * 72)
pr("STEP 4: |+>^2 product state (different Pauli structure)")
pr("=" * 72)
pr()
pr("rho(0) = (1/4)(II + IX + XI + XX)")
pr("  a_II = 1 (w=0), a_IX = 1 (w=1), a_XI = 1 (w=1), a_XX = 1 (w=2)")
pr()
pr("Under Z-dephasing with g = exp(-2*gamma*t):")
pr("  a_II(t) = 1, a_IX(t) = g, a_XI(t) = g, a_XX(t) = g^2")
pr()
pr("Sector weights:")
pr("  W_0 = 1/4           (just II)")
pr("  W_1 = (1/4)(g^2 + g^2) = g^2/2")
pr("  W_2 = (1/4)*g^4 = g^4/4")
pr()
pr("C = (1+g^2)^2/4,  Psi = (2g+g^2)/3")
pr("CPsi = (1+g^2)^2*(2g+g^2)/12")
pr()


def cpsi_plus(g):
    return (1 + g**2)**2 * (2 * g + g**2) / 12


# Find crossing
g_cross = brentq(lambda g: cpsi_plus(g) - 0.25, 0.3, 1.0)
t_cross_plus = -np.log(g_cross) / (2 * gamma)

C_plus = (1 + g_cross**2)**2 / 4
Psi_plus = (2 * g_cross + g_cross**2) / 3
W0_plus = 1 / 4
W1_plus = g_cross**2 / 2
W2_plus = g_cross**4 / 4

pr(f"|+>^2 crossing (CPsi = 1/4):")
pr(f"  g_cross = {g_cross:.10f}")
pr(f"  W_0 = {W0_plus:.10f}")
pr(f"  W_1 = {W1_plus:.10f}")
pr(f"  W_2 = {W2_plus:.10f}")
pr(f"  C   = {C_plus:.10f}")
pr(f"  Psi = {Psi_plus:.10f}")
pr(f"  CPsi  = {C_plus * Psi_plus:.10f}")
pr()

# dCPsi/dt for |+>^2
# dC/dt = -4*gamma*(W_1 + 2*W_2) = -4*gamma*(g^2/2 + g^4/2)
# dPsi/dt = (-4*gamma*g - 4*gamma*g^2)/3 = -4*gamma*g*(1+g)/3
dC_plus = -4 * gamma * (W1_plus + 2 * W2_plus)
dPsi_plus = -4 * gamma * g_cross * (1 + g_cross) / 3
dcpsi_plus = dC_plus * Psi_plus + C_plus * dPsi_plus

# Numerical differentiation to verify
dg = 1e-8
cpsi_p = cpsi_plus(g_cross + dg)
cpsi_m = cpsi_plus(g_cross - dg)
dgdt = -2 * gamma * g_cross
dcpsi_num = (cpsi_p - cpsi_m) / (2 * dg) * dgdt

pr(f"  dCPsi/dt (analytical): {dcpsi_plus:.10f}")
pr(f"  dCPsi/dt (numerical):  {dcpsi_num:.10f}")
pr(f"  Match: {abs(dcpsi_plus - dcpsi_num) < 1e-6}")

prefactor_plus = 2 / abs(dcpsi_plus)
pr()
pr(f"  |+>^2 dwell prefactor: {prefactor_plus:.10f}")
pr(f"  Bell+  dwell prefactor: {prefactor_direct:.10f}")
pr(f"  Ratio: {prefactor_plus / prefactor_direct:.6f}")
pr()

# Can the |+>^2 prefactor be expressed from weights?
# For |+>^2, dCPsi/dt involves dPsi/dt which has sqrt(W_1) terms.
# General: dCPsi/dt = dC/dt * Psi + C * dPsi/dt
#   dC/dt = -4*gamma*(W_1 + 2*W_2)  [weight-only]
#   dPsi/dt = -4*gamma*g*(1+g)/3     [NOT weight-only: needs g, not g^2]
pr("Weight-based analysis for |+>^2:")
pr(f"  dC/dt = -4*gamma*(W_1 + 2*W_2) = {dC_plus:.10f}")
pr(f"    This IS a weight-only quantity.")
pr()
pr(f"  dPsi/dt = -4*gamma*g*(1+g)/3 = {dPsi_plus:.10f}")
pr(f"    This involves g = sqrt(2*W_1), NOT a pure weight function.")
pr(f"    g = {g_cross:.10f}, sqrt(2*W_1) = {np.sqrt(2*W1_plus):.10f}")
pr()
pr("The |+>^2 prefactor requires sqrt(W_1), i.e., the Pauli coefficient")
pr("magnitudes, not just the sector weights (which are sums of squares).")

# Test: does the Bell+ weight formula generalize?
# For Bell+: prefactor = (2+4W_2)/(1+6W_2)
# What would this give for |+>^2?
prefactor_bell_formula = (2 + 4 * W2_plus) / (1 + 6 * W2_plus)
pr()
pr(f"Bell+ weight formula applied to |+>^2 weights:")
pr(f"  (2+4*W_2)/(1+6*W_2) = {prefactor_bell_formula:.10f}")
pr(f"  Actual |+>^2 prefactor: {prefactor_plus:.10f}")
pr(f"  These DIFFER: the Bell+ formula is state-specific (k=1 absent).")


# ============================================================
# Summary
# ============================================================

pr()
pr("=" * 72)
pr("SUMMARY")
pr("=" * 72)
pr()
pr("For Bell+ (only k=0 and k=2 Pauli content):")
pr("  dCPsi/dt = -2*gamma*Psi*(1 + 6*W_2)")
pr(f"  Prefactor = (2+4*W_2)/(1+6*W_2) = {prefactor_weight:.6f}")
pr("  This IS a pure weight function. SUCCESS for Bell+.")
pr()
pr("For |+>^2 (k=0, k=1, and k=2 content):")
pr(f"  Prefactor = {prefactor_plus:.6f}")
pr("  dPsi/dt involves sqrt(W_1), not W_1 itself.")
pr("  The prefactor requires Pauli coefficient magnitudes,")
pr("  not just sector weights. PARTIAL for |+>^2.")
pr()
pr("General verdict: PARTIAL SUCCESS.")
pr("  The static/dynamic connection is algebraic for states with")
pr("  only even-weight Pauli content (Bell+, GHZ). For states with")
pr("  odd-weight content (product states), the L1-coherence Psi")
pr("  introduces square roots of individual coefficients that are")
pr("  not reconstructible from sector weights alone.")
pr()
pr("What the prefactor depends on:")
pr("  1. The sector weights W_k at the crossing (algebraic)")
pr("  2. The relationship between Psi and the W_k (state-specific)")
pr("  The 'extra structure' beyond weights is: the sign pattern of")
pr("  Pauli coefficients determines L1-coherence, which is not")
pr("  reconstructible from W_k = sum |a_P|^2 alone.")

# Write output
out_path = os.path.join(OUT_DIR, "dwell_prefactor_from_weights.txt")
with open(out_path, "w", encoding="utf-8") as f:
    f.write("\n".join(lines))
pr(f"\nOutput written to {out_path}")
print("DONE")
