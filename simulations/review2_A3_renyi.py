"""
Independent from-below verifier for review finding A3 (UNIQUENESS_PROOF.md, Step 5).

Two questions, attacked symbolically (sympy), refute-first:

Q1 (the alleged non-sequitur): does "purity is degree-2 => quadratic recursion =>
    discriminant gives 1/4" actually FORCE 1/4, or only after assuming the specific
    recursion form R = C(Psi+R)^2 with those coefficients?

Q2 (the alleged genuine forcing): does Renyi alpha=2 STATE-INDEPENDENCE single out
    1/4? I derive the fold-bifurcation threshold for the generalized recursion
    R = C(Psi+R)^alpha FROM SCRATCH and check (a) it matches the roadmap's claimed
    CPsi*_alpha = (alpha-1)^(alpha-1)/(alpha^alpha * Psi^(alpha-2)), (b) it equals 1/4
    at alpha=2, (c) it is Psi-independent ONLY at alpha=2.
"""
import sympy as sp

C, Psi, R, x, a = sp.symbols('C Psi R x alpha', positive=True)

print("=" * 70)
print("Q1: the degree-2 -> 1/4 chain")
print("=" * 70)
# The recursion as written in UNIQUENESS_PROOF: R = C*(Psi+R)^2
lhs = C * (Psi + R)**2
quad = sp.expand(lhs - R)                       # = C R^2 + (2C Psi - 1) R + C Psi^2
print("R = C(Psi+R)^2  =>  0 =", quad)
poly = sp.Poly(quad, R)
A_, B_, Cc_ = poly.all_coeffs()
disc = sp.simplify(B_**2 - 4 * A_ * Cc_)
print("coeffs (a,b,c) =", (A_, B_, Cc_))
print("discriminant b^2-4ac =", disc, " => vanishes at CPsi =",
      sp.solve(sp.Eq(disc, 0), C * Psi))
print()
print("Observation: the 1/4 is correct GIVEN this exact recursion + coefficients.")
print("But degree-2-ness of purity does NOT by itself fix the coefficient structure")
print("(2C*Psi-1 linear term, C*Psi^2 constant). A generic degree-2 fixed-point map")
print("a*R^2+b*R+c=0 has its discriminant vanish at b^2=4ac, an arbitrary locus, not 1/4.")
# demonstrate with a generic quadratic whose discriminant-zero locus is NOT 1/4
b1, c1 = sp.symbols('b1 c1')
print("   e.g. R^2 + b1 R + c1 = 0 -> disc 0 at b1^2 = 4 c1 (no 1/4 unless coeffs are the R=C(Psi+R)^2 ones)")

print("\n" + "=" * 70)
print("Q2: derive the fold threshold for R = C(Psi+R)^alpha from scratch")
print("=" * 70)
# Fixed point f(R)=R and fold (tangency) f'(R)=1, with f(R)=C(Psi+R)^alpha.
f = C * (Psi + R)**a
fixed = sp.Eq(f, R)                 # C(Psi+R)^alpha = R
fold = sp.Eq(sp.diff(f, R), 1)      # C*alpha*(Psi+R)^(alpha-1) = 1
print("fixed-point: ", fixed)
print("fold/tangency:", fold)

# From fold: (Psi+R)^(alpha-1) = 1/(C*alpha). Substitute into fixed point.
# Cleanest: let u = Psi+R. fixed: C u^alpha = u - Psi ; fold: C alpha u^(alpha-1)=1.
u = sp.symbols('u', positive=True)
fixed_u = sp.Eq(C * u**a, u - Psi)
fold_u = sp.Eq(C * a * u**(a - 1), 1)
# from fold_u: C = 1/(alpha u^(alpha-1)); plug into fixed_u
C_of_u = sp.solve(fold_u, C)[0]
fixed_sub = sp.simplify(fixed_u.lhs.subs(C, C_of_u) - fixed_u.rhs)
# solve for u in terms of Psi
u_sol = sp.solve(sp.Eq(fixed_sub, 0), u)
print("C(u) from fold =", C_of_u)
print("u at threshold =", u_sol)
u_star = u_sol[0]
C_star = sp.simplify(C_of_u.subs(u, u_star))
CPsi_star = sp.simplify(C_star * Psi)
print("C* =", C_star)
print("CPsi* =", sp.simplify(CPsi_star))

# Compare to the roadmap's claimed formula
claimed = (a - 1)**(a - 1) / (a**a * Psi**(a - 2))
print("\nroadmap claims CPsi*_alpha =", claimed)
# sympy's simplify struggles with symbolic powers; force powsimp + numeric spot-checks
diff_forced = sp.simplify(sp.powsimp(sp.expand_power_base(CPsi_star / claimed, force=True), force=True))
print("CPsi* / claimed  (force-simplified, =1 iff identical) =", diff_forced)
print("numeric spot-checks (CPsi* vs claimed):")
for av, pv in [(2, 0.3), (3, 0.5), (4, 0.7), (sp.Rational(5, 2), 0.4)]:
    lhs_v = float(CPsi_star.subs({a: av, Psi: pv}))
    rhs_v = float(claimed.subs({a: av, Psi: pv}))
    print(f"   alpha={av}, Psi={pv}: CPsi*={lhs_v:.8f}, claimed={rhs_v:.8f}, equal={abs(lhs_v-rhs_v)<1e-12}")

print("\n  alpha=2  -> CPsi* =", sp.simplify(CPsi_star.subs(a, 2)))
print("  alpha=3  -> CPsi* =", sp.simplify(CPsi_star.subs(a, 3)))
print("  alpha=4  -> CPsi* =", sp.simplify(CPsi_star.subs(a, 4)))

# Psi-independence: dCPsi*/dPsi == 0 identically only when alpha=2
dPsi = sp.simplify(sp.diff(CPsi_star, Psi))
print("\nd(CPsi*)/dPsi =", dPsi)
sol_alpha = sp.solve(sp.Eq(dPsi, 0), a)
print("alpha making CPsi* Psi-independent (dCPsi*/dPsi=0):", sol_alpha)
print("=> only alpha=2 gives a state-independent threshold, and there CPsi*=1/4.")
print("\nDONE.")
