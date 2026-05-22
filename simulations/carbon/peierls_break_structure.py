"""Characterize the BROKEN Peierls spectrum (docs/carbon follow-up).

BENZENE_LIOUVILLIAN_PALINDROME.md established THAT the F1 palindrome breaks when
the vibrational bath couples to the C-C bond (Peierls/SSH, jump operator the bond
bilinear B_b = X_aX_b + Y_aY_b) rather than on-site (Holstein, jump operator Z_l).
Residual ~ 11 at C4, ~ 14 at C6 against a ~1e-7 floor for Holstein.

That doc only proved the break. This script asks for its STRUCTURE:

  Q1  gamma-scaling of the break. F1's one known analytic breaker, depolarising
      noise, leaves a residual LINEAR in gamma (F1PalindromeIdentity.cs: "(2/3)
      Sigma_gamma, linear in gamma and N"). D[sqrt(gamma)*B] = gamma*D[B] exactly,
      so the dissipator is gamma-linear at the operator level; the spectrum is not
      obviously so. Scan gamma over decades, fit the residual growth law.

  Q2  structure of the broken spectrum. Total break, or a palindromic residual
      core? A different reflection centre? A different involution (lambda -> -lambda
      with no shift; complex conjugation)? Outliers or spread?

  Q3  framework-classification fit. F87 sorts palindrome behaviour into
      truly / soft / hard (ANALYTICAL_FORMULAS.md F87): truly = operator M=0,
      soft = M!=0 but Spec(L) still pairs under lambda -> -lambda-2*Sigma_gamma,
      hard = spectral pairing fails. Known F1-Brechers: depolarising (gamma-linear),
      T1 / amplitude damping, transverse field. Peierls B_b = XX+YY is itself a
      truly-class bond bilinear. Which Brecher, if any, does it resemble?

Reuses the Liouvillian scaffolding of simulations/carbon/benzene_liouvillian_palindrome.py
(site, bond_op, commutator, dissipator). Investigation only; records nothing.
"""
import sys
import numpy as np

sys.stdout.reconfigure(encoding="utf-8")

# ---------------------------------------------------------------- scaffolding
I2 = np.array([[1, 0], [0, 1]], dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.array([[1, 0], [0, -1]], dtype=complex)


def site(op, l, N):
    m = np.array([[1.0 + 0j]])
    for i in range(N):
        m = np.kron(m, op if i == l else I2)
    return m


def bond_op(a, b, N):                       # XX + YY on one C-C bond
    return site(X, a, N) @ site(X, b, N) + site(Y, a, N) @ site(Y, b, N)


def commutator(A, d):                       # -i[A, .] as a superoperator
    eye = np.eye(d)
    return -1j * (np.kron(A, eye) - np.kron(eye, A.T))


def dissipator(jump, d):                    # D[L] = L.L+ - 1/2 {L+ L, .}
    eye = np.eye(d)
    ld_l = jump.conj().T @ jump
    return (np.kron(jump, jump.conj())
            - 0.5 * np.kron(ld_l, eye)
            - 0.5 * np.kron(eye, ld_l.T))


def ring_pieces(N):
    """Return (d, H, L_H, D_site_unit, D_bond_unit) for the C_N ring.

    D_*_unit are the gamma=1 dissipators; D[sqrt(gamma)*J] = gamma*D[J] so the
    gamma=g Liouvillian is L_H + g * D_unit. Holstein uses Z_l, Peierls uses B_b.
    """
    d = 2 ** N
    bonds = [(l, (l + 1) % N) for l in range(N)]
    H = sum(bond_op(a, b, N) for a, b in bonds)
    L_H = commutator(H, d)
    D_site_unit = sum(dissipator(site(Z, l, N), d) for l in range(N))
    D_bond_unit = sum(dissipator(bond_op(a, b, N), d) for a, b in bonds)
    return d, H, L_H, D_site_unit, D_bond_unit


def palindrome_residual(ev, centre):
    """Largest distance from a reflected eigenvalue (2*centre - lambda) to the
    nearest actual eigenvalue. centre is complex; ev is the spectrum array."""
    reflected = 2.0 * complex(centre) - ev
    # vectorised nearest-neighbour over the full spectrum
    diff = np.abs(ev[None, :] - reflected[:, None])
    return float(diff.min(axis=1).max())


def matched_pairs_residuals(ev, centre):
    """Greedy one-to-one matching of each eigenvalue to its reflected partner;
    returns the per-eigenvalue partner distance array. A palindromic core shows
    up as a cluster of near-zero entries."""
    reflected = 2.0 * complex(centre) - ev
    n = len(ev)
    used = np.zeros(n, dtype=bool)
    out = np.zeros(n)
    for i in range(n):
        if used[i]:
            continue
        d = np.abs(ev - reflected[i])
        d[used] = np.inf
        j = int(np.argmin(d))
        out[i] = out[j] = d[j]
        used[i] = True
        if j != i:
            used[j] = True
    return out


# ============================================================ Q1: gamma scan
print("=" * 78)
print("Q1  gamma-SCALING OF THE PEIERLS BREAK")
print("=" * 78)
print("D[sqrt(g)*B] = g*D[B] exactly; the Liouvillian is L_H + g*D_bond_unit.")
print("As g->0, L -> -i[H,.], palindromic about 0, so residual -> 0.")
print("Holstein centre is -Sigma_gamma = -N*g; Peierls has no F1 centre, so we")
print("test against the spectrum mean (the most generous palindrome test).\n")

gammas = np.array([0.001, 0.003, 0.01, 0.03, 0.1, 0.3, 1.0, 3.0])

for N in (4, 6):
    d, H, L_H, D_site_unit, D_bond_unit = ring_pieces(N)
    print(f"--- C{N} ring  (d={d}, L is {d * d}x{d * d}) ---")
    print(f"{'gamma':>8} {'Holstein resid':>16} {'Peierls resid':>16} "
          f"{'Peierls/gamma':>14} {'Peierls/gamma^2':>16}")
    hol_res, pei_res = [], []
    for g in gammas:
        ev_h = np.linalg.eigvals(L_H + g * D_site_unit)
        ev_p = np.linalg.eigvals(L_H + g * D_bond_unit)
        r_h = palindrome_residual(ev_h, -N * g)            # strict F1 centre
        r_p = palindrome_residual(ev_p, ev_p.mean())       # generous centre
        hol_res.append(r_h)
        pei_res.append(r_p)
        print(f"{g:>8.3f} {r_h:>16.3e} {r_p:>16.3e} "
              f"{r_p / g:>14.4f} {r_p / g**2:>16.4f}")
    # log-log slope fit on the small-gamma half (away from any saturation)
    g_arr = gammas
    pei = np.array(pei_res)
    lo = g_arr <= 0.1
    slope_lo, icpt_lo = np.polyfit(np.log(g_arr[lo]), np.log(pei[lo]), 1)
    slope_all, _ = np.polyfit(np.log(g_arr), np.log(pei), 1)
    print(f"  log-log fit  Peierls residual ~ gamma^p :  "
          f"p = {slope_lo:.4f} (gamma<=0.1),  p = {slope_all:.4f} (all gamma)")
    print(f"  => prefactor C in resid ~ C*gamma^p (small-gamma): "
          f"C = {np.exp(icpt_lo):.4f}\n")

# ====================================================== Q2: spectrum structure
print("=" * 78)
print("Q2  STRUCTURE OF THE BROKEN PEIERLS SPECTRUM")
print("=" * 78)
print("At gamma=1, C6: is the break total, or is a palindromic core left?")
print("Test several involutions and look at the per-eigenvalue partner gaps.\n")

N = 6
g = 1.0
d, H, L_H, D_site_unit, D_bond_unit = ring_pieces(N)
L_pei = L_H + g * D_bond_unit
ev = np.linalg.eigvals(L_pei)
Sigma = N * g

print(f"C{N} Peierls, gamma={g}: {len(ev)} eigenvalues")
print(f"  spectrum mean              = {ev.mean():.6f}")
print(f"  Re range [{ev.real.min():.4f}, {ev.real.max():.4f}]   "
      f"width {ev.real.max() - ev.real.min():.4f}")
print(f"  Im range [{ev.imag.min():.4f}, {ev.imag.max():.4f}]\n")

involutions = [
    ("F1 strict   lambda -> -lambda - 2*Sigma_gamma", -Sigma),
    ("about 0     lambda -> -lambda",                  0.0),
    ("about mean  lambda -> 2*mean - lambda",          ev.mean()),
]
for label, centre in involutions:
    r = palindrome_residual(ev, centre)
    print(f"  {label:48s} residual = {r:.4e}")

# conjugation symmetry: a Lindbladian spectrum is closed under complex conj
# iff L is real-representable; check it as a separate involution
conj_resid = palindrome_residual(ev, None) if False else None
diff_conj = np.abs(ev[None, :] - ev.conj()[:, None])
conj_resid = float(diff_conj.min(axis=1).max())
print(f"  {'conjugation lambda -> conj(lambda)':48s} residual = {conj_resid:.4e}")

# best real centre: minimise the residual over the real axis
from scipy.optimize import minimize_scalar
res_opt = minimize_scalar(lambda c: palindrome_residual(ev, c),
                          bounds=(ev.real.min() - Sigma, ev.real.max()),
                          method="bounded")
print(f"\n  best real reflection centre = {res_opt.x:.6f} "
      f"(residual {res_opt.fun:.4e});  F1 would put it at {-Sigma:.4f}")

# palindromic core: per-eigenvalue partner gaps under the best centre
gaps = matched_pairs_residuals(ev, res_opt.x)
floor = 1e-6
core = int(np.sum(gaps < floor))
print(f"\n  per-eigenvalue partner gaps under the best centre:")
print(f"    eigenvalues with a partner within {floor:.0e} : "
      f"{core} / {len(ev)}  (a palindromic core, if > 0)")
print(f"    gap distribution: min {gaps.min():.3e}, median "
      f"{np.median(gaps):.3e}, max {gaps.max():.3e}")
# how localised is the break: sort gaps, see if a few dominate
gs = np.sort(gaps)[::-1]
print(f"    largest 8 partner gaps: "
      f"{', '.join(f'{x:.3f}' for x in gs[:8])}")
print(f"    fraction of eigenvalues with gap > 1% of spectral width "
      f"({0.01 * (ev.real.max() - ev.real.min()):.4f}): "
      f"{np.mean(gaps > 0.01 * (ev.real.max() - ev.real.min())):.3f}")

# is the break carried by specific eigenvalue regions? bucket by Re(lambda)
print(f"\n  break localisation -- mean partner gap by Re(lambda) decile:")
order = np.argsort(ev.real)
ev_sorted = ev[order]
gaps_sorted = matched_pairs_residuals(ev, res_opt.x)[order]
ndec = 10
for k in range(ndec):
    sl = slice(k * len(ev) // ndec, (k + 1) * len(ev) // ndec)
    re_lo, re_hi = ev_sorted[sl].real.min(), ev_sorted[sl].real.max()
    print(f"    decile {k}: Re in [{re_lo:8.3f}, {re_hi:8.3f}]  "
          f"mean gap {gaps_sorted[sl].mean():.4f}")

# zero-mode / steady-state check: a Lindbladian always has a 0 eigenvalue
zero_idx = np.argmin(np.abs(ev))
print(f"\n  closest-to-zero eigenvalue (steady state): {ev[zero_idx]:.3e}")

# ================================================ Q3: classification fit
print()
print("=" * 78)
print("Q3  FRAMEWORK-CLASSIFICATION FIT  (F87 truly / soft / hard)")
print("=" * 78)
print("F87 (ANALYTICAL_FORMULAS.md): build M = Pi.L.Pi^-1 + L + 2*Sigma*I.")
print("  truly = ||M||_F < eps                   (operator identity holds)")
print("  soft  = ||M|| >= eps but Spec(L) still pairs under -lambda-2*Sigma")
print("  hard  = spectral pairing under -lambda-2*Sigma fails")
print("We do not need Pi explicitly: F87's soft/hard split is a STATEMENT ABOUT")
print("THE SPECTRUM, namely whether Spec(L) closes under lambda->-lambda-2*Sigma.")
print("That is exactly the F1-strict residual already computed above.\n")

op_tol = 1e-10
spec_tol = 1e-6
for N in (4, 6):
    d, H, L_H, D_site_unit, D_bond_unit = ring_pieces(N)
    Sigma = N * 1.0
    for cname, L in [("Holstein D[Z_l]", L_H + 1.0 * D_site_unit),
                     ("Peierls  D[B_b]", L_H + 1.0 * D_bond_unit)]:
        ev_c = np.linalg.eigvals(L)
        r_f1 = palindrome_residual(ev_c, -Sigma)
        # F87 verdict from the spectral-pairing test (the soft/hard discriminator)
        if r_f1 < spec_tol:
            # could be truly or soft; truly needs ||M||=0 which for the
            # Holstein case F1 guarantees. We label by the spectral test and
            # note the operator side separately.
            verdict = "truly/soft (Spec closes under F1 involution)"
        else:
            verdict = "HARD  (Spec does NOT close under F1 involution)"
        print(f"  C{N} {cname}: F1-strict residual = {r_f1:.3e}  -> {verdict}")
print()

# F87-hard's defining trait: contrast with depolarising noise, F1's known
# analytic Brecher, which is hard with a gamma-LINEAR residual. We compute the
# depolarising residual on the SAME ring for a like-for-like gamma-scaling
# comparison: depolarising = (gamma/3) * (D[X_l] + D[Y_l] + D[Z_l]) per site.
print("-" * 78)
print("Like-for-like: depolarising noise (F1's known analytic Brecher) on the")
print("same C4 ring. Depolarising = (g/3) sum_l (D[X_l]+D[Y_l]+D[Z_l]).")
print("F1PalindromeIdentity.cs: depolarising residual is linear in gamma.\n")
N = 4
d, H, L_H, D_site_unit, D_bond_unit = ring_pieces(N)
D_depol_unit = sum(
    (1.0 / 3.0) * (dissipator(site(X, l, N), d)
                   + dissipator(site(Y, l, N), d)
                   + dissipator(site(Z, l, N), d))
    for l in range(N)
)
print(f"{'gamma':>8} {'depol resid':>14} {'depol/gamma':>14} "
      f"{'Peierls resid':>14} {'Peierls/gamma':>14}")
dep_res, pei_res2 = [], []
for g in gammas:
    ev_dep = np.linalg.eigvals(L_H + g * D_depol_unit)
    ev_pei = np.linalg.eigvals(L_H + g * D_bond_unit)
    r_dep = palindrome_residual(ev_dep, -N * g)        # F1 centre
    r_pei = palindrome_residual(ev_pei, ev_pei.mean())
    dep_res.append(r_dep)
    pei_res2.append(r_pei)
    print(f"{g:>8.3f} {r_dep:>14.4e} {r_dep / g:>14.5f} "
          f"{r_pei:>14.4e} {r_pei / g:>14.5f}")
dep = np.array(dep_res)
lo = gammas <= 0.1
sdep, _ = np.polyfit(np.log(gammas[lo]), np.log(dep[lo]), 1)
print(f"\n  depolarising log-log slope (gamma<=0.1): p = {sdep:.4f}  "
      f"(F1 caveat predicts p=1, linear)")
print("  => compare this exponent to the Peierls exponent from Q1.")
