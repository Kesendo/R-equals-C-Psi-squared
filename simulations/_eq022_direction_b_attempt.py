"""Direction (b) doubled-PTF Ansatz exploration: c=2 HWHM_left/Q_peak closed-form attempt.

Empirical anchors (PROOF_F86_QPEAK Statement 2 across N=5..8 at gamma_0=0.05):
  Interior: 0.7506 (mean), Q_peak ~ 1.48..1.60
  Endpoint: 0.7728 (mean), Q_peak ~ 2.50..2.55

Doubled-PTF model:
  - Channel A: HD=1, rate -2 gamma_0
  - Channel B: HD=3, rate -6 gamma_0
  - Inter-channel coupling: sigma_0 (= g_eff)
  - 2-level L_2 = [[-2 gamma_0, +i J g_eff], [+i J g_eff, -6 gamma_0]] with J = Q gamma_0
  - Q_EP = 2/g_eff (eigenvalues coalesce at Q = Q_EP)
  - t_peak = 1/(4 gamma_0)

The bond-class signature comes from how the dL/dJ_b enters the Duhamel formula --
this is the V_b cross-block in the 4-mode picture. We want to derive HWHM_left/Q_peak.
"""

import numpy as np


def compute_hwhm_ratio_from_curve(K, Q_grid):
    """Find Q_peak and HWHM_left from a precomputed |K|(Q) curve."""
    K_abs = np.abs(K)
    i_max = int(np.argmax(K_abs))
    K_max = K_abs[i_max]
    if i_max == 0 or i_max == len(Q_grid) - 1:
        # Peak at boundary; bail
        return Q_grid[i_max], K_max, 0.0, 0.0
    # Parabolic interpolation
    y0, y1, y2 = K_abs[i_max - 1], K_abs[i_max], K_abs[i_max + 1]
    denom = y0 - 2 * y1 + y2
    delta = 0.5 * (y0 - y2) / denom if abs(denom) > 1e-15 else 0.0
    dQ = Q_grid[1] - Q_grid[0]
    Q_peak = Q_grid[i_max] + delta * dQ

    half = K_max / 2.0
    hwhm_left = None
    for j in range(i_max, 0, -1):
        if K_abs[j-1] < half <= K_abs[j]:
            frac = (half - K_abs[j-1]) / (K_abs[j] - K_abs[j-1])
            Q_half = Q_grid[j-1] + frac * (Q_grid[j] - Q_grid[j-1])
            hwhm_left = Q_peak - Q_half
            break
    if hwhm_left is None:
        hwhm_left = Q_peak - Q_grid[0]
    return Q_peak, K_max, hwhm_left, hwhm_left / Q_peak if Q_peak > 0 else 0.0


def k_b_two_level_duhamel(Q, gamma0, g_eff, t):
    """Compute the K_b observable in the simplest 2-level Duhamel model.

    L_2(Q) = [[-2 gamma_0, +i J g_eff], [+i J g_eff, -6 gamma_0]] with J = Q gamma_0.
    Probe rho_0 = [1, 0] (in the slow channel).
    V_b = derivative dL/dJ = [[0, +i g_eff], [+i g_eff, 0]]
    S_kernel = identity (probe-block)
    K_b = 2 Re <rho(t) | S | drho/dJ>
    """
    J = Q * gamma0
    L = np.array([
        [-2.0 * gamma0, 1j * J * g_eff],
        [1j * J * g_eff, -6.0 * gamma0],
    ], dtype=complex)
    Vb = np.array([
        [0.0, 1j * g_eff],
        [1j * g_eff, 0.0],
    ], dtype=complex)
    evals, R = np.linalg.eig(L)
    Rinv = np.linalg.inv(R)
    rho0 = np.array([1.0, 0.0], dtype=complex)

    expLam = np.exp(evals * t)
    expL = R @ np.diag(expLam) @ Rinv
    rho_t = expL @ rho0

    # Duhamel kernel I_jk = (exp(lam_k t) - exp(lam_j t))/(lam_k - lam_j)
    n = 2
    I_mat = np.zeros((n, n), dtype=complex)
    for j in range(n):
        for k in range(n):
            diff = evals[k] - evals[j]
            if abs(diff) < 1e-10:
                I_mat[j, k] = t * expLam[j]
            else:
                I_mat[j, k] = (expLam[k] - expLam[j]) / diff

    X = Rinv @ Vb @ R
    c0 = Rinv @ rho0
    fbC0 = np.zeros(n, dtype=complex)
    for r in range(n):
        for c in range(n):
            fbC0[r] += X[r, c] * I_mat[r, c] * c0[c]
    drho = R @ fbC0
    # S = I (identity, probe-block)
    inner = np.vdot(rho_t, drho)
    return 2.0 * inner.real


def explore_two_level_duhamel_K_b():
    print("=" * 72)
    print("2-level Duhamel K_b observable: HWHM_left/Q_peak")
    print("=" * 72)
    gamma0 = 0.05
    t_peak = 1.0 / (4.0 * gamma0)
    Q_grid = np.linspace(0.05, 5.0, 2000)

    g_eff_values = [
        2.0, 2.0 * np.sqrt(2.0),
        2.765, 2.802, 2.828, 2.839,  # actual sigma_0 at c=2 N=5..8
    ]
    for g in g_eff_values:
        Q_EP = 2.0 / g
        K_curve = np.array([k_b_two_level_duhamel(Q, gamma0, g, t_peak) for Q in Q_grid])
        Q_peak, K_max, hwhm_left, ratio = compute_hwhm_ratio_from_curve(K_curve, Q_grid)
        print(f"  g_eff={g:.4f}, Q_EP={Q_EP:.4f}, Q_peak={Q_peak:.4f}, "
              f"|K|max={K_max:.4f}, HWHM-/Q*={ratio:.4f}, x_peak={Q_peak/Q_EP:.4f}")

    # Universal in dimensionless x = Q/Q_EP?
    print("\n  Same in dimensionless x = Q/Q_EP (gamma_0=1, g_eff=2 -> Q_EP=1):")
    gamma0 = 1.0
    t_peak = 1.0 / (4.0 * gamma0)
    Q_grid = np.linspace(0.005, 8.0, 4000)
    K_curve = np.array([k_b_two_level_duhamel(Q, gamma0, 2.0, t_peak) for Q in Q_grid])
    Q_peak, K_max, hwhm_left, ratio = compute_hwhm_ratio_from_curve(K_curve, Q_grid)
    print(f"  g_eff=2.0, Q_EP=1.0, Q_peak={Q_peak:.6f}, |K|max={K_max:.6f}, "
          f"HWHM-/Q*={ratio:.6f}")
    print(f"  -> x_peak = Q_peak/Q_EP = {Q_peak:.6f}")
    print(f"  -> HWHM_left = {hwhm_left:.6f}")
    print(f"  -> HWHM_left/x_peak = {ratio:.6f}")


def k_four_mode_doubled_ptf(Q, gamma0, g_eff, cb_alpha, cb_beta, t,
                            probe_idx=0, S_pattern='diag1100'):
    """4-mode doubled-PTF with cross-block coupling cb_alpha (top-right) and cb_beta (off).

    Basis: [|c_1>, |c_3>, |u_0>, |v_0>]
    D_eff = diag(-2g0, -6g0, -2g0, -6g0)
    Probe rho_0 = e_probe_idx
    M_H_eff in this basis:
      Probe-block (top-left): zero (channel-uniform basis F73)
      SVD-block (bottom-right): [[0, +i sigma_0], [+i sigma_0, 0]]
      Cross-block (top-right + symmetric): determined by V_b at the bond
    V_b = M_H_per_bond[b] in the 4-mode basis
      Probe-block: zero (single bond doesn't contribute uniformly)
      SVD-block: contributes a portion of sigma_0
      Cross-block: this is what differs by bond class
    """
    # Build M_H_eff at unit J
    M_eff = np.zeros((4, 4), dtype=complex)
    M_eff[2, 3] = 1j * g_eff
    M_eff[3, 2] = 1j * g_eff
    # cross-block (top-right):
    M_eff[0, 2] = 1j * cb_alpha
    M_eff[2, 0] = 1j * cb_alpha
    M_eff[0, 3] = 1j * cb_beta
    M_eff[3, 0] = 1j * cb_beta
    M_eff[1, 2] = 1j * cb_beta
    M_eff[2, 1] = 1j * cb_beta
    M_eff[1, 3] = 1j * cb_alpha
    M_eff[3, 1] = 1j * cb_alpha

    D_eff = np.diag([-2*gamma0, -6*gamma0, -2*gamma0, -6*gamma0]).astype(complex)
    L = D_eff + Q * gamma0 * M_eff

    Vb = M_eff  # using the same matrix as bond coupling for this toy

    evals, R = np.linalg.eig(L)
    Rinv = np.linalg.inv(R)
    rho0 = np.zeros(4, dtype=complex)
    rho0[probe_idx] = 1.0

    expLam = np.exp(evals * t)
    rho_t = R @ np.diag(expLam) @ (Rinv @ rho0)

    n = 4
    I_mat = np.zeros((n, n), dtype=complex)
    for j in range(n):
        for k in range(n):
            diff = evals[k] - evals[j]
            if abs(diff) < 1e-10:
                I_mat[j, k] = t * expLam[j]
            else:
                I_mat[j, k] = (expLam[k] - expLam[j]) / diff

    X = Rinv @ Vb @ R
    c0 = Rinv @ rho0
    fbC0 = np.zeros(n, dtype=complex)
    for r in range(n):
        for c in range(n):
            fbC0[r] += X[r, c] * I_mat[r, c] * c0[c]
    drho = R @ fbC0

    if S_pattern == 'diag1100':
        S = np.diag([1, 1, 0, 0]).astype(complex)
    elif S_pattern == 'identity':
        S = np.eye(4).astype(complex)
    else:
        S = np.diag([1, 1, 0, 0]).astype(complex)

    sDrho = S @ drho
    inner = np.vdot(rho_t, sDrho)
    return 2.0 * inner.real


def explore_four_mode_cross_block():
    print()
    print("=" * 72)
    print("4-mode doubled-PTF: HWHM/Q* vs cross-block strength")
    print("=" * 72)
    gamma0 = 0.05
    t_peak = 1.0 / (4.0 * gamma0)
    Q_grid = np.linspace(0.05, 5.0, 1500)
    g_eff = 2.0 * np.sqrt(2.0)  # asymptotic sigma_0 for c=2

    # Vary cross-block magnitudes
    print(f"  g_eff = sigma_0 = 2*sqrt(2) = {g_eff:.4f}")
    print(f"  Q_EP = 2/g_eff = {2.0/g_eff:.4f}")
    print()
    print("  Format: cb_alpha (top-right [0,2], [1,3]), cb_beta (top-right [0,3], [1,2])")
    print()

    cases = [
        ("trivial (cb=0)", 0.0, 0.0),
        ("alpha=0.05, beta=0", 0.05, 0.0),
        ("alpha=0.1, beta=0", 0.1, 0.0),
        ("alpha=0.2, beta=0", 0.2, 0.0),
        ("alpha=0, beta=0.05", 0.0, 0.05),
        ("alpha=0, beta=0.1", 0.0, 0.1),
        ("alpha=0, beta=0.2", 0.0, 0.2),
        ("alpha=0.1, beta=0.1", 0.1, 0.1),
    ]
    for label, ca, cb in cases:
        K_curve = np.array([k_four_mode_doubled_ptf(Q, gamma0, g_eff, ca, cb, t_peak)
                            for Q in Q_grid])
        Q_peak, K_max, hwhm_left, ratio = compute_hwhm_ratio_from_curve(K_curve, Q_grid)
        print(f"  {label:30s}: Q_peak={Q_peak:.4f}, |K|max={K_max:.4g}, "
              f"HWHM-/Q*={ratio:.4f}")


def search_universal_x_peak():
    """The bare 2-level K_b model gives x_peak = Q_peak/Q_EP = 2.1969 (universal).

    Let's solve dK/dx = 0 numerically with very high precision and see if x_peak is
    a known constant.
    """
    from scipy.optimize import brentq
    print()
    print("=" * 72)
    print("Find universal x_peak and HWHM ratio for bare 2-level Duhamel K_b")
    print("=" * 72)

    def K_universal(x):
        """K(x = Q/Q_EP) in the universal post-EP coordinate, gamma_0=1, g_eff=2."""
        return k_b_two_level_duhamel(Q=x, gamma0=1.0, g_eff=2.0, t=0.25)

    def absK(x):
        return abs(K_universal(x))

    # Very fine grid search starting from near zero
    x_grid = np.linspace(0.001, 4.0, 100000)
    K_curve = np.array([K_universal(x) for x in x_grid])
    K_abs = np.abs(K_curve)
    i_max = int(np.argmax(K_abs))
    print(f"  Grid argmax x = {x_grid[i_max]:.6f}, K = {K_curve[i_max]:.6f}, |K| = {K_abs[i_max]:.6f}")

    # Refine via Brent on d|K|/dx
    def dabsK(x):
        h = 1e-6
        return (absK(x + h) - absK(x - h)) / (2 * h)

    try:
        x_peak = brentq(dabsK, x_grid[i_max] - 0.005, x_grid[i_max] + 0.005)
    except Exception as ex:
        print(f"  Brent failed ({ex}), using grid argmax")
        x_peak = x_grid[i_max]
    print(f"  Refined x_peak = {x_peak:.10f}")
    K_at_peak = absK(x_peak)
    print(f"  |K|(x_peak) = {K_at_peak:.10f}")

    # HWHM_left: |K|(1.001) might be > half; we need to find a region where |K| crosses half.
    # Walk down the grid from x_peak:
    half = K_at_peak / 2.0
    # Find first index where K_abs[i] < half walking left from i_max
    i_left = i_max
    while i_left > 0 and K_abs[i_left] >= half:
        i_left -= 1
    print(f"  Half-max at index {i_left}, x={x_grid[i_left]:.6f}, |K|={K_abs[i_left]:.6f}")
    print(f"  Half-max threshold = {half:.6f}")
    x_half = brentq(lambda x: absK(x) - half, x_grid[i_left], x_grid[i_left+1])
    print(f"  K(x_half_left) = K(x_peak)/2 at x_half_left = {x_half:.10f}")
    hwhm_left = x_peak - x_half
    print(f"  HWHM_left (in x) = {hwhm_left:.10f}")
    print(f"  HWHM_left / x_peak = {hwhm_left / x_peak:.10f}")
    print()

    # Identify candidates
    print(f"  x_peak = {x_peak:.6f}, candidates:")
    candidates = {
        '1+sqrt(3)': 1 + np.sqrt(3),
        '2+1/sqrt(5)': 2 + 1/np.sqrt(5),
        'pi/sqrt(2)': np.pi/np.sqrt(2),
        '7/sqrt(10)': 7/np.sqrt(10),
        '11/5': 11/5,
        '2+pi/16': 2 + np.pi/16,
        'sqrt(5)-pi/4': np.sqrt(5) - np.pi/4,
        '2.197': 2.197,
        # Can we work it from the algebra: K(x) = exp(-1)*[(cos(s/4) + sin(s/4)/(s/2)) ...]
    }
    for label, val in candidates.items():
        print(f"    {label:20s} = {val:.6f}, diff = {x_peak - val:.6f}")

    # Try: is x_peak related to a root of a specific equation?
    # K(x) = 2 Re <rho(t)|drho/dJ> in 2-level. Let's derive the analytic form.
    print()
    return x_peak, hwhm_left / x_peak


def derive_analytic_2level_K_b():
    """Derive K_b(x) analytically for the 2-level Duhamel model.

    L_2 = [[-2, +i*2x], [+i*2x, -6]] with gamma_0 = 1, g_eff = 2, J = x.
    eigenvalues: lambda_+- = -4 +- 2*sqrt(1 - x^2)
    For x > 1: lambda_+- = -4 +- 2*i*sqrt(x^2 - 1)
    """
    print()
    print("=" * 72)
    print("Analytic K_b(x) for bare 2-level Duhamel (universal post-EP)")
    print("=" * 72)
    # Set t = 1/4. Eigenvalues at -4 +- i*2 sqrt(x^2-1) for x > 1.
    # exp(lambda_+ t) = exp(-1) * exp(i*0.5*sqrt(x^2-1))
    # exp(lambda_- t) = exp(-1) * exp(-i*0.5*sqrt(x^2-1))
    # lambda_+ - lambda_- = 4i*sqrt(x^2-1) =: 2is, where s = 2*sqrt(x^2-1)
    # I_jk = (exp(lam_k t) - exp(lam_j t))/(lam_k - lam_j)
    # I_++ = t exp(lam_+ t), I_-- = t exp(lam_- t)
    # I_+- = (exp(lam_- t) - exp(lam_+ t))/(lam_- - lam_+) = ...

    # Eigenvectors of L = [[-2, +i*2x], [+i*2x, -6]]:
    # (L - lambda_+ I) v = 0
    # (-2 - lambda_+) v_0 + 2ix v_1 = 0
    # v = [2ix, 2 + lambda_+] = [2ix, -2 + 2*sqrt(1-x^2)] (for x < 1)
    # For x > 1: v = [2ix, -2 + 2i*sqrt(x^2-1)] (right eigenvector for lambda_+)
    # Normalize later — just use numerical confirmation.
    import sympy as sp
    x = sp.Symbol('x', real=True, positive=True)
    s = 2 * sp.sqrt(x**2 - 1)
    # For x > 1: lam_+ = -4 + i*s, lam_- = -4 - i*s
    lam_p = -4 + sp.I * s
    lam_m = -4 - sp.I * s
    t = sp.Rational(1, 4)
    e_p = sp.exp(lam_p * t)
    e_m = sp.exp(lam_m * t)
    # Eigenvectors (right):
    # (L - lam_+) v_+ = 0 => -(2 + lam_+) v_0 + 2ix v_1 = 0
    # v_+ = [2ix, 2 + lam_+] = [2ix, -2 + i*s]
    v_p = sp.Matrix([2*sp.I*x, -2 + sp.I*s])
    v_m = sp.Matrix([2*sp.I*x, -2 - sp.I*s])
    R = sp.Matrix([[v_p[0], v_m[0]], [v_p[1], v_m[1]]])
    print(f"  R = ")
    sp.pprint(R)
    Rinv = R.inv()
    print(f"  R inv = ")
    sp.pprint(sp.simplify(Rinv))

    # Probe
    rho0 = sp.Matrix([1, 0])
    c0 = sp.simplify(Rinv * rho0)
    print(f"  c0 = ")
    sp.pprint(c0)

    # rho(t) = R * diag(e_p, e_m) * c0
    expL_diag = sp.Matrix([[e_p, 0], [0, e_m]])
    rho_t = sp.simplify(R * expL_diag * c0)
    print(f"  rho_t = ")
    sp.pprint(rho_t)

    # V_b = dL/dJ = [[0, 2i], [2i, 0]] (with J -> J in L = D + J*M_H, here M_H[0,1]=2i)
    # But wait L = D + Q*g_eff... let me re-set.
    # We have L(x) = D + x * M, where M = [[0, 2i], [2i, 0]] (since J*g_eff with g_eff=2, J=x)
    # so dL/dx = M.
    # K_b = 2 Re <rho(t) | S | drho/dx>
    # drho/dx = sum_jk e_j * I_jk(t) * X[r,c] * c0[c] form... but easier:
    # drho/dx = R * [(R^{-1} dL/dx R) (.) I] * R^{-1} * rho0
    # where I is the kernel matrix.
    # In the eigenbasis dL/dx -> X = Rinv * M * R
    M = sp.Matrix([[0, 2*sp.I], [2*sp.I, 0]])
    X = sp.simplify(Rinv * M * R)
    print(f"  X = R^-1 M R = ")
    sp.pprint(X)

    # I_mat[j,k] = (e_k - e_j)/(lam_k - lam_j)
    I00 = t * e_p
    I11 = t * e_m
    I01 = (e_m - e_p) / (lam_m - lam_p)
    I10 = (e_p - e_m) / (lam_p - lam_m)
    Imat = sp.Matrix([[I00, I01], [I10, I11]])
    print(f"  I = ")
    sp.pprint(sp.simplify(Imat))

    # drho/dx in the eigenbasis: fbC0[r] = sum_c X[r,c] * I[r,c] * c0[c]
    fbC0 = sp.zeros(2, 1)
    for r in range(2):
        for c in range(2):
            fbC0[r] = fbC0[r] + X[r, c] * Imat[r, c] * c0[c]
    fbC0 = sp.simplify(fbC0)
    drho = sp.simplify(R * fbC0)
    print(f"  drho = ")
    sp.pprint(drho)

    # K_b = 2 Re <rho_t | drho> with S = identity (probe-block diagonal)
    inner = sp.simplify(sp.conjugate(rho_t[0]) * drho[0] + sp.conjugate(rho_t[1]) * drho[1])
    K_b_expr = 2 * sp.re(inner)
    K_b_simplified = sp.simplify(K_b_expr)
    print(f"  K_b(x) = ")
    sp.pprint(K_b_simplified)

    # Numerically evaluate at x = 2.197 to verify
    K_at_peak = K_b_simplified.subs(x, 2.197).evalf()
    print(f"  K_b(2.197) = {K_at_peak}")


def empirical_x_peak_check():
    """For the c=2 N=5..8 actual system, what is x_peak = Q_peak/Q_EP_naive?

    Q_EP_naive = 2/sigma_0 where sigma_0 ~ 2sqrt(2) asymptotically.
    """
    print()
    print("=" * 72)
    print("Empirical x_peak check: Q_peak/Q_EP_naive for actual c=2 systems")
    print("=" * 72)

    # From PROOF_F86_QPEAK Statement 2 anchor table:
    cases = [
        # N, sigma_0, Q_peak Interior, Q_peak Endpoint
        (5, 2.765, 1.4821, 2.5008),
        (6, 2.802, 1.5801, 2.5470),
        (7, 2.828, 1.5831, 2.5299),
        (8, 2.839, 1.6049, 2.5145),
    ]
    print("  bare doubled-PTF predicts: x_peak = 2.1969 (universal)")
    print()
    print(f"  {'N':>2}  {'sigma_0':>8}  {'Q_EP':>8}  {'Q_peak Int':>11}  {'Q_peak End':>11}"
          f"  {'x_peak Int':>11}  {'x_peak End':>11}")
    for N, sig, Q_int, Q_end in cases:
        Q_EP = 2.0 / sig
        x_int = Q_int / Q_EP
        x_end = Q_end / Q_EP
        print(f"  {N:>2}  {sig:>8.4f}  {Q_EP:>8.4f}  {Q_int:>11.4f}  {Q_end:>11.4f}"
              f"  {x_int:>11.4f}  {x_end:>11.4f}")

    print()
    print("  Interior x_peak ~ 2.05..2.07: very close to bare doubled-PTF 2.1969")
    print("  Endpoint x_peak ~ 3.45..3.55: well above bare doubled-PTF 2.1969")


def explore_probe_mixing():
    """The probe is the projection of the Dicke probe on the 4-mode basis.

    In actual c=2 systems, the probe has nontrivial weights on both |c_1> and |c_3>.
    Let's see how the bare 2-level HWHM ratio depends on the probe mixing.
    """
    print()
    print("=" * 72)
    print("2-level Duhamel K_b: HWHM_left/Q_peak vs probe mixing angle")
    print("=" * 72)
    gamma0 = 1.0  # dimensionless
    g_eff = 2.0   # so Q_EP = 1
    t_peak = 0.25  # = 1/(4 gamma0)
    Q_grid = np.linspace(0.001, 5.0, 50000)

    def k_b_with_mixing(Q, theta):
        """K_b with probe rho_0 = cos(theta)|c_1> + sin(theta)|c_3>."""
        J = Q * gamma0
        L = np.array([
            [-2.0 * gamma0, 1j * J * g_eff],
            [1j * J * g_eff, -6.0 * gamma0],
        ], dtype=complex)
        Vb = np.array([
            [0.0, 1j * g_eff],
            [1j * g_eff, 0.0],
        ], dtype=complex)
        evals, R = np.linalg.eig(L)
        Rinv = np.linalg.inv(R)
        rho0 = np.array([np.cos(theta), np.sin(theta)], dtype=complex)
        expLam = np.exp(evals * t_peak)
        rho_t = R @ np.diag(expLam) @ (Rinv @ rho0)
        n = 2
        I_mat = np.zeros((n, n), dtype=complex)
        for j in range(n):
            for k in range(n):
                diff = evals[k] - evals[j]
                if abs(diff) < 1e-10:
                    I_mat[j, k] = t_peak * expLam[j]
                else:
                    I_mat[j, k] = (expLam[k] - expLam[j]) / diff
        X = Rinv @ Vb @ R
        c0 = Rinv @ rho0
        fbC0 = np.zeros(n, dtype=complex)
        for r in range(n):
            for c in range(n):
                fbC0[r] += X[r, c] * I_mat[r, c] * c0[c]
        drho = R @ fbC0
        # S = identity (probe-block = full)
        return 2.0 * np.vdot(rho_t, drho).real

    # Sweep theta
    print(f"  {'theta':>8}  {'cos':>6}  {'sin':>6}  {'Q_peak':>8}  {'|K|max':>10}  {'HWHM-/Q*':>10}")
    for theta_deg in [0, 15, 30, 45, 60, 75, 90]:
        theta = np.radians(theta_deg)
        K_curve = np.array([k_b_with_mixing(Q, theta) for Q in Q_grid])
        Q_peak, K_max, hwhm_left, ratio = compute_hwhm_ratio_from_curve(K_curve, Q_grid)
        print(f"  {theta_deg:>8}  {np.cos(theta):>6.3f}  {np.sin(theta):>6.3f}  "
              f"{Q_peak:>8.4f}  {K_max:>10.4g}  {ratio:>10.6f}")

    # And also for asymmetric S kernel
    print()
    print("  Same with S_kernel = diag(1, 0) (probe slow channel only):")
    print(f"  {'theta':>8}  {'Q_peak':>8}  {'HWHM-/Q*':>10}")

    def k_b_S_slow(Q, theta):
        """K_b with S_kernel projecting onto |c_1> slow channel only."""
        J = Q * gamma0
        L = np.array([
            [-2.0 * gamma0, 1j * J * g_eff],
            [1j * J * g_eff, -6.0 * gamma0],
        ], dtype=complex)
        Vb = np.array([
            [0.0, 1j * g_eff],
            [1j * g_eff, 0.0],
        ], dtype=complex)
        evals, R = np.linalg.eig(L)
        Rinv = np.linalg.inv(R)
        rho0 = np.array([np.cos(theta), np.sin(theta)], dtype=complex)
        expLam = np.exp(evals * t_peak)
        rho_t = R @ np.diag(expLam) @ (Rinv @ rho0)
        n = 2
        I_mat = np.zeros((n, n), dtype=complex)
        for j in range(n):
            for k in range(n):
                diff = evals[k] - evals[j]
                if abs(diff) < 1e-10:
                    I_mat[j, k] = t_peak * expLam[j]
                else:
                    I_mat[j, k] = (expLam[k] - expLam[j]) / diff
        X = Rinv @ Vb @ R
        c0 = Rinv @ rho0
        fbC0 = np.zeros(n, dtype=complex)
        for r in range(n):
            for c in range(n):
                fbC0[r] += X[r, c] * I_mat[r, c] * c0[c]
        drho = R @ fbC0
        S = np.diag([1, 0]).astype(complex)
        return 2.0 * np.vdot(rho_t, S @ drho).real

    for theta_deg in [0, 15, 30, 45, 60, 75, 90]:
        theta = np.radians(theta_deg)
        K_curve = np.array([k_b_S_slow(Q, theta) for Q in Q_grid])
        Q_peak, K_max, hwhm_left, ratio = compute_hwhm_ratio_from_curve(K_curve, Q_grid)
        print(f"  {theta_deg:>8}  {Q_peak:>8.4f}  {ratio:>10.6f}")


def main():
    np.set_printoptions(precision=6, suppress=True)
    explore_two_level_duhamel_K_b()
    explore_four_mode_cross_block()
    search_universal_x_peak()
    empirical_x_peak_check()
    explore_probe_mixing()
    # Sympy derivation succeeded but is messy; commented out for cleanup
    # try:
    #     derive_analytic_2level_K_b()
    # except Exception as e:
    #     print(f"  sympy derivation failed: {e}")


if __name__ == "__main__":
    main()
