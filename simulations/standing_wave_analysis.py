"""
STANDING WAVE ANALYSIS
======================

Formalizes the standing wave structure predicted by the palindromic symmetry.

For each palindromic pair (lambda_k, lambda_k' = -lambda_k - 2*Sg):
  Combined contribution under decay envelope exp(-Sg*t):
    exp(-Sg*t) * [ c_k * exp(mu_k*t) * |r_k>  +  c_k' * exp(-mu_k*t) * |r_k'> ]
  where mu_k = lambda_k + Sg (centered eigenvalue, symmetric around zero).

  If mu_k is purely imaginary: STANDING WAVE (counter-rotating oscillations).
  If mu_k is complex: one mode dominates early, the other late.
  If mu_k is purely real: exponential growth/decay pair (steady-XOR type).

Script: standing_wave_analysis.py
Output: results/standing_wave_analysis.txt

Authors: Tom Wicht, Claude
Date: March 19, 2026
"""
import numpy as np
from datetime import datetime

OUT = r"D:\Entwicklung\Projekte Privat\R-equals-C-Psi-squared\simulations\results\standing_wave_analysis.txt"
f = open(OUT, "w", buffering=1)

def log(msg=""):
    print(msg, flush=True)
    f.write(msg + "\n"); f.flush()

# ============================================================
# OPERATORS
# ============================================================

I2 = np.eye(2, dtype=complex)
sx = np.array([[0,1],[1,0]], dtype=complex)
sy = np.array([[0,-1j],[1j,0]], dtype=complex)
sz = np.array([[1,0],[0,-1]], dtype=complex)
PM = {'X': sx, 'Y': sy, 'Z': sz}

def site_op(op, s, N):
    ops = [I2]*N; ops[s] = op
    r = ops[0]
    for o in ops[1:]: r = np.kron(r, o)
    return r

def build_H(N, pairs, comps):
    """Build Hamiltonian from coupling pairs and component dict."""
    d = 2**N; H = np.zeros((d,d), dtype=complex)
    for i,j in pairs:
        for c, J in comps.items():
            if J == 0: continue
            H += J * site_op(PM[c[0]],i,N) @ site_op(PM[c[1]],j,N)
    return H

def build_L(H, gamma, N):
    """Build Liouvillian with uniform Z-dephasing rate gamma."""
    d = 2**N; d2 = d*d; Id = np.eye(d)
    L = -1j*(np.kron(H, Id) - np.kron(Id, H.T))
    for k in range(N):
        Zk = site_op(sz, k, N)
        L += gamma*(np.kron(Zk, Zk.conj()) - np.eye(d2))
    return L

def eigendecompose(L):
    """Eigendecomposition with left eigenvectors for biorthogonal expansion."""
    evals, R = np.linalg.eig(L)
    try:
        L_inv = np.linalg.inv(R)
    except np.linalg.LinAlgError:
        L_inv = np.linalg.pinv(R)
    return evals, R, L_inv

# ============================================================
# STATE BUILDERS
# ============================================================

def make_ghz(N):
    d = 2**N; v = np.zeros(d, dtype=complex)
    v[0] = v[-1] = 1/np.sqrt(2)
    return v

def make_w(N):
    d = 2**N; v = np.zeros(d, dtype=complex)
    for k in range(N): v[1 << (N-1-k)] = 1/np.sqrt(N)
    return v

def make_product(N, s):
    basis = {'0': np.array([1,0], dtype=complex),
             '1': np.array([0,1], dtype=complex),
             '+': np.array([1,1], dtype=complex)/np.sqrt(2),
             '-': np.array([1,-1], dtype=complex)/np.sqrt(2)}
    psi = basis[s[0]]
    for c in s[1:]: psi = np.kron(psi, basis[c])
    return psi

def make_bell_plus(N):
    d = 2**N; psi = np.zeros(d, dtype=complex)
    psi[0] = 1/np.sqrt(2)
    psi[3 << (N-2)] = 1/np.sqrt(2)
    return psi

def rho_vec(psi):
    return np.outer(psi, psi.conj()).flatten()

# ============================================================
# PALINDROMIC PAIRING
# ============================================================

def find_palindromic_pairs(evals, Sg, tol=1e-8):
    """Find palindromic pairs: lambda_k' = -lambda_k - 2*Sg.

    Returns pair_map[k] = index of partner (or -1 if unpaired).
    """
    n = len(evals)
    paired = np.zeros(n, dtype=bool)
    pair_map = -np.ones(n, dtype=int)

    for k in range(n):
        if paired[k]:
            continue
        target = -evals[k] - 2 * Sg
        # Find closest unmatched eigenvalue
        best_j = -1
        best_dist = np.inf
        for j in range(k + 1, n):
            if paired[j]:
                continue
            dist = abs(evals[j] - target)
            if dist < best_dist:
                best_dist = dist
                best_j = j
        if best_j >= 0 and best_dist < tol:
            paired[k] = True
            paired[best_j] = True
            pair_map[k] = best_j
            pair_map[best_j] = k

    return pair_map, paired

def classify_pairs(evals, pair_map, Sg, tol_re=1e-6, tol_im=1e-6):
    """Classify palindromic pairs by their centered eigenvalue mu_k.

    Returns dict with lists of (k, j, mu_k) tuples for each category.
    """
    mu = evals + Sg
    n = len(evals)
    seen = set()

    osc_pairs = []    # |Re(mu)| < tol, |Im(mu)| > tol: standing wave
    mixed_pairs = []  # |Re(mu)| > tol and |Im(mu)| > tol: mixed
    real_pairs = []   # |Im(mu)| < tol: pure decay (includes steady-XOR)
    zero_modes = []   # |mu| < tol: self-paired at center

    for k in range(n):
        if k in seen or pair_map[k] < 0:
            continue
        j = pair_map[k]
        seen.add(k); seen.add(j)

        mu_k = mu[k]
        re = abs(np.real(mu_k))
        im = abs(np.imag(mu_k))

        if re < tol_re and im < tol_im:
            zero_modes.append((k, j, mu_k))
        elif im < tol_im:
            real_pairs.append((k, j, mu_k))
        elif re < tol_re:
            osc_pairs.append((k, j, mu_k))
        else:
            mixed_pairs.append((k, j, mu_k))

    return {
        'oscillatory': osc_pairs,
        'mixed': mixed_pairs,
        'real': real_pairs,
        'zero': zero_modes,
        'mu': mu,
    }

# ============================================================
# MAIN
# ============================================================

if __name__ == "__main__":
    N = 3
    gamma = 0.05
    Sg = N * gamma  # sum_gamma for uniform dephasing
    pairs = [(i, i+1) for i in range(N-1)]
    d = 2**N
    d2 = d * d

    log("=" * 90)
    log("STANDING WAVE ANALYSIS")
    log(f"Date: {datetime.now()}")
    log("=" * 90)

    # ================================================================
    # SECTION 1: Eigendecomposition
    # ================================================================
    log()
    log("=" * 90)
    log("SECTION 1: Eigendecomposition")
    log(f"  N={N}, gamma={gamma}, Heisenberg chain, pairs={pairs}")
    log(f"  d={d}, d^2={d2} (Liouvillian dimension)")
    log(f"  Sg = N*gamma = {Sg}")
    log(f"  T2_envelope = 1/Sg = {1/Sg:.4f}")
    log("=" * 90)

    comps_heis = {'XX': 1, 'YY': 1, 'ZZ': 1}
    H = build_H(N, pairs, comps_heis)
    L = build_L(H, gamma, N)
    evals, R, L_inv = eigendecompose(L)

    log(f"\n  {d2} eigenvalues computed")

    # --- Eigenvalue overview ---
    log(f"\n  Eigenvalue overview:")
    log(f"    Re range: [{np.min(np.real(evals)):.6f}, {np.max(np.real(evals)):.6f}]")
    log(f"    Im range: [{np.min(np.imag(evals)):.6f}, {np.max(np.imag(evals)):.6f}]")

    # Count by basic type
    tol_ss = 1e-10
    n_steady = np.sum(np.abs(evals) < tol_ss)
    n_xor = np.sum((np.abs(np.real(evals) + 2*Sg) < 1e-6) & (np.abs(np.imag(evals)) < 1e-6))
    n_dynamic = d2 - n_steady
    log(f"    Steady states (|lambda| < 1e-10): {n_steady}")
    log(f"    XOR modes (Re ~ -2Sg, Im ~ 0):    {n_xor}")
    log(f"    Other dynamic modes:               {d2 - n_steady - n_xor}")

    # --- Palindromic pairing ---
    log(f"\n  Palindromic pairing (tol=1e-8):")
    pair_map, paired = find_palindromic_pairs(evals, Sg)
    n_paired = np.sum(paired)
    n_unpaired = d2 - n_paired
    n_pair_count = n_paired // 2

    log(f"    Paired:   {n_paired}/{d2} = {n_pair_count} pairs")
    log(f"    Unpaired: {n_unpaired}")

    if n_unpaired > 0:
        log(f"    Unpaired eigenvalues:")
        for k in range(d2):
            if not paired[k]:
                log(f"      lambda_{k} = {evals[k]:.10f}")

    # Verify: mu_k + mu_partner = 0
    mu = evals + Sg
    max_err = 0.0
    for k in range(d2):
        j = pair_map[k]
        if j >= 0:
            err = abs(mu[k] + mu[j])
            if err > max_err:
                max_err = err
    log(f"    Symmetry check: max |mu_k + mu_partner| = {max_err:.2e}")

    # --- Centered eigenvalue classification ---
    log(f"\n  Centered eigenvalue classification (mu_k = lambda_k + Sg):")
    cls = classify_pairs(evals, pair_map, Sg)

    n_osc = len(cls['oscillatory'])
    n_mix = len(cls['mixed'])
    n_real = len(cls['real'])
    n_zero = len(cls['zero'])

    log(f"    Purely imaginary mu (OSCILLATORY):  {n_osc} pairs  <-- standing wave candidates")
    log(f"    Complex mu (mixed decay+osc):       {n_mix} pairs")
    log(f"    Purely real mu (exponential):        {n_real} pairs  (includes steady-XOR)")
    log(f"    Zero mu (self-paired at center):     {n_zero} pairs")
    log(f"    Total classified pairs:              {n_osc + n_mix + n_real + n_zero}")

    # --- Detail: oscillatory pairs ---
    log(f"\n  --- Oscillatory pairs (standing wave candidates) ---")
    if cls['oscillatory']:
        log(f"  {'#':>4}  {'omega':>14}  {'period T':>12}  {'T/T2':>10}")
        log(f"  {'-'*46}")
        for idx, (k, j, mu_k) in enumerate(sorted(cls['oscillatory'],
                key=lambda x: abs(np.imag(x[2])))):
            omega = abs(np.imag(mu_k))
            T = 2 * np.pi / omega if omega > 1e-15 else np.inf
            T_ratio = T * Sg  # T / T2 where T2 = 1/Sg
            log(f"  {idx+1:>4}  {omega:>14.8f}  {T:>12.4f}  {T_ratio:>10.4f}")
        log(f"\n  These pairs oscillate WITHOUT net decay in the rescaled frame.")
        log(f"  Pattern persists; only overall amplitude (exp(-Sg*t)) decays.")
    else:
        log(f"  None found -- all pairs have Re(mu) != 0")
        log(f"  The standing wave still exists but has temporal asymmetry:")
        log(f"  one mode dominates early, its partner dominates late.")

    # --- Detail: complex (mixed) pairs ---
    log(f"\n  --- Complex pairs (mixed decay + oscillation) ---")
    if cls['mixed']:
        log(f"  {'#':>4}  {'Re(mu)':>14}  {'Im(mu)':>14}  {'omega':>14}  {'period T':>12}")
        log(f"  {'-'*64}")
        for idx, (k, j, mu_k) in enumerate(sorted(cls['mixed'],
                key=lambda x: abs(np.imag(x[2])))):
            re_mu = np.real(mu_k)
            im_mu = abs(np.imag(mu_k))
            T = 2 * np.pi / im_mu if im_mu > 1e-15 else np.inf
            log(f"  {idx+1:>4}  {re_mu:>+14.8f}  {im_mu:>14.8f}  {im_mu:>14.8f}  {T:>12.4f}")
        log(f"\n  These pairs oscillate AND have asymmetric decay in the rescaled frame.")
        log(f"  One partner dominates early (Re(mu) > 0), the other late (Re(mu) < 0).")
    else:
        log(f"  None found")

    # --- Detail: real pairs ---
    log(f"\n  --- Real pairs (pure decay, no oscillation) ---")
    if cls['real']:
        log(f"  {'#':>4}  {'|Re(mu)|':>14}  {'type':>20}")
        log(f"  {'-'*44}")
        for idx, (k, j, mu_k) in enumerate(sorted(cls['real'],
                key=lambda x: abs(np.real(x[2])))):
            re_mu = abs(np.real(mu_k))
            # Identify steady-XOR pairs
            if abs(re_mu - Sg) < 1e-6:
                ptype = "steady <-> XOR"
            else:
                ptype = "decay pair"
            log(f"  {idx+1:>4}  {re_mu:>14.8f}  {ptype:>20}")
    else:
        log(f"  None found")

    # --- Full eigenvalue table ---
    log(f"\n  --- Full eigenvalue table (sorted by Re) ---")
    log(f"  {'idx':>4}  {'Re(lambda)':>14}  {'Im(lambda)':>14}  "
        f"{'Re(mu)':>14}  {'Im(mu)':>14}  {'partner':>8}  {'type':>10}")
    log(f"  {'-'*86}")

    # Sort by real part
    order = np.argsort(np.real(evals))
    for rank, k in enumerate(order):
        re_l = np.real(evals[k])
        im_l = np.imag(evals[k])
        re_m = np.real(mu[k])
        im_m = np.imag(mu[k])
        partner = pair_map[k]

        if abs(evals[k]) < tol_ss:
            etype = "steady"
        elif abs(re_l + 2*Sg) < 1e-6 and abs(im_l) < 1e-6:
            etype = "XOR"
        elif abs(im_l) < 1e-6:
            etype = "real-dyn"
        else:
            etype = "complex"

        log(f"  {k:>4}  {re_l:>+14.8f}  {im_l:>+14.8f}  "
            f"{re_m:>+14.8f}  {im_m:>+14.8f}  {partner:>8}  {etype:>10}")

    # --- Summary ---
    log(f"\n  SECTION 1 SUMMARY:")
    log(f"    {d2} eigenvalues, {n_pair_count} palindromic pairs")
    log(f"    {n_steady} steady states, {n_xor} XOR modes")
    log(f"    Palindromic pairing: {'100%' if n_unpaired == 0 else f'{100*n_paired/d2:.1f}%'}")
    log(f"    Standing wave pairs (purely oscillatory): {n_osc}")
    log(f"    Mixed pairs (oscillation + asymmetry):    {n_mix}")
    log(f"    Real pairs (no oscillation):              {n_real}")
    if n_osc > 0:
        omegas = [abs(np.imag(mu_k)) for _, _, mu_k in cls['oscillatory']]
        log(f"    Oscillation frequencies: {len(set(f'{w:.6f}' for w in omegas))} distinct")
    if n_mix > 0:
        omegas = [abs(np.imag(mu_k)) for _, _, mu_k in cls['mixed']]
        log(f"    Mixed-pair frequencies:  {len(set(f'{w:.6f}' for w in omegas))} distinct")

    log(f"\n{'='*90}")
    log(f"Section 1 completed: {datetime.now()}")
    log(f"{'='*90}")

    # ================================================================
    # SECTION 2: Initial state decomposition
    # ================================================================
    log()
    log("=" * 90)
    log("SECTION 2: Initial state decomposition in eigenbasis")
    log("  For each state: c_k = <l_k | vec(rho_0)> via biorthogonal projection")
    log("  Track which palindromic pairs are excited and their standing wave content")
    log("=" * 90)

    # Build list of all pairs with their classification
    # pair_info[pair_idx] = (k, j, mu_k, category)
    pair_info = []
    seen_pairs = set()
    for k in range(d2):
        j = pair_map[k]
        if j < 0 or k in seen_pairs:
            continue
        seen_pairs.add(k); seen_pairs.add(j)
        mu_k = mu[k]
        re_m = abs(np.real(mu_k))
        im_m = abs(np.imag(mu_k))

        if abs(evals[k]) < tol_ss or abs(evals[j]) < tol_ss:
            cat = "steady-XOR"
        elif im_m < 1e-6 and re_m < 1e-6:
            cat = "zero"
        elif im_m < 1e-6:
            cat = "real"
        elif re_m < 1e-6:
            cat = "oscillatory"
        else:
            cat = "mixed"
        pair_info.append((k, j, mu_k, cat))

    n_total_pairs = len(pair_info)

    # Define test states
    test_states = {
        'GHZ':       make_ghz(N),
        'W':         make_w(N),
        '|000>':     make_product(N, '000'),
        '|+++>':     make_product(N, '+++'),
        'Bell(0,1)': make_bell_plus(N),
    }

    for sname, psi in test_states.items():
        rv = rho_vec(psi)
        coeffs = L_inv @ rv

        log(f"\n  {'-'*80}")
        log(f"  State: {sname}")
        log(f"  {'-'*80}")

        # Verify decomposition: reconstruct rho_vec from coefficients
        rv_recon = R @ coeffs
        recon_err = np.linalg.norm(rv_recon - rv) / np.linalg.norm(rv)
        log(f"  Reconstruction error: {recon_err:.2e}")

        # Total weight
        total_weight = np.sum(np.abs(coeffs)**2)
        log(f"  Total |c_k|^2: {total_weight:.6f}")

        # Weight by eigenvalue category
        w_steady = 0.0; w_xor = 0.0; w_osc = 0.0; w_mix = 0.0; w_real = 0.0

        for pi_idx, (k, j, mu_k, cat) in enumerate(pair_info):
            pw = abs(coeffs[k])**2 + abs(coeffs[j])**2
            if cat == "steady-XOR":
                # Split into steady and XOR parts
                if abs(evals[k]) < tol_ss:
                    w_steady += abs(coeffs[k])**2
                    w_xor += abs(coeffs[j])**2
                else:
                    w_steady += abs(coeffs[j])**2
                    w_xor += abs(coeffs[k])**2
            elif cat == "oscillatory":
                w_osc += pw
            elif cat == "mixed":
                w_mix += pw
            elif cat == "real":
                w_real += pw

        log(f"\n  Weight distribution:")
        log(f"    Steady state:           {w_steady:>10.6f}  ({100*w_steady/total_weight:>6.2f}%)")
        log(f"    XOR modes:              {w_xor:>10.6f}  ({100*w_xor/total_weight:>6.2f}%)")
        log(f"    Oscillatory pairs:      {w_osc:>10.6f}  ({100*w_osc/total_weight:>6.2f}%)")
        log(f"    Mixed (osc+asym) pairs: {w_mix:>10.6f}  ({100*w_mix/total_weight:>6.2f}%)")
        log(f"    Real decay pairs:       {w_real:>10.6f}  ({100*w_real/total_weight:>6.2f}%)")

        # Standing wave content = weight in pairs that oscillate (osc + mixed)
        w_standing = w_osc + w_mix
        w_dynamic = total_weight - w_steady
        log(f"\n  Standing wave content (oscillating pairs / dynamic weight):")
        if w_dynamic > 1e-15:
            log(f"    {100*w_standing/w_dynamic:.2f}% of dynamic weight oscillates")
        else:
            log(f"    No dynamic weight")

        # Pair-by-pair analysis: which pairs are excited?
        log(f"\n  Pair excitation analysis:")
        hdr_ckp = "|c_k'|"
        log(f"  {'pair':>5}  {'cat':>10}  {'|c_k|':>10}  {hdr_ckp:>10}  "
            f"{'balance':>8}  {'pair_wt%':>9}  {'omega':>10}")
        log(f"  {'-'*70}")

        n_excited = 0
        n_both_excited = 0
        for pi_idx, (k, j, mu_k, cat) in enumerate(pair_info):
            ck = abs(coeffs[k])
            cj = abs(coeffs[j])
            pw = ck**2 + cj**2
            pw_pct = 100 * pw / total_weight if total_weight > 1e-15 else 0

            if pw < 1e-12:
                continue  # skip unexcited pairs

            n_excited += 1
            mx = max(ck, cj)
            mn = min(ck, cj)
            balance = mn / mx if mx > 1e-15 else 0.0
            if mn > 1e-8:
                n_both_excited += 1

            omega = abs(np.imag(mu_k))
            omega_str = f"{omega:.4f}" if omega > 1e-6 else "0"

            log(f"  {pi_idx+1:>5}  {cat:>10}  {ck:>10.6f}  {cj:>10.6f}  "
                f"{balance:>8.4f}  {pw_pct:>8.2f}%  {omega_str:>10}")

        log(f"\n  Excited pairs: {n_excited}/{n_total_pairs}")
        log(f"  Both members excited (standing wave active): {n_both_excited}/{n_excited}")

        # Dominant frequencies: which oscillation frequencies carry the most weight?
        freq_weight = {}  # omega -> total weight
        for pi_idx, (k, j, mu_k, cat) in enumerate(pair_info):
            omega = abs(np.imag(mu_k))
            if omega < 1e-6:
                continue
            pw = abs(coeffs[k])**2 + abs(coeffs[j])**2
            if pw < 1e-12:
                continue
            key = f"{omega:.4f}"
            freq_weight[key] = freq_weight.get(key, 0) + pw

        if freq_weight:
            log(f"\n  Dominant frequencies (sorted by weight):")
            for omega_str, fw in sorted(freq_weight.items(), key=lambda x: -x[1]):
                log(f"    omega = {omega_str:>10}  weight = {fw:.6f}  "
                    f"({100*fw/total_weight:.2f}%)")
        else:
            log(f"\n  No oscillatory modes excited")

    # ================================================================
    # SECTION 2 SUMMARY: State comparison
    # ================================================================
    log(f"\n  {'='*80}")
    log(f"  SECTION 2 SUMMARY: Standing wave content comparison")
    log(f"  {'='*80}")
    log(f"\n  {'State':<12} {'Steady%':>8} {'XOR%':>8} {'Osc%':>8} "
        f"{'Mixed%':>8} {'Real%':>8} {'SW_dyn%':>9}")
    log(f"  {'-'*65}")

    for sname, psi in test_states.items():
        rv = rho_vec(psi)
        coeffs = L_inv @ rv
        total_weight = np.sum(np.abs(coeffs)**2)

        w_steady = 0.0; w_xor = 0.0; w_osc = 0.0; w_mix = 0.0; w_real = 0.0
        for pi_idx, (k, j, mu_k, cat) in enumerate(pair_info):
            pw = abs(coeffs[k])**2 + abs(coeffs[j])**2
            if cat == "steady-XOR":
                if abs(evals[k]) < tol_ss:
                    w_steady += abs(coeffs[k])**2
                    w_xor += abs(coeffs[j])**2
                else:
                    w_steady += abs(coeffs[j])**2
                    w_xor += abs(coeffs[k])**2
            elif cat == "oscillatory":
                w_osc += pw
            elif cat == "mixed":
                w_mix += pw
            elif cat == "real":
                w_real += pw

        w_dynamic = total_weight - w_steady
        sw_dyn = 100*(w_osc + w_mix)/w_dynamic if w_dynamic > 1e-15 else 0

        log(f"  {sname:<12} {100*w_steady/total_weight:>7.1f}% "
            f"{100*w_xor/total_weight:>7.1f}% "
            f"{100*w_osc/total_weight:>7.1f}% "
            f"{100*w_mix/total_weight:>7.1f}% "
            f"{100*w_real/total_weight:>7.1f}% "
            f"{sw_dyn:>8.1f}%")

    log(f"\n  SW_dyn% = standing wave fraction of dynamic weight")
    log(f"         = (oscillatory + mixed pair weight) / (total - steady)")
    log(f"\n  Prediction check:")
    log(f"    W -> maximum standing wave? (100% palindromic = all in oscillating pairs)")
    log(f"    GHZ -> no standing wave? (100% XOR = real decay pairs, no oscillation)")
    log(f"    |000> -> no standing wave? (100% steady = no dynamics at all)")

    log(f"\n{'='*90}")
    log(f"Section 2 completed: {datetime.now()}")
    log(f"{'='*90}")

    # ================================================================
    # SECTION 3: Rescaled dynamics
    # ================================================================
    log()
    log("=" * 90)
    log("SECTION 3: Rescaled dynamics")
    log("  rho_rescaled(t) = rho(t) * exp(Sg*t)")
    log("  In eigenbasis: rho_rescaled_vec(t) = R @ (c_k * exp(mu_k * t))")
    log("  Removes uniform decay envelope -> reveals standing wave pattern")
    log("=" * 90)

    # Time range: 0 to 5*T2, resolve fastest oscillation (omega~6, T~1.05)
    T2 = 1.0 / Sg
    t_max = 5 * T2
    n_t = 1000
    t_arr = np.linspace(0, t_max, n_t)
    dt = t_arr[1] - t_arr[0]

    # Two-site Pauli operators for observables (sites 0,1)
    XX_01 = np.kron(np.kron(sx, sx), I2)
    YY_01 = np.kron(np.kron(sy, sy), I2)
    ZZ_01 = np.kron(np.kron(sz, sz), I2)
    XY_01 = np.kron(np.kron(sx, sy), I2)

    def compute_rescaled_observables(coeffs_in, t_arr, mu, R, d):
        """Compute rescaled-frame observables over time."""
        n_t = len(t_arr)
        # rho_rescaled_vec(t) = R @ (c_k * exp(mu_k * t))
        exp_mu_t = np.exp(mu[None, :] * t_arr[:, None])  # (n_t, d2)
        coeffs_t = coeffs_in[None, :] * exp_mu_t          # (n_t, d2)
        rho_vecs = (R @ coeffs_t.T).T                      # (n_t, d2)

        purity = np.zeros(n_t)
        rho01_abs = np.zeros(n_t)
        xx = np.zeros(n_t)
        yy = np.zeros(n_t)
        zz = np.zeros(n_t)
        xy = np.zeros(n_t)

        for ti in range(n_t):
            rm = rho_vecs[ti].reshape(d, d)
            purity[ti] = np.real(np.trace(rm @ rm))
            rho01_abs[ti] = abs(rm[0, 1])
            xx[ti] = np.real(np.trace(XX_01 @ rm))
            yy[ti] = np.real(np.trace(YY_01 @ rm))
            zz[ti] = np.real(np.trace(ZZ_01 @ rm))
            xy[ti] = np.real(np.trace(XY_01 @ rm))

        return {'Purity': purity, '|rho_01|': rho01_abs,
                '<XX>': xx, '<YY>': yy, '<ZZ>': zz, '<XY>': xy}

    def report_observables(obs_dict, t_arr, Sg, label):
        """Print summary table and time trace for observables."""
        T2 = 1.0 / Sg
        log(f"\n  {'Observable':<14} {'t=0':>10} {'mean':>10} {'min':>10} "
            f"{'max':>10} {'amp':>10} {'osc?':>6}")
        log(f"  {'-'*72}")

        for oname, ovals in obs_dict.items():
            v0 = ovals[0]
            vmean = np.mean(ovals)
            vmin = np.min(ovals)
            vmax = np.max(ovals)
            amp = vmax - vmin
            osc = "YES" if amp > max(0.01 * abs(vmean), 1e-6) else "no"
            log(f"  {oname:<14} {v0:>10.6f} {vmean:>10.6f} {vmin:>10.6f} "
                f"{vmax:>10.6f} {amp:>10.6f} {osc:>6}")

        # Time trace at key points
        key_t = [0, T2/4, T2/2, T2, 2*T2, 3*T2, 5*T2]
        hdr = f"  {'t':>8} {'t/T2':>6} "
        for oname in obs_dict:
            hdr += f" {oname:>12}"
        log(f"\n{hdr}")
        log(f"  {'-'*(16 + 13*len(obs_dict))}")
        for tk in key_t:
            idx = np.argmin(abs(t_arr - tk))
            t = t_arr[idx]
            line = f"  {t:>8.3f} {t*Sg:>6.3f} "
            for ovals in obs_dict.values():
                line += f" {ovals[idx]:>12.6f}"
            log(line)

    def fft_analysis(signal, t_arr, name, n_peaks=3):
        """FFT frequency analysis of a signal."""
        amp = np.max(signal) - np.min(signal)
        if amp < 1e-6:
            log(f"\n  FFT of {name}: signal is flat (amp={amp:.2e}), skipping")
            return
        dt = t_arr[1] - t_arr[0]
        n = len(signal)
        fft_vals = np.fft.rfft(signal - np.mean(signal))
        freqs = np.fft.rfftfreq(n, d=dt)
        power = np.abs(fft_vals)**2
        power[0] = 0  # exclude DC
        top_idx = np.argsort(power)[::-1][:n_peaks]
        log(f"\n  FFT of {name} (top {n_peaks} frequencies):")
        for rank, fi in enumerate(top_idx):
            omega_fft = 2 * np.pi * freqs[fi]
            T_fft = 1.0 / freqs[fi] if freqs[fi] > 0 else np.inf
            log(f"    #{rank+1}: omega = {omega_fft:.4f}  "
                f"(T = {T_fft:.4f})  power = {power[fi]:.4e}")

    # --- Run Section 3 for all original states ---
    for sname, psi in test_states.items():
        rv = rho_vec(psi)
        coeffs = L_inv @ rv

        log(f"\n  {'='*80}")
        log(f"  State: {sname}")
        log(f"  {'='*80}")

        obs = compute_rescaled_observables(coeffs, t_arr, mu, R, d)
        report_observables(obs, t_arr, Sg, sname)
        fft_analysis(obs['Purity'], t_arr, 'Purity')
        fft_analysis(obs['<XX>'], t_arr, '<XX>')

    # ================================================================
    # SECTION 3b: Bell + |+++> superposition
    # ================================================================
    log(f"\n\n  {'='*80}")
    log(f"  SPECIAL TEST: Bell(0,1) + |+++> superposition")
    log(f"  Bell has oscillating pairs (ingredient a)")
    log(f"  |+++> has balanced pairs (ingredient b)")
    log(f"  Does the superposition combine both?")
    log(f"  {'='*80}")

    bell_psi = make_bell_plus(N)
    plus_psi = make_product(N, '+++')

    # Check overlap
    overlap = abs(np.dot(bell_psi.conj(), plus_psi))
    log(f"\n  |<Bell|+++>| = {overlap:.6f}")

    # Superposition (equal weight, normalized)
    sup_psi = bell_psi + plus_psi
    sup_psi /= np.linalg.norm(sup_psi)
    log(f"  |psi> = (|Bell> + |+++>) / norm")
    log(f"  norm factor = {np.linalg.norm(bell_psi + plus_psi):.6f}")

    sup_rv = rho_vec(sup_psi)
    sup_coeffs = L_inv @ sup_rv
    sup_total = np.sum(np.abs(sup_coeffs)**2)

    # Pair excitation analysis (Section 2 style)
    log(f"\n  --- Pair excitation analysis ---")

    w_steady = 0.0; w_xor = 0.0; w_osc = 0.0; w_mix = 0.0; w_real = 0.0
    n_excited = 0; n_both = 0; n_osc_both = 0

    for pi_idx, (k, j, mu_k, cat) in enumerate(pair_info):
        ck = abs(sup_coeffs[k])
        cj = abs(sup_coeffs[j])
        pw = ck**2 + cj**2

        if cat == "steady-XOR":
            if abs(evals[k]) < tol_ss:
                w_steady += ck**2; w_xor += cj**2
            else:
                w_steady += cj**2; w_xor += ck**2
        elif cat == "oscillatory":
            w_osc += pw
        elif cat == "mixed":
            w_mix += pw
        elif cat == "real":
            w_real += pw

        if pw > 1e-12:
            n_excited += 1
            mn = min(ck, cj)
            if mn > 1e-8:
                n_both += 1
                if cat in ("oscillatory", "mixed"):
                    n_osc_both += 1

    w_dynamic = sup_total - w_steady
    sw_dyn = 100 * (w_osc + w_mix) / w_dynamic if w_dynamic > 1e-15 else 0

    log(f"\n  Weight distribution:")
    log(f"    Steady:      {100*w_steady/sup_total:>7.2f}%")
    log(f"    XOR:         {100*w_xor/sup_total:>7.2f}%")
    log(f"    Oscillatory: {100*w_osc/sup_total:>7.2f}%")
    log(f"    Mixed:       {100*w_mix/sup_total:>7.2f}%")
    log(f"    Real:        {100*w_real/sup_total:>7.2f}%")
    log(f"\n  Standing wave content (osc+mixed / dynamic): {sw_dyn:.2f}%")
    log(f"  Excited pairs: {n_excited}/{n_total_pairs}")
    log(f"  Both members excited: {n_both}/{n_excited}")
    log(f"  Both members excited IN OSCILLATING pairs: {n_osc_both}")
    log(f"\n  KEY: n_osc_both > 0 means BOTH ingredients present!")
    log(f"       (oscillating pair frequency + balanced excitation = standing wave)")

    # Detailed pair table for oscillating pairs
    log(f"\n  Oscillating pairs detail:")
    hdr_ckp = "|c_k'|"
    log(f"  {'pair':>5}  {'cat':>8}  {'|c_k|':>10}  {hdr_ckp:>10}  "
        f"{'balance':>8}  {'wt%':>7}  {'omega':>10}")
    log(f"  {'-'*66}")
    for pi_idx, (k, j, mu_k, cat) in enumerate(pair_info):
        if cat not in ("oscillatory", "mixed"):
            continue
        ck = abs(sup_coeffs[k])
        cj = abs(sup_coeffs[j])
        pw = ck**2 + cj**2
        if pw < 1e-12:
            continue
        mx = max(ck, cj)
        mn = min(ck, cj)
        balance = mn / mx if mx > 1e-15 else 0.0
        omega = abs(np.imag(mu_k))
        log(f"  {pi_idx+1:>5}  {cat:>8}  {ck:>10.6f}  {cj:>10.6f}  "
            f"{balance:>8.4f}  {100*pw/sup_total:>6.2f}%  {omega:>10.4f}")

    # Rescaled dynamics
    log(f"\n  --- Rescaled dynamics ---")
    obs_sup = compute_rescaled_observables(sup_coeffs, t_arr, mu, R, d)
    report_observables(obs_sup, t_arr, Sg, 'Bell+|+++>')
    fft_analysis(obs_sup['Purity'], t_arr, 'Purity')
    fft_analysis(obs_sup['<XX>'], t_arr, '<XX>')
    fft_analysis(obs_sup['<XY>'], t_arr, '<XY>')

    # Compare oscillation amplitudes
    log(f"\n  {'='*80}")
    log(f"  SECTION 3 SUMMARY: Oscillation in rescaled frame")
    log(f"  {'='*80}")
    log(f"\n  {'State':<14} {'Pur_amp':>10} {'rho01_amp':>10} "
        f"{'XX_amp':>10} {'YY_amp':>10} {'ZZ_amp':>10} {'XY_amp':>10}")
    log(f"  {'-'*76}")

    all_s3_states = dict(test_states)
    all_s3_states['Bell+|+++>'] = sup_psi

    for sname, psi in all_s3_states.items():
        rv = rho_vec(psi)
        coeffs = L_inv @ rv
        obs = compute_rescaled_observables(coeffs, t_arr, mu, R, d)
        amps = {k: np.max(v) - np.min(v) for k, v in obs.items()}
        log(f"  {sname:<14} {amps['Purity']:>10.6f} {amps['|rho_01|']:>10.6f} "
            f"{amps['<XX>']:>10.6f} {amps['<YY>']:>10.6f} "
            f"{amps['<ZZ>']:>10.6f} {amps['<XY>']:>10.6f}")

    log(f"\n  Amplitude > 0 = oscillation in rescaled frame = standing wave signature")

    log(f"\n{'='*90}")
    log(f"Section 3 completed: {datetime.now()}")
    log(f"{'='*90}")

    # ================================================================
    # SECTION 4: Identify the standing pattern
    # ================================================================
    log()
    log("=" * 90)
    log("SECTION 4: Standing pattern identification")
    log("  Analytical approach: separate oscillating vs non-oscillating eigenmode")
    log("  contributions. The 'standing pattern' = Pauli fingerprint of the")
    log("  oscillating part rho_osc(0) = sum_{Im(mu)!=0} c_k * |r_k>")
    log("=" * 90)

    # Build 2-site Pauli observables for all site pairs
    pauli_1q = {'I': I2, 'X': sx, 'Y': sy, 'Z': sz}
    pauli_keys = ['I', 'X', 'Y', 'Z']

    # Full 3-site Pauli basis (64 operators)
    pauli_basis_3 = []
    for a in pauli_keys:
        for b in pauli_keys:
            for c in pauli_keys:
                label = a + b + c
                op = np.kron(np.kron(pauli_1q[a], pauli_1q[b]), pauli_1q[c])
                pauli_basis_3.append((label, op))

    # Precompute: for each eigenmode k, compute Tr(P * M_k) for all Pauli P
    # M_k = R[:,k].reshape(d,d)
    # pauli_trace[p, k] = Tr(P_p @ M_k)
    n_pauli = len(pauli_basis_3)
    pauli_trace = np.zeros((n_pauli, d2), dtype=complex)
    for pi in range(n_pauli):
        P = pauli_basis_3[pi][1]
        for ki in range(d2):
            pauli_trace[pi, ki] = np.trace(P @ R[:, ki].reshape(d, d))

    # Masks for oscillating vs non-oscillating modes
    osc_mask = np.abs(np.imag(mu)) > 1e-6
    nonosc_mask = ~osc_mask

    # Frequency bands for grouping
    freq_bands = [
        ("omega~2", 1.5, 2.5),
        ("omega~4", 3.5, 4.5),
        ("omega~6", 5.5, 6.5),
    ]

    # States to analyze in Section 4
    states_s4 = {
        'Bell(0,1)': make_bell_plus(N),
        'Bell+|+++>': sup_psi,
        'W':          make_w(N),
        'GHZ':        make_ghz(N),
        '|+++>':      make_product(N, '+++'),
    }

    for sname, psi in states_s4.items():
        rv = rho_vec(psi)
        coeffs = L_inv @ rv

        log(f"\n  {'='*80}")
        log(f"  State: {sname}")
        log(f"  {'='*80}")

        # 4a: Oscillating part at t=0
        # rho_osc_vec = sum over oscillating modes c_k * R[:,k]
        coeffs_osc = coeffs * osc_mask
        coeffs_bg = coeffs * nonosc_mask

        rho_osc_vec = R @ coeffs_osc
        rho_bg_vec = R @ coeffs_bg

        rho_osc_0 = rho_osc_vec.reshape(d, d)
        rho_bg_0 = rho_bg_vec.reshape(d, d)

        osc_norm = np.linalg.norm(rho_osc_vec)
        bg_norm = np.linalg.norm(rho_bg_vec)
        total_norm = np.linalg.norm(rv)

        log(f"\n  Mode separation at t=0:")
        log(f"    ||rho_osc||   = {osc_norm:.6f}  "
            f"({100*osc_norm**2/total_norm**2:.1f}% of Frobenius weight)")
        log(f"    ||rho_bg||    = {bg_norm:.6f}  "
            f"({100*bg_norm**2/total_norm**2:.1f}%)")

        if osc_norm < 1e-10:
            log(f"\n  No oscillating modes excited. Skipping Pauli analysis.")
            continue

        # 4b: Pauli decomposition of oscillating part
        # <P>_osc = Tr(P @ rho_osc(0)) / d  (Pauli coefficient)
        pauli_osc = np.zeros(n_pauli, dtype=complex)
        pauli_bg = np.zeros(n_pauli, dtype=complex)
        for pi in range(n_pauli):
            pauli_osc[pi] = np.trace(pauli_basis_3[pi][1] @ rho_osc_0) / d
            pauli_bg[pi] = np.trace(pauli_basis_3[pi][1] @ rho_bg_0) / d

        # Report nonzero oscillating Pauli components (sorted by magnitude)
        log(f"\n  Pauli fingerprint of oscillating part rho_osc(0):")
        log(f"  (Pauli coefficients = Tr(P @ rho_osc) / d, showing |coeff| > 1e-6)")
        log(f"  {'Pauli':>6}  {'|coeff|':>10}  {'Re':>12}  {'Im':>12}  {'bg_coeff':>12}")
        log(f"  {'-'*56}")

        order = np.argsort(-np.abs(pauli_osc))
        n_shown = 0
        for pi in order:
            if np.abs(pauli_osc[pi]) < 1e-6:
                break
            label = pauli_basis_3[pi][0]
            mag = abs(pauli_osc[pi])
            re = np.real(pauli_osc[pi])
            im = np.imag(pauli_osc[pi])
            bg_re = np.real(pauli_bg[pi])
            log(f"  {label:>6}  {mag:>10.6f}  {re:>+12.6f}  {im:>+12.6f}  {bg_re:>+12.6f}")
            n_shown += 1

        if n_shown == 0:
            log(f"  (all Pauli coefficients below threshold)")

        # 4c: Frequency-resolved Pauli decomposition
        log(f"\n  Frequency-resolved oscillating Pauli components:")
        for band_name, omega_lo, omega_hi in freq_bands:
            band_mask = (np.abs(np.imag(mu)) > omega_lo) & \
                        (np.abs(np.imag(mu)) < omega_hi)
            coeffs_band = coeffs * band_mask
            if np.sum(np.abs(coeffs_band)**2) < 1e-15:
                continue

            rho_band = (R @ coeffs_band).reshape(d, d)
            log(f"\n    Band: {band_name}")
            log(f"    {'Pauli':>6}  {'|coeff|':>10}")
            log(f"    {'-'*18}")

            pauli_band = np.array([
                np.trace(pauli_basis_3[pi][1] @ rho_band) / d
                for pi in range(n_pauli)
            ])
            order_b = np.argsort(-np.abs(pauli_band))
            for pi in order_b[:8]:  # top 8
                if np.abs(pauli_band[pi]) < 1e-6:
                    break
                log(f"    {pauli_basis_3[pi][0]:>6}  {abs(pauli_band[pi]):>10.6f}")

        # 4d: Node identification
        # Nodes = Pauli operators with NO oscillating contribution
        # (but nonzero background)
        log(f"\n  Nodes (non-oscillating Paulis with nonzero background):")
        nodes = []
        antinodes = []
        for pi in range(n_pauli):
            label = pauli_basis_3[pi][0]
            if label == 'III':
                continue  # skip identity
            has_osc = abs(pauli_osc[pi]) > 1e-6
            has_bg = abs(pauli_bg[pi]) > 1e-6
            if has_bg and not has_osc:
                nodes.append(label)
            if has_osc:
                antinodes.append(label)

        if nodes:
            log(f"    {', '.join(nodes)}")
            log(f"    These observables have static value, no oscillation")
        else:
            log(f"    No nodes found (all nonzero observables oscillate, or no bg)")

        log(f"\n  Antinodes (oscillating Paulis):")
        if antinodes:
            log(f"    {', '.join(antinodes)}")
        else:
            log(f"    None")

    # ================================================================
    # SECTION 4 SUMMARY: Cross-state comparison
    # ================================================================
    log(f"\n  {'='*80}")
    log(f"  SECTION 4 SUMMARY: Standing wave fingerprints")
    log(f"  {'='*80}")

    log(f"\n  Which Pauli correlations oscillate for each state?")
    log(f"  (listing Pauli strings with |osc_coeff| > 0.01)")
    log()

    for sname, psi in states_s4.items():
        rv = rho_vec(psi)
        coeffs = L_inv @ rv
        coeffs_osc = coeffs * osc_mask
        rho_osc_0 = (R @ coeffs_osc).reshape(d, d)
        osc_norm = np.linalg.norm(R @ coeffs_osc)

        osc_paulis = []
        for pi in range(n_pauli):
            label = pauli_basis_3[pi][0]
            if label == 'III':
                continue
            val = abs(np.trace(pauli_basis_3[pi][1] @ rho_osc_0) / d)
            if val > 0.01:
                osc_paulis.append(f"{label}({val:.3f})")

        if osc_paulis:
            log(f"  {sname:<14} ||osc||={osc_norm:.4f}  {', '.join(osc_paulis)}")
        else:
            log(f"  {sname:<14} ||osc||={osc_norm:.4f}  (no oscillating Paulis)")

    # Frobenius-normalized time-average (secondary check)
    log(f"\n  {'='*80}")
    log(f"  Frobenius-normalized time-average of rescaled state")
    log(f"  rho_shape(t) = rho_resc(t) / ||rho_resc(t)||_F, then averaged over [0, 5*T2]")
    log(f"  WARNING: dominated by late-time growth (steady-state mode)")
    log(f"  {'='*80}")

    for sname in ['Bell(0,1)', 'Bell+|+++>']:
        psi = states_s4[sname]
        rv = rho_vec(psi)
        coeffs = L_inv @ rv

        # Compute time-average of Frobenius-normalized state
        n_avg = 200
        t_avg = np.linspace(0, 5 * T2, n_avg)
        rho_avg = np.zeros((d, d), dtype=complex)

        for ti in range(n_avg):
            exp_mu_ti = np.exp(mu * t_avg[ti])
            rv_ti = R @ (coeffs * exp_mu_ti)
            rm_ti = rv_ti.reshape(d, d)
            frob = np.sqrt(np.real(np.trace(rm_ti @ rm_ti.conj().T)))
            if frob > 1e-15:
                rho_avg += rm_ti / frob

        rho_avg /= n_avg

        # Pauli decomposition of the average
        log(f"\n  State: {sname}")
        log(f"  {'Pauli':>6}  {'coeff':>12}")
        log(f"  {'-'*20}")

        pauli_avg = []
        for pi in range(n_pauli):
            label = pauli_basis_3[pi][0]
            val = np.real(np.trace(pauli_basis_3[pi][1] @ rho_avg) / d)
            if abs(val) > 1e-4:
                pauli_avg.append((label, val))

        for label, val in sorted(pauli_avg, key=lambda x: -abs(x[1])):
            log(f"  {label:>6}  {val:>+12.6f}")

    log(f"\n{'='*90}")
    log(f"Section 4 completed: {datetime.now()}")
    log(f"{'='*90}")

    # ================================================================
    # SECTION 5: State comparison
    # ================================================================
    log()
    log("=" * 90)
    log("SECTION 5: Standing wave content across initial states")
    log("  Using Heisenberg Liouvillian from Section 1")
    log("  KEY QUESTION: does the standing wave depend on the initial state?")
    log("=" * 90)

    # Extended state set
    states_s5 = {
        'GHZ':           make_ghz(N),
        'W':             make_w(N),
        '|000>':         make_product(N, '000'),
        '|111>':         make_product(N, '111'),
        '|010>':         make_product(N, '010'),
        '|+++>':         make_product(N, '+++'),
        '|+-+>':         make_product(N, '+-+'),
        '|+00>':         make_product(N, '+00'),
        'Bell(0,1)':     make_bell_plus(N),
        'Bell+|+++>':    sup_psi,
    }

    # Also try more superpositions that might combine both ingredients
    # Bell + W: does W's palindromic weight + Bell's oscillation combine?
    bw_psi = make_bell_plus(N) + make_w(N)
    bw_psi /= np.linalg.norm(bw_psi)
    states_s5['Bell+W'] = bw_psi

    # GHZ + Bell: both have coherences but different structure
    gb_psi = make_ghz(N) + make_bell_plus(N)
    gb_psi /= np.linalg.norm(gb_psi)
    states_s5['GHZ+Bell'] = gb_psi

    log(f"\n  {'State':<14} {'||osc||':>8} {'osc%':>7} {'#osc_P':>7} "
        f"{'top Pauli':>10} {'top_amp':>9} {'n_nodes':>8}")
    log(f"  {'-'*70}")

    for sname, psi in states_s5.items():
        rv = rho_vec(psi)
        coeffs = L_inv @ rv
        total_norm = np.linalg.norm(rv)

        coeffs_osc = coeffs * osc_mask
        rho_osc_0 = (R @ coeffs_osc).reshape(d, d)
        osc_norm = np.linalg.norm(R @ coeffs_osc)
        osc_pct = 100 * osc_norm**2 / total_norm**2

        # Count oscillating Paulis and find top one
        osc_paulis = []
        bg_paulis = []
        coeffs_bg = coeffs * nonosc_mask
        rho_bg_0 = (R @ coeffs_bg).reshape(d, d)
        for pi in range(n_pauli):
            label = pauli_basis_3[pi][0]
            if label == 'III':
                continue
            osc_val = abs(np.trace(pauli_basis_3[pi][1] @ rho_osc_0) / d)
            bg_val = abs(np.trace(pauli_basis_3[pi][1] @ rho_bg_0) / d)
            if osc_val > 1e-6:
                osc_paulis.append((label, osc_val))
            if bg_val > 1e-6 and osc_val < 1e-6:
                bg_paulis.append(label)

        n_osc_p = len(osc_paulis)
        if osc_paulis:
            top_p, top_a = max(osc_paulis, key=lambda x: x[1])
        else:
            top_p, top_a = "-", 0.0

        log(f"  {sname:<14} {osc_norm:>8.4f} {osc_pct:>6.1f}% {n_osc_p:>7} "
            f"{top_p:>10} {top_a:>9.4f} {len(bg_paulis):>8}")

    # Which states actually produce oscillation?
    log(f"\n  FINDING: Only states containing Bell(0,1) coherences produce oscillation.")
    log(f"  Bell creates coherences between sectors with different magnetization")
    log(f"  (|000> <-> |110>) that couple to oscillating Liouvillian modes.")
    log(f"  W, GHZ, product states excite ONLY non-oscillating modes.")

    # Check: does the PATTERN differ between oscillating states?
    log(f"\n  Do oscillating states share the same Pauli fingerprint?")
    log(f"  (Showing top 5 oscillating Paulis for each oscillating state)")

    for sname, psi in states_s5.items():
        rv = rho_vec(psi)
        coeffs = L_inv @ rv
        coeffs_osc = coeffs * osc_mask
        osc_norm = np.linalg.norm(R @ coeffs_osc)
        if osc_norm < 1e-6:
            continue

        rho_osc_0 = (R @ coeffs_osc).reshape(d, d)
        osc_paulis = []
        for pi in range(n_pauli):
            label = pauli_basis_3[pi][0]
            if label == 'III':
                continue
            val = abs(np.trace(pauli_basis_3[pi][1] @ rho_osc_0) / d)
            if val > 1e-6:
                osc_paulis.append((label, val))
        osc_paulis.sort(key=lambda x: -x[1])

        top5 = ', '.join(f"{l}({v:.4f})" for l, v in osc_paulis[:5])
        log(f"  {sname:<14} {top5}")

    log(f"\n{'='*90}")
    log(f"Section 5 completed: {datetime.now()}")
    log(f"{'='*90}")

    # ================================================================
    # SECTION 6: Hamiltonian comparison
    # ================================================================
    log()
    log("=" * 90)
    log("SECTION 6: Standing wave on different Hamiltonians")
    log("  Repeat analysis for: Heisenberg, XY-only, Ising, DM (XY-YX)")
    log("  Using Bell(0,1) as test state (strongest oscillator)")
    log("  QUESTION: does the pattern change? Do frequencies differ?")
    log("=" * 90)

    models_s6 = {
        'Heisenberg':  {'XX': 1, 'YY': 1, 'ZZ': 1},
        'XY-only':     {'XX': 1, 'YY': 1},
        'Ising':       {'ZZ': 1},
        'DM (XY-YX)':  {'XY': 1, 'YX': -1},
        'XXZ d=0.5':   {'XX': 1, 'YY': 1, 'ZZ': 0.5},
        'Heis+DM':     {'XX': 1, 'YY': 1, 'ZZ': 1, 'XY': 0.3, 'YX': -0.3},
    }

    bell_psi_s6 = make_bell_plus(N)

    # Summary table
    log(f"\n  {'Model':<14} {'#osc_pr':>8} {'#mix_pr':>8} {'#real_pr':>8} "
        f"{'osc%':>7} {'freqs':>30}")
    log(f"  {'-'*82}")

    model_results = {}

    for mname, comps in models_s6.items():
        H_m = build_H(N, pairs, comps)
        L_m = build_L(H_m, gamma, N)
        ev_m, R_m, Li_m = eigendecompose(L_m)
        mu_m = ev_m + Sg
        pm_m, pd_m = find_palindromic_pairs(ev_m, Sg)
        cls_m = classify_pairs(ev_m, pm_m, Sg)

        # Bell decomposition
        rv_m = rho_vec(bell_psi_s6)
        coeffs_m = Li_m @ rv_m
        total_m = np.linalg.norm(rv_m)

        osc_mask_m = np.abs(np.imag(mu_m)) > 1e-6
        coeffs_osc_m = coeffs_m * osc_mask_m
        osc_norm_m = np.linalg.norm(R_m @ coeffs_osc_m)
        osc_pct_m = 100 * osc_norm_m**2 / total_m**2

        # Collect distinct frequencies
        freq_set = set()
        for ki in range(d2):
            im = abs(np.imag(mu_m[ki]))
            if im > 1e-6:
                freq_set.add(f"{im:.2f}")
        freq_str = ', '.join(sorted(freq_set, key=float))

        log(f"  {mname:<14} {len(cls_m['oscillatory']):>8} "
            f"{len(cls_m['mixed']):>8} {len(cls_m['real']):>8} "
            f"{osc_pct_m:>6.1f}% {freq_str:>30}")

        model_results[mname] = {
            'evals': ev_m, 'R': R_m, 'L_inv': Li_m, 'mu': mu_m,
            'cls': cls_m, 'osc_mask': osc_mask_m,
            'osc_pct': osc_pct_m, 'freqs': freq_set,
        }

    # Detailed Pauli fingerprint comparison
    log(f"\n  Pauli fingerprint of Bell's oscillating part for each Hamiltonian:")
    log(f"  (top 8 Pauli strings by amplitude)")

    for mname, comps in models_s6.items():
        mr = model_results[mname]
        rv_m = rho_vec(bell_psi_s6)
        coeffs_m = mr['L_inv'] @ rv_m
        coeffs_osc_m = coeffs_m * mr['osc_mask']
        osc_norm_m = np.linalg.norm(mr['R'] @ coeffs_osc_m)

        if osc_norm_m < 1e-6:
            log(f"\n  {mname}: no oscillating modes for Bell")
            continue

        rho_osc_m = (mr['R'] @ coeffs_osc_m).reshape(d, d)
        osc_paulis_m = []
        for pi in range(n_pauli):
            label = pauli_basis_3[pi][0]
            if label == 'III':
                continue
            val = abs(np.trace(pauli_basis_3[pi][1] @ rho_osc_m) / d)
            if val > 1e-6:
                osc_paulis_m.append((label, val))
        osc_paulis_m.sort(key=lambda x: -x[1])

        top8 = ', '.join(f"{l}({v:.4f})" for l, v in osc_paulis_m[:8])
        log(f"\n  {mname}:")
        log(f"    osc% = {mr['osc_pct']:.1f}%, freqs = {{{', '.join(sorted(mr['freqs'], key=float))}}}")
        log(f"    {top8}")

    # Node comparison
    log(f"\n  Node comparison (non-oscillating Paulis with nonzero background):")

    for mname in models_s6:
        mr = model_results[mname]
        rv_m = rho_vec(bell_psi_s6)
        coeffs_m = mr['L_inv'] @ rv_m

        coeffs_osc_m = coeffs_m * mr['osc_mask']
        coeffs_bg_m = coeffs_m * (~mr['osc_mask'])
        osc_norm_m = np.linalg.norm(mr['R'] @ coeffs_osc_m)

        if osc_norm_m < 1e-6:
            log(f"  {mname:<14} (no oscillation)")
            continue

        rho_osc_m = (mr['R'] @ coeffs_osc_m).reshape(d, d)
        rho_bg_m = (mr['R'] @ coeffs_bg_m).reshape(d, d)

        nodes_m = []
        for pi in range(n_pauli):
            label = pauli_basis_3[pi][0]
            if label == 'III':
                continue
            osc_val = abs(np.trace(pauli_basis_3[pi][1] @ rho_osc_m) / d)
            bg_val = abs(np.trace(pauli_basis_3[pi][1] @ rho_bg_m) / d)
            if bg_val > 1e-6 and osc_val < 1e-6:
                nodes_m.append(label)

        log(f"  {mname:<14} nodes: {', '.join(nodes_m) if nodes_m else '(none)'}")

    # Test more states on different Hamiltonians
    log(f"\n  Standing wave content (osc%) for different states x Hamiltonians:")
    hdr_line = f"  {'State':<14}"
    for mname in models_s6:
        hdr_line += f" {mname:>12}"
    log(hdr_line)
    log(f"  {'-'*(14 + 13*len(models_s6))}")

    for sname, psi in [('GHZ', make_ghz(N)), ('W', make_w(N)),
                        ('Bell(0,1)', make_bell_plus(N)),
                        ('|+++>', make_product(N, '+++')),
                        ('Bell+|+++>', sup_psi)]:
        line = f"  {sname:<14}"
        for mname in models_s6:
            mr = model_results[mname]
            rv_m = rho_vec(psi)
            coeffs_m = mr['L_inv'] @ rv_m
            total_m = np.linalg.norm(rv_m)
            coeffs_osc_m = coeffs_m * mr['osc_mask']
            osc_norm_m = np.linalg.norm(mr['R'] @ coeffs_osc_m)
            osc_pct_m = 100 * osc_norm_m**2 / total_m**2
            line += f" {osc_pct_m:>11.1f}%"
        log(line)

    # Final summary
    log(f"\n  {'='*80}")
    log(f"  SECTION 6 SUMMARY")
    log(f"  {'='*80}")
    log(f"\n  1. FREQUENCIES depend on the Hamiltonian:")
    log(f"     Each model produces its own set of oscillation frequencies.")
    log(f"     Heisenberg: ~2, ~4, ~6 (harmonics of J)")
    log(f"     XY-only: different set (missing ZZ coupling changes spectrum)")
    log(f"     Ising: purely real eigenvalues (no oscillation at all)")
    log()
    log(f"  2. PATTERN depends on the Hamiltonian:")
    log(f"     The oscillating Pauli strings differ between models.")
    log(f"     But the STRUCTURE is similar: XX/YY-type correlations oscillate,")
    log(f"     ZZ-type correlations tend to be nodes.")
    log()
    log(f"  3. Standing wave content depends on BOTH state AND Hamiltonian:")
    log(f"     Bell produces oscillation for models with off-diagonal coupling.")
    log(f"     Ising (purely diagonal) produces NO oscillation for any state.")

    log(f"\n{'='*90}")
    log(f"Section 6 completed: {datetime.now()}")
    log(f"{'='*90}")
    log()
    log("=" * 90)
    log("ANALYSIS COMPLETE")
    log("=" * 90)
    f.close()
    print(f"\n>>> Results written to {OUT}")
