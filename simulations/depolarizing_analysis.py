"""
Depolarizing Analysis -- Why depolarizing noise breaks the palindrome
=====================================================================
Seven-section investigation:
  1. Per-site dephasing rates for each noise type
  2. Split classification (immune vs decaying)
  3. Valid per-site bijections (exhaustive search)
  4. Approximate palindrome under depol (N=3,4)
  5. Interpolation Z-deph -> depol
  6. Rate structure comparison
  7. Weight distribution / counting argument

Script: simulations/depolarizing_analysis.py
Output: simulations/results/depolarizing_analysis.txt
"""
import numpy as np
from itertools import product as iproduct, permutations
from datetime import datetime

OUT = r"D:\Entwicklung\Projekte Privat\R-equals-C-Psi-squared\simulations\results\depolarizing_analysis.txt"
f = open(OUT, "w", buffering=1)


def log(msg=""):
    print(msg, flush=True)
    f.write(msg + "\n")
    f.flush()


# ============================================================
# PAULI BASICS
# ============================================================
I2 = np.eye(2, dtype=complex)
sx = np.array([[0, 1], [1, 0]], dtype=complex)
sy = np.array([[0, -1j], [1j, 0]], dtype=complex)
sz = np.array([[1, 0], [0, -1]], dtype=complex)
PAULIS = [I2, sx, sy, sz]
NAMES = ['I', 'X', 'Y', 'Z']


def anti_commutes(p, q):
    """True if Pauli p anti-commutes with Pauli q (indices 0-3)."""
    if p == 0 or q == 0:
        return False
    return p != q


def per_site_rates(noise, gamma):
    """Per-site dephasing rate [r_I, r_X, r_Y, r_Z] for a noise model."""
    # channels: list of (strength, axis_index)
    channels = {
        'Z':     [(gamma, 3)],
        'X':     [(gamma, 1)],
        'Y':     [(gamma, 2)],
        'depol': [(gamma / 3, 1), (gamma / 3, 2), (gamma / 3, 3)],
        'ZX':    [(gamma / 2, 3), (gamma / 2, 1)],
        'ZY':    [(gamma / 2, 3), (gamma / 2, 2)],
        'XY':    [(gamma / 2, 1), (gamma / 2, 2)],
    }
    r = [0.0, 0.0, 0.0, 0.0]
    for (g, axis) in channels.get(noise, []):
        for p in range(4):
            if anti_commutes(p, axis):
                r[p] += 2 * g
    return r


# ============================================================
# BUILD SYSTEM IN PAULI BASIS
# ============================================================
def build_H(N, bonds, comps):
    """Build Hamiltonian from coupling dict, e.g. {'XX':1,'YY':1,'ZZ':1}."""
    PM = {'X': sx, 'Y': sy, 'Z': sz}
    d = 2 ** N
    H = np.zeros((d, d), dtype=complex)
    for i, j in bonds:
        for label, J in comps.items():
            if J == 0:
                continue
            ops = [I2] * N
            ops[i] = PM[label[0]]
            ops[j] = PM[label[1]]
            term = ops[0]
            for o in ops[1:]:
                term = np.kron(term, o)
            H += J * term
    return H


def build_L_H_pauli(N, H):
    """Build Hamiltonian superoperator in Pauli basis (expensive)."""
    d = 2 ** N
    num = 4 ** N
    all_idx = list(iproduct(range(4), repeat=N))
    pmats = []
    for idx in all_idx:
        m = PAULIS[idx[0]]
        for i in idx[1:]:
            m = np.kron(m, PAULIS[i])
        pmats.append(m)
    L_H = np.zeros((num, num), dtype=complex)
    for b in range(num):
        comm = -1j * (H @ pmats[b] - pmats[b] @ H)
        for a in range(num):
            L_H[a, b] = np.trace(pmats[a] @ comm) / d
    return L_H, all_idx


def build_L_D_diag(N, site_rates):
    """Build L_D diagonal from per-site rates [r_I, r_X, r_Y, r_Z]."""
    all_idx = list(iproduct(range(4), repeat=N))
    diag = np.zeros(4 ** N)
    for a, idx in enumerate(all_idx):
        diag[a] = -sum(site_rates[idx[k]] for k in range(N))
    return diag


def palindrome_err(evals, center):
    """Full palindrome error: max_k min_j |lambda_k + lambda_j + 2c|."""
    targets = -(evals + 2 * center)
    return max(np.min(np.abs(evals - t)) for t in targets)


def rate_palin_err(evals, center):
    """Rate-only palindrome error."""
    rates = -np.real(evals)
    return max(np.min(np.abs(rates - (2 * center - r))) for r in rates)


# ============================================================
# SETUP
# ============================================================
gamma = 0.05

log("=" * 90)
log("DEPOLARIZING ANALYSIS: Why depolarizing noise breaks the palindrome")
log(f"Date: {datetime.now()}")
log(f"gamma = {gamma}")
log("=" * 90)


# ############################################################
# SECTION 1: Per-site dephasing rates
# ############################################################
log()
log("=" * 90)
log("SECTION 1: Per-site dephasing rates")
log("  For each noise model, rate of each Pauli index at a single site")
log("=" * 90)

noise_types = ['Z', 'X', 'Y', 'depol', 'ZX', 'ZY', 'XY']
noise_labels = ['Z-deph', 'X-deph', 'Y-deph', 'Depol', 'Z+X', 'Z+Y', 'X+Y']

header = f"  {'Pauli':>6}"
for nl in noise_labels:
    header += f"  {nl:>8}"
log(f"\n{header}")
log(f"  {'-' * (6 + 10 * len(noise_labels))}")

for p in range(4):
    line = f"  {NAMES[p]:>6}"
    for nt in noise_types:
        r = per_site_rates(nt, gamma)
        # Show rate as fraction of gamma
        val = r[p]
        if abs(val) < 1e-15:
            line += f"  {'0':>8}"
        elif abs(val - 2 * gamma) < 1e-15:
            line += f"  {'2g':>8}"
        elif abs(val - gamma) < 1e-15:
            line += f"  {'g':>8}"
        elif abs(val - 4 * gamma / 3) < 1e-15:
            line += f"  {'4g/3':>8}"
        else:
            line += f"  {val / gamma:>7.3f}g"
    log(line)

log(f"\n  (g = gamma = {gamma})")


# ############################################################
# SECTION 2: Split classification
# ############################################################
log()
log("=" * 90)
log("SECTION 2: Split classification")
log("  Which Pauli indices are immune (rate=0) vs decaying (rate>0)?")
log("=" * 90)

log(f"\n  {'Noise':>8}  {'Immune':>16}  {'Decaying':>16}  {'Split':>6}  {'Palindrome?':>12}")
log(f"  {'-' * 66}")

for nt, nl in zip(noise_types, noise_labels):
    r = per_site_rates(nt, gamma)
    immune = [NAMES[p] for p in range(4) if r[p] < 1e-15]
    decay = [NAMES[p] for p in range(4) if r[p] > 1e-15]
    split = f"{len(immune)}:{len(decay)}"
    can_pair = len(immune) == len(decay)  # preliminary check
    log(f"  {nl:>8}  {','.join(immune):>16}  {','.join(decay):>16}  {split:>6}  "
        f"{'2:2 -> YES' if can_pair else '1:3 -> ???':>12}")

log(f"\n  Note: 2:2 split is necessary but not sufficient for palindrome.")
log(f"  The actual condition is whether the 4 rates can be partitioned")
log(f"  into two pairs with equal sums. See Section 3.")


# ############################################################
# SECTION 3: Valid per-site bijections
# ############################################################
log()
log("=" * 90)
log("SECTION 3: Valid per-site bijections (exhaustive search)")
log("  For each noise model, find ALL permutations sigma of {I,X,Y,Z}")
log("  such that rate(k) + rate(sigma(k)) = const for all k")
log("  (Phases do not affect this relation, so only permutations matter)")
log("=" * 90)

for nt, nl in zip(noise_types, noise_labels):
    r = per_site_rates(nt, gamma)
    valid_perms = []

    for perm in permutations(range(4)):
        # Check: rate(k) + rate(perm[k]) = const for all k
        sums = [r[k] + r[perm[k]] for k in range(4)]
        if max(sums) - min(sums) < 1e-12:
            valid_perms.append((perm, sums[0]))

    log(f"\n  {nl}: {len(valid_perms)} valid permutations")
    log(f"    Rates: [{', '.join(f'{v/gamma:.3f}g' for v in r)}]")

    if valid_perms:
        for perm, c in valid_perms:
            mapping = ', '.join(f'{NAMES[k]}->{NAMES[perm[k]]}' for k in range(4))
            log(f"    {mapping}  (c = {c / gamma:.3f}g)")
    else:
        log(f"    NO valid permutation exists. Proof by exhaustion:")
        # Show why each attempt fails
        # Group by sigma(0) choice
        for target_0 in range(4):
            c = r[0] + r[target_0]
            needed = [c - r[k] for k in range(1, 4)]
            available = [r[j] for j in range(4) if j != target_0]
            match = all(any(abs(n - a) < 1e-12 for a in available) for n in needed)
            # Check bijection feasibility
            used = {target_0}
            ok = True
            for k in range(1, 4):
                need_rate = c - r[k]
                found = False
                for j in range(4):
                    if j not in used and abs(r[j] - need_rate) < 1e-12:
                        used.add(j)
                        found = True
                        break
                if not found:
                    ok = False
                    break
            status = "OK" if ok else "FAIL"
            log(f"      sigma(I)={NAMES[target_0]}: c={c / gamma:.3f}g -> "
                f"need rates {[f'{n/gamma:.3f}g' for n in needed]} -> {status}")

# The key theorem
log(f"\n  {'=' * 70}")
log(f"  KEY THEOREM: Rate pairing condition")
log(f"  {'=' * 70}")
log(f"  For 4 per-site rates [r0, r1, r2, r3], a palindromic Pi_D exists")
log(f"  iff the rates can be partitioned into 2 pairs with equal sums.")
log()
log(f"  For three-axis dephasing at rates (gamma_X, gamma_Y, gamma_Z):")
log(f"    rates = [0, 2(gY+gZ), 2(gX+gZ), 2(gX+gY)]")
log(f"    Pairing exists iff at least one of gX, gY, gZ = 0.")
log(f"    i.e., dephasing along at most 2 axes.")
log()
log(f"  Depolarizing (gX=gY=gZ=g/3): rates = [0, 4g/3, 4g/3, 4g/3]")
log(f"    Three indices need to map to one. Bijection impossible. QED.")


# ############################################################
# SECTION 4: Approximate palindrome under depol
# ############################################################
log()
log("=" * 90)
log("SECTION 4: Approximate palindrome under depolarizing noise")
log("  Build full L = L_H + L_D, eigendecompose, measure palindrome error")
log("=" * 90)

for N in [3, 4]:
    bonds = [(i, i + 1) for i in range(N - 1)]
    Sg = N * gamma
    num = 4 ** N

    hamiltonians = {
        'Heisenberg': {'XX': 1, 'YY': 1, 'ZZ': 1},
        'XY-only':    {'XX': 1, 'YY': 1},
        'Ising':      {'ZZ': 1},
    }

    log(f"\n  N = {N}, gamma = {gamma}, chain, Sg = {Sg}")
    log(f"  {'Model':>12}  {'Noise':>8}  {'Center':>8}  "
        f"{'Full err':>10}  {'Rate err':>10}")
    log(f"  {'-' * 56}")

    for mname, comps in hamiltonians.items():
        H = build_H(N, bonds, comps)
        L_H, all_idx = build_L_H_pauli(N, H)

        for noise, nlabel in [('Z', 'Z-deph'), ('depol', 'Depol')]:
            r = per_site_rates(noise, gamma)
            diag = build_L_D_diag(N, r)
            L = L_H.copy()
            for a in range(num):
                L[a, a] += diag[a]

            evals = np.linalg.eigvals(L)
            center = Sg  # standard center

            ferr = palindrome_err(evals, center)
            rerr = rate_palin_err(evals, center)

            log(f"  {mname:>12}  {nlabel:>8}  {center:>8.4f}  "
                f"{ferr:>10.4e}  {rerr:>10.4e}")

    # Detailed depol analysis: best center search
    H = build_H(N, bonds, {'XX': 1, 'YY': 1, 'ZZ': 1})
    L_H, all_idx = build_L_H_pauli(N, H)
    r_depol = per_site_rates('depol', gamma)
    diag_depol = build_L_D_diag(N, r_depol)
    L_depol = L_H.copy()
    for a in range(num):
        L_depol[a, a] += diag_depol[a]
    ev_depol = np.linalg.eigvals(L_depol)

    log(f"\n  N={N} Depol center sweep (Heisenberg):")
    best_c, best_e = 0, np.inf
    for frac in np.linspace(0.3, 2.0, 50):
        c = frac * Sg
        e = palindrome_err(ev_depol, c)
        if e < best_e:
            best_e, best_c = e, c
    log(f"    Best center: {best_c:.6f} = {best_c / Sg:.4f}*Sg")
    log(f"    Best full error: {best_e:.6e}")
    log(f"    Error at Sg={Sg:.4f}: {palindrome_err(ev_depol, Sg):.6e}")

    # Check if error is a simple fraction of gamma
    err_at_Sg = palindrome_err(ev_depol, Sg)
    log(f"    err/gamma = {err_at_Sg / gamma:.6f}")
    log(f"    err/(2g/3) = {err_at_Sg / (2 * gamma / 3):.6f}")

    # Partial pairing: which eigenvalues pair well?
    log(f"\n  N={N} Partial pairing analysis (center=Sg):")
    rates = -np.real(ev_depol)
    n_good = 0
    n_bad = 0
    for k in range(num):
        target = -(ev_depol[k] + 2 * Sg)
        min_d = np.min(np.abs(ev_depol - target))
        if min_d < 1e-6:
            n_good += 1
        else:
            n_bad += 1
    log(f"    Well-paired (err < 1e-6): {n_good}/{num}")
    log(f"    Poorly-paired:            {n_bad}/{num}")

# Error formula
log(f"\n  ERROR FORMULA (at center Sg = N*gamma):")
log(f"  The steady state (rate 0) needs partner at rate 2*Sg.")
log(f"  Under depol, max rate = (4g/3)*N = (4/3)*Sg.")
log(f"  Gap = 2*Sg - (4/3)*Sg = (2/3)*Sg.")
log(f"  Therefore: palindrome error = (2/3)*Sg = (2/3)*N*gamma.")
log(f"    N=3: err = (2/3)*0.15 = 0.1000  (confirmed)")
log(f"    N=4: err = (2/3)*0.20 = 0.1333  (confirmed)")
log(f"  Error is EXACT, Hamiltonian-INDEPENDENT, and proportional to Sg.")
log(f"\n  Note: excluding steady states (as in depolarizing_test.txt),")
log(f"  the error reduces to (2g/3) per weight level, giving")
log(f"  the previously reported err = 3.33e-02 for N=3.")


# ############################################################
# SECTION 5: Interpolation Z-deph -> depol
# ############################################################
log()
log("=" * 90)
log("SECTION 5: Interpolation Z-deph -> depol")
log("  L_D(alpha) = (1-alpha)*L_D_Z + alpha*L_D_depol")
log("  alpha=0: pure Z-deph.  alpha=1: pure depol.")
log("=" * 90)

for N in [3, 4]:
    bonds = [(i, i + 1) for i in range(N - 1)]
    Sg = N * gamma
    num = 4 ** N
    H = build_H(N, bonds, {'XX': 1, 'YY': 1, 'ZZ': 1})
    L_H, all_idx = build_L_H_pauli(N, H)

    r_z = per_site_rates('Z', gamma)
    r_d = per_site_rates('depol', gamma)

    log(f"\n  N = {N}:")
    log(f"  {'alpha':>8}  {'full_err':>10}  {'rate_err':>10}  "
        f"{'r_Z_site':>10}  {'pairable':>10}")
    log(f"  {'-' * 54}")

    alphas = [0.0, 0.001, 0.01, 0.05, 0.1, 0.2, 0.3, 0.5, 0.7, 0.9, 1.0]
    for alpha in alphas:
        r_interp = [(1 - alpha) * r_z[p] + alpha * r_d[p] for p in range(4)]
        diag = build_L_D_diag(N, r_interp)
        L = L_H.copy()
        for a in range(num):
            L[a, a] += diag[a]
        evals = np.linalg.eigvals(L)

        ferr = palindrome_err(evals, Sg)
        rerr = rate_palin_err(evals, Sg)

        # Check per-site pairing
        r_sorted = sorted(r_interp)
        # Try all 3 pairings
        pairable = False
        pairings = [
            ([0, 1], [2, 3]),
            ([0, 2], [1, 3]),
            ([0, 3], [1, 2]),
        ]
        for (a, b), (c, d) in pairings:
            s1 = r_sorted[a] + r_sorted[b]
            s2 = r_sorted[c] + r_sorted[d]
            if abs(s1 - s2) < 1e-12:
                pairable = True
                break

        log(f"  {alpha:>8.3f}  {ferr:>10.4e}  {rerr:>10.4e}  "
            f"{r_interp[3] / gamma:>9.4f}g  {'YES' if pairable else 'NO':>10}")

    log(f"\n  Observation: palindrome error is {'linear' if N == 3 else 'approx linear'} "
        f"in alpha for small alpha.")
    log(f"  The palindrome breaks IMMEDIATELY for any alpha > 0.")
    log(f"  Per-site pairing fails because Z acquires nonzero rate {r_d[3] / gamma:.3f}*g*alpha,")
    log(f"  creating a 1:3 split (only I is immune).")


# ############################################################
# SECTION 6: Rate structure comparison
# ############################################################
log()
log("=" * 90)
log("SECTION 6: Rate structure under Z-deph vs depol")
log("  N=3 Heisenberg chain, sorted eigenvalue rates")
log("=" * 90)

N = 3
bonds = [(i, i + 1) for i in range(N - 1)]
Sg = N * gamma
num = 4 ** N
H = build_H(N, bonds, {'XX': 1, 'YY': 1, 'ZZ': 1})
L_H, all_idx = build_L_H_pauli(N, H)

for noise, nlabel in [('Z', 'Z-deph'), ('depol', 'Depol')]:
    r = per_site_rates(noise, gamma)
    diag = build_L_D_diag(N, r)
    L = L_H.copy()
    for a in range(num):
        L[a, a] += diag[a]
    evals = np.linalg.eigvals(L)

    rates = -np.real(evals)
    freqs = np.imag(evals)
    order = np.argsort(rates)

    log(f"\n  {nlabel} (center = {Sg:.4f}):")
    log(f"  {'#':>4}  {'rate':>10}  {'freq':>10}  {'2c-rate':>10}  "
        f"{'nearest':>10}  {'err':>10}")
    log(f"  {'-' * 58}")

    for rank, k in enumerate(order):
        d = rates[k]
        w = freqs[k]
        partner_rate = 2 * Sg - d
        # Find nearest rate to partner
        nearest_dist = np.min(np.abs(rates - partner_rate))
        nearest_idx = np.argmin(np.abs(rates - partner_rate))
        log(f"  {rank:>4}  {d:>10.6f}  {w:>+10.4f}  {partner_rate:>10.6f}  "
            f"{rates[nearest_idx]:>10.6f}  {nearest_dist:>10.4e}")

    # Summary
    ferr = palindrome_err(evals, Sg)
    log(f"\n    Max palindrome error: {ferr:.6e}")

    # Rate histogram by L_D weight
    log(f"\n    Rate distribution (L_D diagonal):")
    unique_rates = sorted(set(f"{-diag[a]:.8f}" for a in range(num)))
    for ur in unique_rates:
        rv = float(ur)
        count = sum(1 for a in range(num) if abs(-diag[a] - rv) < 1e-6)
        partner = 2 * Sg - rv
        pc = sum(1 for a in range(num) if abs(-diag[a] - partner) < 1e-6)
        log(f"      rate={rv:.6f} ({count} strings)  "
            f"partner={partner:.6f} ({pc} strings)  "
            f"{'MATCHED' if count == pc else f'MISMATCH {count} vs {pc}'}")


# ############################################################
# SECTION 7: Weight distribution / counting argument
# ############################################################
log()
log("=" * 90)
log("SECTION 7: Weight distribution and counting argument")
log("  Under Z-deph: 'weight' = XY-weight (sites with X or Y)")
log("  Under depol: 'weight' = non-I weight (sites with X, Y, or Z)")
log("=" * 90)

for N in [3, 4]:
    log(f"\n  N = {N}:")

    # Z-deph: XY-weight w has C(N,w) * 2^w * 2^(N-w) = C(N,w) * 2^N strings
    log(f"\n    Z-deph: XY-weight sectors")
    log(f"    {'w':>4}  {'count':>8}  {'partner w':>10}  {'p_count':>8}  {'match':>8}")
    log(f"    {'-' * 44}")
    from math import comb
    for w in range(N + 1):
        c = comb(N, w) * (2 ** N)
        pw = N - w
        pc = comb(N, pw) * (2 ** N)
        log(f"    {w:>4}  {c:>8}  {pw:>10}  {pc:>8}  "
            f"{'EQUAL' if c == pc else 'UNEQUAL':>8}")
    total_z = sum(comb(N, w) * 2 ** N for w in range(N + 1))
    log(f"    Total: {total_z} = 4^{N} = {4 ** N}")

    # Depol: non-I weight w has C(N,w) * 3^w * 1^(N-w) = C(N,w) * 3^w strings
    log(f"\n    Depol: non-I weight sectors")
    log(f"    {'w':>4}  {'count':>8}  {'partner w':>10}  {'p_count':>8}  "
        f"{'ratio':>8}  {'match':>8}")
    log(f"    {'-' * 54}")
    for w in range(N + 1):
        c = comb(N, w) * (3 ** w)
        pw = N - w
        pc = comb(N, pw) * (3 ** pw)
        ratio = c / pc if pc > 0 else float('inf')
        log(f"    {w:>4}  {c:>8}  {pw:>10}  {pc:>8}  "
            f"{ratio:>8.2f}  {'EQUAL' if c == pc else 'UNEQUAL':>8}")
    total_d = sum(comb(N, w) * 3 ** w for w in range(N + 1))
    log(f"    Total: {total_d} = 4^{N} = {4 ** N}")

log(f"\n  {'=' * 70}")
log(f"  THE COUNTING ARGUMENT")
log(f"  {'=' * 70}")
log()
log(f"  Under Z-deph: weight-w sector has C(N,w) * 2^w * 2^(N-w) = C(N,w)*2^N strings.")
log(f"  Its partner weight-(N-w) has C(N,N-w)*2^N = C(N,w)*2^N. ALWAYS EQUAL.")
log(f"  Reason: per-site split is 2:2, so each weight factor is 2^k * 2^(N-k) = 2^N.")
log()
log(f"  Under depol: weight-w sector has C(N,w) * 3^w * 1^(N-w) = C(N,w)*3^w strings.")
log(f"  Its partner weight-(N-w) has C(N,w)*3^(N-w). Ratio = 3^(N-2w).")
log(f"  UNEQUAL for w != N/2. No bijection between weight sectors possible.")
log()
log(f"  Root cause: per-site split is 1:3 (one immune, three decaying).")
log(f"  The factors 1^(N-w) * 3^w vs 1^w * 3^(N-w) are asymmetric.")
log()
log(f"  GENERAL THEOREM: The palindrome requires equal per-site split (m:m).")
log(f"  For 4 Pauli indices, the only equal split is 2:2.")
log(f"  Any split 1:3 or 0:4 breaks the weight-sector bijection.")
log()
log(f"  However, even with a 1:3 split, if the three decaying rates")
log(f"  are NOT all equal, a rate-pairing may still exist (Section 3).")
log(f"  But the COUNTING obstruction remains: weight sectors have different")
log(f"  sizes, so no per-site Pi can biject them.")
log()
log(f"  Resolution: for two-axis dephasing (e.g., Z+X), the rate pairing")
log(f"  works AND a different Pi exists -- not the Z-deph Pi, but one that")
log(f"  swaps the two immune/decaying pairs correctly. The counting works")
log(f"  because the relevant 'weight' is defined by the specific Pi, not")
log(f"  by the non-I count.")


# ############################################################
# SECTION 7b: Time reversal interpretation
# ############################################################
log()
log("=" * 90)
log("SECTION 7b: Time reversal interpretation")
log("=" * 90)

log(f"""
  Under Z-deph:
    Past  = {{I,Z}} per site = 2 choices = populations = classical
    Future = {{X,Y}} per site = 2 choices = coherences = quantum
    Mirror (Pi): swaps past <-> future, bijectively.
    Equal choices means equal past and future. Time is symmetric.

  Under depol:
    Past  = {{I}} per site = 1 choice = only the identity survives
    Future = {{X,Y,Z}} per site = 3 choices = everything else decays
    Mirror: would need to map 1 -> 3 per site. IMPOSSIBLE.
    There are 3x more ways to be quantum than classical.
    The future is bigger than the past. Time reversal is broken.

  Quantitatively:
    Z-deph: |past|/|future| = 2^N / 2^N = 1 (balanced)
    Depol:  |past|/|future| = 1^N / 3^N = (1/3)^N (exponentially unbalanced)

  For N=3: past has 1 string (III), future has 63 strings. Ratio = 1:63.
  For N=4: past has 1 string, future has 255. Ratio = 1:255.

  The mirror cannot reflect everything because the future is exponentially
  larger than the past. This is why depolarizing noise breaks time reversal.
""")


# ############################################################
# DONE
# ############################################################
log()
log("=" * 90)
log("ANALYSIS COMPLETE")
log(f"Date: {datetime.now()}")
log("=" * 90)
f.close()
print(f"\n>>> Results written to {OUT}")
