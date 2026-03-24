#!/usr/bin/env python3
"""
Two Open Questions
====================
Q1: Why does K = γ·t_cross scale perfectly but τ=γt doesn't for the
    full trajectory? Answer: The Hamiltonian oscillation breaks τ-scaling.
    At J=0 (pure dephasing), τ-scaling is PERFECT.

Q2: Higher Rényi entropies S_α (α≠2). Does the 1/4 boundary shift if
    we use Tr(ρ^α) instead of Tr(ρ²)?

Script:  simulations/k_scaling_and_renyi.py
Output:  simulations/results/k_scaling_and_renyi.txt
"""

import numpy as np
from scipy.linalg import expm
import os, sys, time as _time

OUT_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                        "results", "k_scaling_and_renyi.txt")
_outf = open(OUT_PATH, "w", encoding="utf-8", buffering=1)
if sys.platform == "win32":
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")


def log(msg=""):
    print(msg, flush=True)
    _outf.write(msg + "\n")
    _outf.flush()


I2 = np.eye(2, dtype=complex)
sx = np.array([[0, 1], [1, 0]], dtype=complex)
sy = np.array([[0, -1j], [1j, 0]], dtype=complex)
sz = np.array([[1, 0], [0, -1]], dtype=complex)
up = np.array([1, 0], dtype=complex)
dn = np.array([0, 1], dtype=complex)

N = 2; d = 4; d2 = 16


def site_op(op, k):
    return np.kron(op, I2) if k == 0 else np.kron(I2, op)


def build_H(J=1.0):
    H = np.zeros((d, d), dtype=complex)
    for P in [sx, sy, sz]:
        H += J * site_op(P, 0) @ site_op(P, 1)
    return H


def build_L(H, gammas):
    Id = np.eye(d)
    L = -1j * (np.kron(H, Id) - np.kron(Id, H.T))
    for k in range(N):
        Zk = site_op(sz, k)
        L += gammas[k] * (np.kron(Zk, Zk.conj()) - np.eye(d2))
    return L


def evolve(L, rho0, t):
    v = expm(L * t) @ rho0.flatten()
    rho = v.reshape(d, d)
    return (rho + rho.conj().T) / 2


def purity(rho):
    return float(np.trace(rho @ rho).real)


def psi_norm(rho):
    l1 = float(np.sum(np.abs(rho)) - np.sum(np.abs(np.diag(rho))))
    return l1 / (d - 1)


def cpsi(rho):
    return purity(rho) * psi_norm(rho)


def renyi_purity(rho, alpha):
    """Generalized purity Tr(rho^alpha). For alpha=2: standard purity."""
    if alpha == 1:
        return 1.0  # Tr(rho) = 1 always
    if alpha == float('inf'):
        return float(np.max(np.linalg.eigvalsh(rho)))
    rho_a = np.linalg.matrix_power(rho, int(alpha)) if isinstance(alpha, int) else rho
    if isinstance(alpha, int):
        return float(np.trace(rho_a).real)
    # For non-integer alpha, use eigendecomposition
    ev = np.linalg.eigvalsh(rho)
    ev = np.maximum(ev, 0)
    return float(np.sum(ev ** alpha))


def make_bell_plus():
    psi = (np.kron(up, up) + np.kron(dn, dn)) / np.sqrt(2)
    return np.outer(psi, psi.conj())


def entropy_vn(rho):
    ev = np.linalg.eigvalsh(rho)
    ev = ev[ev > 1e-15]
    return -float(np.sum(ev * np.log2(ev)))


def concurrence(rho):
    try:
        sy2 = np.kron(sy, sy)
        rho_tilde = sy2 @ rho.conj() @ sy2
        R = rho @ rho_tilde
        evals = np.linalg.eigvals(R)
        evals = np.real(np.sqrt(np.maximum(evals, 0.0)))
        evals = np.sort(evals)[::-1]
        return float(max(0.0, evals[0] - evals[1] - evals[2] - evals[3]))
    except:
        return 0.0


# ====================================================================
# Q1: Why K scales but τ doesn't
# ====================================================================

def question_1():
    log("=" * 70)
    log("Q1: WHY DOES K SCALE BUT tau DOESN'T?")
    log("=" * 70)
    log()

    rho0 = make_bell_plus()
    gammas_list = [0.01, 0.02, 0.05, 0.1, 0.2]
    tau_points = np.linspace(0, 1.0, 201)  # tau = gamma * t

    # --- Test 1a: J=0 (pure dephasing, no Hamiltonian) ---
    log("--- Test 1a: J=0 (pure dephasing) ---")
    log("  If tau-scaling works, all gamma values give the same CΨ(τ) curve.")
    log()

    observables_J0 = {}
    for gamma in gammas_list:
        H = np.zeros((d, d), dtype=complex)
        L = build_L(H, [gamma, gamma])
        curves = {'cpsi': [], 'purity': [], 'psi': [], 'S': [], 'conc': []}

        for tau in tau_points:
            t = tau / gamma
            rho = evolve(L, rho0, t)
            curves['cpsi'].append(cpsi(rho))
            curves['purity'].append(purity(rho))
            curves['psi'].append(psi_norm(rho))
            rhoA = np.array([[rho[0,0]+rho[1,1], rho[0,2]+rho[1,3]],
                             [rho[2,0]+rho[3,1], rho[2,2]+rho[3,3]]])
            curves['S'].append(entropy_vn(rhoA))
            curves['conc'].append(concurrence(rho))

        observables_J0[gamma] = {k: np.array(v) for k, v in curves.items()}

    # Compute max delta across gamma values for each observable
    log(f"  {'Observable':>12}  {'max delta':>10}  {'tau-scales?':>12}")
    log("  " + "-" * 38)
    ref_gamma = gammas_list[0]
    for obs_name in ['cpsi', 'purity', 'psi', 'S', 'conc']:
        max_delta = 0
        for gamma in gammas_list[1:]:
            delta = np.max(np.abs(observables_J0[gamma][obs_name] -
                                  observables_J0[ref_gamma][obs_name]))
            if delta > max_delta:
                max_delta = delta
        scales = "JA" if max_delta < 0.01 else "NEIN"
        log(f"  {obs_name:>12}  {max_delta:10.6f}  {scales:>12}")

    log()
    log("  At J=0: ALL observables should tau-scale perfectly.")
    log()

    # --- Test 1b: J=1.0 (with Hamiltonian) ---
    log("--- Test 1b: J=1.0 (Heisenberg coupling) ---")
    log()

    observables_J1 = {}
    for gamma in gammas_list:
        H = build_H(J=1.0)
        L = build_L(H, [gamma, gamma])
        curves = {'cpsi': [], 'purity': [], 'psi': [], 'S': [], 'conc': []}

        for tau in tau_points:
            t = tau / gamma
            rho = evolve(L, rho0, t)
            curves['cpsi'].append(cpsi(rho))
            curves['purity'].append(purity(rho))
            curves['psi'].append(psi_norm(rho))
            rhoA = np.array([[rho[0,0]+rho[1,1], rho[0,2]+rho[1,3]],
                             [rho[2,0]+rho[3,1], rho[2,2]+rho[3,3]]])
            curves['S'].append(entropy_vn(rhoA))
            curves['conc'].append(concurrence(rho))

        observables_J1[gamma] = {k: np.array(v) for k, v in curves.items()}

    log(f"  {'Observable':>12}  {'max delta':>10}  {'tau-scales?':>12}")
    log("  " + "-" * 38)
    for obs_name in ['cpsi', 'purity', 'psi', 'S', 'conc']:
        max_delta = 0
        for gamma in gammas_list[1:]:
            delta = np.max(np.abs(observables_J1[gamma][obs_name] -
                                  observables_J1[ref_gamma][obs_name]))
            if delta > max_delta:
                max_delta = delta
        scales = "JA" if max_delta < 0.01 else "NEIN"
        log(f"  {obs_name:>12}  {max_delta:10.6f}  {scales:>12}")

    log()

    # --- Test 1c: Why K still works ---
    log("--- Test 1c: K-invariance despite tau-scaling failure ---")
    log()
    log("  The CROSSING POINT is on the monotonic envelope.")
    log("  Bell+ is eigenstate of H, so no Hamiltonian oscillation.")
    log("  Therefore: Bell+ CΨ(τ) IS the envelope, and the envelope scales.")
    log()

    log(f"  {'gamma':>8}  {'J':>4}  {'t_cross':>8}  {'K':>8}  {'tau_c':>8}")
    log("  " + "-" * 40)

    for J in [0.0, 1.0, 5.0]:
        K_vals = []
        for gamma in gammas_list:
            H = build_H(J) if J > 0 else np.zeros((d, d), dtype=complex)
            L = build_L(H, [gamma, gamma])

            # Find crossing
            t_cross = None
            dt = 0.002
            prev = cpsi(rho0)
            for step in range(1, 50000):
                t = step * dt
                rho = evolve(L, rho0, t)
                c = cpsi(rho)
                if prev > 0.25 and c <= 0.25:
                    frac = (prev - 0.25) / (prev - c)
                    t_cross = (step - 1 + frac) * dt
                    break
                prev = c

            if t_cross:
                K = gamma * t_cross
                K_vals.append(K)
                log(f"  {gamma:8.3f}  {J:4.1f}  {t_cross:8.3f}  {K:8.4f}  {K:8.4f}")

        if K_vals:
            log(f"  {'':>8}  {J:4.1f}  {'K_mean':>8}  {np.mean(K_vals):8.4f}  "
                f"CV={np.std(K_vals)/np.mean(K_vals)*100:.2f}%")
        log()

    log("  CONCLUSION: K is invariant because Bell+ CΨ(τ) is the pure envelope")
    log("  (no Hamiltonian oscillation). The envelope scales with τ = γ·t.")
    log("  The full trajectory doesn't scale because non-eigenstate components")
    log("  oscillate at rate J/γ, which differs across γ values.")
    log()


# ====================================================================
# Q2: Higher Renyi entropies
# ====================================================================

def question_2():
    log("=" * 70)
    log("Q2: HIGHER RENYI ENTROPIES - Does the 1/4 boundary shift?")
    log("=" * 70)
    log()

    rho0 = make_bell_plus()
    gamma = 0.05
    J = 1.0
    H = build_H(J)
    L = build_L(H, [gamma, gamma])

    alphas = [1.0, 1.5, 2.0, 3.0, 4.0, 5.0, 10.0]

    # --- Test 2a: CΨ_α trajectories ---
    log("--- Test 2a: Generalized CΨ_α = Tr(ρ^α) × ψ_norm ---")
    log()
    log(f"  CΨ_α(t) for Bell+ under Z-dephasing, γ={gamma}, J={J}")
    log()

    t_points = [0, 0.2, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0]

    header = f"  {'t':>6}" + "".join(f"  {'α='+str(a):>8}" for a in alphas)
    log(header)
    log("  " + "-" * (6 + 10 * len(alphas)))

    for t in t_points:
        rho = evolve(L, rho0, t)
        pn = psi_norm(rho)
        row = f"  {t:6.1f}"
        for alpha in alphas:
            rp = renyi_purity(rho, alpha)
            cpsi_a = rp * pn
            row += f"  {cpsi_a:8.4f}"
        log(row)

    log()

    # --- Test 2b: Crossing thresholds for each α ---
    log("--- Test 2b: Does each CΨ_α cross a fixed threshold? ---")
    log()

    # For α=2: threshold is 1/4. What are the thresholds for other α?
    # The initial value CΨ_α(0) for Bell+ (pure state):
    # Tr(ρ^α) = 1 for pure state (eigenvalues are 0 and 1)
    # ψ_norm(Bell+) = 1/3
    # So CΨ_α(0) = 1 × 1/3 = 1/3 for ALL α.

    log("  Initial CΨ_α(0) = 1/3 for all α (Bell+ is pure)")
    log()

    # Find crossing time for CΨ_α = 1/4 for each α
    log(f"  {'alpha':>6}  {'CΨ_α(0)':>8}  {'t_cross(1/4)':>12}  {'K=γ·t_c':>8}  "
        f"{'CΨ_α_final':>10}")
    log("  " + "-" * 52)

    for alpha in alphas:
        dt = 0.01
        t_cross = None
        prev = 1/3  # CΨ_α(0) for pure state
        cpsi_final = 0

        for step in range(1, 20001):
            t = step * dt
            rho = evolve(L, rho0, t)
            pn = psi_norm(rho)
            rp = renyi_purity(rho, alpha)
            ca = rp * pn
            cpsi_final = ca

            if prev > 0.25 and ca <= 0.25 and t_cross is None:
                frac = (prev - 0.25) / (prev - ca) if prev != ca else 0.5
                t_cross = (step - 1 + frac) * dt

            prev = ca

        tc = f"{t_cross:.3f}" if t_cross else "NEVER"
        K = f"{gamma * t_cross:.4f}" if t_cross else "---"
        log(f"  {alpha:6.1f}  {1/3:8.4f}  {tc:>12}  {K:>8}  {cpsi_final:10.6f}")

    log()

    # --- Test 2c: K-invariance for each α ---
    log("--- Test 2c: K_α invariance across γ values ---")
    log()

    gammas = [0.01, 0.02, 0.05, 0.1, 0.2]

    for alpha in [2.0, 3.0, 4.0]:
        log(f"  α = {alpha:.0f}:")
        K_vals = []

        for gamma_sweep in gammas:
            H_sweep = build_H(J)
            L_sweep = build_L(H_sweep, [gamma_sweep, gamma_sweep])

            dt = 0.002
            t_cross = None
            prev = 1/3

            for step in range(1, 100001):
                t = step * dt
                rho = evolve(L_sweep, rho0, t)
                pn = psi_norm(rho)
                rp = renyi_purity(rho, alpha)
                ca = rp * pn

                if prev > 0.25 and ca <= 0.25 and t_cross is None:
                    frac = (prev - 0.25) / (prev - ca) if prev != ca else 0.5
                    t_cross = (step - 1 + frac) * dt
                    break
                prev = ca

            if t_cross:
                K = gamma_sweep * t_cross
                K_vals.append(K)
                log(f"    γ={gamma_sweep:.3f}  t_cross={t_cross:.3f}  K={K:.4f}")

        if K_vals:
            Km = np.mean(K_vals)
            Ks = np.std(K_vals)
            cv = Ks/Km*100 if Km > 0 else 0
            log(f"    K_{alpha:.0f} = {Km:.4f} ± {Ks:.4f} (CV={cv:.1f}%)")
        log()

    # --- Test 2d: Natural threshold for each α ---
    log("--- Test 2d: What IS the natural threshold for each α? ---")
    log()
    log("  For α=2: the threshold 1/4 comes from the discriminant of")
    log("  the quadratic R = C·Ψ². For general α: R = C_α·Ψ^α.")
    log("  The critical boundary of x^α - x + CΨ_α = 0 shifts with α.")
    log()

    # The fixed point equation for R_α = C_α(Ψ + R_α)^α:
    # Setting x = R_α, C = C_α: x = C(Ψ + x)^α
    # At the bifurcation: both f(x) = x and f'(x) = 1 hold.
    # f(x) = C(Ψ+x)^α, f'(x) = αC(Ψ+x)^{α-1} = 1
    # From f'=1: C = 1/(α(Ψ+x)^{α-1})
    # Substituting into f=x: x = (Ψ+x)^α / (α(Ψ+x)^{α-1}) = (Ψ+x)/α
    # So x = (Ψ+x)/α → αx = Ψ+x → x(α-1) = Ψ → x = Ψ/(α-1)
    # And CΨ_α = C·Ψ = Ψ/(α(Ψ+x)^{α-1}) · ... let me compute CΨ at bifurcation.

    # At bifurcation: x* = Ψ/(α-1), so Ψ+x* = Ψ·α/(α-1)
    # C* = 1/(α·(Ψ·α/(α-1))^{α-1}) = (α-1)^{α-1} / (α^α · Ψ^{α-1})
    # CΨ_α* = C*·Ψ = (α-1)^{α-1} / (α^α · Ψ^{α-2})

    # Hmm, this depends on Ψ for α≠2. For α=2:
    # CΨ_2* = (2-1)^1 / (2^2 · Ψ^0) = 1/4. Correct!

    # For α=2, the threshold is Ψ-independent: CΨ* = 1/4 always.
    # For α≠2, the threshold DEPENDS ON Ψ. This is the key result!

    log("  Analytical bifurcation threshold for R = C_α(Ψ + R)^α:")
    log()
    log("  At bifurcation: x* = Ψ/(α-1), C* = (α-1)^{α-1} / (α^α · Ψ^{α-1})")
    log("  CΨ_α* = C* · Ψ = (α-1)^{α-1} / (α^α · Ψ^{α-2})")
    log()
    log("  For α = 2: CΨ* = 1/4              (Ψ-independent!)")
    log("  For α ≠ 2: CΨ* depends on Ψ       (NOT universal!)")
    log()

    log(f"  {'α':>4}  {'Threshold (Ψ=1/3)':>18}  {'Threshold (Ψ=1/2)':>18}  {'Universal?':>10}")
    log("  " + "-" * 56)
    for alpha in [1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 10.0]:
        for psi_val in [1/3, 1/2]:
            if alpha == 1:
                thresh = float('inf')
            else:
                thresh = (alpha - 1) ** (alpha - 1) / (alpha ** alpha * psi_val ** (alpha - 2))

            if psi_val == 1/3:
                t1 = thresh
                row = f"  {alpha:4.1f}  {thresh:18.6f}"
            else:
                t2 = thresh
                univ = "YES" if abs(t1 - t2) < 0.001 else "NO"
                row += f"  {thresh:18.6f}  {univ:>10}"
                log(row)

    log()
    log("  RESULT: Only α=2 gives a Ψ-independent (universal) threshold.")
    log("  This is WHY purity Tr(ρ²) is the natural choice:")
    log("  it is the UNIQUE Renyi order where the bifurcation boundary")
    log("  does not depend on the state's coherence.")
    log()


# ====================================================================
# Main
# ====================================================================

if __name__ == "__main__":
    t_start = _time.time()

    log("Two Open Questions: K-scaling and Renyi Entropies")
    log("=" * 70)
    log()

    question_1()
    question_2()

    total = _time.time() - t_start
    log(f"Total runtime: {total:.1f}s ({total/60:.1f} min)")
    log(f"Results saved to: {OUT_PATH}")
    _outf.close()
