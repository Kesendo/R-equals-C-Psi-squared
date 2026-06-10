# Verifier for JDefectLightMigrationClaim (compute/RCPsiSquared.Diagnostics).
#
# Along the jdefect axis (uniform XY chain + Z-dephasing, J-bond defect δJ on
# bond (0,1); the palindrome holds at every δJ but the spectrum moves), the
# per-mode absorption identity tracks the motion exactly:
#
#   (a) for EVERY eigenmode v(δJ):  Re λ(δJ) = −2γ·light(v(δJ))  at every δJ
#       (δJ-pointwise; the Re-drift IS light migration of the eigenvectors,
#       the Hamiltonian contributes zero throughout),
#   (b) the N+1 kernel modes stay at light ≡ 0 (U(1) protection, no migration),
#   (c) palindrome partners migrate oppositely: complete pairing
#       λ_s + λ_f = −2Nγ with light_s + light_f = N at every δJ,
#   (d) the content of the axis: the per-mode light migration
#       |light(v(δJ_max)) − light(v(0))| is O(δJ) and nonzero (marks move
#       exactly as much as the eigenvectors' light migrates).
#
# Conventions match the C# layer (see F8PartnerLightComplementarityTests and
# simulations/f8_partner_light_sanity.py): row-major vec (x = a·d + b), site 0
# = most significant bit, H = Σ_bonds (J_b/2)(XX+YY), dissipator diagonal
# −2γ·popcount(a⊕b); light(v) = Σ_x popcount-weight(x)·|v_x|² / ‖v‖².
import sys

import numpy as np

if sys.stdout.encoding and sys.stdout.encoding.lower() != 'utf-8':
    sys.stdout.reconfigure(encoding='utf-8')

X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
I2 = np.eye(2, dtype=complex)

TOL_IDENTITY = 1e-9
TOL_KERNEL = 1e-10
TOL_PAIR = 1e-8


def site_pair(N, i, j, P, Q):
    ops = [I2] * N
    ops[i] = P
    ops[j] = Q
    out = np.array([[1.0 + 0j]])
    for o in ops:
        out = np.kron(out, o)
    return out


def liouvillian(N, delta_j, gamma, J=1.0, defect_bond=0):
    d = 1 << N
    H = np.zeros((d, d), dtype=complex)
    for b in range(N - 1):
        Jb = J + (delta_j if b == defect_bond else 0.0)
        H += Jb / 2.0 * (site_pair(N, b, b + 1, X, X) + site_pair(N, b, b + 1, Y, Y))
    Id = np.eye(d)
    L = -1j * (np.kron(H, Id) - np.kron(Id, H.T))
    w = np.array([bin((x // d) ^ (x % d)).count("1") for x in range(d * d)], dtype=float)
    L += np.diag(-2.0 * gamma * w)
    return L, w


def eig_with_light(L, w, gamma):
    vals, R = np.linalg.eig(L)
    prob = np.abs(R) ** 2
    light = (w @ prob) / prob.sum(axis=0)
    identity_resid = np.abs(vals.real + 2.0 * gamma * light)
    return vals, R, light, identity_resid


def palindrome_pairing(vals, light, N, gamma):
    """Greedy multiset pairing λ_s + λ_f ≈ −2Nγ; returns the worst pairing
    residual and the worst per-pair complementarity defect |light_s+light_f−N|."""
    target = -2.0 * N * gamma
    unpaired = list(range(len(vals)))
    worst_pair = 0.0
    worst_comp = 0.0
    while unpaired:
        i = unpaired.pop(0)
        resid = [abs(vals[i] + vals[j] - target) for j in unpaired]
        k = int(np.argmin(resid))
        worst_pair = max(worst_pair, resid[k])
        j = unpaired.pop(k)
        worst_comp = max(worst_comp, abs(light[i] + light[j] - N))
    return worst_pair, worst_comp


def match_by_overlap(R0, R1):
    """Greedy mode matching δJ=0 → δJ=max by maximal eigenvector overlap.
    Inside degenerate clusters at δJ=0 the matched basis is eig's choice, so
    the migration there reads basis rotation too; the reported maximum is the
    descriptive content of the axis, not an asserted identity."""
    A = np.abs(R0.conj().T @ R1)  # columns are unit-norm from np.linalg.eig
    n = A.shape[0]
    order = np.argsort(A, axis=None)[::-1]
    used_i = np.zeros(n, dtype=bool)
    used_j = np.zeros(n, dtype=bool)
    match = np.empty(n, dtype=int)
    found = 0
    for flat in order:
        i, j = divmod(int(flat), n)
        if used_i[i] or used_j[j]:
            continue
        match[i] = j
        used_i[i] = used_j[j] = True
        found += 1
        if found == n:
            break
    return match


def run(N, gamma=0.05, deltas=(0.0, 0.05, 0.1)):
    d2 = (1 << N) ** 2
    print(f"N={N} (dim {d2}), XY chain J=1, γ={gamma} uniform, defect bond (0,1), δJ ∈ {list(deltas)}")
    lights = {}
    vecs = {}
    for dj in deltas:
        L, w = liouvillian(N, dj, gamma)
        vals, R, light, resid = eig_with_light(L, w, gamma)

        # (a) the absorption identity, δJ-pointwise, every mode
        assert resid.max() < TOL_IDENTITY, f"identity broken at δJ={dj}: {resid.max():.3e}"

        # (b) kernel: exactly N+1 modes at λ ≈ 0, all with light ≡ 0
        kernel = np.where(np.abs(vals) < TOL_KERNEL)[0]
        assert len(kernel) == N + 1, f"kernel count {len(kernel)} ≠ {N + 1} at δJ={dj}"
        kernel_light = light[kernel].max()
        assert kernel_light < TOL_IDENTITY, f"kernel light {kernel_light:.3e} at δJ={dj}"

        # (c) complete palindrome pairing with per-pair light complementarity
        worst_pair, worst_comp = palindrome_pairing(vals, light, N, gamma)
        assert worst_pair < TOL_PAIR, f"pairing incomplete at δJ={dj}: {worst_pair:.3e}"
        assert worst_comp < TOL_IDENTITY, f"complementarity broken at δJ={dj}: {worst_comp:.3e}"

        lights[dj] = light
        vecs[dj] = R
        print(f"  δJ={dj:>4}: max |Re λ + 2γ·light| = {resid.max():.3e}; "
              f"kernel {len(kernel)} modes, max light {kernel_light:.3e}; "
              f"pairing complete (worst |λ_s+λ_f+2Nγ| = {worst_pair:.3e}), "
              f"max |light_s+light_f−N| = {worst_comp:.3e}")

    # (d) the content: per-mode light migration δJ=0 → δJ=max
    dj0, djm = deltas[0], deltas[-1]
    match = match_by_overlap(vecs[dj0], vecs[djm])
    migration = np.abs(lights[djm][match] - lights[dj0])
    print(f"  light migration δJ={dj0}→{djm}: max = {migration.max():.6f}, "
          f"mean = {migration.mean():.6f}, modes moved > 1e-6: {(migration > 1e-6).sum()}/{d2}")
    print()


if __name__ == "__main__":
    run(4)
    run(5)
    print("ALL ASSERTIONS PASSED: the identity is δJ-pointwise, the kernel is dark, "
          "partner migrations cancel; the marks move as the light migrates.")
