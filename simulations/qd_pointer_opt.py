"""qd_pointer_opt.py — optimizing the pointers.

Continuation of qd_scout.py (the pointer door). Question: MAXIMIZE the record
quality / redundancy over bond strengths (per-bond Delta), topology (leaves on S),
and gamma — and find what survives J > 0.

The closed form this scout gates (J=0, arbitrary graph, per-bond Delta_ab,
per-site gamma; |+>^N start; z in {+1,-1}, bit 0 -> +1):

  rho_{ab}[(za,zb),(za',zb')] = 1/4
      * exp(-i t Delta_ab (za zb - za' zb'))
      * exp(-2 t (gamma_a [za!=za'] + gamma_b [zb!=zb']))
      * prod_{k not in {a,b}} cos(t [Delta_ak (za-za') + Delta_bk (zb-zb')])

(Derivation: at J=0 both -i[H,.] for H = sum Delta_ab Z_a Z_b and the dephasing
mask are diagonal in the |z><z'| basis; the partial trace sets z_k = z'_k on
traced sites, so traced-site gammas drop out and each traced z_k = +-1 averages
the S/j-facing phases into the cos product.)

Laws read off the closed form (all gated numerically below):
  LAW A (write/forgive/blind parities): a witness j adjacent to S holds a
    PERFECT record iff the write angle 2 Delta_Sj t = pi/2 (mod pi) AND every
    watcher angle 2 Delta_jk t = 0 (mod pi). At the readout t* = pi/(4 Delta_Sj)
    this is arithmetic: watcher ratio r = Delta_jk / Delta_Sj even integer ->
    forgiven (record survives), odd integer -> exactly blind (leaf law is the
    special case "no watcher"), non-integer -> generic contrast |prod cos|.
  LAW B (alignment beats topology): R_perfect = #{j ~ S : all watcher ratios
    even} <= deg(S), and the bound is achievable on ANY graph by choosing
    even-ratio env bonds. Uniform coupling makes every watcher ratio 1 (odd),
    so uniformly only leaves record — the leaf law as parity, not as geometry.
  LAW C (the gamma race, exact): record distinguishability (trace distance of the
    two conditional pointer states)
    D_j(t) = exp(-2 gamma_j t) * |sin(2 Delta_Sj t)| * prod_k |cos(2 Delta_jk t)|;
    for an unwatched witness its argmax is t_opt = arctan(Delta_Sj/gamma_j)
    / (2 Delta_Sj). t_opt maximizes D, NOT the mutual information (MI is a
    different functional; away from t* the pair is not classical-quantum).
  OPEN (J > 0): hopping breaks the commuting structure; the J-sweep at the end
    measures how fast the aligned records degrade.

Conventions as in qd_scout.py: site 0 = MSB, Pauli couplings,
D[rho] = sum gamma_l (Z rho Z - rho) as the entrywise mask, delta_QD = 0.1,
H_ref = 1 bit (X^N-pinned pointer entropy).
"""

import numpy as np
from itertools import combinations
import time as walltime

np.set_printoptions(linewidth=200, precision=6, suppress=True)

I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.array([[1, 0], [0, -1]], dtype=complex)

DELTA_QD = 0.10
N = 8
DT = 0.005
TSTAR = np.pi / 4          # readout for Delta_S-bond = 1


def op_at(n, site, P):
    M = np.array([[1.0]], dtype=complex)
    for s in range(n):
        M = np.kron(M, P if s == site else I2)
    return M


def build_H(n, wbonds, J):
    """wbonds: list of (i, j, Delta_ij). Hopping J uniform on the same bonds."""
    d = 2 ** n
    H = np.zeros((d, d), dtype=complex)
    for (i, j, Delta) in wbonds:
        if J != 0.0:
            H += J * (op_at(n, i, X) @ op_at(n, j, X) + op_at(n, i, Y) @ op_at(n, j, Y))
        if Delta != 0.0:
            Zi = op_at(n, i, Z)
            Zj = op_at(n, j, Z)
            H += Delta * (Zi @ Zj)
    return H


def deph_mask(n, gammas):
    d = 2 ** n
    idx = np.arange(d)
    mask = np.zeros((d, d))
    for l in range(n):
        bit = (idx >> (n - 1 - l)) & 1
        diff = (bit[:, None] != bit[None, :])
        mask -= 2.0 * gammas[l] * diff
    return mask


def rk4_evolve(rho, H, mask, dt, steps):
    def rhs(r):
        return -1j * (H @ r - r @ H) + mask * r
    for _ in range(steps):
        k1 = rhs(rho)
        k2 = rhs(rho + 0.5 * dt * k1)
        k3 = rhs(rho + 0.5 * dt * k2)
        k4 = rhs(rho + dt * k3)
        rho = rho + (dt / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4)
    return 0.5 * (rho + rho.conj().T)


def ptrace(rho, n, keep):
    keep = sorted(keep)
    t = rho.reshape((2,) * (2 * n))
    row = list(range(n))
    col = [(n + i if i in keep else i) for i in range(n)]
    out = keep + [n + i for i in keep]
    red = np.einsum(t, row + col, out)
    k = len(keep)
    return red.reshape((2 ** k, 2 ** k))


def vn_entropy(rho):
    w = np.linalg.eigvalsh(rho)
    w = w[w > 1e-14]
    return float(-(w * np.log2(w)).sum())


def mutual_info(rho, n, A, B):
    SA = vn_entropy(ptrace(rho, n, A))
    SB = vn_entropy(ptrace(rho, n, B))
    SAB = vn_entropy(ptrace(rho, n, sorted(A + B)))
    return SA + SB - SAB


def h2(p):
    s = 0.0
    for v in (p, 1.0 - p):
        if v > 1e-15:
            s -= v * np.log2(v)
    return s


def mi_of_pair(rho4):
    S1 = vn_entropy(np.array([[rho4[0, 0] + rho4[1, 1], rho4[0, 2] + rho4[1, 3]],
                              [rho4[2, 0] + rho4[3, 1], rho4[2, 2] + rho4[3, 3]]]))
    S2 = vn_entropy(np.array([[rho4[0, 0] + rho4[2, 2], rho4[0, 1] + rho4[2, 3]],
                              [rho4[1, 0] + rho4[3, 2], rho4[1, 1] + rho4[3, 3]]]))
    return S1 + S2 - vn_entropy(rho4)


def closed_pair(t, wbonds, gammas, a, b):
    """The exact J=0 pair state rho_{ab}(t); basis order (z_a, z_b), a < b."""
    assert a < b
    n = N
    D = np.zeros((n, n))
    for (i, j, Delta) in wbonds:
        D[i, j] = D[j, i] = Delta
    zval = [1.0, -1.0]                     # bit 0 -> z=+1
    rho = np.zeros((4, 4), dtype=complex)
    for ra in range(2):
        for rb in range(2):
            for ca in range(2):
                for cb in range(2):
                    za, zb = zval[ra], zval[rb]
                    za_, zb_ = zval[ca], zval[cb]
                    val = 0.25 * np.exp(-1j * t * D[a, b] * (za * zb - za_ * zb_))
                    val *= np.exp(-2 * t * (gammas[a] * (ra != ca) + gammas[b] * (rb != cb)))
                    for k in range(n):
                        if k in (a, b):
                            continue
                        arg = t * (D[a, k] * (za - za_) + D[b, k] * (zb - zb_))
                        if arg != 0.0:
                            val *= np.cos(arg)
                    rho[2 * ra + rb, 2 * ca + cb] = val
    return rho


def record_quality(t, wbonds, gammas, S, j):
    """LAW C: D_j(t) = e^{-2 gamma_j t} |sin(2 Delta_Sj t)| prod |cos(2 Delta_jk t)|."""
    D = np.zeros((N, N))
    for (i, jj, Delta) in wbonds:
        D[i, jj] = D[jj, i] = Delta
    q = np.exp(-2 * gammas[j] * t) * abs(np.sin(2 * D[S, j] * t))
    for k in range(N):
        if k in (S, j) or D[j, k] == 0.0:
            continue
        q *= abs(np.cos(2 * D[j, k] * t))
    return q


def evolve_to(wbonds, J, gammas, t_target):
    # dt adjusted so the endpoint lands EXACTLY on t_target (the parity laws
    # live at exact angles; a grid-rounded readout misses pi by ~2e-3 and
    # costs ~1e-5 MI, which a 1e-6 gate sees).
    H = build_H(N, wbonds, J)
    mask = deph_mask(N, gammas)
    plus = np.array([1, 1], dtype=complex) / np.sqrt(2)
    psi = np.array([1.0], dtype=complex)
    for _ in range(N):
        psi = np.kron(psi, plus)
    rho = np.outer(psi, psi.conj())
    steps = max(1, int(np.ceil(t_target / DT)))
    dt = t_target / steps
    return rk4_evolve(rho, H, mask, dt, steps), t_target


def fragment_sweep(rho, S, env, tag):
    by_size = {}
    for m in range(1, len(env) + 1):
        vals = [mutual_info(rho, N, [S], list(F)) for F in combinations(env, m)]
        by_size[m] = float(np.mean(vals))
    thr = (1 - DELTA_QD) * 1.0
    m_delta = next((m for m in range(1, len(env) + 1) if by_size[m] >= thr), None)
    R = (len(env) / m_delta) if m_delta else 0.0
    curve = "  ".join(f"m={m}:{by_size[m]:.4f}" for m in range(1, len(env) + 1))
    print(f"    [{tag}] Ibar(m): {curve}")
    print(f"    [{tag}] m_delta={m_delta}  R_delta={R:.3f}")
    return R


GATES = []


def gate(name, ok, detail):
    GATES.append((name, ok))
    print(f"  [GATE {name}] {detail}  {'OK' if ok else 'FAIL'}")


def chain_wbonds(delta_S_bonds, delta_env, S):
    """Chain 0-1-...-7; bonds touching S get delta_S_bonds, the rest delta_env."""
    wb = []
    for i in range(N - 1):
        d = delta_S_bonds if S in (i, i + 1) else delta_env
        wb.append((i, i + 1, d))
    return wb


def broom_wbonds(k, delta_spoke, delta_tail):
    """S=0 with k leaf spokes (1..k) and a tail 0-(k+1)-(k+2)-...-7."""
    wb = [(0, s, delta_spoke) for s in range(1, k + 1)]
    if k < N - 1:
        wb.append((0, k + 1, delta_spoke))
        for i in range(k + 1, N - 1):
            wb.append((i, i + 1, delta_tail))
    return wb


def main():
    t0 = walltime.time()
    print(f"QD pointer optimization  N={N}  dt={DT}  delta={DELTA_QD}  t*=pi/4")
    print("Closed form gated: J=0 pair reduction with per-bond Delta, per-site gamma.")

    S = 3
    env = [j for j in range(N) if j != S]

    # ------------------------------------------------------------------
    print()
    print("=" * 100)
    print("GATE 1 — closed form vs RK4 (three configs x with/without gamma, t grid)")
    print("=" * 100)
    configs = [
        ("uniform chain", chain_wbonds(1.0, 1.0, S), [0.0] * N),
        ("aligned chain (env Delta=2)", chain_wbonds(1.0, 2.0, S), [0.0] * N),
        ("aligned chain + gamma_E=0.05", chain_wbonds(1.0, 2.0, S),
         [0.05 if j != S else 0.0 for j in range(N)]),
        ("aligned broom k=4", broom_wbonds(4, 1.0, 2.0), [0.0] * N),
    ]
    for (name, wb, gams) in configs:
        Sc = 0 if "broom" in name else S
        worst = 0.0
        for t_target in [0.3, TSTAR, 1.1]:
            rho, t_now = evolve_to(wb, 0.0, gams, t_target)
            for j in range(N):
                if j == Sc:
                    continue
                a, b = min(Sc, j), max(Sc, j)
                num = ptrace(rho, N, [a, b])
                ana = closed_pair(t_now, wb, gams, a, b)
                worst = max(worst, float(np.abs(num - ana).max()))
        gate(f"closed-form {name}", worst < 1e-7,
             f"max|numeric - closed| over all pairs, t in {{0.3, pi/4, 1.1}}: {worst:.2e}")

    # ------------------------------------------------------------------
    print()
    print("=" * 100)
    print("LAW A — write/forgive/blind parities at t* = pi/4 (chain, S=3 interior)")
    print("=" * 100)
    print("  watcher ratio r = Delta_env / Delta_S-bond:  even -> forgiven, odd -> blind,")
    print("  non-integer -> generic contrast |cos(r pi/2)|.")
    for (ratio, expect) in [(1.0, "odd -> blind (leaf law)"), (2.0, "even -> PERFECT"),
                            (3.0, "odd -> blind"), (1.5, "generic")]:
        wb = chain_wbonds(1.0, ratio, S)
        rho, t_now = evolve_to(wb, 0.0, [0.0] * N, TSTAR)
        i2 = mutual_info(rho, N, [S], [2])
        i4 = mutual_info(rho, N, [S], [4])
        pred = mi_of_pair(closed_pair(t_now, wb, [0.0] * N, 2, 3))
        print(f"  ratio={ratio}: I(S:E_2)={i2:.6f}  I(S:E_4)={i4:.6f}  closed-form pred={pred:.6f}   [{expect}]")
        if ratio == 2.0:
            gate("LAW-A forgiven", abs(i2 - 1.0) < 1e-6 and abs(i4 - 1.0) < 1e-6
                 and abs(pred - 1.0) < 1e-10,
                 f"even watcher ratio: both neighbor records PERFECT "
                 f"(RK4 I={i2:.7f}, closed form 1-I={1.0 - pred:.1e})")
        if ratio == 3.0:
            gate("LAW-A blind", i2 < 1e-6 and i4 < 1e-6 and pred < 1e-12,
                 f"odd watcher ratio: both neighbors exactly blind "
                 f"(RK4 I={i2:.2e}, closed form I={pred:.1e})")
        if ratio == 1.5:
            gate("LAW-A generic", abs(i2 - pred) < 1e-6,
                 f"generic ratio: RK4 matches closed form ({i2:.6f} vs {pred:.6f})")
            # the explicit MI formula: I = 1 - h2((1+beta)/2), beta = |prod cos|
            beta = abs(np.cos(ratio * np.pi / 2))
            mi_formula = 1.0 - h2((1.0 + beta) / 2.0)
            gate("LAW-A MI formula", abs(mi_formula - pred) < 1e-9,
                 f"1 - h2((1+beta)/2) with beta=|cos(r pi/2)|={beta:.6f}: "
                 f"{mi_formula:.6f} vs closed form {pred:.6f}")

    # distance-2 in the aligned chain: shared-dresser channel closed by alignment
    wb = chain_wbonds(1.0, 2.0, S)
    rho, _ = evolve_to(wb, 0.0, [0.0] * N, TSTAR)
    i_d2 = max(mutual_info(rho, N, [S], [1]), mutual_info(rho, N, [S], [5]))
    gate("LAW-A distance-2 closed", i_d2 < 1e-8,
         f"aligned chain: distance-2 MI at t* = {i_d2:.2e} (alignment closes the dresser channel)")

    # ------------------------------------------------------------------
    print()
    print("=" * 100)
    print("LAW B — R_perfect = #even-aligned neighbors (topology x alignment)")
    print("=" * 100)
    cases = [
        ("chain S=1 uniform (end leaf)", chain_wbonds(1.0, 1.0, 1), 1, 1),
        ("chain S=3 uniform", chain_wbonds(1.0, 1.0, S), S, 0),
        ("chain S=3 aligned", chain_wbonds(1.0, 2.0, S), S, 2),
        ("broom k=2 uniform", broom_wbonds(2, 1.0, 1.0), 0, 2),
        ("broom k=2 aligned tail", broom_wbonds(2, 1.0, 2.0), 0, 3),
        ("broom k=4 uniform", broom_wbonds(4, 1.0, 1.0), 0, 4),
        ("broom k=4 aligned tail", broom_wbonds(4, 1.0, 2.0), 0, 5),
        ("star (k=7)", broom_wbonds(7, 1.0, 1.0), 0, 7),
    ]
    for (name, wb, Sc, expected) in cases:
        rho, _ = evolve_to(wb, 0.0, [0.0] * N, TSTAR)
        envc = [j for j in range(N) if j != Sc]
        singles = [mutual_info(rho, N, [Sc], [j]) for j in envc]
        n_perfect = sum(1 for v in singles if v > 1.0 - 1e-6)
        prof = "  ".join(f"{v:.4f}" for v in singles)
        print(f"  {name:26s} singles: {prof}")
        gate(f"LAW-B {name}", n_perfect == expected,
             f"R_perfect={n_perfect} (predicted {expected})")
        if "aligned" in name or "star" in name:
            fragment_sweep(rho, Sc, envc, name)

    # ------------------------------------------------------------------
    print()
    print("=" * 100)
    print("LAW B' — the plaquette record (4-cycle, S and j opposite corners)")
    print("=" * 100)
    print("  Predicted BEFORE running (closed form): on a uniform 4-cycle both dressers are")
    print("  shared, the sum angle 4*Delta*t* = pi and the difference angle 0 both have |cos|=1,")
    print("  so the distance-2 witness holds a PERFECT 1-bit record in the Bell channel")
    print("  (both double coherences 1/4, eigenvalues {1/2,1/2,0,0}), while both neighbors")
    print("  stay blind (watcher ratio 1, odd). Law B's neighbor-only reading is TREE-scoped.")
    print("  NOTE the channel: <X_S X_j> = 1, <Z_S Z_j> = 0, the Z-conditional states of j")
    print("  are identical -- this is a record of the ANTI-pointer X_S, counted by MI,")
    print("  invisible to the pointer Z_S and contributing nothing to pointer redundancy.")
    sq = [(0, 1, 1.0), (1, 2, 1.0), (2, 3, 1.0), (3, 0, 1.0)]   # sites 4-7 unbonded
    rho, _ = evolve_to(sq, 0.0, [0.0] * N, TSTAR)
    i_opp = mutual_info(rho, N, [0], [2])
    i_nb = max(mutual_info(rho, N, [0], [1]), mutual_info(rho, N, [0], [3]))
    gate("LAW-B' plaquette record", abs(i_opp - 1.0) < 1e-6 and i_nb < 1e-6,
         f"I(S:E_opposite)={i_opp:.7f} (predicted 1), neighbors max I={i_nb:.2e} (predicted 0)")
    pair = ptrace(rho, N, [0, 2])              # basis (z_S, z_j)
    exx = float(np.real(np.trace(pair @ np.kron(X, X))))
    ezz = float(np.real(np.trace(pair @ np.kron(Z, Z))))
    dcond = float(np.abs(2 * pair[:2, :2] - 2 * pair[2:, 2:]).max())
    gate("LAW-B' anti-pointer channel",
         abs(exx - 1.0) < 1e-7 and abs(ezz) < 1e-7 and dcond < 1e-7,
         f"<X_S X_j>={exx:.7f} (predicted 1), <Z_S Z_j>={ezz:.1e} (predicted 0), "
         f"Z-conditional states of j identical to {dcond:.1e}")

    # ------------------------------------------------------------------
    print()
    print("=" * 100)
    print("LAW C — the gamma race, exact (aligned chain, gamma on witnesses)")
    print("=" * 100)
    gamsE = [0.05 if j != S else 0.0 for j in range(N)]
    wb = chain_wbonds(1.0, 2.0, S)
    rho, t_now = evolve_to(wb, 0.0, gamsE, TSTAR)
    i2 = mutual_info(rho, N, [S], [2])
    pred = mi_of_pair(closed_pair(t_now, wb, gamsE, 2, 3))
    gate("LAW-C gamma-dressed MI", abs(i2 - pred) < 1e-6,
         f"I(S:E_2) at t* with gamma_E=0.05: RK4 {i2:.6f} vs closed form {pred:.6f}")
    beta = np.exp(-2 * 0.05 * TSTAR)          # aligned: |prod cos| = 1, radius = e^{-2 gamma_j t*}
    mi_formula = 1.0 - h2((1.0 + beta) / 2.0)
    gate("LAW-C MI formula", abs(mi_formula - pred) < 1e-9,
         f"1 - h2((1+beta)/2) with beta=e^(-2 gamma t*)={beta:.6f}: "
         f"{mi_formula:.6f} vs closed form {pred:.6f}")
    # t_opt of the record quality D(t) (closed form, fine grid) vs the arctan law
    # for the UNWATCHED witness case: leaf E_0 next to S=1 (chain), gamma_0 = 0.05.
    wb_leaf = [(i, i + 1, 1.0) for i in range(N - 1)]
    gam_leaf = [0.0] * N
    gam_leaf[0] = 0.05
    ts = np.linspace(0.01, np.pi / 2, 4000)
    Dvals = [record_quality(t, wb_leaf, gam_leaf, 1, 0) for t in ts]
    t_num = float(ts[int(np.argmax(Dvals))])
    t_ana = float(np.arctan(1.0 / 0.05) / 2.0)
    gate("LAW-C t_opt", abs(t_num - t_ana) < 2e-3,
         f"argmax D(t) leaf record: numeric {t_num:.4f} vs arctan(Delta/gamma)/(2Delta) = {t_ana:.4f}")
    print(f"    (t* = pi/4 = {TSTAR:.4f}; the collector pulls the distinguishability peak EARLIER by {TSTAR - t_ana:.4f};")
    print(f"     t_opt maximizes D = trace distance, not the MI)")

    # ------------------------------------------------------------------
    print()
    print("=" * 100)
    print("OPEN — J > 0: how fast does hopping degrade the aligned records?")
    print("=" * 100)
    wb = chain_wbonds(1.0, 2.0, S)
    print("  aligned chain, gamma=0, readout at t* = pi/4:")
    js = [0.0, 0.05, 0.1, 0.2, 0.4]
    deficits = []
    for J in js:
        rho, _ = evolve_to(wb, J, [0.0] * N, TSTAR)
        i2 = mutual_info(rho, N, [S], [2])
        i4 = mutual_info(rho, N, [S], [4])
        deficits.append(1.0 - 0.5 * (i2 + i4))
        print(f"    J={J:4.2f}:  I(S:E_2)={i2:.6f}  I(S:E_4)={i4:.6f}  deficit={deficits[-1]:.6f}")
    fit = [d for (J, d) in zip(js, deficits) if J > 0]
    if all(d > 0 for d in fit):
        lj = np.log([J for J in js if J > 0])
        ld = np.log(fit)
        slope = float(np.polyfit(lj, ld, 1)[0])
        print(f"  power-law fit: deficit ~ J^{slope:.2f}  (2.0 = perturbatively quadratic)")
    rho, _ = evolve_to(wb, 0.1, [0.0] * N, TSTAR)
    fragment_sweep(rho, S, env, "aligned chain J=0.1 at t*")

    # ------------------------------------------------------------------
    print()
    print("=" * 100)
    nfail = sum(1 for (_, ok) in GATES if not ok)
    for (name, ok) in GATES:
        print(f"  {'OK  ' if ok else 'FAIL'} {name}")
    print(f"GATES: {len(GATES) - nfail}/{len(GATES)} OK")
    print(f"total wall time: {walltime.time() - t0:.1f} s")


if __name__ == "__main__":
    main()
