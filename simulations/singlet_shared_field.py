"""
The singlet under a shared field: exact closed forms against exact propagation
==============================================================================
experiments/SINGLET_UNDER_A_SHARED_FIELD.md. Two qubits, the Heisenberg bond
H = J (XX + YY + ZZ) (J = 1, the Pauli book). For every noise axis a in a set S (S = {z} or
S = {x, y, z}) each site feels sigma_a at rate gamma, and eta in [-1, 1] is the correlation
coefficient of the two sites' fields: per axis the Kossakowski matrix is gamma [[1, eta], [eta, 1]],
realised by the shared jump sqrt(|eta| gamma) (s_a x 1 + sign(eta) 1 x s_a) and the local jumps
sqrt((1 - |eta|) gamma) s_a on each site.

The script checks, by exact propagation (expm of the Liouvillian) against closed forms:
  - every site alone dephases at rate gamma for every eta, on random mixed states;
  - the general law for a pure state: initial fidelity loss gamma sum_a [Var A_a + Var B_a + 2 eta Cov_a],
    A_a = s_a x 1, B_a = 1 x s_a, Cov_a = <A_a B_a> - <A_a><B_a>, on random pure states; at t = 0
    correlation lowers the loss of a state whose summed covariance is negative and raises it where positive;
  - for the Bell states (Var = 1, Cov_a = s_a, the parity under s_a x s_a) the rule 2 gamma sum_a (1 + s_a eta);
  - which Bell states each shared jump annihilates;
  - the singlet's trajectory is the eta = 0 one at the dose (1 - eta) gamma t, with its closed forms
    (z only: f = exp(-4 (1 - eta) gamma t), C Psi = f^2/3; x, y, z: Werner p = exp(-8 (1 - eta) gamma t),
    C Psi = p (3p - 1)/6; C Psi in the concurrence book, concurrence times l1/3);
  - the crossing of C Psi = 1/4 at t_x = K/((1 - eta) gamma), K = ln(4/3)/8 (z, F14's concurrence
    constant) and K = ln(6/(1 + sqrt19))/8 (x, y, z), root-found on exact propagation;
  - under z only, Psi+ crosses with the singlet; under x, y, z the triplets cross sooner as eta grows;
  - at eta = 1 the singlet is a fixed point (exactly), and under x, y, z the kernel is span{1, |s><s|}.

Tolerances: expm's error is of order eps ||L t||, about 1e-14 at the times evaluated (t <= 7; the
crossing scan stops at the first sign change, before t = 4); trajectories and crossing times are gated at
1e-12, and every measured deviation is printed (they come out near 1e-15). The triplets' common crossing
is read through the concurrence's square roots of near-degenerate eigenvalues, an error of order
sqrt(eps) ~ 1.5e-8, and gated there.

Import-inert; prints only, and exits with status 1 if any check fails.
"""

import sys

import numpy as np
from scipy.linalg import expm, null_space
from scipy.optimize import brentq

I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.diag([1.0, -1.0]).astype(complex)
AX = {"x": X, "y": Y, "z": Z}
H = sum(np.kron(P, P) for P in (X, Y, Z))
ID = np.eye(4)
EPS = np.finfo(float).eps
GAMMA = 0.1
STATES = {"singlet": [0, 1, -1, 0], "Psi+": [0, 1, 1, 0], "Phi+": [1, 0, 0, 1], "Phi-": [1, 0, 0, -1]}
TOL = 1e-12
CHECKS, CLAIMS = [], []


def ket(name):
    return np.array(STATES[name], dtype=complex) / np.sqrt(2)


def bell_proj(name):
    """|v><v| with exact entries 0, +-1/2."""
    v = np.array(STATES[name], dtype=float)
    return (np.outer(v, v) / 2).astype(complex)


def proj(v):
    return np.outer(v, v.conj())


def jumps(eta, axes, g=GAMMA):
    out = []
    for a in axes:
        P = AX[a]
        if eta != 0:
            out.append(np.sqrt(abs(eta) * g) * (np.kron(P, I2) + np.sign(eta) * np.kron(I2, P)))
        if abs(eta) < 1:
            out += [np.sqrt((1 - abs(eta)) * g) * np.kron(P, I2), np.sqrt((1 - abs(eta)) * g) * np.kron(I2, P)]
    return out


def liouvillian(js, with_h=True):
    L = -1j * (np.kron(H, ID) - np.kron(ID, H.T)) if with_h else np.zeros((16, 16), dtype=complex)
    for J in js:
        JdJ = J.conj().T @ J
        L += np.kron(J, J.conj()) - 0.5 * (np.kron(JdJ, ID) + np.kron(ID, JdJ.T))
    return L


def apply(L, r):
    return (L @ r.reshape(-1)).reshape(4, 4)


def evolve(L, r, t):
    return (expm(L * t) @ r.reshape(-1)).reshape(4, 4)


def cpsi(r):
    yy = np.kron(Y, Y)
    lam = np.sqrt(np.abs(np.sort(np.real(np.linalg.eigvals(r @ yy @ r.conj() @ yy)))[::-1]))
    c = max(0.0, lam[0] - lam[1] - lam[2] - lam[3])
    return c * float(np.abs(r).sum() - np.abs(np.diag(r)).sum()) / 3


def loss_rate(L, v):
    return -np.real(v.conj() @ apply(L, proj(v)) @ v)


def check(label, page, value, digits=3):
    ok = value is not None and round(value, digits) == round(page, digits)
    CHECKS.append(ok)
    shown = "none" if value is None else f"{value:.{digits + 3}f}"
    print(f"  {label:<66s} page {page:.{digits}f}  computed {shown}  {'match' if ok else 'DIFFERS'}")


def claim(label, ok, detail=""):
    CLAIMS.append(bool(ok))
    print(f"  {label} {detail}  {'holds' if ok else 'FAILS'}")


def crossing(L, r0, hi=200.0):
    f = lambda t: cpsi(evolve(L, r0, t)) - 0.25
    ts = np.linspace(0, hi, 8001)
    prev = f(ts[0])
    for k in range(1, len(ts)):
        cur = f(ts[k])
        if cur < 0 <= prev:
            return brentq(f, ts[k - 1], ts[k], xtol=1e-14)
        prev = cur
    return None


def main():
    rng = np.random.default_rng(11)
    etas = (-1.0, -0.5, 0.0, 0.5, 0.9, 1.0)
    print("Each site alone dephases at gamma for every eta (reduced dissipator, 20 random mixed states):")
    worst = 0.0
    for _ in range(20):
        A = rng.normal(size=(4, 4)) + 1j * rng.normal(size=(4, 4))
        r = A @ A.conj().T
        r /= np.trace(r)
        rA = np.einsum("ijkj->ik", r.reshape(2, 2, 2, 2))
        for eta in etas:
            for a in AX:
                drA = np.einsum("ijkj->ik", apply(liouvillian(jumps(eta, [a]), with_h=False), r).reshape(2, 2, 2, 2))
                worst = max(worst, np.abs(drA - GAMMA * (AX[a] @ rA @ AX[a] - rA)).max())
    claim("the reduced dissipator of site A is gamma (s_a rho s_a - rho), every eta and axis", worst <= TOL, f"({worst:.1e})")

    print("\nThe general law for a pure state, 200 random states, S = {x, y, z} and {z}:")
    dev, spared, sped = 0.0, 0, 0
    for _ in range(200):
        v = rng.normal(size=4) + 1j * rng.normal(size=4)
        v /= np.linalg.norm(v)
        for axes in (list(AX), ["z"]):
            cov = 0.0
            var = 0.0
            for a in axes:
                Aa, Ba = np.kron(AX[a], I2), np.kron(I2, AX[a])
                ea, eb = np.real(v.conj() @ Aa @ v), np.real(v.conj() @ Ba @ v)
                var += (1 - ea ** 2) + (1 - eb ** 2)
                cov += np.real(v.conj() @ Aa @ Ba @ v) - ea * eb
            for eta in (-0.5, 0.0, 0.7):
                dev = max(dev, abs(loss_rate(liouvillian(jumps(eta, axes)), v) - GAMMA * (var + 2 * eta * cov)))
            if len(axes) == 3:
                spared += cov < 0
                sped += cov > 0
    claim("initial fidelity loss = gamma sum_a [Var A + Var B + 2 eta Cov]", dev <= TOL, f"({dev:.1e})")
    claim("both signs of the summed covariance occur among random pure states under x, y, z",
          spared > 0 and sped > 0, f"(negative {spared}, positive {sped} of 200)")
    v = np.array([0, np.cos(0.3), -np.sin(0.3), 0], dtype=complex)
    r3 = [loss_rate(liouvillian(jumps(eta, list(AX))), v) for eta in (0.0, 0.5, 1.0)]
    for eta, page, val in zip((0.0, 0.5, 1.0), (0.464, 0.319, 0.174), r3):
        check(f"cos0.3|01> - sin0.3|10>, x, y, z, eta {eta}: not dark, but spared", page, val)
    claim("... its loss falls with eta and its shared jump does not annihilate it", r3[0] > r3[1] > r3[2] > 0)
    v = np.array([1, 0, 0, 0], dtype=complex)
    r00 = [loss_rate(liouvillian(jumps(eta, list(AX))), v) for eta in (0.0, 1.0)]
    claim("a product state (|00>): its initial loss is the same at eta 0 and 1 (<= 1e-12)",
          abs(r00[0] - r00[1]) <= TOL, f"({r00[0]:.6f}, {r00[1]:.6f})")

    print("  a finite-time reversal: a state with negative initial covariance, lower fidelity under correlation")
    found = None
    Lz0, Lz1 = liouvillian(jumps(0.0, ["z"])), liouvillian(jumps(1.0, ["z"]))
    for _ in range(300):
        v = rng.normal(size=4) + 1j * rng.normal(size=4)
        v /= np.linalg.norm(v)
        Az, Bz = np.kron(Z, I2), np.kron(I2, Z)
        cz = np.real(v.conj() @ Az @ Bz @ v) - np.real(v.conj() @ Az @ v) * np.real(v.conj() @ Bz @ v)
        if cz >= 0:
            continue
        r0 = proj(v)
        for t in (0.5, 1.0, 2.0, 5.0):
            f0 = np.real(v.conj() @ evolve(Lz0, r0, t) @ v)
            f1 = np.real(v.conj() @ evolve(Lz1, r0, t) @ v)
            if f1 < f0 - 1e-6:
                found = (cz, t, f0, f1)
                break
        if found:
            break
    claim("z only: some state with Cov_z < 0 has lower fidelity at eta 1 than at eta 0 at a finite t",
          found is not None, "" if found is None else f"(Cov_z {found[0]:.3f}, t {found[1]}: {found[3]:.4f} < {found[2]:.4f})")

    print("\nThe Bell states: Var = 1, Cov_a = s_a, so 2 gamma sum_a (1 + s_a eta):")
    parity = {name: {a: int(round(np.real(ket(name).conj() @ np.kron(AX[a], AX[a]) @ ket(name)))) for a in AX}
              for name in STATES}
    claim("parities: singlet (-, -, -); Psi+ (+, +, -), Phi+ (+, -, +), Phi- (-, +, +)",
          parity["singlet"] == {"x": -1, "y": -1, "z": -1} and parity["Psi+"] == {"x": 1, "y": 1, "z": -1}
          and parity["Phi+"] == {"x": 1, "y": -1, "z": 1} and parity["Phi-"] == {"x": -1, "y": 1, "z": 1})
    dev = 0.0
    for axes in (["z"], list(AX)):
        for eta in etas:
            L = liouvillian(jumps(eta, axes))
            for name in STATES:
                dev = max(dev, abs(loss_rate(L, ket(name)) - 2 * GAMMA * sum(1 + parity[name][a] * eta for a in axes)))
    claim("rule holds for all four Bell states, eta in {-1, -0.5, 0, 0.5, 0.9, 1}, S = {z} and {x, y, z}",
          dev <= TOL, f"({dev:.1e})")
    for eta, page in ((0.5, 0.3), (0.9, 0.06)):
        check(f"x, y, z, eta {eta}: singlet, 6 gamma (1 - eta)", page, loss_rate(liouvillian(jumps(eta, list(AX))), ket("singlet")), 6)
    for eta, page in ((0.5, 0.7), (1.0, 0.8)):
        check(f"x, y, z, eta {eta}: a triplet, (6 + 2 eta) gamma", page, loss_rate(liouvillian(jumps(eta, list(AX))), ket("Psi+")), 6)

    print("\nWhich Bell states each shared jump annihilates:")
    dark = {}
    for a in AX:
        for sgn in (1, -1):
            S = np.kron(AX[a], I2) + sgn * np.kron(I2, AX[a])
            dark[(a, sgn)] = sorted(n for n in STATES if np.abs(S @ ket(n)).max() == 0.0)
    expect = {(a, sgn): sorted(n for n in STATES if parity[n][a] == -sgn) for a in AX for sgn in (1, -1)}
    claim("the sum annihilates exactly the Bell states with s_a = -1, the difference those with s_a = +1 (exact)",
          dark == expect, "(" + "; ".join(f"{a}{'+' if s > 0 else '-'}: {', '.join(v)}" for (a, s), v in dark.items()) + ")")
    claim("only the singlet is annihilated by all three sums", set.intersection(*(set(dark[(a, 1)]) for a in AX)) == {"singlet"})

    print("\nThe singlet's trajectory is the eta = 0 one at the dose (1 - eta) gamma t (t = 0.5, 2, 7):")
    s0 = bell_proj("singlet")
    for eta in (-0.5, 0.0, 0.5, 0.9):
        Lz, Li = liouvillian(jumps(eta, ["z"])), liouvillian(jumps(eta, list(AX)))
        Lz0, Li0 = liouvillian(jumps(0.0, ["z"])), liouvillian(jumps(0.0, list(AX)))
        dz = di = ds = 0.0
        for t in (0.5, 2.0, 7.0):
            f = np.exp(-4 * (1 - eta) * GAMMA * t)
            rz = np.zeros((4, 4), dtype=complex)
            rz[1, 1] = rz[2, 2] = 0.5
            rz[1, 2] = rz[2, 1] = -0.5 * f
            dz = max(dz, np.abs(evolve(Lz, s0, t) - rz).max())
            p = np.exp(-8 * (1 - eta) * GAMMA * t)
            di = max(di, np.abs(evolve(Li, s0, t) - (p * s0 + (1 - p) * ID / 4)).max())
            ds = max(ds, np.abs(evolve(Lz, s0, t) - evolve(Lz0, s0, (1 - eta) * t)).max(),
                     np.abs(evolve(Li, s0, t) - evolve(Li0, s0, (1 - eta) * t)).max())
        claim(f"eta {eta}: z form with f, Werner with p, and rho_eta(t) = rho_0((1 - eta) t)",
              max(dz, di, ds) <= TOL, f"({dz:.1e}, {di:.1e}, {ds:.1e})")

    print("\nThe crossing of C Psi = 1/4 (concurrence book), root-found on exact propagation:")
    Kz, Ki = np.log(4 / 3) / 8, np.log(6 / (1 + np.sqrt(19))) / 8
    for eta, pz, pi in ((0.0, 0.360, 0.141), (0.5, 0.719, 0.283), (0.9, 3.596, 1.413)):
        tz = crossing(liouvillian(jumps(eta, ["z"])), s0)
        ti = crossing(liouvillian(jumps(eta, list(AX))), s0)
        check(f"z only, eta {eta}", pz, tz)
        check(f"x, y, z, eta {eta}", pi, ti)
        d = max(abs(tz - Kz / ((1 - eta) * GAMMA)), abs(ti - Ki / ((1 - eta) * GAMMA)))
        claim(f"eta {eta}: both equal K/((1 - eta) gamma)", d <= TOL, f"({d:.1e})")
    tp = {eta: crossing(liouvillian(jumps(eta, ["z"])), bell_proj("Psi+")) for eta in (0.0, 0.5)}
    claim("z only: Psi+ crosses with the singlet (eta 0, 0.5)",
          abs(tp[0.0] - Kz / GAMMA) <= TOL and abs(tp[0.5] - Kz / (0.5 * GAMMA)) <= TOL, f"({tp[0.0]:.6f}, {tp[0.5]:.6f})")
    tt = [crossing(liouvillian(jumps(eta, list(AX))), bell_proj("Psi+")) for eta in (0.0, 0.5, 1.0)]
    claim("x, y, z: a triplet crosses sooner as eta grows (eta 0, 0.5, 1)", tt[0] > tt[1] > tt[2], "(" + ", ".join(f"{t:.4f}" for t in tt) + ")")
    check("x, y, z, eta 0: a triplet crosses at", 0.141, tt[0])
    check("x, y, z, eta 0.5: a triplet crosses at", 0.129, tt[1])

    print("\neta = 1:")
    for axes in (["z"], list(AX)):
        dev = np.abs(apply(liouvillian(jumps(1.0, axes)), s0)).max()
        claim(f"S = {{{', '.join(axes)}}}: L(|s><s|) = 0 (exactly; |s><s| with exact entries +-1/2)", dev == 0.0, f"({dev:.1e})")
    L1 = liouvillian(jumps(1.0, list(AX)))
    ker = null_space(L1, rcond=1e-10)
    basis = np.column_stack([ID.reshape(-1) / 2, s0.reshape(-1)])
    resid = np.abs(ker - basis @ np.linalg.lstsq(basis, ker, rcond=None)[0]).max()
    claim("S = {x, y, z}: ker L = span{1, |s><s|}, every Werner state stationary", ker.shape[1] == 2 and resid <= TOL,
          f"(dim {ker.shape[1]}, {resid:.1e})")
    tts = {tr: crossing(L1, bell_proj(tr)) for tr in ("Psi+", "Phi+", "Phi-")}
    check("S = {x, y, z}, eta 1: the triplets cross at", 0.119, tts["Psi+"])
    # each triplet's trajectory stays Bell-diagonal with its own weight the largest, so concurrence and
    # l1 = |w_Phi+ - w_Phi-| + |w_Psi+ - w_s| run alike for the three; the root-found values differ through
    # the concurrence's square roots of near-degenerate eigenvalues, an error of order sqrt(eps) ~ 1.5e-8
    claim("all three triplets cross at the same time (<= 1.5e-8, see comment)", max(tts.values()) - min(tts.values()) <= 1.5e-8,
          "(" + ", ".join(f"{k} {v:.6f}" for k, v in tts.items()) + ")")

    held, matched = sum(CLAIMS), sum(CHECKS)
    print(f"\n{matched} of {len(CHECKS)} printed values match at the page's precision; {held} of {len(CLAIMS)} claims hold.")
    if held < len(CLAIMS) or matched < len(CHECKS):
        sys.exit(1)


if __name__ == "__main__":
    main()
