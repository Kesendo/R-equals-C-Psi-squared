"""PXP's scars under Z dephasing, against Ma, Guo, Gao, Papić and Ying (arXiv:2504.12291).

Model, as their Eq. 3 and Eq. 8 (Fig. 4): a ring of N sites in the Rydberg blockade basis (bit 1 = excited,
Z = -1; no two adjacent 1s, the ring closed), H = sum_j P_(j-1) X_j P_(j+1) with coefficient 1 and P = (1 + Z)/2,
and the dissipator gamma * sum_j (Z_j rho Z_j - rho). They write P = (1 - Z)/2; X on every site carries one onto
the other and leaves the dissipator unchanged. In the basis |x><y| the dissipator is diagonal, -2 gamma times the
Hamming distance of x and y.

Part 1: C = prod_j Z_j anticommutes with H, entry by entry; U(rho) = C rho^T C then commutes with L at every gamma
        (H real; derived, not computed here: its consequence B = A is checked in Part 5). The two shifts: the palindrome centre of the full space, -N gamma, and the trace centre of the
        blockaded space, -2<Hamming> gamma = -N(1 - <z>^2) gamma, the zero of Ma et al.'s Re(lambda') as their Fig. 4(a)
        is read by eye.
Part 2: the tower. S_0 is the projection of Z2 onto the zero-energy eigenspace; S_1, S_2, ... are walked outward,
        each the projection of Z2 onto the eigenspace of largest Z2 weight between 0.5 and 2.0 above the last;
        S_-k is the projection onto the eigenspace at -E_k, its overlap with C S_k checked to have modulus 1. <S|T|S>, T the translation by
        one site, gives each level's momentum (+1: 0, -1: pi).
Part 3: first order, route 1. For each coherence |S_a><S_b| of the tower in layer 1 and for the partners that straddle
        the middle in the odd layers, |S_(m+1)><S_-m| and |S_m><S_-(m+1)| (partners under U), the block
        M[(a',b'),(a,b)] = sum_l <a'|Z_l|a><b|Z_l|b'> - N delta over all pairs of eigenvectors at that frequency
        (its eigenvalues are <= 0, the rates -gamma times them), and the 2 x 2 in the coherence and its partner with
        the U-parities of its modes. Then every coherence of the tower against its partner |S_-b><S_-a|: the half
        splitting of the 2 x 2 against the coupling sum_l <S_-b|Z_l|S_a>^2 / N, the largest by distance |a + b|, the
        U-even mode checked to be the slower wherever the half splitting exceeds 1e-9, and the coherences one and two
        steps through the middle in their full
        blocks with the U-parities of the block's modes. Rates in units of gamma N.
Part 4: route 2, the exact L. At N = 8 the dense L at gamma = 1e-5 and 1e-4: every block's eigenvalues against the
        exact ones at that frequency (the real parts, the frequency matched to 1e-4), the difference over (gamma N)^2
        at both rates; the blocks are those of layer 1, of the straddling partners and of |S_2><S_0|. At N = 10
        (--route2-n10, about twenty-five minutes) the middle partners by shift-invert.
Part 5: the revival from Z2 at N = 10, by expm_multiply: A(t) = <S_1|rho|S_0>, B(t) = <S_0|rho|S_-1>,
        O(t) = <S_2|rho|S_1>. At gamma = 0 their moduli are constant, so their decay is the dephasing's alone. The
        block's first-order prediction beside the exact amplitudes, and at the smallest rate the beat of ln|A| about
        its fitted line against the detuning of the neighbouring layer-1 coherence.
Part 6 (--fig4, about 30 minutes, 4 GB per matrix): the full spectrum at N = 10, gamma = 0.0004, in Ma et al.'s
        frame: each mode of the straddling partners as the exact eigenvalue nearest its first-order position, with its
        position, its momentum and its overlap with rho_Z2 = |Z2><Z2|; the two-end test for a symmetry reflecting the
        spectrum about the frame's zero; and a figure in the frame of their Fig. 4(a), drawn like theirs at momentum 0
        and pi only.
Runtime without options: about two minutes.
"""
import os
import sys
import time
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8", newline="\n")

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results", "pxp_scars_under_dephasing")


class Ring:
    """The blockaded ring: basis, H, the site Z values, C = prod Z, the eigenbasis of H."""

    def __init__(self, n):
        self.n = n
        m = (1 << n) - 1
        self.basis = [s for s in range(1 << n) if (s & (((s << 1) | (s >> (n - 1))) & m)) == 0]
        self.d = d = len(self.basis)
        self.idx = {s: i for i, s in enumerate(self.basis)}
        rows, cols = [], []
        for i, s in enumerate(self.basis):
            for l in range(n):
                if (s >> ((l - 1) % n)) & 1 or (s >> ((l + 1) % n)) & 1:
                    continue
                rows.append(self.idx[s ^ (1 << l)])
                cols.append(i)
        self.H = sp.csr_matrix((np.ones(len(rows)), (rows, cols)), shape=(d, d))
        self.Hd = self.H.toarray()
        self.Z = np.array([[1 - 2 * ((s >> j) & 1) for j in range(n)] for s in self.basis], dtype=float)
        self.C = np.prod(self.Z, axis=1)
        self.ham = (self.Z[:, None, :] != self.Z[None, :, :]).sum(2).astype(float)
        self.E, self.V = np.linalg.eigh(self.Hd)
        self.z2 = self.idx[int("01" * (n // 2), 2)]
        self.psi = np.zeros(d)
        self.psi[self.z2] = 1.0
        # eigenspaces (clusters of equal energy)
        order = np.argsort(self.E)
        cl, cur = [], [order[0]]
        for a, b in zip(order, order[1:]):
            if self.E[b] - self.E[a] < 1e-8:
                cur.append(b)
            else:
                cl.append(cur)
                cur = [b]
        cl.append(cur)
        self.clusters = [np.array(c) for c in cl]
        self.cE = np.array([self.E[c[0]] for c in self.clusters])
        self.cW = np.array([np.linalg.norm(self.V[:, c].T @ self.psi) ** 2 for c in self.clusters])

    _A = None

    @property
    def Zeig(self):
        """<i|Z_l|j> in the eigenbasis of H, stacked over l."""
        if self._A is None:
            self._A = np.stack([self.V.T @ (self.Z[:, l][:, None] * self.V) for l in range(self.n)])
        return self._A

    def proj(self, k):
        P = self.V[:, self.clusters[k]]
        v = P @ (P.T @ self.psi)
        return v / np.linalg.norm(v)

    def dissipator(self, X):
        return sum(self.Z[:, l][:, None] * X * self.Z[:, l][None, :] for l in range(self.n)) - self.n * X

    def U(self, X):
        return self.C[:, None] * X.T * self.C[None, :]

    def liouvillian(self, gamma):
        I = sp.identity(self.d, format="csr")
        LH = -1j * (sp.kron(self.H, I) - sp.kron(I, self.H.T))
        return (LH + sp.diags(-2.0 * gamma * self.ham.reshape(-1))).tocsr()


def tower(r):
    """S_k for k = -N/2 .. N/2, walked outward from E = 0; returns dict k -> (E, vector, weight, multiplicity)."""
    k0 = int(np.argmin(np.abs(r.cE)))
    assert abs(r.cE[k0]) < 1e-8
    S = {0: (0.0, r.proj(k0), r.cW[k0], len(r.clusters[k0]))}
    e, chosen = 0.0, [k0]
    for k in range(1, r.n // 2 + 1):
        cand = [c for c in range(len(r.cE)) if e + 0.5 < r.cE[c] < e + 2.0]
        if not cand:
            break
        c = max(cand, key=lambda c: r.cW[c])
        e = r.cE[c]
        S[k] = (e, r.proj(c), r.cW[c], len(r.clusters[c]))
        cm = int(np.argmin(np.abs(r.cE + e)))
        assert abs(r.cE[cm] + e) < 1e-8
        vm = r.proj(cm)
        ov = vm @ (r.C * S[k][1])
        assert abs(abs(ov) - 1) < 1e-9, ov
        S[-k] = (-e, vm, r.cW[cm], len(r.clusters[cm]))
    return S


def block_matrix(r, omega):
    """The first-order block at frequency omega: the pairs (i, j) of eigenvectors with E_i - E_j = omega,
    M[(i',j'),(i,j)] = sum_l <i'|Z_l|i><j|Z_l|j'> - N delta (real symmetric), and its eigensystem."""
    E, n = r.E, r.n
    pa, pb = np.nonzero(np.abs(E[:, None] - E[None, :] - omega) < 1e-8)
    A = r.Zeig
    X = A[:, pa[:, None], pa[None, :]]
    Y = A[:, pb[None, :], pb[:, None]]
    M = (X * Y).sum(0) - n * np.eye(len(pa))
    assert np.abs(M - M.T).max() < 1e-12
    w, W = np.linalg.eigh(M)
    return pa, pb, w, W


def block(r, a, b, omega):
    """First-order block at the frequency omega of |a><b| (a, b vectors): eigenvalues in gamma N, and the
    coherence's weights on the block's eigenvectors."""
    pa, pb, w, W = block_matrix(r, omega)
    c = (r.V[:, pa].T @ a) * (r.V[:, pb].T @ b)
    share = (W.T @ c) ** 2
    share /= share.sum()
    return -w / r.n, share, len(pa)


def block_modes(r, a, b, omega, ntop=2):
    """The block modes carrying the largest shares of |a><b|: (rate in gamma N, share, U-parity of the mode)."""
    pa, pb, w, W = block_matrix(r, omega)
    c = (r.V[:, pa].T @ a) * (r.V[:, pb].T @ b)
    share = (W.T @ c) ** 2
    share /= share.sum()
    out = []
    for k in np.argsort(share)[::-1][:ntop]:
        X = r.V[:, pa] @ (W[:, k][:, None] * r.V[:, pb].T)
        out.append((-w[k] / r.n, share[k], np.sum(X * r.U(X)) / np.sum(X * X)))
    return out, len(pa)


def pair2x2(r, X1, X2):
    """The 2 x 2 of the dissipator in span{X1, X2} (X2 = +-U(X1)): rates in gamma N, U-parities of the modes."""
    G = np.array([[np.sum(P * r.dissipator(Q)) for Q in (X1, X2)] for P in (X1, X2)])
    w, W = np.linalg.eigh(G)
    out = []
    for k in range(2):
        X = W[0, k] * X1 + W[1, k] * X2
        out.append((-w[k] / r.n, np.sum(X * r.U(X)) / np.sum(X * X)))
    return sorted(out), -G[0, 0] / r.n


def part1(r):
    print(f"\n== Part 1, N = {r.n}: the model ==")
    acom = np.abs(r.C[:, None] * r.Hd * r.C[None, :] + r.Hd).max()
    print(f"  blockaded dimension {r.d}; max |C H C + H| = {acom}")
    zbar = r.Z.mean(0)
    print(f"  palindrome centre of the full space: -N gamma = {-r.n} gamma")
    print(f"  trace centre of the blockaded space: -2<Hamming> gamma = {-2 * r.ham.mean():.6f} gamma;"
          f" N(1 - <z>^2) = {r.n * (1 - zbar.mean() ** 2):.6f}, <z> = {zbar.mean():.6f};"
          f" the mean rate over all cells |x><y|, populations included, {2 * r.ham.mean() / r.n:.4f} gamma N")


def part2(r, S):
    print(f"\n== Part 2, N = {r.n}: the tower ==")
    rot = np.array([r.idx[((s << 1) | (s >> (r.n - 1))) & ((1 << r.n) - 1)] for s in r.basis])
    for k in sorted(S):
        e, v, w, m = S[k]
        zj = (v ** 2) @ r.Z
        Tv = np.zeros_like(v)
        Tv[rot] = v                                    # translation by one site
        print(f"  S_{k:+d}: E = {e: .4f}, Z2 weight {w:.4f}, eigenspace dimension {m}, <Z_j> = {zj.mean():.4f}"
              f" (spread over j {zj.max() - zj.min():.0e}), <S|T|S> = {v @ Tv:+.6f}")
    print(f"  Z2 weight on the tower: {sum(S[k][2] for k in S):.4f}")


def part3(r, S):
    print(f"\n== Part 3, N = {r.n}: first order (route 1), rates in units of gamma N ==")
    res = {}
    ks = sorted(k for k in S if k + 1 in S)
    print("  layer 1, |S_(k+1)><S_k|:")
    for k in ks:
        (ea, va, _, _), (eb, vb, _, _) = S[k + 1], S[k]
        rates, share, size = block(r, va, vb, ea - eb)
        top = np.argsort(share)[::-1][:2]
        diag = -np.sum(np.outer(va, vb) * r.dissipator(np.outer(va, vb))) / r.n
        print(f"    k = {k:+d}: omega = {ea - eb:.4f}, block {size}, diagonal {diag:.4f}, coherence's largest shares "
              + ", ".join(f"{share[t]:.3f} at {rates[t]:.4f}" for t in top))
        res[("l1", k)] = (ea - eb, rates, share)
    print("  odd layers l = 2m + 1, the straddling U-partners |S_(m+1)><S_-m| and |S_m><S_-(m+1)|:")
    for m in range(0, r.n // 2):
        if m + 1 not in S or -(m + 1) not in S:
            break
        X1 = np.outer(S[m + 1][1], S[-m][1])
        X2 = np.outer(S[m][1], S[-(m + 1)][1])
        partner = np.sum(r.U(X1) * X2)
        modes, diag = pair2x2(r, X1, X2)
        omega = S[m + 1][0] - S[-m][0]
        rates, share, size = block(r, S[m + 1][1], S[-m][1], omega)
        top = np.argsort(share)[::-1][:2]
        print(f"    l = {2 * m + 1}: omega = {omega:.4f}, <U(X1), X2> = {partner:+.6f}, diagonal {diag:.4f}; 2 x 2 modes "
              + ", ".join(f"{rt:.4f} (U {u:+.3f})" for rt, u in modes)
              + f"; block {size}, largest shares " + ", ".join(f"{share[t]:.3f} at {rates[t]:.4f}" for t in top))
        res[("odd", m)] = (omega, modes, rates, share)
    # every coherence of the tower: |S_a><S_b| (a > b) and its U-partner |S_-b><S_-a| share the frequency
    # exactly; the first-order coupling of the two is +-sum_l <S_-b|Z_l|S_a>^2, a sum of squares, read here
    # against the distance |a + b| of the two levels it connects (a + b = 0: its own partner)
    print("  every coherence |S_a><S_b|, a > b, a + b >= 0, against its U-partner |S_-b><S_-a|:")
    print("    layer a-b, distance a+b: half the 2 x 2 splitting (gamma N), coupling sum_l <S_-b|Z_l|S_a>^2 / N")
    rows, bydist = {}, {}
    for a in S:
        for b in S:
            if a <= b or a + b < 0:
                continue
            Sa, Sb, Sma, Smb = S[a][1], S[b][1], S[-a][1], S[-b][1]
            coup = sum((Smb @ (r.Z[:, l] * Sa)) ** 2 for l in range(r.n)) / r.n
            if a + b == 0:
                X1 = np.outer(Sa, Sb)
                u = np.sum(r.U(X1) * X1)
                rows.setdefault(a - b, []).append(f"({a:+d},{b:+d}) self, U {u:+.3f}")
            else:
                modes, diag = pair2x2(r, np.outer(Sa, Sb), np.outer(Smb, Sma))
                half = abs(modes[1][0] - modes[0][0]) / 2
                assert modes[0][1] > 0 or half < 1e-9         # the U-even mode is the slower
                rows.setdefault(a - b, []).append(f"({a:+d},{b:+d}) d={a + b}: {half:.4f} [{coup:.4f}]")
                bydist.setdefault(min(a + b, 4), []).append((half, a, b))
    for l in sorted(rows):
        print(f"    l = {l:2d}: " + "; ".join(rows[l]))
    print("    largest half splitting by distance: " + "; ".join(
        f"d {'>= 4' if d == 4 else d}: {max(bydist[d])[0]:.2e} at ({max(bydist[d])[1]:+d},{max(bydist[d])[2]:+d})"
        for d in sorted(bydist)))
    # the pairs at distance 1 and 2 through the middle, in their full blocks
    for a, b in ((1, 0), (2, 0)):
        modes, size = block_modes(r, S[a][1], S[b][1], S[a][0] - S[b][0])
        print(f"    full block of |S_{a:+d}><S_{b:+d}| ({size} pairs of eigenvectors): " + ", ".join(
            f"{rt:.4f} (share {sh:.3f}, U {u:+.3f})" for rt, sh, u in modes))
    omega2 = S[2][0] - S[0][0]
    rates2, share2, _ = block(r, S[2][1], S[0][1], omega2)
    res[("l2", 0)] = (omega2, rates2, share2)
    return res


def part4_dense(r, S, res, gammas=(1e-5, 1e-4)):
    print(f"\n== Part 4, N = {r.n}: route 2, the exact L (dense) at gamma = {', '.join(str(g) for g in gammas)} ==")
    print("  the real parts against the block's, in units of gamma N: largest difference, and that over (gamma N)^2")
    devs = {}
    for g in gammas:
        w = np.linalg.eigvals(r.liouvillian(g).toarray())
        for key, val in res.items():
            omega, rates = val[0], (val[1] if key[0] in ("l1", "l2") else val[2])
            ex = np.sort(-w[np.abs(w.imag + omega) < 1e-4].real / g / r.n)
            r1 = np.sort(rates)
            devs.setdefault(key, []).append((np.abs(ex - r1).max() if len(ex) == len(r1) else np.inf, len(ex), len(r1)))
    for key, dl in devs.items():
        omega = res[key][0]
        print(f"  {key[0]} {key[1]:+d}: omega = {omega:.4f}, block {dl[0][2]}, exact {dl[0][1]}, {dl[1][1]}: "
              + "; ".join(f"gamma = {g}: {d:.2e}, / (gamma N)^2 = {d / (g * r.n) ** 2:.3f}" for g, (d, _, _) in zip(gammas, dl)))
    worst = max(dl[0][0] for dl in devs.values())
    print(f"  largest difference over all blocks at gamma = {gammas[0]}: {worst:.2e} gamma N")


def part4_sparse(r, S, res, gamma=1e-5):
    print(f"\n== Part 4, N = {r.n}: route 2 for the middle partners, shift-invert at gamma = {gamma} ==")
    omega, modes, rates, share = res[("odd", 0)]
    L = r.liouvillian(gamma).tocsc()
    for rt, _ in modes:
        rb = rates[np.argmin(np.abs(rates - rt))]
        vals = spla.eigs(L, k=6, sigma=-1j * omega - gamma * r.n * rb, return_eigenvectors=False)
        best = vals[np.argmin(np.abs(vals - (-1j * omega - gamma * r.n * rb)))]
        print(f"  first order (full block) {rb:.4f}, in the 2 x 2 {rt:.4f}: nearest exact eigenvalue"
              f" {-best.real / gamma / r.n:.4f} gamma N at Im {best.imag:.6f}")


def part5(r, S, gN_list=(0.0, 0.0005, 0.002, 0.005, 0.05), window=1.5, num=601):
    print(f"\n== Part 5, N = {r.n}: the revival from Z2, read on the tower coherences ==")
    d = r.d
    A_op = (S[1][1], S[0][1])
    B_op = (S[0][1], S[-1][1])
    O_op = (S[2][1], S[1][1])
    modes, diag = pair2x2(r, np.outer(*A_op), np.outer(*B_op))
    oX = np.outer(*O_op)
    rO = -np.sum(oX * r.dissipator(oX)) / r.n
    (rs, us), (rf, uf) = modes
    print(f"  first order in the 2 x 2: middle partners diagonal {diag:.4f}, U-even {[rs, rf][[us, uf].index(max(us, uf))]:.4f},"
          f" U-odd {[rs, rf][[us, uf].index(min(us, uf))]:.4f}; next coherence diagonal {rO:.4f}")
    # the first-order prediction of the full block: rho(0) = |Z2><Z2| on the block's pairs, evolved by
    # exp(gamma M t), read on <S_a|.|S_b>; and the shares of the U-even vector A + B on the block's modes
    def first_order(Sa, Sb, omega):
        pa, pb, w, W = block_matrix(r, omega)
        c0 = r.V[r.z2, pa] * r.V[r.z2, pb]
        read = (r.V[:, pa].T @ Sa) * (r.V[:, pb].T @ Sb)
        amp = (W.T @ c0) * (W.T @ read)
        return -w / r.n, amp, W, pa, pb
    rA, ampA, WA, paA, pbA = first_order(*A_op, S[1][0] - S[0][0])
    rO_, ampO, _, _, _ = first_order(*O_op, S[2][0] - S[1][0])
    ev = (r.V[:, paA].T @ A_op[0]) * (r.V[:, pbA].T @ A_op[1]) + (r.V[:, paA].T @ B_op[0]) * (r.V[:, pbA].T @ B_op[1])
    sh = (WA.T @ ev) ** 2
    sh /= sh.sum()
    top = np.argsort(sh)[::-1][:3]
    print("  the U-even operator |S_1><S_0| + |S_0><S_-1| on the middle block's modes (the full block): "
          + ", ".join(f"{sh[t]:.4f} at {rA[t]:.4f}" for t in top))
    rho0 = np.zeros(d * d, dtype=complex)
    rho0[r.z2 * d + r.z2] = 1
    stag = np.array([(-1) ** j for j in range(r.n)])
    imb_x = (r.Z @ stag) / r.n
    curves = {}
    for gN in gN_list:
        t0 = time.time()
        T = window / gN if gN > 0 else 300.0
        R = spla.expm_multiply(r.liouvillian(gN / r.n), rho0, start=0, stop=T, num=num, endpoint=True)
        t = np.linspace(0, T, num)
        rho = R.reshape(num, d, d)
        A = np.einsum("i,tij,j->t", A_op[0], rho, A_op[1])
        B = np.einsum("i,tij,j->t", B_op[0], rho, B_op[1])
        O = np.einsum("i,tij,j->t", O_op[0], rho, O_op[1])
        imb = np.einsum("tii,i->t", rho, imb_x).real
        sel = np.abs(B - A).max() / np.abs(A).max()
        line = f"  gamma N = {gN}: max|B - A| / max|A| = {sel:.0e}"
        if gN > 0:
            kA = -np.polyfit(t, np.log(np.abs(A)), 1)[0] / gN
            kO = -np.polyfit(t, np.log(np.abs(O)), 1)[0] / gN
            predA = np.abs(np.exp(-np.outer(gN * t, rA)) @ ampA) / abs(ampA.sum())
            predO = np.abs(np.exp(-np.outer(gN * t, rO_)) @ ampO) / abs(ampO.sum())
            devA = np.abs(np.abs(A) / abs(A[0]) - predA).max()
            devO = np.abs(np.abs(O) / abs(O[0]) - predO).max()
            line += (f"; fitted slope of -ln|A| {kA:.4f} gamma N, of -ln|O| {kO:.4f} gamma N;"
                     f" |A|/|A(0)| at gamma N t = {window}: {abs(A[-1] / A[0]):.4f} (block {predA[-1]:.4f});"
                     f" largest |exact - block| over the window: A {devA:.2e} ({devA / gN:.2f} gamma N), O {devO:.2e}")
        else:
            line += f"; |A(T)|/|A(0)| = {abs(A[-1] / A[0]):.10f}, |O(T)|/|O(0)| = {abs(O[-1] / O[0]):.10f} at T = {T:g}"
        print(line + f"  [{time.time() - t0:.0f} s]")
        curves[gN] = (t, A, O, imb)
    rU = [rs, rf][[us, uf].index(max(us, uf))]
    print(f"  e^(-{rU:.4f} * {window}) = {np.exp(-rU * window):.4f}; e^(-{diag:.4f} * {window}) = {np.exp(-diag * window):.4f}")
    # the distance from first order at the smallest rate: the beat of ln|A| about its fitted line, against the
    # detuning of the neighbouring layer-1 coherence |S_2><S_1| from the middle one
    g0 = min(g for g in curves if g > 0)
    t, A = curves[g0][0], curves[g0][1]
    la = np.log(np.abs(A))
    resid = la - np.polyval(np.polyfit(t, la, 1), t)
    f = np.fft.rfftfreq(len(t), t[1] - t[0]) * 2 * np.pi
    beat = f[1 + np.argmax(np.abs(np.fft.rfft(resid))[1:])]
    print(f"  at gamma N = {g0}: ln|A| beats about its fitted line at angular frequency {beat:.4f} (resolution {f[1]:.4f});"
          f" the detuning of |S_2><S_1| from |S_1><S_0| is {(S[1][0] - S[0][0]) - (S[2][0] - S[1][0]):.4f}")
    np.savez(os.path.join(OUT, "revival_N10.npz"), **{f"{k}_{g}": v for g, c in curves.items() for k, v in zip("tAOI", c)})
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(1, 2, figsize=(11, 4.3))
        cols = {0.0005: "#1f4e9c", 0.002: "#3a8ee6", 0.005: "#e69a3a", 0.05: "#d62020"}
        for gN, c in cols.items():
            if gN not in curves:
                continue
            t, A, O, _ = curves[gN]
            ax[0].semilogy(gN * t, np.abs(A) / abs(A[0]), c=c, lw=1.6, label=f"γN = {gN}")
            ax[1].semilogy(gN * t, np.abs(O) / abs(O[0]), c=c, lw=1.6, label=f"γN = {gN}")
        xx = np.linspace(0, window, 60)
        rU = [rs, rf][[us, uf].index(max(us, uf))]
        rD = rA[np.argmin(np.abs(rA - [rs, rf][[us, uf].index(min(us, uf))]))]
        ax[0].semilogy(xx, np.abs(np.exp(-np.outer(xx, rA)) @ ampA) / abs(ampA.sum()), "k--", lw=1, label="first order, full block")
        ax[0].semilogy(xx, np.exp(-diag * xx), "k:", lw=1, label=f"e^(−{diag:.3f}·γNt), its own diagonal")
        ax[0].semilogy(xx, np.exp(-rD * xx), "k-.", lw=1, label=f"e^(−{rD:.3f}·γNt), U-odd mode of the block")
        ax[1].semilogy(xx, np.abs(np.exp(-np.outer(xx, rO_)) @ ampO) / abs(ampO.sum()), "k--", lw=1, label="first order, full block")
        for a_, ttl in [(ax[0], f"⟨S₁|ρ(t)|S₀⟩ from Z₂, ring of {r.n}"), (ax[1], f"⟨S₂|ρ(t)|S₁⟩ from Z₂, ring of {r.n}")]:
            a_.set_xlim(0, window)
            a_.set_ylim(0.15, 1.05)
            a_.set_xlabel("γN·t")
            a_.set_title(ttl, fontsize=10)
            a_.legend(fontsize=7.5)
        ax[0].set_ylabel("|amplitude| / |amplitude at t = 0|")
        plt.tight_layout()
        plt.savefig(os.path.join(OUT, "revival_N10.png"), dpi=120)
    except ImportError:
        print("  (matplotlib missing, no figure)")
    return curves


def part6(r, S, res, gamma=0.0004):
    print(f"\n== Part 6, N = {r.n}: the full spectrum at gamma = {gamma}, in Ma et al.'s frame ==")
    d = r.d
    t0 = time.time()
    L = r.liouvillian(gamma).toarray()
    w, Vr = np.linalg.eig(L)
    del L
    Vr /= np.linalg.norm(Vr, axis=0)
    ov = np.abs(Vr[r.z2 * d + r.z2, :])
    rot = [r.idx[((s << 1) | (s >> (r.n - 1))) & ((1 << r.n) - 1)] for s in r.basis]
    T = np.array(rot)
    perm = (T[:, None] * d + T[None, :]).reshape(-1)
    tv = np.einsum("ik,ik->k", Vr.conj(), Vr[np.argsort(perm), :]).real
    del Vr
    shift = 2 * r.ham.mean()
    x = w.real / gamma + shift
    print(f"  {len(w)} eigenvalues in {time.time() - t0:.0f} s; frame Re(lambda')/gamma = Re(lambda)/gamma + {shift:.4f}")
    # each mode of the straddling partners: the exact eigenvalue nearest its first-order position i*omega - gamma*N*r
    # (the conjugate of the partners' -i*omega, the spectrum being closed under conjugation), r the block's
    # eigenvalue nearest the 2 x 2 mode
    marks = []
    for m in range(0, r.n // 2):
        if ("odd", m) not in res:
            break
        omega, modes, rates, share = res[("odd", m)]
        out = []
        for rt, u in modes:
            rb = rates[np.argmin(np.abs(rates - rt))]
            target = 1j * omega - gamma * r.n * rb
            k = int(np.argmin(np.abs(w - target)))
            marks.append((k, u))
            out.append(f"U {u:+.0f}: x = {x[k]:.3f} (first order {r.n * (1 - rb) - (r.n - shift):.3f}), overlap with rho_Z2 {ov[k]:.1e},"
                       f" q = {'pi' if tv[k] < -0.9 else ('0' if tv[k] > 0.9 else 'other')}, |Im - omega| = {abs(w[k].imag - omega):.0e}")
        print(f"  l = {2 * m + 1}, Im = {omega:.4f}: " + "; ".join(out))
    # a symmetry T with T L' T^-1 = -L', L' = L + shift*gamma, linear or antilinear, fixes its centre by the trace
    # (the zero of the frame) and carries the stationary eigenvalue lambda = 0 (x = +shift) onto lambda = -2*shift*gamma
    # (x = -shift) with the same multiplicity: the necessary half of the two-end count
    stat = int(np.sum(np.abs(w) < 1e-9))
    near = np.min(np.abs(w - (-2 * shift * gamma))) / gamma
    print(f"  two ends: {stat} stationary eigenvalue(s) at x = +{shift:.4f}; the leftmost eigenvalue at x = {x.min():.3f};"
          f" the nearest eigenvalue to lambda = -2*{shift:.4f}*gamma (x = -{shift:.4f}) lies {near:.3f} gamma away")
    q0n, qpin = int(np.sum(tv > 0.9)), int(np.sum(tv < -0.9))
    print(f"  drawn: {q0n} eigenvalues at q = 0 and {qpin} at q = pi, {q0n + qpin} of {len(w)}")
    np.savez(os.path.join(OUT, "spectrum_N10_g0.0004.npz"), w=w, ov=ov, tv=tv, shift=shift)
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(5.6, 5.2))
        q0, qpi = tv > 0.9, tv < -0.9
        ax.scatter(x[q0], w.imag[q0], s=5, c="#f5c542", label="q = 0")
        ax.scatter(x[qpi], w.imag[qpi], s=5, c="#bbbbbb", label="q = π")
        for k, u in marks:
            for sgn in (1, -1):
                ax.scatter(x[k], sgn * w.imag[k], s=36, facecolors="#3a8ee6", edgecolors="k" if u > 0 else "#d62020",
                           linewidths=1.3, zorder=3)
        ax.scatter([], [], s=36, facecolors="#3a8ee6", edgecolors="k", label="straddling partners, U-even mode")
        ax.scatter([], [], s=36, facecolors="#3a8ee6", edgecolors="#d62020", label="straddling partners, U-odd mode")
        ax.set_xlim(-6, 6)
        ax.set_ylim(-13, 13)
        ax.set_xlabel("Re(λ′/γ)")
        ax.set_ylabel("Im(λ)")
        ax.set_title(f"dephased PXP, ring of {r.n}, γ = {gamma}", fontsize=10)
        ax.legend(fontsize=7.5, loc="lower left")
        plt.tight_layout()
        plt.savefig(os.path.join(OUT, "spectrum_N10_g0.0004.png"), dpi=120)
    except ImportError:
        print("  (matplotlib missing, no figure)")


def main():
    os.makedirs(OUT, exist_ok=True)
    args = set(sys.argv[1:])
    print("PXP's scars under Z dephasing (simulations/pxp_scars_under_dephasing.py)")
    for n in (8, 10, 12, 14):
        r = Ring(n)
        S = tower(r)
        part1(r)
        part2(r, S)
        res = part3(r, S)
        if n == 8:
            part4_dense(r, S, res)
        if n == 10:
            if "--route2-n10" in args:
                part4_sparse(r, S, res)
            part5(r, S)
            if "--fig4" in args:
                part6(r, S, res)


if __name__ == "__main__":
    main()
