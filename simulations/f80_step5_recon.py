"""F80 Step 5 reconnaissance: make visible how the residual

    M = Π·L_H·Π⁻¹ + L_H

acquires the structure −2i·(H ⊗ I_bra) for the chain Π²-odd 2-body Hamiltonian.

Tom 2026-05-22: at knotty junctures, get the data and look, do not armchair it.

N (argv, default 3), H = Σ_l X_l Y_{l+1} (chain, Π²-odd 2-body). Π is built from the
verified PiOperator.cs definition (Z-dephasing, per site: I→X ×1, X→I ×1,
Z→Y ×i, Y→Z ×i) and self-checked. Outputs: a self-check report and the
H-eigenbasis structure of Π, L_H, M to stdout, plus heatmaps under
simulations/results/f80_step5_recon/.
"""
import os
import sys
import itertools
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# Windows consoles default to cp1252; force UTF-8 so Σ/Π/σ print.
sys.stdout.reconfigure(encoding="utf-8")
sys.stderr.reconfigure(encoding="utf-8")
np.set_printoptions(precision=4, suppress=True, linewidth=220)

N = int(sys.argv[1]) if len(sys.argv) > 1 else 3
d = 2 ** N
dd = d * d

# Pauli matrices, indexed I=0, X=1, Z=2, Y=3 (PiOperator.cs enum order)
I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.array([[1, 0], [0, -1]], dtype=complex)
PAULI = [I2, X, Z, Y]


def kron_all(mats):
    out = np.array([[1.0 + 0j]])
    for m in mats:
        out = np.kron(out, m)
    return out


def site(op, l):
    return kron_all([op if i == l else I2 for i in range(N)])


# H = Σ_l X_l Y_{l+1}: chain Π²-odd 2-body Hamiltonian
H = np.zeros((d, d), dtype=complex)
for l in range(N - 1):
    H = H + site(X, l) @ site(Y, l + 1)
assert np.allclose(H, H.conj().T), "H must be Hermitian"

# L_H = −i[H, ·] as a dd×dd superoperator on vec(ρ) with vec(ρ)[i*d+j] = ρ[i,j];
# vec(Hρ) = (H⊗I)vec(ρ), vec(ρH) = (I⊗Hᵀ)vec(ρ).
LH = -1j * (np.kron(H, np.eye(d)) - np.kron(np.eye(d), H.T))

# Π as a dd×dd superoperator, from PiOperator.cs (Z-dephasing).
PI_LETTER = {0: (1, 1.0 + 0j), 1: (0, 1.0 + 0j), 2: (3, 1.0j), 3: (2, 1.0j)}
Pi = np.zeros((dd, dd), dtype=complex)
for letters in itertools.product(range(4), repeat=N):
    sigma = kron_all([PAULI[L] for L in letters])
    new = [PI_LETTER[L][0] for L in letters]
    phase = np.prod([PI_LETTER[L][1] for L in letters])
    sigma_new = kron_all([PAULI[L] for L in new])
    # Pauli strings are HS-orthogonal with Tr(σ†σ) = d, so the dual of vec(σ) is vec(σ)/d.
    Pi += phase * np.outer(sigma_new.reshape(-1), sigma.reshape(-1).conj()) / d

# M = Π·L_H·Π⁻¹ + L_H
M = Pi @ LH @ np.linalg.inv(Pi) + LH

# Self-checks: if any fails, the recon below is not trustworthy.
chk_order4 = np.allclose(np.linalg.matrix_power(Pi, 4), np.eye(dd))
XN = kron_all([X] * N)
XN_conj = np.kron(XN, XN.T)            # vec(XN ρ XN) = (XN⊗XNᵀ)vec(ρ)
chk_pi2 = np.allclose(Pi @ Pi, XN_conj)
H_eigs = np.linalg.eigvalsh(H)
M_eigs = np.linalg.eigvals(M)


def multiset(z, tol=6):
    return sorted(np.round(z, tol).tolist(), key=lambda c: (c.real, c.imag))


pred = multiset(np.repeat(2j * H_eigs, d))   # 2i·E_a, each with bra multiplicity d
chk_spec = np.allclose(multiset(M_eigs), pred, atol=1e-6)

print("=== F80 Step 5 reconnaissance, N=" + str(N) + ", H = Σ X_l Y_{l+1} ===\n")
print("self-checks (must all be True):")
print(f"  Π is order 4 (Π⁴ = I) ............ {chk_order4}")
print(f"  Π² = X⊗N conjugation ............ {chk_pi2}")
print(f"  Spec(M) = 2i·Spec(H), mult ×2^N .. {chk_spec}")
if not (chk_order4 and chk_pi2 and chk_spec):
    print("\n  A self-check FAILED -- recon data below is NOT trustworthy.")
print()

# Change to the H-eigen-operator basis σ_(a,b) = |a⟩⟨b|.
E, V = np.linalg.eigh(H)               # H|a⟩ = E[a]|a⟩, |a⟩ = V[:,a], E ascending
U = np.kron(V, V.conj())               # column a*d+b = vec(|a⟩⟨b|) = V[:,a]⊗conj(V[:,b])
LH_e = U.conj().T @ LH @ U
Pi_e = U.conj().T @ Pi @ U
M_e = U.conj().T @ M @ U


def offdiag_norm(A):
    return np.linalg.norm(A - np.diag(np.diag(A)))


# Norm outside the d×d ket-energy blocks (block a = fixed ket index a).
ketblock_mask = np.ones((dd, dd), dtype=bool)
for a in range(d):
    ketblock_mask[a * d:(a + 1) * d, a * d:(a + 1) * d] = False


def ketblock_offnorm(A):
    return np.linalg.norm(A[ketblock_mask])


print("H many-body spectrum E_a:", np.round(E, 4), "\n")
print("L_H in the σ_(a,b) basis:")
print(f"  off-diagonal norm ............... {offdiag_norm(LH_e):.2e}   (expect ~0; L_H diagonal as −i(E_a−E_b))\n")
print("M in the σ_(a,b) basis  (index = a·d + b, a = ket, b = bra):")
print(f"  off-diagonal norm ............... {offdiag_norm(M_e):.2e}")
print(f"  norm outside ket-energy blocks .. {ketblock_offnorm(M_e):.2e}   (if ~0: M block-diagonal in the ket index)")
print(f"  Frobenius norm ‖M‖_F ............ {np.linalg.norm(M_e):.4f}")
Mdiag = np.diag(M_e).reshape(d, d)     # row a = ket, col b = bra
print("  Im part of diag(M) reshaped [ket a, bra b]:")
print(Mdiag.imag.round(3))
print("  compare each row to  −2·E_a =", np.round(-2 * E, 3), "\n")
print("Π in the σ_(a,b) basis:")
print(f"  off-diagonal norm ............... {offdiag_norm(Pi_e):.2e}")
print(f"  norm outside ket-energy blocks .. {ketblock_offnorm(Pi_e):.2e}")

# Tom's "Blindenschrift" read of |Π|: is Π_e a signed permutation of the σ_(a,b)
# basis?  The data said no -- and the recon below says why.  H = Σ X_l Y_{l+1} has a
# degenerate spectrum, so np.linalg.eigh picks an ARBITRARY orthonormal basis inside
# each H-eigenspace.  The σ_(a,b) basis, hence the fine structure of Π_e, is gauge;
# only sector-level facts survive.  The checks below split gauge noise from structure.
tol = 1e-6
nnz_col = (np.abs(Pi_e) > tol).sum(axis=0)
print("\nΠ in the σ_(a,b) basis -- is it a signed permutation?  (Tom's Blindenschrift)")
print(f"  nonzeros per column: min {nnz_col.min()}, max {nnz_col.max()}")
print(f"  Π_e is a signed permutation of σ_(a,b): "
      f"{bool(nnz_col.min() == 1 and nnz_col.max() == 1)}")

# Gauge test: re-orient the arbitrary basis inside every degenerate H-eigenspace.
rng = np.random.default_rng(0)
uniqE, grp = np.unique(np.round(E, 6), return_inverse=True)
Vg = V.copy()
for g in range(len(uniqE)):
    cols = np.where(grp == g)[0]
    if len(cols) > 1:
        A = (rng.standard_normal((len(cols),) * 2)
             + 1j * rng.standard_normal((len(cols),) * 2))
        Vg[:, cols] = Vg[:, cols] @ np.linalg.qr(A)[0]
Ug = np.kron(Vg, Vg.conj())
Pi_g = Ug.conj().T @ Pi @ Ug
nnz_g = (np.abs(Pi_g) > tol).sum(axis=0)
M_g = Ug.conj().T @ M @ Ug
gauge_moved = (int(nnz_col.min()), int(nnz_col.max())) != (int(nnz_g.min()), int(nnz_g.max()))
print(f"  H-eigenspace degeneracies ......... {[int((grp == g).sum()) for g in range(len(uniqE))]}")
print(f"  |Π| nonzeros/col, random basis .... min {nnz_g.min()}, max {nnz_g.max()}")
print(f"  Π_e count changed under regauge ... {gauge_moved}"
      "   (True => Π_e fine structure is gauge)")
print(f"  M off-diagonal norm, random basis . {offdiag_norm(M_g):.2e}"
      "   (~0 => M diagonal stays diagonal, gauge-invariant)")

# The real Step-5 object: Π·L_H·Π⁻¹.  M = L_H + Π·L_H·Π⁻¹ is diagonal iff this is.
a_idx, b_idx = np.arange(dd) // d, np.arange(dd) % d
diff_val = np.round(E[a_idx] - E[b_idx], 6)        # L_H eigenvalue on σ_(a,b) is -i*diff
sum_val = np.round(E[a_idx] + E[b_idx], 6)         # {H,.} eigenvalue on σ_(a,b) is -i*sum
D = Pi_e @ LH_e @ np.linalg.inv(Pi_e)
off_blk = np.sqrt(sum(np.linalg.norm(Pi_e[~np.isclose(sum_val, diff_val[c]), c]) ** 2
                      for c in range(dd)))
print("\nΠ·L_H·Π⁻¹  (the real Step-5 object, replacing the dead permutation guess):")
print(f"  off-diagonal norm ................. {offdiag_norm(D):.2e}   (diagonal => M diagonal)")
print(f"  diag(Π·L_H·Π⁻¹) = -i*(E_a+E_b) .... {bool(np.allclose(np.diag(D), -1j * sum_val))}"
      "   (so Π·[H,.]·Π⁻¹ = {H,.})")
print(f"  Π maps diff-sector onto sum-sector, off-block norm .. {off_blk:.2e}")
print("  the discrete 'Blindenschrift' is this sector map; the per-block analysis")
print("  below separates what is gauge from what is real inside each block.")

# Tom's sharper read: each "Viereck" of the heatmap carries its own pattern, and at
# N=5 the patterns persist -- only the scale shifts.  Regauging is Π_e -> W†·Π_e·W,
# W block-diagonal over the (ε_ket, ε_bra) sectors, so every sector block B obeys
# B -> W'†·B·W: its entries (the visible pattern) move, its singular values do not.
# The singular values are therefore the gauge-invariant fingerprint of each block.
sec = grp[a_idx] * len(uniqE) + grp[b_idx]      # (ε_ket, ε_bra) sector of each σ index
sec_ids = sorted(set(sec.tolist()))
sec_dim = int((sec == sec_ids[0]).sum())


def block_fingerprints(P):
    fp = {}
    for s in sec_ids:
        cols = np.where(sec == s)[0]
        for t in sec_ids:
            rows = np.where(sec == t)[0]
            B = P[np.ix_(rows, cols)]
            if np.linalg.norm(B) > 1e-6:
                fp[(s, t)] = np.linalg.svd(B, compute_uv=False)
    return fp


def svtag(v):
    vals, cnts = np.unique(np.round(v, 3), return_counts=True)
    return "  ".join(f"{c}x{val:.3f}" for val, c in zip(vals, cnts))


fp_eigh = block_fingerprints(Pi_e)
fp_rand = block_fingerprints(Pi_g)
gauge_inv = (set(fp_eigh) == set(fp_rand)
             and all(np.allclose(fp_eigh[k], fp_rand[k]) for k in fp_eigh))
groups = {}
for v in fp_eigh.values():
    groups[svtag(v)] = groups.get(svtag(v), 0) + 1
print("\nper-block fingerprint (Tom: each Viereck its own pattern):")
print(f"  (ε_ket,ε_bra) sectors .......... {len(sec_ids)}, each {sec_dim}-dim"
      "   (sector dim = the 'Skala' that shifts with N)")
print(f"  nonzero sector-blocks .......... {len(fp_eigh)}")
print(f"  distinct singular-value spectra  {len(groups)}   (>1 => the Vierecke really differ)")
print(f"  fingerprints gauge-invariant ... {gauge_inv}   (same spectra in eigh and random basis)")
for tag, cnt in sorted(groups.items(), key=lambda kv: -kv[1]):
    print(f"    [{cnt:2d} block(s)]  SV = {tag}")

# The nonzero blocks form a permutation π of the sectors (forced: each is a full
# unitary and Π is globally unitary).  π must convert sum-value into diff-value.
K = len(uniqE)
pi_perm = {s: t for (s, t) in fp_eigh}
is_bijection = (sorted(pi_perm) == sec_ids and sorted(pi_perm.values()) == sec_ids)
sum_to_diff = all(np.isclose(uniqE[t // K] + uniqE[t % K], uniqE[s // K] - uniqE[s % K])
                  for s, t in pi_perm.items())
print(f"  the blocks form a sector permutation π: bijection = {is_bijection}")
print(f"  π converts sum to diff  (sum(π(s)) = diff(s)) ...... {sum_to_diff}")

# --- Step 5 proof verification (2026-05-22). ---
# Π·[H,·]·Π⁻¹ = {H,·} reduces to a per-site Pauli computation. For each bond
# bond_l = X_l Y_{l+1} and every Pauli string Q,  Π(bond_l·Q) = −σ(l,Q)·bond_l·Π(Q),
# σ = +1/−1 as bond_l commutes/anticommutes with Q. On the anticommuting bonds (those
# [H,·] keeps) this is Π(bond_l·Q) = bond_l·Π(Q); Π also flips commute<->anticommute,
# so Π[H,Q] = {H,Π Q}. The per-site core is μ(P·a) = ε_P·c_P(a)·P·μ(a) with
# ε_X = ε_Z = +1, ε_Y = −1, c = commute sign; (I)/(II)/(III) checked directly below.
Pi1 = np.zeros((4, 4), dtype=complex)
for _k in range(4):
    _nk, _ph = PI_LETTER[_k]
    Pi1 += _ph * np.outer(PAULI[_nk].reshape(-1), PAULI[_k].reshape(-1).conj()) / 2


def mu1(A):
    return (Pi1 @ A.reshape(-1)).reshape(2, 2)


def apply_pi(O):
    return (Pi @ O.reshape(-1)).reshape(d, d)


def csign(P, A):
    return 1 if np.allclose(P @ A, A @ P) else -1


id_I = all(np.allclose(mu1(X @ a), csign(X, a) * X @ mu1(a)) for a in PAULI)
id_II = all(np.allclose(mu1(Y @ a), -csign(Y, a) * Y @ mu1(a)) for a in PAULI)
bonds = [site(X, l) @ site(Y, l + 1) for l in range(N - 1)]
bond_lemma = True
flip_ok = True
for _letters in itertools.product(range(4), repeat=N):
    Q = kron_all([PAULI[L] for L in _letters])
    PiQ = apply_pi(Q)
    for bond in bonds:
        if not np.allclose(apply_pi(bond @ Q), -csign(bond, Q) * bond @ PiQ):
            bond_lemma = False
        if csign(bond, PiQ) != -csign(bond, Q):
            flip_ok = False
id_III = all(np.allclose(mu1(Z @ a), csign(Z, a) * Z @ mu1(a)) for a in PAULI)
eye = np.eye(d)
Pi_inv = np.linalg.inv(Pi)
comm = np.kron(H, eye) - np.kron(eye, H.T)
anti = np.kron(H, eye) + np.kron(eye, H.T)
pi_comm = np.allclose(Pi @ comm @ Pi_inv, anti)
print("\nStep 5 proof verification (per-site identities + bond lemma):")
print(f"  (I)   μ(X·a) = c_X(a)·X·μ(a) ................. {id_I}")
print(f"  (II)  μ(Y·a) = −c_Y(a)·Y·μ(a) ................ {id_II}")
print(f"  (III) μ(Z·a) = c_Z(a)·Z·μ(a) ................. {id_III}")
print(f"  bond lemma Π(bond·Q) = −σ·bond·Π(Q), all {4 ** N} Q  {bond_lemma}")
print(f"  Π flips every bond relation, σ(l,ΠQ) = −σ(l,Q) ... {flip_ok}")
print(f"  => Π·[H,·]·Π⁻¹ = {{H,·}} for H = Σ X_l Y_(l+1) . {pi_comm}")
print("  all four Π²-odd pairs, sign s where Π·[H,·]·Π⁻¹ = s·{H,·}:")
for pn, P, qn, Qm in [("X", X, "Y", Y), ("X", X, "Z", Z),
                      ("Y", Y, "X", X), ("Z", Z, "X", X)]:
    Hp = sum(site(P, l) @ site(Qm, l + 1) for l in range(N - 1))
    cp = np.kron(Hp, eye) - np.kron(eye, Hp.T)
    ap = np.kron(Hp, eye) + np.kron(eye, Hp.T)
    conj = Pi @ cp @ Pi_inv
    if np.allclose(conj, ap):
        s = "+1   M = -2i(H⊗I_bra)"
    elif np.allclose(conj, -ap):
        s = "-1   M = +2i(I_ket⊗Hᵀ)"
    else:
        s = "??   not ±{H,·}"
    print(f"    ({pn},{qn}): s = {s}")

# Heatmaps.
outdir = "simulations/results/f80_step5_recon"
os.makedirs(outdir, exist_ok=True)
fig, axes = plt.subplots(1, 3, figsize=(19, 6))
for ax, A, t in [(axes[0], LH_e, "L_H"), (axes[1], Pi_e, "Π"), (axes[2], M_e, "M")]:
    im = ax.imshow(np.abs(A), cmap="magma")
    ax.set_title(f"|{t}|  in σ_(a,b) basis (index a·d+b)")
    for g in range(d, dd, d):
        ax.axhline(g - 0.5, color="cyan", lw=0.6)
        ax.axvline(g - 0.5, color="cyan", lw=0.6)
    fig.colorbar(im, ax=ax, fraction=0.046)
fig.suptitle("F80 Step 5 recon, N=" + str(N) + ": |L_H|, |Π|, |M| in the H-eigen-operator basis "
             "(cyan = ket-energy blocks)")
plt.tight_layout()
path = os.path.join(outdir, f"heatmaps_N{N}.png")
plt.savefig(path, dpi=110)
print(f"\nheatmaps written to {path}")

# Sector-ordered |Π|: reindex σ so each (ε_ket,ε_bra) sector is contiguous, then the
# sector permutation is visible directly (cyan grid = sector boundaries).
order = np.argsort(sec, kind="stable")
fig2, ax2 = plt.subplots(figsize=(7.5, 6.5))
im2 = ax2.imshow(np.abs(Pi_e[np.ix_(order, order)]), cmap="magma")
ax2.set_title(f"|Π| N={N} reordered by (ε_ket,ε_bra) sector "
              f"({len(sec_ids)} sectors x {sec_dim}-dim, cyan = sector boundaries)")
for g in range(sec_dim, dd, sec_dim):
    ax2.axhline(g - 0.5, color="cyan", lw=0.7)
    ax2.axvline(g - 0.5, color="cyan", lw=0.7)
fig2.colorbar(im2, ax=ax2, fraction=0.046)
plt.tight_layout()
path2 = os.path.join(outdir, f"pi_sector_order_N{N}.png")
plt.savefig(path2, dpi=110)
print(f"sector-ordered |Π| written to {path2}")
