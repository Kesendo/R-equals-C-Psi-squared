"""Verify the Q<->Delta handover UNIFICATION (the xxz_axis_handover arc conjecture).

Claim: the Delta-handover (XXZ anisotropy axis, experiments/XXZ_AXIS_BANDEDGE_TO_LEBENSADER.md)
is the SAME band-edge-floor crossing as the Q-handover (simulations/carbon/handover_q.py,
HandoverFloorClaim). Both = the interior/Lebensader dressed mode's darkness <n_XY> crossing 1
(= the Absorption-Theorem band-edge floor 2*gamma). Three falsifiable predictions:

  (1) THE FLOOR is exact and model-independent: the bright band-edge single-magnon coherence
      |vac><magnon| sits at rate = 2*gamma (darkness <n_XY> = 1) for ALL Delta - it is an
      eigenoperator of [H_XXZ, .] (lambda = -2g + i*E_k; the ZZ shifts only the frequency E_k,
      not Re). So even though F50's 2N-degeneracy COUNT breaks for Delta!=1, the floor persists.
  (2) Delta* = where the Lebensader rate crosses that floor: at Delta*, the slowest mode's
      darkness = 1 exactly (the Lebensader rate = 2*gamma). Same condition as the Q-handover.
  (3) It is a LEVEL CROSSING (the ring family), not an EP: at Delta* the frozen Lebensader
      (Im = 0) meets the oscillating band edge (Im >> 0) - two distinct branches crossing in Re,
      not a coalescence (the chain Q-handover was the EP).

Also checks the Absorption Theorem holds in XXZ (rate = 2*gamma * Sum_k k*weight_k, both modes).
Convention matches xxz_axis_bandedge_lebensader.py: J=1, gamma=0.05 (Q=20), open chain, Z-dephasing.
"""
import sys
import itertools
import warnings

sys.path.insert(0, "simulations")
import numpy as np
import framework as fw

warnings.filterwarnings("ignore")
GAMMA = 0.05
FLOOR = 2.0 * GAMMA  # the band-edge Absorption-Theorem rate (darkness <n_XY> = 1)
I = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.array([[1, 0], [0, -1]], dtype=complex)
PAULI = {"I": I, "X": X, "Y": Y, "Z": Z}


def site_op(M, i, N):
    mats = [I] * N
    mats[i] = M
    out = mats[0]
    for m in mats[1:]:
        out = np.kron(out, m)
    return out


def H_xxz(N, Delta, J=1.0):
    H = np.zeros((2 ** N, 2 ** N), dtype=complex)
    for i in range(N - 1):
        H += J * (site_op(X, i, N) @ site_op(X, i + 1, N) + site_op(Y, i, N) @ site_op(Y, i + 1, N))
        H += Delta * (site_op(Z, i, N) @ site_op(Z, i + 1, N))
    return H


def pauli_basis(N):
    d = 2 ** N
    mats = np.empty((4 ** N, d, d), dtype=complex)
    nxy = np.empty(4 ** N, dtype=int)
    for idx, letters in enumerate(itertools.product("IXYZ", repeat=N)):
        S = PAULI[letters[0]]
        for l in letters[1:]:
            S = np.kron(S, PAULI[l])
        mats[idx] = S
        nxy[idx] = sum(1 for c in letters if c in "XY")
    return mats, nxy


def branches(N, Delta, basis_mats, basis_nxy, gamma=GAMMA):
    """The slowest LEBENSADER (n_XY=0-dominated) and slowest BAND-EDGE (else) modes.
    Each: dict(rate, im, darkness=rate/2g, nxy_mean). darkness and nxy_mean are two readings of
    <n_XY> (rate/2g vs Sum_k k*weight_k) - the Absorption Theorem says they're equal."""
    d = 2 ** N
    L = fw.lindbladian_z_dephasing(H_xxz(N, Delta), [gamma] * N)
    ev, evec = np.linalg.eig(L)
    rate = -ev.real
    out = {"leb": None, "band": None}
    nz = np.where(rate > 1e-9)[0]
    for s in nz[np.argsort(rate[nz])]:
        Op = evec[:, s].reshape(d, d)
        coeffs = np.einsum("sij,ij->s", basis_mats.conj(), Op) / d
        w = np.abs(coeffs) ** 2
        by = {k: float(w[basis_nxy == k].sum()) for k in range(N + 1)}
        tot = sum(by.values())
        by = {k: v / tot for k, v in by.items()}
        kind = "leb" if by.get(0, 0.0) > by.get(1, 0.0) else "band"
        if out[kind] is None:
            out[kind] = dict(rate=float(rate[s]), im=float(abs(ev[s].imag)),
                             darkness=float(rate[s] / (2 * gamma)),
                             nxy_mean=float(sum(k * v for k, v in by.items())))
        if out["leb"] and out["band"]:
            break
    return out


def delta_star_rate(N, basis_mats, basis_nxy, lo=0.5, hi=2.5):
    """Delta* by the RATE crossing: where the Lebensader rate = the band-edge floor 2*gamma
    (i.e. where the Lebensader darkness crosses 1). Bisection on (leb_rate - FLOOR)."""
    for _ in range(40):
        mid = 0.5 * (lo + hi)
        b = branches(N, mid, basis_mats, basis_nxy)
        leb = b["leb"]
        # below Delta*: leb is faster than the floor (rate > 2g); above: slower (rate < 2g).
        if leb is None or leb["rate"] > FLOOR:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


if __name__ == "__main__":
    for N in (4, 5):
        basis_mats, basis_nxy = pauli_basis(N)
        print(f"=== N={N}  (J=1, gamma={GAMMA}, floor 2*gamma={FLOOR}) ===")
        dstar = delta_star_rate(N, basis_mats, basis_nxy)

        # (1) THE FLOOR: the band edge sits at exactly 2*gamma (darkness 1) for all Delta.
        floor_ok = True
        print("  Delta | band-edge (rate, darkness, |Im|) | Lebensader (rate, darkness, |Im|)")
        for Delta in [0.0, 1.0, dstar - 0.1, dstar, dstar + 0.1, 2.0]:
            b = branches(N, Delta, basis_mats, basis_nxy)
            be, le = b["band"], b["leb"]
            floor_ok &= be is not None and abs(be["darkness"] - 1.0) < 1e-9
            bs = f"rate={be['rate']:.5f} dark={be['darkness']:.6f} Im={be['im']:.3f}" if be else "—"
            ls = f"rate={le['rate']:.5f} dark={le['darkness']:.6f} Im={le['im']:.2e}" if le else "—"
            tag = "  <- Delta*" if abs(Delta - dstar) < 1e-6 else ""
            print(f"  {Delta:5.3f} | {bs:42s} | {ls}{tag}")

        # (2) at Delta*, the slowest mode's darkness = 1 (the Lebensader rate crosses the floor).
        at = branches(N, dstar, basis_mats, basis_nxy)
        leb_at, band_at = at["leb"], at["band"]
        # (3) level crossing: Lebensader frozen (Im~0), band edge oscillating (Im>>0).
        print(f"\n  [1] FLOOR exact (band edge darkness = 1 for all Delta):        {floor_ok}")
        print(f"  [2] Delta* (rate crossing, leb darkness->1) = {dstar:.4f}  "
              f"(doc: N=5~1.525, N=4~1.618=phi)")
        print(f"      at Delta*: Lebensader darkness = {leb_at['darkness']:.5f}  (predicted 1.000)")
        print(f"  [3] level crossing: Lebensader |Im|={leb_at['im']:.2e} (frozen), "
              f"band-edge |Im|={band_at['im']:.3f} (oscillating)")

        assert floor_ok, f"N={N}: the band-edge floor is NOT exactly 2*gamma for all Delta"
        assert abs(leb_at["darkness"] - 1.0) < 5e-3, \
            f"N={N}: at Delta* the slowest darkness {leb_at['darkness']} != 1 (the unification fails)"
        assert leb_at["im"] < 1e-6, f"N={N}: Lebensader not frozen at Delta* (Im={leb_at['im']}) - would be an EP"
        assert band_at["im"] > 0.1, f"N={N}: band edge not oscillating at Delta* (Im={band_at['im']})"
        # Absorption Theorem in XXZ: rate/2g == Sum_k k*weight_k for both modes.
        assert abs(leb_at["darkness"] - leb_at["nxy_mean"]) < 1e-6, "Absorption Theorem fails for the Lebensader"
        assert abs(band_at["darkness"] - band_at["nxy_mean"]) < 1e-6, "Absorption Theorem fails for the band edge"
        print(f"  [AT] Absorption Theorem holds in XXZ: rate/2g == Sum_k k*w_k "
              f"(leb {leb_at['nxy_mean']:.5f}, band {band_at['nxy_mean']:.5f})")
        print()

    print("All asserts passed: the Delta-handover IS the band-edge-floor crossing (darkness=1), a LEVEL")
    print("CROSSING (frozen Lebensader meets oscillating band edge) - the same condition as the Q-handover,")
    print("driven along the anisotropy axis. The unification holds.")
