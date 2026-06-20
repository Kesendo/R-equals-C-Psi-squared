"""GATE-FIRST Stage 1: is F124's object guard the DefectDecoder's missing observable?

Stage 0 showed F124's full-matrix conditioning kappa(M)~N^2 is the bond-recovery inverse
problem. But the live DefectDecoder operates on the LOCATION dictionary k=2..N (the f-profiles
span the PT slow-mode weights), which is EXACTLY F124's object-guard object: rank N-2, the
K-partner k=N a null column, lambda_min=0. So the decoder inherits a rank deficiency, and its
known ambiguity (residual ratio ~1.5 at N=5, edge bond 3 weakened vs interior bond 1 strengthened)
is that deficiency, NOT kappa(M) (sqrt(kappa(5))=2.30 != 1.5, the naive gate correctly fires).

F124 proves the STRENGTH column k=1 (the carrier-to-carrier response M[b,1]=2 c_a c_{a+1}, which a
bond defect DOES produce) lifts the floor from 0 to E. PREDICTION (F124 as the decoder's fix):
re-including the strength channel resolves the sign-location ambiguity the location-only decoder
cannot. GATE: location-only -> the two mirror-pair signatures are anti-collinear (cos ~ -1,
ambiguity); with the strength column -> they separate (cos away from -1, residual ratio >> 1).
If the strength column does NOT resolve it, the gate fires.
"""
import sys
import numpy as np

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
np.set_printoptions(precision=4, suppress=True, linewidth=130)


def carrier(N, k=1):
    j = np.arange(N)
    return np.sqrt(2.0 / (N + 1)) * np.sin(np.pi * k * (j + 1) / (N + 1))


def build_M(N):
    """M[b,k] = <psi_k|V_b|psi_1>, bonds b=0..N-2, modes k=1..N (column 0 = strength channel k=1)."""
    c = carrier(N, 1)
    P = np.array([carrier(N, k) for k in range(1, N + 1)])
    M = np.zeros((N - 1, N))
    for a in range(N - 1):
        for ki in range(N):
            M[a, ki] = P[ki][a] * c[a + 1] + P[ki][a + 1] * c[a]
    return M


def decode_residual_ratio(dictionary, truth_bond, truth_dj, rng=None, noise=0.0):
    """Least-squares decode mirroring DefectDecoder: signature r = truth_dj * row(truth_bond);
    for each candidate bond, best-fit dj and residual ||r - dj*row(b)||; return (winner, runner-up,
    ratio = runnerup_residual / winner_residual). A ratio near 1 = ambiguous."""
    rows = dictionary
    r = truth_dj * rows[truth_bond].copy()
    if noise > 0 and rng is not None:
        r = r + rng.normal(0, noise, size=r.shape)
    resids = []
    for b in range(rows.shape[0]):
        f = rows[b]
        dj = (f @ r) / (f @ f)            # best-fit strength
        resids.append(np.linalg.norm(r - dj * f))
    order = np.argsort(resids)
    winner, runner = order[0], order[1]
    ratio = resids[runner] / resids[winner] if resids[winner] > 1e-15 else np.inf
    return winner, runner, ratio, resids


def main():
    GATE = {"fired": []}

    def gate(name, cond, detail=""):
        flag = "ok " if cond else "GATE-FIRE"
        if not cond:
            GATE["fired"].append(name)
        print(f"   [{flag}] {name}" + (f"   {detail}" if detail else ""))

    print("=" * 96)
    print("STAGE 1 - is F124's strength column the fix for the decoder's sign-location ambiguity?")
    print("=" * 96)

    # --- The strength channel separates edge from interior bonds by magnitude ---
    print("\nThe strength channel M[b,1] = 2 c_a c_{a+1} (carrier-to-carrier response) per bond, N=5:")
    M5 = build_M(5)
    for b in range(4):
        print(f"   bond {b}=({b},{b+1}): M[b,1] = {M5[b,0]:.4f}  ({'edge' if b in (0,3) else 'interior'})")
    edge = abs(M5[0, 0] - M5[3, 0]) < 1e-9 and abs(M5[1, 0] - M5[2, 0]) < 1e-9
    sep = abs(M5[0, 0] - M5[1, 0]) > 0.1
    gate("strength channel separates edge {0,3} from interior {1,2} by magnitude",
         edge and sep, f"edge={M5[0,0]:.4f}, interior={M5[1,0]:.4f}")

    # --- The mirror-pair anti-collinearity: location-only vs with-strength ---
    print("\nThe N=5 mirror pair: bond 3 weakened (dJ=-0.02) vs bond 1 strengthened (dJ=+0.02)")
    M = build_M(5)
    loc = M[:, 1:]                 # location dictionary k=2..N (drop strength column k=1)
    full = M                       # full dictionary k=1..N (with strength)
    sig3 = -0.02 * np.array([row for row in M])[3]    # truth: bond 3 weakened
    sig1 = +0.02 * np.array([row for row in M])[1]    # confuser: bond 1 strengthened
    def cosang(u, v, cols):
        u, v = u[cols], v[cols]
        return (u @ v) / (np.linalg.norm(u) * np.linalg.norm(v))
    cos_loc = cosang(sig3, sig1, slice(1, None))      # location modes k=2..N
    cos_full = cosang(sig3, sig1, slice(0, None))     # full k=1..N
    print(f"   cos(bond3-weakened, bond1-strengthened):  location-only = {cos_loc:+.4f},  with strength = {cos_full:+.4f}")
    gate("location-only: the two signatures are nearly ANTI-collinear (cos < -0.9, the ambiguity)",
         cos_loc < -0.9, f"cos_loc={cos_loc:+.4f}")
    gate("with strength column: the signatures SEPARATE (cos moves well away from -1)",
         cos_full > cos_loc + 0.2, f"cos_full={cos_full:+.4f} vs cos_loc={cos_loc:+.4f}")

    # --- The decode residual ratio: location-only ~1.5 (ambiguous) vs full >> 1 (resolved) ---
    print("\nDecode residual ratio (runner-up / winner; ~1 = ambiguous, >>1 = clean), truth = bond 3 weakened:")
    rng = np.random.default_rng(20260620)
    for N in (4, 5, 6, 7):
        M = build_M(N)
        loc = M[:, 1:]
        full = M
        truth = N - 2                 # the edge bond (3 at N=5); use the last interior-adjacent edge
        # weaken the edge bond
        w_loc, r_loc, ratio_loc, _ = decode_residual_ratio(loc, truth, -0.02)
        w_full, r_full, ratio_full, _ = decode_residual_ratio(full, truth, -0.02)
        print(f"   N={N}: location-only winner=bond{w_loc} ratio={ratio_loc:.3f} {'(AMBIGUOUS)' if ratio_loc < 3 else ''}"
              f"  |  with strength winner=bond{w_full} ratio={'inf' if np.isinf(ratio_full) else f'{ratio_full:.3f}'}"
              f" {'(RESOLVED, picks truth)' if w_full == truth else '(still wrong!)'}")
    # gate on N=5 specifically (the documented 1.5 case)
    M = build_M(5)
    w_loc, r_loc, ratio_loc, _ = decode_residual_ratio(M[:, 1:], 3, -0.02)
    w_full, r_full, ratio_full, _ = decode_residual_ratio(M, 3, -0.02)
    gate("N=5 location-only is ambiguous (winner != truth bond 3, OR ratio < 3)",
         w_loc != 3 or ratio_loc < 3.0, f"winner=bond{w_loc}, ratio={ratio_loc:.3f}")
    gate("N=5 WITH strength column resolves it (winner = truth bond 3)",
         w_full == 3, f"winner=bond{w_full}")

    print("\n" + "=" * 96)
    if GATE["fired"]:
        print(f"STAGE 1: {len(GATE['fired'])} GATE(S) FIRED -> {GATE['fired']}")
        print("Diagnose, do NOT loosen: is the strength channel NOT the fix, or is our setup off?")
    else:
        print("STAGE 1: ALL GATES PASS -> F124's object guard IS the decoder's missing observable.")
        print("The location-only decoder (k=2..N) inherits the K-partner rank deficiency: edge-weakened")
        print("and interior-strengthened bonds are anti-collinear (the ~1.5 ambiguity). The STRENGTH")
        print("channel k=1 (the carrier's own response, which F124 proves lifts lambda_min from 0 to E)")
        print("separates edge from interior by magnitude and RESOLVES the sign-location ambiguity.")
        print("F124 hands the live handshake_decoder its fix: read the carrier strength channel.")
    print("=" * 96)


if __name__ == "__main__":
    main()
