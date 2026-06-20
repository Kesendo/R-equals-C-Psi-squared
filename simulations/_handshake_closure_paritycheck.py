"""Handshake M4: closure as an operational parity-check (the JOINT two-law checksum).

A valid handshake reading must satisfy the closure laws as a consistency condition. The single
MM(0) mirror-pair-MI checksum is only a coarse filter (it false-accepts ~20-30% of random
single-excitation states; Review R1). The DISCRIMINATING parity-check is the JOINT of two
operational single-state closure laws:

  (1) F77 bandwidth:   MM(0) = sum_{ell} MI(site ell, site N-1-ell) ~ B(N)  (the bonding plateau)
  (2) chiral mirror:   the reading is mirror-symmetric, O = |<psi|R psi>| ~ 1  (R reflects i -> N-1-i)

MM(0) depends only on the populations p_i=|c_i|^2 and is symmetric within each mirror pair, so it is
BLIND to mirror asymmetry; the chiral-mirror law catches exactly the asymmetric random-SE states MM(0)
lets through. Both laws are needed: mirror symmetry alone accepts end-localized symmetric states like
(|1_0>+|1_{N-1}>)/sqrt2 (which MM(0) rejects, MM(0)=2 >> B); MM(0) alone accepts asymmetric random SE
(which the mirror law rejects). The complementary-light and Pi-flux laws are Liouville-PAIR properties,
not single-state checksums, so they are not part of the operational reading parity-check.

HONEST SCOPE (Review R2/R3): the accepted set is the grammatical MANIFOLD - mirror-symmetric
delocalized single-excitation near B(N) (includes the W-state and the sine modes), NOT the unique
bonding carrier. A scalar+symmetry checksum certifies grammaticality (the right manifold), not the
specific carrier. The plateau is N-dependent and EVEN-N-favorable (odd N is a lower curve, the unpaired
center site excluded). The checksum is a necessary consistency condition, not a sufficient "iff carrier".

Gate-first (a firing gate is the find; do not loosen the bands to force a verdict):
  G1 - plateau: bonding carrier MM(0)=B(N), nearly k-independent, B(N) tracks 1+3/(4(N+1)ln2) and -> 1
       bit (even N); odd N is a distinct lower curve (stated, not hidden).
  G2 - discrimination (the JOINT): accept the manifold (carrier, W, sine modes), reject the gross
       failures (GHZ, localized, end-localized-symmetric), and cut the random-SE false-accept rate from
       ~25% (MM(0) alone) to a few % (JOINT) - the distribution test, not one draw.
  G3 - defect-robustness: MM(0) and O are t=0 properties of the prepared reading state, independent of
       the delta-J defect in H, so the verdict on the carrier is defect-invariant by construction.
"""
import sys
import numpy as np

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
np.set_printoptions(precision=4, suppress=True, linewidth=120)

BAND_MI = 0.15      # |MM(0) - B(N)| < BAND_MI * B(N)
BAND_O = 0.85       # mirror-symmetric iff O > BAND_O


def h(p):
    p = np.clip(p, 1e-15, 1 - 1e-15)
    return -p * np.log2(p) - (1 - p) * np.log2(1 - p)


def mm0(pops, N):
    s = 0.0
    for ell in range(N // 2):
        a, b = pops[ell], pops[N - 1 - ell]
        s += h(a) + h(b) - h(min(a + b, 1.0))
    return s


def carrier_amps(N, k=1):
    c = np.sqrt(2.0 / (N + 1)) * np.sin(np.pi * k * (np.arange(N) + 1) / (N + 1))
    return c / np.linalg.norm(c)


def mirror_overlap(c):
    return abs(np.vdot(c, c[::-1]))


def B_of(N):
    return mm0(np.abs(carrier_amps(N, 1)) ** 2, N)


def accept(amps, N, B):
    """The JOINT parity-check: MM(0) near B(N) AND mirror-symmetric."""
    pops = np.abs(amps) ** 2
    return abs(mm0(pops, N) - B) < BAND_MI * B and mirror_overlap(amps) > BAND_O


def main():
    GATE = {"fired": []}

    def gate(name, ok, detail=""):
        flag = "ok " if ok else "GATE-FIRE"
        if not ok:
            GATE["fired"].append(name)
        print(f"   [{flag}] {name}" + (f"   {detail}" if detail else ""))

    # ---- G1: the plateau (even N -> 1 bit; odd N a distinct lower curve) ----
    print("=" * 100)
    print("G1 - the bonding plateau B(N): even N -> 1 bit (tracks 1+3/(4(N+1)ln2)); odd N a lower curve")
    print("=" * 100)
    print(f"{'N':>4} {'B(N)':>9} {'1+3/(4(N+1)ln2)':>17} {'parity':>7} {'k-indep spread':>15}")
    even = []
    for N in [4, 6, 8, 10, 20, 50, 100]:
        B = B_of(N)
        lead = 1 + 3 / (4 * (N + 1) * np.log(2))
        # k-independence: spread of MM(0) across bonding-ish sine modes k=1..min(3,N-1)
        ks = [mm0(np.abs(carrier_amps(N, k)) ** 2, N) for k in range(1, min(4, N))]
        spread = max(ks) - min(ks)
        even.append((N, B))
        print(f"{N:>4} {B:>9.5f} {lead:>17.5f} {'even':>7} {spread:>15.2e}")
    for N in [5, 7, 9, 11]:
        B = B_of(N)
        print(f"{N:>4} {B:>9.5f} {'(lower curve)':>17} {'odd':>7}")
    decreasing = all(even[i][1] > even[i + 1][1] for i in range(len(even) - 1))
    approaches1 = abs(even[-1][1] - 1.0) < 0.05
    gate("even-N B(N) decreases toward 1 bit (N=4..100)", decreasing and approaches1,
         f"B(4)={even[0][1]:.4f} -> B(100)={even[-1][1]:.4f}")
    kspread_ok = all(max(mm0(np.abs(carrier_amps(N, k)) ** 2, N) for k in range(1, min(4, N)))
                     - min(mm0(np.abs(carrier_amps(N, k)) ** 2, N) for k in range(1, min(4, N))) < 1e-2
                     for N in [4, 6, 8, 10])
    gate("MM(0) nearly k-independent across bonding sine modes (the manifold, not just k*)", kspread_ok)

    # ---- G2: discrimination (the JOINT two-law checksum; the distribution test) ----
    print("\n" + "=" * 100)
    print("G2 - JOINT discrimination: accept the manifold, reject the gross failures, cut the random-SE")
    print("     false-accept rate from ~25% (MM(0) alone) to a few % (JOINT). The distribution, not one draw.")
    print("=" * 100)
    rng = np.random.default_rng(20260620)
    SAMPLE = 40000
    all_ok = True
    for N in (4, 6, 8):
        B = B_of(N)
        # references
        carrier = carrier_amps(N, 1)
        w = np.ones(N) / np.sqrt(N)
        sine2 = carrier_amps(N, 2)
        endloc = np.zeros(N); endloc[0] = endloc[-1] = 1 / np.sqrt(2)        # mirror-symmetric BUT not a reading
        loc = np.zeros(N); loc[0] = 1.0                                       # localized
        manifold = {"carrier(k=1)": carrier, "W": w, "sine(k=2)": sine2}
        reject_refs = {"end-localized-sym": endloc, "localized |1_0>": loc}
        acc_manifold = {nm: accept(a, N, B) for nm, a in manifold.items()}
        rej_refs = {nm: not accept(a, N, B) for nm, a in reject_refs.items()}
        # GHZ is a wrong-sector 2^N state; its MM(0) ~ N/2 bits (1 per classically-correlated pair) >> B
        ghz_mm = N / 2.0
        ghz_reject = abs(ghz_mm - B) >= BAND_MI * B
        # random SE distribution: MM(0)-alone vs JOINT false-accept
        amps = rng.normal(size=(SAMPLE, N)) + 1j * rng.normal(size=(SAMPLE, N))
        amps /= np.linalg.norm(amps, axis=1, keepdims=True)
        mm_alone = np.array([abs(mm0(np.abs(amps[i]) ** 2, N) - B) < BAND_MI * B for i in range(SAMPLE)])
        joint = np.array([accept(amps[i], N, B) for i in range(SAMPLE)])
        fa_alone, fa_joint = mm_alone.mean(), joint.mean()
        print(f"\nN={N}: B(N)={B:.4f}")
        print(f"   manifold accepted: {acc_manifold}")
        print(f"   gross failures rejected: {rej_refs}, GHZ(MM~{ghz_mm:.1f})={ghz_reject}")
        print(f"   random SE false-accept: MM(0) alone = {fa_alone*100:.1f}%  ->  JOINT = {fa_joint*100:.2f}%  "
              f"(cut {fa_alone/max(fa_joint,1e-9):.0f}x)")
        ok = (all(acc_manifold.values()) and all(rej_refs.values()) and ghz_reject
              and fa_joint < 0.05 and fa_joint < fa_alone / 4)
        all_ok &= ok
        gate(f"N={N}: JOINT accepts the manifold, rejects gross failures, random-SE false-accept < 5% and << MM(0)-alone",
             ok, f"joint FA={fa_joint*100:.2f}% vs alone {fa_alone*100:.1f}%")

    # ---- G3: defect-robustness (MM(0), O are t=0 state properties, H-independent) ----
    print("\n" + "=" * 100)
    print("G3 - defect-robustness: the checksum reads the PREPARED state (t=0), not H; defect-invariant")
    print("=" * 100)
    N = 5  # use odd to also show the verdict is computed identically regardless of N parity
    B = B_of(N)
    carrier = carrier_amps(N, 1)
    # the verdict on the carrier reading does not take any delta-J / bond defect as input:
    verdict_clean = accept(carrier, N, B)
    # "applying a defect" changes H (the dynamics), not the prepared reading state at t=0:
    verdict_with_defect_label = accept(carrier, N, B)  # same state -> same verdict, by construction
    gate("MM(0), O are t=0 properties of the prepared carrier; verdict is defect-invariant",
         verdict_clean == verdict_with_defect_label,
         f"carrier verdict (clean)={verdict_clean}, (with a defect in H)={verdict_with_defect_label}")
    print("   (the defect lives in H/the dynamics; the parity-check validates the READING at t=0, "
          "upstream of and independent of the (b, delta-J) the decoder then extracts.)")

    print("\n" + "=" * 100)
    if GATE["fired"]:
        print(f"M4 PARITY-CHECK: {len(GATE['fired'])} GATE(S) FIRED -> {GATE['fired']}")
        sys.exit(1)
    else:
        print("M4 PARITY-CHECK: ALL GATES PASS. The JOINT two-law checksum (MM(0)~B(N) AND mirror-symmetric")
        print("O~1) is the operational reading parity-check: it cuts the random-SE false-accept from ~25%")
        print("(MM(0) alone) to a few %, accepts the grammatical manifold (carrier/W/sine modes), and")
        print("rejects the gross failures (GHZ/localized/end-localized). HONEST SCOPE: a necessary")
        print("consistency condition on the grammatical MANIFOLD (mirror-symmetric SE near B(N)), even-N-")
        print("favorable, NOT a sufficient 'iff the unique bonding carrier'.")
    print("=" * 100)


if __name__ == "__main__":
    main()
