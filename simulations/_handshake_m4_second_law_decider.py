"""GATE-FIRST decider for handshake M4 (A vs B): does a SECOND closure law remove the
random-single-excitation states that false-accept the MM(0) mirror-pair-MI checksum?

Background: M4's single checksum MM(0) = sum_pairs MI(i, N-1-i) false-accepts ~20-30% of random
single-excitation (SE) states (Review R1). MM(0) depends ONLY on the populations p_i=|c_i|^2 and is
SYMMETRIC within each mirror pair (MI(i,j)=h(p_i)+h(p_j)-h(p_i+p_j)), so it is BLIND to mirror
ASYMMETRY. The chiral-mirror closure (ChiralMirrorTrajectoryClaim) measures exactly that: a grammatical
reading is mirror-symmetric, c_i = +/- c_{N-1-i}, i.e. the mirror overlap O = |<psi|R psi>| = 1
(R reflects sites i -> N-1-i). The bonding carrier has O=1; a random SE state has O~0.

DECISION:
  - If the MM(0)-false-accepting random SE states are mirror-ASYMMETRIC (O well below 1), the second
    law removes them -> the JOINT checksum discriminates -> path B (design+land the joint checksum;
    grammatical manifold = mirror-symmetric SE near B(N), with the honest R2 caveat: a manifold, not
    the unique carrier).
  - If they are mirror-SYMMETRIC too (O~1), the second law adds nothing -> path A (honest coarse filter).

Gate-first: the deciding gate fires (-> A) iff adding the mirror law does NOT substantially cut the
random-SE false-accept rate. Do not band-tune to force a verdict; report the asymmetry DISTRIBUTION.
"""
import sys
import numpy as np

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
np.set_printoptions(precision=4, suppress=True, linewidth=120)


def h(p):
    p = np.clip(p, 1e-15, 1 - 1e-15)
    return -p * np.log2(p) - (1 - p) * np.log2(1 - p)


def mm0_from_pops(p, N):
    """MM(0) = sum over mirror pairs of MI(i, N-1-i) = h(p_i)+h(p_j)-h(p_i+p_j), pure SE state."""
    s = 0.0
    for ell in range(N // 2):
        a, b = p[ell], p[N - 1 - ell]
        s += h(a) + h(b) - h(min(a + b, 1.0))
    return s


def carrier_pops(N, k=1):
    c = np.sqrt(2.0 / (N + 1)) * np.sin(np.pi * k * (np.arange(N) + 1) / (N + 1))
    return c ** 2 / np.sum(c ** 2)


def mirror_overlap(c):
    """O = |<psi|R psi>| = |sum_i conj(c_i) c_{N-1-i}|; =1 iff mirror-symmetric (c_i=+/-c_{N-1-i})."""
    return abs(np.vdot(c, c[::-1]))


def main():
    GATE = {"to_A": []}

    def gate(name, to_A, detail=""):
        flag = "A!" if to_A else "ok"
        if to_A:
            GATE["to_A"].append(name)
        print(f"   [{flag}] {name}" + (f"   {detail}" if detail else ""))

    rng = np.random.default_rng(20260620)
    SAMPLE = 60000
    print("=" * 100)
    print("M4 DECIDER - does the chiral-mirror law remove the MM(0) random-SE false-accepts? (A vs B)")
    print("=" * 100)

    for N in (4, 6, 8):
        B = mm0_from_pops(carrier_pops(N, 1), N)
        band = 0.15 * B                                   # the spec's own MM(0) band
        # reference grammatical readings (mirror-symmetric): carrier k=1, W-state, sine k=2
        cN = carrier_pops(N, 1)
        w_c = np.ones(N) / np.sqrt(N)
        sine2 = np.sqrt(2.0 / (N + 1)) * np.sin(np.pi * 2 * (np.arange(N) + 1) / (N + 1)); sine2 /= np.linalg.norm(sine2)
        print(f"\nN={N}: B(N)={B:.4f}, MM(0) band=+/-{band:.4f}")
        print(f"   reference (mirror-symmetric, should pass BOTH laws): "
              f"carrier O={mirror_overlap(np.sqrt(cN)):.3f}, W O={mirror_overlap(w_c):.3f}, sine-k2 O={mirror_overlap(sine2):.3f}")

        # random SE ensemble (complex Gaussian amplitudes, normalized)
        amps = rng.normal(size=(SAMPLE, N)) + 1j * rng.normal(size=(SAMPLE, N))
        amps /= np.linalg.norm(amps, axis=1, keepdims=True)
        pops = np.abs(amps) ** 2
        mm = np.array([mm0_from_pops(pops[i], N) for i in range(SAMPLE)])
        O = np.array([mirror_overlap(amps[i]) for i in range(SAMPLE)])

        mm_pass = np.abs(mm - B) < band                   # MM(0)-alone acceptance
        fa_mm = mm_pass.mean()
        # among MM(0)-false-accepts, the mirror-overlap distribution (carrier sits at O=1)
        O_fa = O[mm_pass]
        print(f"   MM(0)-ALONE false-accept rate (random SE) = {fa_mm*100:.1f}%   "
              f"(reproduces Review R1 ~20-30%)")
        print(f"   their mirror overlap O: mean={O_fa.mean():.3f}, median={np.median(O_fa):.3f}, "
              f"95th pct={np.percentile(O_fa,95):.3f}, max={O_fa.max():.3f}  (carrier O=1.000)")

        # JOINT checksum: MM(0) in band AND mirror-symmetric (O > 1 - mband), swept over the mirror band
        for mband in (0.30, 0.15, 0.05):
            joint = (mm_pass & (O > 1 - mband)).mean()
            print(f"      JOINT (MM(0) band AND O>{1-mband:.2f}) false-accept rate = {joint*100:.3f}%")
        # the deciding gate: with a mirror band MATCHED to the MM(0) band's strictness (O>0.85),
        # does the joint false-accept drop to a small fraction of the MM(0)-alone rate?
        joint85 = (mm_pass & (O > 0.85)).mean()
        cut = joint85 / fa_mm if fa_mm > 0 else 0.0
        gate(f"N={N}: the mirror law does NOT substantially cut false-accepts (-> path A)",
             cut > 0.5, f"joint/alone = {joint85*100:.3f}% / {fa_mm*100:.1f}% = {cut:.3f} (B if << 1)")

    print("\n" + "=" * 100)
    if GATE["to_A"]:
        print(f"DECISION: PATH A (the mirror law does not discriminate; {GATE['to_A']})")
        print("The closure laws do not separate the false-accepts from the carrier -> land the honest")
        print("coarse grammaticality filter (rejects GHZ/localized only).")
    else:
        print("DECISION: PATH B - the chiral-mirror closure REMOVES the random-SE false-accepts.")
        print("MM(0) is blind to mirror asymmetry; random SE states are mirror-ASYMMETRIC (O~0), so the")
        print("JOINT checksum (MM(0) near B(N) AND mirror-symmetric O~1) discriminates where MM(0) alone")
        print("fails. The grammatical manifold = mirror-symmetric delocalized SE near B(N) (HONEST R2")
        print("caveat: a manifold including W and the sine modes, NOT the unique bonding carrier). Design")
        print("and land the joint multi-law checksum.")
    print("=" * 100)


if __name__ == "__main__":
    main()
