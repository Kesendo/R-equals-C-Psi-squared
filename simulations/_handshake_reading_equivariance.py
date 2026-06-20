"""handshake M3 (2026-06-20, gate-first): the defect-reading map is spatial-reflection equivariant.

The location dictionary M[b,k] = <psi_k|V_b|psi_1> (modes k=2..N, carrier psi_1, V_b the bond
hopping) commutes with the geometric chain mirror R (i -> N-1-i):

    M[N-2-b, k] = (-1)^(k-1) M[b, k]          (exact, single-excitation algebra)

because R V_b R = V_{N-2-b}, R psi_1 = +psi_1, R psi_k = (-1)^(k-1) psi_k. The sign-location
confusability is then a closed-form parity-weighted mode sum:

    cos(b, N-2-b) = sum_k (-1)^(k-1) w_k,  w_k = M[b,k]^2 / ||M[b,.]||^2

NOTE (honest scope): this is the BARE matrix-element dictionary; the bare cosine is -0.33 at N=5,
NOT -0.97. The -0.97 is the PAINTED (propagated alpha-profile) instance (_handshake_gram_metric.py,
Q=20), where the readout concentrates location weight on the R-odd seesaw k=2. R-equivariance is
preserved by the painting (R commutes with the Liouvillian for reflection-symmetric couplings);
only the weights w_k change. This is the cross-language ECHO of the C# witness
(DefectReadingEquivarianceClaim), not the validation.

GATES (any firing IS the find):
  G1  R-parity: |psi_k(N-1-j) - (-1)^(k-1) psi_k(j)| ~ 0 for all k, j.
  G2  equivariance: max_{b,k} |M[N-2-b,k] - (-1)^(k-1) M[b,k]| ~ 0.
  G3  cosine formula: measured cos(b,N-2-b) == parity-weighted sum, exactly.
  G4  anti-collinear: every mirror-pair cosine < 0 (R-odd seesaw carries the location), N=4,5,6.
"""
import sys

import numpy as np

sys.stdout.reconfigure(encoding="utf-8", errors="replace")


def psi(N, k):
    # unit-normalized via the norm; equivalent to BondingMode's baked-in sqrt(2/(N+1))
    # prefactor (every tested quantity -- parity, equivariance, cosine -- is scale-invariant).
    v = np.array([np.sin(np.pi * k * (i + 1) / (N + 1)) for i in range(N)])
    return v / np.linalg.norm(v)


def Vb(N, b):
    M = np.zeros((N, N))
    M[b, b + 1] = 1.0
    M[b + 1, b] = 1.0
    return M


def main():
    print("=== handshake M3: defect-reading spatial-reflection equivariance (gate-first) ===\n", flush=True)
    g1 = g2 = g3 = g4 = True
    for N in (4, 5, 6):
        p1 = psi(N, 1)
        # dictionary M[b][kIdx], b=0..N-2, k=2..N
        Md = np.array([[psi(N, k) @ Vb(N, b) @ p1 for k in range(2, N + 1)] for b in range(N - 1)])
        # G1 R-parity
        rp = max(abs(psi(N, k)[::-1] - ((-1) ** (k - 1)) * psi(N, k)).max() for k in range(1, N + 1))
        # G2 equivariance
        eq = 0.0
        for b in range(N - 1):
            for ki, k in enumerate(range(2, N + 1)):
                eq = max(eq, abs(Md[N - 2 - b, ki] - ((-1) ** (k - 1)) * Md[b, ki]))
        # G3/G4 cosine
        cos_err = 0.0
        cosines = []
        for b in range(N - 1):
            bm = N - 2 - b
            if bm <= b:
                continue
            u, v = Md[b], Md[bm]
            cm = u @ v / (np.linalg.norm(u) * np.linalg.norm(v))
            w = u ** 2 / (u @ u)
            pred = sum(((-1) ** (k - 1)) * w[ki] for ki, k in enumerate(range(2, N + 1)))
            cos_err = max(cos_err, abs(cm - pred))
            cosines.append(round(cm, 4))
        print(f"N={N}: R-parity={rp:.1e}  equivariance={eq:.1e}  cos-formula-err={cos_err:.1e}  "
              f"mirror-pair cos={cosines}", flush=True)
        g1 &= rp < 1e-9
        g2 &= eq < 1e-9
        g3 &= cos_err < 1e-9
        g4 &= all(c < 0 for c in cosines)
    print(f"\nG1 R-parity [{'ok' if g1 else 'FIRED'}]  G2 equivariance [{'ok' if g2 else 'FIRED'}]  "
          f"G3 cosine formula [{'ok' if g3 else 'FIRED'}]  G4 anti-collinear [{'ok' if g4 else 'FIRED'}]",
          flush=True)
    ok = g1 and g2 and g3 and g4
    print("\nVERDICT:", "all gates green; spatial-reflection equivariance confirmed (bare, exact)."
          if ok else "a gate FIRED -- diagnose, do not loosen.", flush=True)
    return ok


if __name__ == "__main__":
    sys.exit(0 if main() else 1)
