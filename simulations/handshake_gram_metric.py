"""handshake_decoder: the location-metric, CORRECTED after the math+physics review convergence (2026-06-19).

Both lenses (mathematical-review + physics-first-review) converged: the earlier "abstract well-conditioned
vs painted ill-conditioned" contrast was an INDEXING ARTIFACT (k=1..N-1 silently drops the forbidden k=N
column and substitutes the strength k=1), and "near-null eigenvalue = ambiguity" was a CATEGORY ERROR
(a global rank property identified with a pairwise event). The corrected picture, three DISTINCT objects:

1. THE DERIVABLE NULL (identifiability). The location dictionary M[b,k]=<psi_k|V_b|psi_1>, k=2..N, has rank
   EXACTLY N-2 -- one exact-zero singular value -- already in the BARE couplings, BEFORE any painter, with
   ZERO Q dependence. That null IS the K-partner selection rule <psi_N|V_b|psi_1>=0 (KPartnerSelectionRule
   Claim, Tier1Derived, Pi-mirror-protected). The painted reading INHERITS this null; it does not create it.
2. THE SECTORS. The painted bond-Gram's eigenvectors split by bond-MIRROR PARITY, not magnitude:
   antisymmetric {dominant seesaw, the K-partner location null} vs symmetric {the closure/F123 strength
   direction, a dim symmetric mode}. The painter adds a GRADED (E1-Ek) tail on top of the one exact null.
3. PAIRWISE CONFUSABILITY (the DefectDecoder's actual AmbiguityFactor). This is the COSINE matrix of the
   normalized letters and its largest off-diagonal (edge-vs-interior, near-ANTI-collinear), NOT a Gram
   eigenvalue. cos ~ -1 = "+dJ here reads like -dJ there", the sign-location ambiguity.

Identifiability (rank of a noiseless dictionary) is NOT FI(Q) (statistical precision, Cramer-Rao, ~Q): at
Q->inf, FI->inf but the null stays exactly null. They are complementary diagnostics, not two channels of one.
The IsDeadEnd dead-end lives in the letter VALUES (R_k); the metric's null STRUCTURE is derivable (this file).
"""
import importlib.util
import sys

import numpy as np

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
sys.path.insert(0, "simulations")
sys.path.insert(0, "simulations/carbon")
from incompleteness_survivor import bonds


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


blk = _load("blk", "simulations/handshake_rk_block.py")
psi = blk.psi


def location_dict(N):
    """The location dictionary M[b,k]=<psi_k|V_b|psi_1>, bonds b=0..N-2, channels k=2..N (INCLUDING the
    forbidden k=N so the rank deficit is VISIBLE, not indexed away). Closed form, no painter, no Q."""
    return np.array([[psi(k, b, N) * psi(1, b + 1, N) + psi(k, b + 1, N) * psi(1, b, N)
                      for k in range(2, N + 1)] for b in range(N - 1)])


def parity(v):
    """+1 if the bond-vector is symmetric under b<->N-2-b, -1 if antisymmetric (the chain mirror)."""
    r = v[::-1]
    return float(np.dot(v, r) / max(np.dot(v, v), 1e-30))


def main():
    J, dJ, g = 1.0, 0.02, 0.05            # canonical DefectDecoder protocol: gamma=0.05 => Q=20
    print("=== handshake location-metric (corrected): the K-partner null, the sectors, the cosine ===\n",
          flush=True)

    # GATE 1: the DERIVABLE null is in the BARE couplings, rank N-2, Q-free, = the K-partner rule.
    print("GATE 1 -- the null is in the BARE couplings (k=2..N), rank N-2 = the K-partner rule:", flush=True)
    print(f"{'N':>3} {'rank(M)':>8} {'N-2':>4} {'min sv (the K-partner null)':>28} {'|c_N(b)| max':>13}", flush=True)
    g1 = True
    for N in range(3, 9):
        M = location_dict(N)
        sv = np.linalg.svd(M, compute_uv=False)
        rank = int(np.sum(sv > 1e-9 * sv[0]))
        cN = max(abs(psi(N, b, N) * psi(1, b + 1, N) + psi(N, b + 1, N) * psi(1, b, N)) for b in range(N - 1))
        ok = rank == N - 2
        g1 &= ok
        print(f"{N:>3} {rank:>8} {N-2:>4} {sv[-1]:>28.2e} {cN:>13.2e}  {'' if ok else 'CHECK'}", flush=True)
    print(f"  => {'PASS: the bare location couplings already carry the rank-(N-2) null (no painter, no Q)' if g1 else 'FAIL'}\n",
          flush=True)

    # GATE 2 + 3: the PAINTED letters -- sectors (parity) and pairwise confusability (cosine).
    for N in (4, 5):
        bnds = bonds(N, "chain")
        f, rel, _ = blk.paint_block(N, J, g, dJ, bnds)
        G = f @ f.T
        w, V = np.linalg.eigh(G)
        order = np.argsort(w)[::-1]
        w, V = w[order], V[:, order]
        print(f"--- N={N} painted bond-Gram (canonical Q=20) ---", flush=True)
        print(f"  spectrum (norm) = [{', '.join(f'{x/w[0]:.4f}' for x in w)}]", flush=True)
        print("  GATE 2 -- eigenvectors split by MIRROR PARITY (sym=closure/strength, antisym=seesaw+K-partner):",
              flush=True)
        with np.printoptions(precision=3, suppress=True):
            for i in range(len(w)):
                p = parity(V[:, i])
                lab = "ANTISYM" if p < -0.5 else ("SYM" if p > 0.5 else "mixed")
                role = ("seesaw" if i == 0 and p < -0.5 else
                        "K-partner null" if p < -0.5 else "closure/strength")
                print(f"    eig {w[i]/w[0]:>7.4f}  parity {p:+.2f} [{lab:>7}]  {role:>16}  vec={V[:, i]}", flush=True)
        # GATE 3: cosine matrix of the normalized letters; the worst off-diagonal = the AmbiguityFactor pair
        norm = f / np.linalg.norm(f, axis=1, keepdims=True)
        C = norm @ norm.T
        off = C - np.diag(np.diag(C))
        ij = np.unravel_index(np.argmax(np.abs(off)), off.shape)
        print(f"  GATE 3 -- pairwise confusability = the COSINE matrix (NOT an eigenvalue):", flush=True)
        print(f"    worst pair bonds {bnds[ij[0]]}<->{bnds[ij[1]]}: cos = {off[ij]:+.3f} "
              f"({'anti-collinear: +dJ here ~ -dJ there (sign-location ambiguity)' if off[ij] < 0 else 'collinear'})",
              flush=True)
        print(flush=True)

    print("CORRECTED STATEMENT: the location-metric's null is the K-partner selection rule (rank N-2, in the",
          flush=True)
    print("bare couplings, Q-free, Tier1Derived) -- the painter INHERITS it and adds a graded (E1-Ek) tail.",
          flush=True)
    print("Pairwise confusability is the cosine matrix, a different object. Identifiability != FI(Q). The",
          flush=True)
    print("dead-end (IsDeadEnd) is the letter VALUES (R_k); the metric's null STRUCTURE is derivable.", flush=True)
    return g1


if __name__ == "__main__":
    sys.exit(0 if main() else 1)
