"""Lift DEGENERACY_PALINDROME.md from chain to ALL topologies (its Open Question 2, systematic).

The document proves, for the HEISENBERG chain under Z-dephasing, that the Liouvillian degeneracy
sequence d(k) (eigenvalue multiplicity at grid position Re = -k*gamma_q, gamma_q = 2*gamma_phys = 0.1)
is palindromic d(k)=d(N-k) from Pi, with boundary closed forms d_real(0)=N+1, d_real(1)=2N. It is
CHAIN-ONLY. The C# rmt export uses the SAME isotropic Heisenberg H for every topology (Topology.cs:
"isotropic Heisenberg (X+Y+Z)"), so ring/star/complete are the same model -- only the bond graph
changes. This feeds those spectra through the same machinery.

Gate-first questions (an assert that CAN fire is the find):
  G1. Palindrome d(k)=d(N-k): must hold for EVERY topology (Pi is topology-independent). A failure
      would mean a data/model problem, not physics.
  G2. d_real(0) = N+1 universal? (the N+1 magnetization-sector conserved quantities; U(1) holds for
      any Heisenberg graph -> expected universal.)
  G3. d_real(1) = 2N universal? (the document's SWAP-invariance proof; SWAP-invariance is exact only
      on the complete graph -> this is the one most likely to BREAK off-chain. If it breaks, the 2N
      is chain/graph-specific, a real finding for Open Question 2.)
  Then: compare the inner d_real(k>=2) sequences across topologies (the document says these are
  topology-dependent; here is the systematic table).

Run:  python simulations/degeneracy_palindrome_topology.py
"""
from pathlib import Path
import numpy as np

RESULTS = Path(__file__).parent / "results"
GAMMA_Q = 0.1     # the grid quantum = 2 * gamma_phys (0.05); matches the chain verifier
TOL = 1e-8
TOPOS = ("chain", "ring", "star", "complete")


def load(topo, N):
    name = f"rmt_eigenvalues_N{N}.csv" if topo == "chain" else f"rmt_eigenvalues_{topo}_N{N}.csv"
    p = RESULTS / name
    if not p.exists():
        return None
    data = np.loadtxt(p, delimiter="\t", skiprows=1)
    return data[:, 0] + 1j * data[:, 1]


def grid_seq(re_vals, N):
    """Multiplicity at each grid position Re = -k*GAMMA_Q, k=0..N."""
    return [int(np.sum(np.abs(re_vals - (-k * GAMMA_Q)) < TOL)) for k in range(N + 1)]


def sequences(topo, N):
    ev = load(topo, N)
    if ev is None:
        return None, None
    real = ev[np.abs(ev.imag) < TOL].real
    return grid_seq(real, N), grid_seq(ev.real, N)   # d_real(k), d_total(k)


# ---- chain reference gate (port fidelity: must reproduce the document) ----
print("=" * 92)
print("REFERENCE GATE: chain d_real(k) must reproduce DEGENERACY_PALINDROME.md exactly")
print("=" * 92)
chain_expected = {2: [3, 4, 3], 3: [4, 6, 6, 4], 4: [5, 8, 14, 8, 5],
                  5: [6, 10, 14, 14, 10, 6], 6: [7, 12, 19, 16, 19, 12, 7]}
for N in range(3, 7):
    dr, _ = sequences("chain", N)
    ok = dr == chain_expected[N]
    print(f"  chain N={N}: {dr}  {'OK' if ok else 'MISMATCH vs ' + str(chain_expected[N])}")
    assert ok, f"REFERENCE GATE FIRED: chain N={N} {dr} != document {chain_expected[N]}"
print("REFERENCE GATE PASS: the machinery reproduces the chain document.\n")

# ---- the topology sweep ----
print("=" * 92)
print("TOPOLOGY SWEEP: d_real(k) and the gates G1 (palindrome), G2 (N+1), G3 (2N)")
print("=" * 92)
g1_all, g2_all, g3_all = True, True, True
inner = {t: {} for t in TOPOS}
for N in range(3, 7):
    print(f"\n  N={N}:")
    for topo in TOPOS:
        dr, dt = sequences(topo, N)
        if dr is None:
            print(f"    {topo:9} (no CSV)")
            continue
        pal = dr == dr[::-1]
        d0_ok = dr[0] == N + 1
        d1_ok = dr[1] == 2 * N
        g1_all &= pal
        g2_all &= d0_ok
        g3_all &= d1_ok
        inner[topo][N] = dr
        flags = f"pal={'Y' if pal else 'N'} d0={dr[0]}{'=N+1' if d0_ok else '!=N+1'} " \
                f"d1={dr[1]}{'=2N' if d1_ok else '!=2N'}"
        print(f"    {topo:9} d_real={dr}  [{flags}]")

print("\n" + "-" * 92)
print(f"  G1 palindrome d(k)=d(N-k), all topologies:  {'HOLDS' if g1_all else 'BROKEN'}")
print(f"  G2 d_real(0) = N+1, all topologies:         {'HOLDS' if g2_all else 'BROKEN'}")
print(f"  G3 d_real(1) = 2N,  all topologies:         {'HOLDS' if g3_all else 'BROKEN (chain/graph-specific!)'}")

# ---- inner-position table (Open Question 2: topology-dependent inner degeneracies) ----
print("\n" + "=" * 92)
print("INNER DEGENERACY d_real(2) across topologies (Open Question 2 of the document)")
print("=" * 92)
print(f"{'N':>2} | " + " | ".join(f"{t:>9}" for t in TOPOS))
for N in range(4, 7):
    cells = []
    for topo in TOPOS:
        seq = inner.get(topo, {}).get(N)
        cells.append(f"{seq[2]:>9}" if seq and len(seq) > 2 else f"{'-':>9}")
    print(f"{N:>2} | " + " | ".join(cells))
print("\n(The document's data point: d_real(2) at N=4 is Chain=14, Star=16, Complete=36. "
      "G1 must hold everywhere (Pi); G3 is the interesting one.)")
