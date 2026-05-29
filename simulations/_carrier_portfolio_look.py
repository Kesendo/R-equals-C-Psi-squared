"""_carrier_portfolio_look.py - point the carrier-vector portfolio at our own systems.

The C# tool (CarrierVectorPortfolio) holds a carrier vector + every mode's per-channel
difference-portfolio + the law. Here we just LOOK: aim it at a few Heisenberg chains and
read off the slow (memory) modes as portfolios, to see what the perspective surfaces when
we point it at ourselves rather than at Si:P.

Three systems, same machinery:
  A. uniform γ chain         - the popcount baseline (degenerate rungs).
  B. one good qubit, COUPLED - a slow qubit (small γ) Heisenberg-coupled to fast ones.
  C. one good qubit, DECOUPLED - the same slow qubit, left off the coupling graph.

The question the tool answers at a glance: where does the memory live, and does a good
qubit stay good once it is coupled?

Run: python simulations/_carrier_portfolio_look.py
"""
from __future__ import annotations

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

if sys.platform == "win32":
    try:
        sys.stdout.reconfigure(encoding="utf-8", errors="replace")
        sys.stderr.reconfigure(encoding="utf-8", errors="replace")
    except Exception:
        pass

import numpy as np

X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.array([[1, 0], [0, -1]], dtype=complex)
I2 = np.eye(2, dtype=complex)


def site_op(N, site, P):
    op = P if site == 0 else I2
    for s in range(1, N):
        op = np.kron(op, P if s == site else I2)
    return op


def heisenberg(N, bonds, J=1.0):
    d = 2 ** N
    H = np.zeros((d, d), dtype=complex)
    for (a, b) in bonds:
        for P in (X, Y, Z):
            H = H + J * (site_op(N, a, P) @ site_op(N, b, P))
    return H


def lindbladian(H, c_list, gammas):
    d = H.shape[0]
    Id = np.eye(d, dtype=complex)
    L = -1j * (np.kron(H, Id) - np.kron(Id, H.T))
    for c, g in zip(c_list, gammas):
        cdc = c.conj().T @ c
        L = L + g * (np.kron(c, c.conj()) - 0.5 * (np.kron(cdc, Id) + np.kron(Id, cdc.T)))
    return L


def look(label, N, bonds, gammas, n_show=4):
    d = 2 ** N
    Zc = [site_op(N, l, Z) for l in range(N)]
    Id2 = np.eye(d * d, dtype=complex)
    N_ops = [(Id2 - np.kron(Zc[l], Zc[l])) / 2.0 for l in range(N)]
    H = heisenberg(N, bonds)
    L = lindbladian(H, Zc, gammas)
    w, V = np.linalg.eig(L)

    print(f"\n{'='*84}\n{label}")
    print(f"  carrier γ = {gammas},  bonds = {bonds}")
    print(f"  {'rate −Re(λ)':>12}   portfolio ⟨Δ⟩ per qubit            memory read")
    order = np.argsort(w.real)[::-1]  # slowest (largest Re, least negative) first
    shown = 0
    for k in order:
        rate = -w[k].real
        if rate < 1e-9:
            continue  # skip the stored steady state (kernel)
        v = V[:, k]
        nrm = np.vdot(v, v).real
        port = [np.vdot(v, N_ops[l] @ v).real / nrm for l in range(N)]
        pstr = "  ".join(f"q{l} {100*port[l]:3.0f}%" for l in range(N))
        # where does it live: the channel carrying most weighted decay
        dom = int(np.argmax([gammas[l] * port[l] for l in range(N)]))
        slow = int(np.argmin(gammas))
        read = (f"pure q{slow} (protected)" if port[slow] > 0.99 and sum(port) - port[slow] < 0.01
                else f"leaks into q{dom}" if dom != slow and port[dom] > 0.01
                else f"on q{dom}")
        print(f"  {rate:>12.4f}   {pstr}     {read}")
        shown += 1
        if shown >= n_show:
            break


def main():
    print("=" * 84)
    print("POINTING THE CARRIER-VECTOR PORTFOLIO AT OUR OWN CHAINS")
    print("=" * 84)

    # A. uniform γ: the popcount baseline
    look("A. uniform γ Heisenberg chain (the popcount baseline)",
         N=3, bonds=[(0, 1), (1, 2)], gammas=[0.1, 0.1, 0.1])

    # B. one good qubit (q0, small γ) COUPLED into the chain
    look("B. one good qubit q0 (γ=0.01) COUPLED to fast q1,q2 (γ=1.0)",
         N=3, bonds=[(0, 1), (1, 2)], gammas=[0.01, 1.0, 1.0])

    # C. the SAME good qubit, left OFF the coupling graph (decoupled)
    look("C. the same good qubit q0, DECOUPLED (bond only 1-2)",
         N=3, bonds=[(1, 2)], gammas=[0.01, 1.0, 1.0])

    print(f"\n{'='*84}")
    print("""What the tool shows:
  A. Uniform γ: every slow mode's portfolio is a flat popcount mix; the rungs are
     degenerate, nothing distinguishes the qubits. The vector lies flat.
  B. A good qubit COUPLED is not a good qubit. Its would-be protected coherence picks up
     a fat slice of the fast neighbours' difference (q1/q2 %), and the rate jumps decades
     above 2γ_q0. The portfolio shows the leak directly.
  C. The same qubit DECOUPLED keeps a clean ~100%-q0 portfolio at rate 2γ_q0: protected.
  The engineering reading falls straight out of the portfolio: protection is not having a
  slow qubit, it is keeping its bra-ket agreement in the fast channels (⟨Δ_fast⟩ → 0). The
  same lesson Si:P teaches with hyperfine, now read off any chain at a glance.
""")


if __name__ == "__main__":
    main()
