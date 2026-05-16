#!/usr/bin/env python3
"""Multi-lens data tour on |+−+⟩ N=3.

Locker take: pick one familiar setup, run it through every framework lens we
have, lay everything out side-by-side. No interpretation forced — just look.

Lenses (one column per anchor's reading):
  - F87 trichotomy           classify_pauli_pair → 'truly' / 'soft' / 'hard'
  - F49 Frobenius scaling    ‖M‖² (including γ_T1 contribution)
  - F81 Π-decomposition      |M_anti|²/|M|² ∈ [0, 1]
  - F80 Bloch sign-walk      Π²-odd 2-body M-spectrum (max |imag|; — if N/A)
  - polarity                  per-site Bloch axis (state-only — same for all H)
  - d=0 axis                 fraction of ρ on d=0 substrate (state-only)
  - cockpit Lebensader       cusp pattern + skeleton/trace rating
"""
from __future__ import annotations

import math
import sys
from pathlib import Path

import numpy as np

if sys.platform == "win32":
    try:
        sys.stdout.reconfigure(encoding="utf-8")
    except Exception:
        pass

SCRIPT_DIR = Path(__file__).parent
sys.path.insert(0, str(SCRIPT_DIR))
import framework as fw
from framework.lebensader import cockpit_panel as lebensader_cockpit_panel


def fmt_num(x, w=8, prec=4):
    if x is None:
        return f"{'—':>{w}s}"
    if isinstance(x, (int, np.integer)):
        return f"{x:>{w}d}"
    try:
        return f"{x:>{w}.{prec}f}"
    except Exception:
        return f"{str(x):>{w}s}"


def main():
    N = 3
    GAMMA_DEPH = 0.1
    GAMMA_T1 = 0.01
    J = 1.0
    bonds = [(i, i + 1) for i in range(N - 1)]
    chain = fw.ChainSystem(N=N, gamma_0=GAMMA_DEPH, J=J, topology='chain', H_type='xy')

    plus = np.array([1, 1], dtype=complex) / math.sqrt(2)
    minus = np.array([1, -1], dtype=complex) / math.sqrt(2)
    psi = np.kron(plus, np.kron(minus, plus))
    rho_0 = np.outer(psi, psi.conj())

    cases = [
        ('truly XX+YY', [('X', 'X'), ('Y', 'Y')]),
        ('XY+YX',       [('X', 'Y'), ('Y', 'X')]),
        ('IY+YI',       [('I', 'Y'), ('Y', 'I')]),
        ('YZ+ZY',       [('Y', 'Z'), ('Z', 'Y')]),
        ('XZ+ZX',       [('X', 'Z'), ('Z', 'X')]),
        ('XZ+XZ',       [('X', 'Z'), ('X', 'Z')]),
    ]

    # State-only lenses (computed once on ρ_0):
    pol = fw.polarity_diagnostic(rho_0, N)
    d0 = fw.d_zero_decomposition(rho_0, chain)
    state_pol_str = ','.join('+' if x > 0 else ('-' if x < 0 else '0')
                              for x in pol['polarity_axis'])
    state_pol_agg = pol['aggregate_polarity']
    state_on_axis = pol['on_axis_fraction']
    state_d0 = d0['d0_weight']
    state_d2 = d0['d2_norm']
    state_kerdim = d0['kernel_dimension']

    print(f"Multi-lens data tour")
    print(f"  N={N}, |+−+⟩, γ_deph={GAMMA_DEPH}, γ_T1={GAMMA_T1}, J={J}, topology=chain")
    print()
    print(f"State-only lenses (depend on ρ_0, not on H):")
    print(f"  polarity axis      = ({state_pol_str})   aggregate = {state_pol_agg:.4f}")
    print(f"  on-axis fraction   = {state_on_axis:.4f}   (1 = pure X-basis polarity state)")
    print(f"  d=0 weight         = {state_d0:.4f}   d2 norm = {state_d2:.4f}   kernel dim = {state_kerdim}")
    print()

    # Header for operator-side lens table
    cols = [
        ('case',           14),
        ('F87 class',      10),
        ('‖M‖²',           12),
        ('anti %',          8),
        ('F80 |spec|max',  14),
        ('cusp pattern',   12),
        ('Lebensader',     45),
    ]
    print(''.join(f"{name:<{w}s}" for name, w in cols))
    print('-' * sum(w for _, w in cols))

    for label, terms in cases:
        # F87
        f87 = fw.classify_pauli_pair(chain, terms)

        # F49 ‖M‖² (including T1 contribution)
        try:
            f49 = fw.predict_residual_norm_squared_from_terms(
                chain, terms, gamma_t1=GAMMA_T1,
            )
        except Exception:
            f49 = None

        # F81 M_anti fraction (Π²-odd content as fraction of total)
        try:
            anti = fw.predict_pi_decomposition_anti_fraction(chain, terms)
        except Exception:
            anti = None

        # F80 Π²-odd spectrum: max |imag| of spectrum (zero if pure Π²-even)
        try:
            spec_dict = fw.predict_M_spectrum_pi2_odd(chain, terms)
            # spec_dict: {eigenvalue: multiplicity}; pull out max |imag|
            vals = [abs(complex(k).imag) for k in spec_dict.keys()]
            f80_max = max(vals) if vals else None
        except Exception:
            f80_max = None  # N/A (single-body or Π²-even)

        # Cockpit panel
        terms_j = [(t[0], t[1], J) for t in terms]
        H = fw._build_bilinear(N, bonds, terms_j)
        panel = lebensader_cockpit_panel(
            H, [GAMMA_DEPH] * N, rho_0, N,
            gamma_t1_l=[GAMMA_T1] * N, t_max=8.0, dt=0.005,
        )
        cusp_pat = panel['cusp']['pattern']
        rating = panel['lebensader']['rating']

        anti_pct = f"{anti*100:>7.2f}%" if anti is not None else "      —"

        row = [
            f"{label:<14s}",
            f"{f87:<10s}",
            fmt_num(f49, 12, 4),
            anti_pct + ' ',
            fmt_num(f80_max, 14, 4),
            f"{cusp_pat:<12s}",
            f"{rating:<45s}",
        ]
        print(''.join(row))

    print()
    print("Reading guide (anchor-perspective):")
    print("  ‖M‖²       = F49 closed form. T1 adds a fixed γ_T1 contribution.")
    print("               truly XX+YY has nearly zero (Π²-even bilinears).")
    print("               Z-containing pairs lift to O(1000) at N=3.")
    print("  anti %     = F81: |M_anti|²/|M|². 0 = pure Π²-even or single-body;")
    print("               50% = 2-body Π²-odd in clean balance (XY+YX); >50% = Π²-odd-dominated.")
    print("  F80 |spec|max = magnitude of largest M-eigenvalue (always imaginary on Π²-odd 2-body).")
    print("               '—' = single-body (F78 territory) or pure Π²-even.")
    print("  cusp pattern = monotonic (1 crossing of K=γt threshold) | heartbeat (>1).")
    print("  Lebensader = composite: skeleton (Π-protected count drop) +")
    print("               trace (θ-tail) + cusp + rating.")
    print()
    print("Observation prompt: which lens differentiates {truly, soft-Z-free, soft-Z-containing}?")


if __name__ == "__main__":
    main()
