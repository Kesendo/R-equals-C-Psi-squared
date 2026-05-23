"""Check whether the 8 cluster-doublet residual degeneracies match free-fermion
accidental energy degeneracies in the N=6 open chain.

Uses existing framework: XX+YY chain Jordan-Wigner maps to free fermions with
energies eps_k = 2J*cos(pi*k/(N+1)) for k=1..N (cf XyJordanWignerModes.cs,
PROOF_F80_BLOCH_SIGNWALK.md).

For Sz_L=0 (bra side has total Sz = 0): 3 fermions out of 6 sites.
For Sz_R=±1 (ket side total Sz=±1): 4 or 2 fermions.

The cluster's L-eigenvalue Im part = sum of (bra-fermion energies) - sum of
(ket-fermion energies). If two distinct (bra-config, ket-config) pairs give
the same Im AND have the same Sz_L, Sz_R, that's an accidental free-fermion
degeneracy and explains the residual cluster-doublet.

Investigation only.
"""
import sys
import math
import itertools
from collections import Counter, defaultdict

import numpy as np

sys.stdout.reconfigure(encoding="utf-8")

N = 6
J = 0.075


def epsilon_k(k, N=N, J=J):
    """Open-chain XX+YY single-particle energy."""
    return 2 * J * math.cos(math.pi * k / (N + 1))


def enumerate_configs(N, n_fermions):
    """All subsets of {1,...,N} of size n_fermions, returned as sorted tuples."""
    return list(itertools.combinations(range(1, N + 1), n_fermions))


def config_energy(config, J=J, N=N):
    return sum(epsilon_k(k, N, J) for k in config)


def main():
    print(f"N={N}, J={J}")
    print(f"single-particle energies eps_k = 2J*cos(pi*k/{N+1}):")
    for k in range(1, N + 1):
        print(f"  eps_{k} = {epsilon_k(k):+.6f}")
    print(f"  (note: eps_k + eps_{{N+1-k}} = 0 by mirror symmetry)\n")

    # For Sz_L=0 (3 fermions on bra): enumerate and group by energy
    print(f"=== 3-fermion configs (Sz_L=0 or Sz_R=0, sums of 3 eps_k) ===")
    configs_3f = enumerate_configs(N, 3)
    energies_3f = [(c, config_energy(c)) for c in configs_3f]
    # Group by rounded energy
    energy_groups = defaultdict(list)
    for c, e in energies_3f:
        energy_groups[round(e, 6)].append(c)
    # Report each group
    print(f"3-fermion energy spectrum ({len(configs_3f)} configs, "
          f"{len(energy_groups)} distinct energies):")
    for e in sorted(energy_groups.keys(), reverse=True):
        configs = energy_groups[e]
        marker = ' <- DEGENERATE' if len(configs) > 1 else ''
        print(f"  E = {e:+.6f}: {len(configs)} configs {configs}{marker}")

    n_doubled_3f = sum(1 for e, configs in energy_groups.items() if len(configs) >= 2)
    n_modes_doubled_3f = sum(len(configs) for e, configs in energy_groups.items() if len(configs) >= 2)
    print(f"\n  -> {n_doubled_3f} energy levels with multiplicity >= 2; "
          f"{n_modes_doubled_3f} configs total participating")

    # Enumerate 4-fermion configs (Sz_L = +1 or Sz_R = +1)
    print(f"\n=== 4-fermion configs (Sz_L=+1 or Sz_R=+1) ===")
    configs_4f = enumerate_configs(N, 4)
    energies_4f = [(c, config_energy(c)) for c in configs_4f]
    energy_groups_4f = defaultdict(list)
    for c, e in energies_4f:
        energy_groups_4f[round(e, 6)].append(c)
    print(f"4-fermion energy spectrum ({len(configs_4f)} configs, "
          f"{len(energy_groups_4f)} distinct energies):")
    for e in sorted(energy_groups_4f.keys(), reverse=True):
        configs = energy_groups_4f[e]
        marker = ' <- DEGENERATE' if len(configs) > 1 else ''
        print(f"  E = {e:+.6f}: {len(configs)} configs {configs}{marker}")

    # Now: L-eigenvalue Im part = E_bra - E_ket
    # We need to find (bra_config, ket_config) pairs with given Sz_L (bra fermion count)
    # and Sz_R (ket fermion count) that give specific Im.

    # Cluster |Im| = 0.02338. Look for (E_bra - E_ket) = +/-0.02338 across various
    # (n_bra, n_ket) compatible with cluster sector.

    target_Im = 0.02338
    print(f"\n=== Finding (E_bra - E_ket) = +/-{target_Im} across Sz sectors ===")
    Im_tol = 0.001
    matches_doublets = defaultdict(list)  # (n_bra, n_ket, e_bra, e_ket) -> [(c_bra, c_ket), ...]

    # Iterate over sectors: bra fermion count n_b, ket fermion count n_k
    # Cluster Sz_L from {-2..+2} via Sz_L = (n_up_bra - n_down_bra)/2 = n_b - N/2 = n_b - 3 (for N=6)
    # So Sz_L = -2 means n_b = 1, Sz_L = -1 means n_b = 2, Sz_L = 0 means n_b = 3, etc.
    sz_to_nb = {-2: 1, -1: 2, 0: 3, 1: 4, 2: 5}
    for sz_L, n_b in sz_to_nb.items():
        for sz_R, n_k in sz_to_nb.items():
            cfg_b = enumerate_configs(N, n_b)
            cfg_k = enumerate_configs(N, n_k)
            # Enumerate (c_b, c_k) pairs with |E_b - E_k| ≈ target_Im
            for c_b in cfg_b:
                for c_k in cfg_k:
                    diff = config_energy(c_b) - config_energy(c_k)
                    if abs(abs(diff) - target_Im) < Im_tol:
                        key = (sz_L, sz_R, round(diff, 4))
                        matches_doublets[key].append((c_b, c_k))

    print(f"\nFound matches for |E_b - E_k| ≈ {target_Im} (tol {Im_tol}):")
    print(f"per (Sz_L, Sz_R) sector — counting compatible (c_b, c_k) pairs:")
    for (sz_L, sz_R, diff), pairs in sorted(matches_doublets.items()):
        if (sz_L, sz_R) in [(0, 1), (0, -1), (1, 0), (-1, 0)]:
            print(f"\n  (Sz_L={sz_L}, Sz_R={sz_R}, Im={diff:+.4f}): {len(pairs)} compatible config-pairs")
            for c_b, c_k in pairs[:6]:
                e_b = config_energy(c_b)
                e_k = config_energy(c_k)
                print(f"    bra {c_b} E={e_b:+.4f}, ket {c_k} E={e_k:+.4f}, diff={e_b - e_k:+.4f}")
            if len(pairs) > 6:
                print(f"    ... and {len(pairs) - 6} more")

    # Check the doublet structure: in sector (Sz_L=0, Sz_R=+1), how many DISTINCT energies
    # give rise to multiple (c_b, c_k) combinations?
    print(f"\n=== Multiplicity of (Sz_L, Sz_R, Im) tuples in doublet sectors ===")
    for sec in [(0, 1), (0, -1), (1, 0), (-1, 0)]:
        sz_L, sz_R = sec
        keys_in_sec = [k for k in matches_doublets if k[:2] == sec]
        for k in sorted(keys_in_sec):
            print(f"  Sz_L={k[0]}, Sz_R={k[1]}, Im={k[2]:+.4f}: "
                  f"{len(matches_doublets[k])} compatible pairs")


if __name__ == "__main__":
    main()
