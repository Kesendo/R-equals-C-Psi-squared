"""
Pairing Structure of the V-Effect: With whom are the 109 new frequencies paired?

Two N=2 resonators have 2 frequencies each (4 total).
Coupled through a mediator (N=5), they have 109 frequencies, all new.
And every oscillating mode is palindromically paired (proven).

Question: are the new frequencies paired with each other (NEW-NEW),
with the old ones (OLD-NEW), or mixed?

The answer tells us what the V-Effect does to the palindrome.
"""
import numpy as np

I2 = np.eye(2, dtype=complex)
sx = np.array([[0, 1], [1, 0]], dtype=complex)
sy = np.array([[0, -1j], [1j, 0]], dtype=complex)
sz = np.array([[1, 0], [0, -1]], dtype=complex)


def site_op(op, site, N):
    ops = [I2] * N
    ops[site] = op
    r = ops[0]
    for o in ops[1:]:
        r = np.kron(r, o)
    return r


def heisenberg_H(N, J=1.0):
    d = 2 ** N
    H = np.zeros((d, d), dtype=complex)
    for i in range(N - 1):
        for P in [sx, sy, sz]:
            H += J * site_op(P, i, N) @ site_op(P, i + 1, N)
    return H


def build_liouvillian(H, gamma, N):
    d = 2 ** N
    Id = np.eye(d, dtype=complex)
    L = -1j * (np.kron(H, Id) - np.kron(Id, H.T))
    for k in range(N):
        Zk = site_op(sz, k, N)
        L += gamma * (np.kron(Zk, Zk.conj()) - np.eye(d * d, dtype=complex))
    return L


def find_pairs(evals, Sg, tol=1e-4):
    """Find palindromic partner for each eigenvalue."""
    n = len(evals)
    partners = [-1] * n
    errors = [float('inf')] * n

    for k in range(n):
        target = -(evals[k] + 2 * Sg)
        dists = np.abs(evals - target)
        best = np.argmin(dists)
        if dists[best] < tol * max(1, abs(target)):
            partners[k] = int(best)
            errors[k] = float(dists[best])

    return partners, errors


# Step 1: N=2 spectrum
print("=" * 60)
print("STEP 1: N=2 SINGLE RESONATOR")
print("=" * 60)

N2 = 2
gamma = 0.05
Sg2 = N2 * gamma  # 0.10

H2 = heisenberg_H(N2)
L2 = build_liouvillian(H2, gamma, N2)
evals2 = np.linalg.eigvals(L2)

freqs2 = sorted(set(round(abs(ev.imag), 6) for ev in evals2 if abs(ev.imag) > 0.01))
print(f"Eigenvalues: {len(evals2)}")
print(f"Distinct frequencies (|Im|>0.01): {len(freqs2)}")
print(f"Frequencies: {freqs2}")
print(f"Sg = {Sg2}")

partners2, errors2 = find_pairs(evals2, Sg2)
n_paired2 = sum(1 for p in partners2 if p >= 0)
print(f"Paired: {n_paired2}/{len(evals2)}")


# Step 2: N=5 spectrum
print("\n" + "=" * 60)
print("STEP 2: N=5 COUPLED SYSTEM (MediatorBridge)")
print("=" * 60)

N5 = 5
Sg5 = N5 * gamma  # 0.25

H5 = heisenberg_H(N5)
L5 = build_liouvillian(H5, gamma, N5)
evals5 = np.linalg.eigvals(L5)

freqs5 = sorted(set(round(abs(ev.imag), 6) for ev in evals5 if abs(ev.imag) > 0.01))
print(f"Eigenvalues: {len(evals5)}")
print(f"Distinct frequencies (|Im|>0.01): {len(freqs5)}")
print(f"Sg = {Sg5}")

partners5, errors5 = find_pairs(evals5, Sg5)
n_paired5 = sum(1 for p in partners5 if p >= 0)
print(f"Paired: {n_paired5}/{len(evals5)}")


# Step 3: Classify each N=5 frequency as OLD or NEW
print("\n" + "=" * 60)
print("STEP 3: OLD vs NEW FREQUENCY CLASSIFICATION")
print("=" * 60)

freq_tol = 1e-4  # broader tolerance for shifted modes

old_freqs = []
new_freqs = []
for f in freqs5:
    is_old = any(abs(f - f2) < freq_tol for f2 in freqs2)
    if is_old:
        old_freqs.append(f)
    else:
        new_freqs.append(f)

print(f"OLD frequencies (exist in N=2): {len(old_freqs)}")
print(f"NEW frequencies (coupling only): {len(new_freqs)}")


# Step 4: THE PAIRING MAP
print("\n" + "=" * 60)
print("STEP 4: PAIRING MAP")
print("=" * 60)

def freq_is_old(f, old_set, tol=1e-4):
    return any(abs(f - fo) < tol for fo in old_set)

old_set = set(round(f, 6) for f in freqs2)

counts = {"OLD-OLD": 0, "NEW-NEW": 0, "OLD-NEW": 0,
           "DECAY-DECAY": 0, "UNPAIRED": 0}

pair_details = {"OLD-OLD": [], "NEW-NEW": [], "OLD-NEW": []}

# Palindromic pairing via optimal assignment. The spectrum is heavily
# degenerate, so per-mode "best partner" maps are ill-defined: the original
# asymmetric count let several modes share one partner (556 "pairs" out of
# only 904 oscillating modes, arithmetically impossible, disjoint maximum
# 452), and mutual-best matching starves on ties. Optimal assignment
# handles degeneracy exactly and pairs the full spectrum.
from scipy.optimize import linear_sum_assignment

mirror5 = -evals5.conj() - 2 * Sg5
cost = np.abs(evals5[:, None] - mirror5[None, :])
ri, ci = linear_sum_assignment(cost)
max_pair_err = cost[ri, ci].max()
print(f"Optimal-assignment palindromy: max pairing error {max_pair_err:.2e} "
      f"over {len(evals5)} modes (exact palindrome => ~0)")

osc_mask = np.abs(evals5.imag) > 0.01
n_osc = int(osc_mask.sum())
n_dec = len(evals5) - n_osc
# The mirror partner of an oscillating mode is oscillating (|Im| is
# preserved by lam -> -conj(lam) - 2*Sg), so oscillating modes pair among
# themselves: n_osc/2 disjoint pairs, and pure-decay modes likewise.
n_old_osc = sum(1 for k in range(len(evals5)) if osc_mask[k]
                and freq_is_old(round(abs(evals5[k].imag), 6), old_set))
counts["DECAY-DECAY"] = n_dec // 2
if n_old_osc == 0:
    counts["NEW-NEW"] = n_osc // 2
else:
    print(f"WARNING: {n_old_osc} oscillating modes carry OLD frequencies; "
          f"classify pairs individually before trusting the NEW-NEW count.")
counts["UNPAIRED"] = len(evals5) - 2 * (counts["NEW-NEW"] + counts["DECAY-DECAY"])
print(f"Oscillating modes: {n_osc} ({n_osc // 2} pairs), "
      f"pure-decay modes: {n_dec} ({n_dec // 2} pairs), "
      f"oscillating modes at OLD frequencies: {n_old_osc}")

print(f"Total palindromic pairs (oscillating + decay):")
for cat, count in counts.items():
    print(f"  {cat}: {count}")

print(f"\nOscillating pairs only:")
for cat in ["OLD-OLD", "NEW-NEW", "OLD-NEW"]:
    n = counts[cat]
    print(f"  {cat}: {n}")
    if pair_details[cat]:
        for fk, fj, rk, rj in pair_details[cat][:5]:
            print(f"    freq {fk:.4f}/{fj:.4f}  rate {rk:.4f}/{rj:.4f}")
        if len(pair_details[cat]) > 5:
            print(f"    ... and {len(pair_details[cat])-5} more")

# Summary
print("\n" + "=" * 60)
print("SUMMARY")
print("=" * 60)
total_osc = counts["OLD-OLD"] + counts["NEW-NEW"] + counts["OLD-NEW"]
if total_osc > 0:
    print(f"OLD-OLD: {counts['OLD-OLD']} ({counts['OLD-OLD']/total_osc*100:.1f}% of oscillating)")
    print(f"NEW-NEW: {counts['NEW-NEW']} ({counts['NEW-NEW']/total_osc*100:.1f}%)")
    print(f"OLD-NEW: {counts['OLD-NEW']} ({counts['OLD-NEW']/total_osc*100:.1f}%)")
    print(f"DECAY-DECAY: {counts['DECAY-DECAY']}")
    print(f"UNPAIRED: {counts['UNPAIRED']}")

if counts["OLD-NEW"] > 0:
    print("\n>>> OLD-NEW PAIRS EXIST: coupling re-wires the palindrome across old/new boundary <<<")
elif counts["NEW-NEW"] == total_osc:
    print("\n>>> ALL oscillating pairs are NEW-NEW: coupling creates a self-contained palindrome <<<")


# Step 5: XY-weight distribution of the oscillating modes (the balance histogram)
# For each oscillating right eigenvector (|Im| > 0.01), project onto the
# orthonormal Pauli basis vec(P)/sqrt(d) and accumulate the probability mass
# per XY-weight (number of X/Y letters). The histogram is a mass SHARE over
# modes, not an integer count.
print("\n" + "=" * 60)
print("STEP 5: XY-WEIGHT DISTRIBUTION (oscillating modes, N=5)")
print("=" * 60)

from itertools import product as _product

_evals5v, _evecs5 = np.linalg.eig(L5)
_osc_idx = [k for k in range(len(_evals5v)) if abs(_evals5v[k].imag) > 0.01]
print(f"Oscillating modes: {len(_osc_idx)}/{len(_evals5v)}")

_d5 = 2 ** N5
_letters = [I2, sx, sy, sz]
_weights = np.zeros(N5 + 1)
_P_cache = []
for _combo in _product(range(4), repeat=N5):
    _P = _letters[_combo[0]]
    for _c in _combo[1:]:
        _P = np.kron(_P, _letters[_c])
    _w = sum(1 for _c in _combo if _c in (1, 2))  # X or Y letters
    _P_cache.append((_P.reshape(-1) / np.sqrt(_d5), _w))

for _k in _osc_idx:
    _v = _evecs5[:, _k]
    _v = _v / np.linalg.norm(_v)
    for _pvec, _w in _P_cache:
        _weights[_w] += abs(np.vdot(_pvec, _v)) ** 2

_weights_pct = 100 * _weights / _weights.sum()
for _w in range(N5 + 1):
    print(f"  w={_w}: {_weights_pct[_w]:5.1f}%  (mass {_weights[_w]:.3f})")
print(f"  interior (w=2,3): {_weights_pct[2]+_weights_pct[3]:.1f}%")
print(f"  boundary (w=1,4): {_weights_pct[1]+_weights_pct[4]:.1f}%")
print(f"  extremes (w=0,5): {_weights_pct[0]+_weights_pct[5]:.1f}%")
