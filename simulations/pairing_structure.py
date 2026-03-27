"""
Pairing Structure of the V-Effect: With whom are the 100 new frequencies paired?

Two N=2 resonators have 2 frequencies each (4 total).
Coupled through a mediator (N=5), they have 104 frequencies.
100 are new. But every oscillating mode is palindromically paired (proven).

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

seen = set()
for k in range(len(evals5)):
    if k in seen:
        continue

    j = partners5[k]
    if j < 0 or j == k:
        counts["UNPAIRED"] += 1
        continue

    seen.add(k)
    seen.add(j)

    freq_k = abs(evals5[k].imag)
    freq_j = abs(evals5[j].imag)
    rate_k = -evals5[k].real
    rate_j = -evals5[j].real

    is_osc_k = freq_k > 0.01
    is_osc_j = freq_j > 0.01

    if not is_osc_k and not is_osc_j:
        counts["DECAY-DECAY"] += 1
        continue

    k_old = freq_is_old(round(freq_k, 6), old_set)
    j_old = freq_is_old(round(freq_j, 6), old_set)

    if k_old and j_old:
        counts["OLD-OLD"] += 1
        pair_details["OLD-OLD"].append((freq_k, freq_j, rate_k, rate_j))
    elif not k_old and not j_old:
        counts["NEW-NEW"] += 1
        pair_details["NEW-NEW"].append((freq_k, freq_j, rate_k, rate_j))
    else:
        counts["OLD-NEW"] += 1
        pair_details["OLD-NEW"].append((freq_k, freq_j, rate_k, rate_j))

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
