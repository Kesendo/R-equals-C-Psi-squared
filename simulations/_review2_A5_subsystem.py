"""
Independent from-below verifier for review finding A5 (PROOF_SUBSYSTEM_CROSSING.md, Case C).

Conjecture 2.1 / F28 claim under attack:
  "For any PRIMITIVE CPTP map on a 2-qubit system, every state with CPsi>1/4
   eventually has CPsi<1/4. The 1/4 boundary is an eventual absorber for ALL
   primitive quantum channels."

CPsi metric is the proof's OWN definition (PROOF_SUBSYSTEM_CROSSING.md:157):
   CPsi(rho) = Tr(rho^2) * L1(rho) / (d-1),   d = 4,   L1 = sum_{i!=j} |rho_ij|.

This script is refute-first: it tries to BREAK the claim three ways, and also
tries to break the proposed REPLACEMENT (scope to physical noise) by checking the
standard physical channels actually stay <= 1/4.
"""
import numpy as np

np.set_printoptions(precision=5, suppress=True)
d = 4
I4 = np.eye(d)

# ---- the proof's own CPsi metric ----
def L1(rho):
    return float(np.sum(np.abs(rho)) - np.sum(np.abs(np.diag(rho))))  # off-diagonal |.| sum
def purity(rho):
    return float(np.real(np.trace(rho @ rho)))
def CPsi(rho):
    return purity(rho) * L1(rho) / (d - 1)

# ---- Bell+ and the counterexample fixed point sigma ----
phi = np.zeros(d); phi[0] = 1/np.sqrt(2); phi[3] = 1/np.sqrt(2)
P = np.outer(phi, phi.conj())           # |Phi+><Phi+|
sigma = 0.95 * P + 0.05 * I4 / 4

print("=" * 70)
print("PILLAR 1: is sigma a legitimate full-rank state, and is CPsi(sigma)>1/4")
print("          in the proof's OWN metric?")
print("=" * 70)
ev = np.linalg.eigvalsh(sigma)
print("sigma eigenvalues:", ev, " (all > 0 => full-rank:", bool(np.all(ev > 1e-12)), ")")
print(f"purity(sigma) = {purity(sigma):.6f}")
print(f"L1(sigma)     = {L1(sigma):.6f}")
print(f"CPsi(sigma)   = {CPsi(sigma):.6f}   (>1/4: {CPsi(sigma) > 0.25})")
print(f"CPsi(Bell+ pure) = {CPsi(P):.6f}")

# ---- the channel eps(rho) = (1-p) rho + p Tr(rho) sigma ----
p = 0.30
def eps(rho):
    return (1 - p) * rho + p * np.trace(rho) * sigma

def superop(channel):
    D2 = d * d
    M = np.zeros((D2, D2), dtype=complex)
    for k in range(D2):
        e = np.zeros(D2, dtype=complex); e[k] = 1
        M[:, k] = channel(e.reshape(d, d)).reshape(-1)
    return M

print("\n" + "=" * 70)
print("PILLAR 2: is eps genuinely PRIMITIVE (unique fixed point, aperiodic)?")
print("=" * 70)
M = superop(eps)
evals = np.linalg.eigvals(M)
order = np.argsort(-np.abs(evals))
print("top-6 channel eigenvalues by |.|:")
for z in evals[order][:6]:
    print(f"   {z:.5f}   |.|={abs(z):.5f}")
peripheral = np.sum(np.abs(np.abs(evals) - 1.0) < 1e-8)
print(f"# eigenvalues on the unit circle = {peripheral}  "
      f"(==1 simple => primitive/aperiodic, unique fixed point)")

# ---- fixed point + convergence from random states ----
rng = np.random.default_rng(7)
def random_state():
    A = rng.normal(size=(d, d)) + 1j * rng.normal(size=(d, d))
    r = A @ A.conj().T
    return r / np.trace(r)

print("\n" + "=" * 70)
print("PILLAR 2b: does every state converge to sigma, and does the boundary")
print("           1/4 FAIL to absorb (states starting >1/4 stay >1/4)?")
print("=" * 70)
worst_min = 1.0
for trial in range(20):
    r = random_state()
    for _ in range(4000):
        r = eps(r)
    err = np.linalg.norm(r - sigma)
    if trial < 3:
        print(f"  random start -> ||rho_inf - sigma|| = {err:.2e}, CPsi_inf = {CPsi(r):.5f}")
# from Bell+ (CPsi=1/3 > 1/4): does it ever cross below 1/4?
r = P.copy(); ever_below = False; cmin = CPsi(r); traj = []
for _ in range(5000):
    r = eps(r); c = CPsi(r); cmin = min(cmin, c); traj.append(c)
    if c < 0.25:
        ever_below = True
print(f"  from Bell+ (CPsi0={CPsi(P):.4f}): ever CPsi<1/4? {ever_below};  "
      f"min CPsi={cmin:.5f};  final CPsi={CPsi(r):.5f}")
print(f"  => boundary 1/4 is {'NOT ' if not ever_below else ''}an absorber for this primitive channel")

print("\n" + "=" * 70)
print("PILLAR 3: is 'max 0.138 over random maps' a SAMPLING artifact?")
print("          (the prior session claims env-dim-2 / near-identity ensembles")
print("           VIOLATE; I test the violation FRACTION among valid+primitive fps)")
print("=" * 70)
def random_cptp(nk, rng):
    G = rng.normal(size=(d * nk, d)) + 1j * rng.normal(size=(d * nk, d))
    Q, _ = np.linalg.qr(G)                # Haar isometry: Q^dag Q = I_d => sum K_i^dag K_i = I
    Ks = [Q[i * d:(i + 1) * d, :] for i in range(nk)]
    return lambda rho: sum(K @ rho @ K.conj().T for K in Ks)
def fixed_point_valid(channel):
    """Return (fp, is_valid_state, is_primitive). Only count genuine PSD trace-1
    fixed points of channels with a SIMPLE peripheral spectrum (true primitivity)."""
    M = superop(channel)
    ev, evec = np.linalg.eig(M)
    on_circle = int(np.sum(np.abs(np.abs(ev) - 1.0) < 1e-6))
    v = evec[:, np.argmin(np.abs(ev - 1.0))].reshape(d, d)
    v = (v + v.conj().T) / 2
    tr = np.real(np.trace(v))
    if abs(tr) < 1e-9:
        return None, False, False
    v = v / tr
    w = np.linalg.eigvalsh(v)
    is_state = bool(w[0] > -1e-7 and abs(np.real(np.trace(v)) - 1.0) < 1e-6)
    return v, is_state, (on_circle == 1)

def report_ensemble(name, make_channel, trials=400):
    n_valid = n_viol = 0; maxc = 0.0
    for _ in range(trials):
        fp, ok, prim = fixed_point_valid(make_channel())
        if fp is None or not (ok and prim):
            continue
        n_valid += 1; c = CPsi(fp); maxc = max(maxc, c)
        if c > 0.25 + 1e-9:
            n_viol += 1
    frac = 100.0 * n_viol / max(n_valid, 1)
    print(f"  {name:34s}: valid&primitive={n_valid:3d}/{trials}, "
          f"max CPsi={maxc:.4f}, >1/4: {n_viol} ({frac:.1f}%)")

# Forward sampling at several Kraus ranks (nk=4 ~ the proof's likely ensemble)
for nk in [2, 3, 4, 6, 8, 16]:
    report_ensemble(f"Haar/Ginibre n_kraus={nk}", lambda nk=nk: random_cptp(nk, rng))

# Reverse construction: sample a random TARGET state, build the primitive
# trace-and-replace channel eps=(1-p)Id + p*R_tau (fixed point = tau, always primitive).
# This samples the corner of channel-space the forward Ginibre sampler never reaches.
def random_target_state(rank, rng):
    A = rng.normal(size=(d, rank)) + 1j * rng.normal(size=(d, rank))
    r = A @ A.conj().T
    return r / np.trace(r)
def trace_replace_channel(tau, p=0.3):
    return lambda rho: (1 - p) * rho + p * np.trace(rho) * tau
for rank in [1, 2, 3]:
    report_ensemble(f"trace-replace, target rank={rank}",
                    lambda rank=rank: trace_replace_channel(random_target_state(rank, rng)))

print("\n  Constructed counterexample FAMILY (vary mixing p and entanglement q):")
for q in [0.80, 0.90, 0.95, 0.99]:
    s = q * P + (1 - q) * I4 / 4
    print(f"     sigma = {q:.2f}|Phi+><Phi+| + {1-q:.2f} I/4 : "
          f"full-rank={bool(np.all(np.linalg.eigvalsh(s) > 1e-12))}, CPsi={CPsi(s):.4f}")

print("\n" + "=" * 70)
print("REPLACEMENT CHECK: do the STANDARD PHYSICAL channels keep fixed point <= 1/4?")
print("  (this is what the scope-retraction claims SURVIVES)")
print("=" * 70)
X = np.array([[0, 1], [1, 0]], complex); Y = np.array([[0, -1j], [1j, 0]]); Z = np.array([[1, 0], [0, -1]], complex)
def local_pauli_dephase(g):
    # Z-dephasing on both qubits
    Kset = []
    import itertools
    a = np.sqrt(1 - g); b = np.sqrt(g)
    Z1 = np.kron(Z, np.eye(2)); Z2 = np.kron(np.eye(2), Z)
    return [a * a * np.eye(4), a * b * Z1, b * a * Z2, b * b * Z1 @ Z2]
def amp_damp(gamma):
    K0 = np.array([[1, 0], [0, np.sqrt(1 - gamma)]], complex)
    K1 = np.array([[0, np.sqrt(gamma)], [0, 0]], complex)
    # local AD on both qubits
    Ks = []
    for A in (K0, K1):
        for B in (K0, K1):
            Ks.append(np.kron(A, B))
    return Ks
chans = {
    "depolarizing(0.3)": lambda rho: 0.7 * rho + 0.3 * np.trace(rho) * I4 / 4,
    "Z-dephasing(0.3)":  (lambda Ks: (lambda rho: sum(K @ rho @ K.conj().T for K in Ks)))(local_pauli_dephase(0.3)),
    "amp-damp(0.3)":     (lambda Ks: (lambda rho: sum(K @ rho @ K.conj().T for K in Ks)))(amp_damp(0.3)),
}
for name, ch in chans.items():
    fp, ok, prim = fixed_point_valid(ch)
    print(f"  {name:20s}: fixed-point CPsi = {CPsi(fp):.5f}  (<=1/4: {CPsi(fp) <= 0.25 + 1e-9}); "
          f"valid_state={ok}, primitive(unique-fp)={prim}")

print("\nDONE.")
