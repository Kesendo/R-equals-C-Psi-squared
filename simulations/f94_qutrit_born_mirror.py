"""Gate-first: does F94's 4/3 = a_-1/3 become d^2/3 = 3 at the qutrit? (the deepest family-A test)

CONTEXT. F94: for |0+0+> on an N=4 Heisenberg ring + Z-dephasing, pair (0,2), |00> outcome,
  Delta_|00>(Q,K) = (4/3) Q^2 K^3,   c = sym3_element / (6 * P_u0) = 8/6 = 4/3,
where sym3 = L_H^2 L'_dis + L_H L'_dis L_H + L'_dis L_H^2 is the gamma^1 piece of L^3 in the time Taylor.
The discriminant anchor reads 4/3 = a_-1/3 with a_-1 = d^2 = 4 (the (J/4)^2 Heisenberg-coupling factor),
i.e. it classifies F94's 4 in FAMILY A (squared dimension). The qutrit prism tests it: if family A, then
c(d=3) -> d^2/3 = 9/3 = 3.

THE MODEL (pre-registered, faithful: reduces to F94 EXACTLY at d=2).
  * H = (J/4) * sum_bonds sum_a lambda^a_i lambda^a_j  (generalized Gell-Mann, ring N=4). At d=2 the
    lambda^a are the 3 Paulis, so H = (J/4)(XX+YY+ZZ) = F94's H. NOTE the Fierz identity
    sum_a lambda^a (x) lambda^a = 2*SWAP - (2/d) I, so the DYNAMICS L_H = -i[H,.] = -i(J/2)[SWAP,.] for
    EVERY d (the identity part drops from the commutator). So the physical coupling is /2 (the SWAP), NOT
    /4 and NOT /d^2 -- the "/4 = d^2" reading is the d=2 coincidence (1/2)^2 = 1/d^2. (Cross-checked
    against an explicit (J/2)*SWAP build below: identical c.)
  * L'_dis[rho]_{a,b} = -2 * Hamming(a,b) * rho_{a,b}  (full-Cartan equidistant dephasing, gamma=1). At
    d=2 this equals sum_l (Z_l rho Z_l - rho) exactly. The F121 qudit dephasing convention.
  * |+> = (|0> + ... + |d-1>) / sqrt(d)  (equal superposition; at d=2 = (|0>+|1>)/sqrt2). rho_0 = |0+0+>.
  * pair (0,2), |00> outcome, P_u0 = <00|Tr_{1,3}[rho_0]|00> = 1 (sites 0,2 are |0>).

  STAGE 0 (port fidelity): d=2 reproduces c = 4/3 bit-close. MUST pass.
  STAGE 1 (the probe): d=3 -> c(3). Gate asserts the FAMILY-A hypothesis c = d^2/3 = 3. A firing gate IS
    the find: diagnose whether c stays 4/3 (the 4 is a setup-count, d-independent) or is something else
    (setup-specific, as the proof's own universality remark warns). Decompose to show what carries d.

Run:  python simulations/f94_qutrit_born_mirror.py
"""
import itertools
import numpy as np

# ----------------------------------------------------------------------------------------------------
# generalized Gell-Mann generators (d^2 - 1 of them, Tr(lambda^a lambda^b) = 2 delta; Paulis at d=2)
# ----------------------------------------------------------------------------------------------------
def gell_mann(d):
    gens = []
    for j in range(d):
        for k in range(j + 1, d):
            s = np.zeros((d, d), dtype=complex); s[j, k] = 1; s[k, j] = 1          # symmetric
            gens.append(s)
            a = np.zeros((d, d), dtype=complex); a[j, k] = -1j; a[k, j] = 1j        # antisymmetric
            gens.append(a)
    for l in range(1, d):                                                          # diagonal
        diag = np.zeros(d, dtype=complex)
        for j in range(l):
            diag[j] = 1.0
        diag[l] = -l
        h = np.sqrt(2.0 / (l * (l + 1))) * np.diag(diag)
        gens.append(h)
    return gens   # length d^2 - 1


def kron_at(d, N, ops_by_site):
    out = np.array([[1]], dtype=complex)
    for k in range(N):
        out = np.kron(out, ops_by_site.get(k, np.eye(d, dtype=complex)))
    return out


def swap_op(d):
    S = np.zeros((d * d, d * d), dtype=complex)
    for i in range(d):
        for j in range(d):
            S[i * d + j, j * d + i] = 1.0
    return S


def heisenberg_ring(d, N, J=1.0, use_swap=False):
    """H = (J/4) sum_bonds sum_a lambda^a_i lambda^a_j  (Gell-Mann), or (J/2)*SWAP per bond if use_swap."""
    bonds = [(i, (i + 1) % N) for i in range(N)]
    D = d ** N
    H = np.zeros((D, D), dtype=complex)
    if use_swap:
        S = swap_op(d)
        for (i, j) in bonds:
            # embed the 2-site SWAP on sites (i,j); build via permutation of basis labels
            H += (J / 2.0) * embed_two_site(d, N, i, j, S)
        return H
    gens = gell_mann(d)
    for (i, j) in bonds:
        for lam in gens:
            H += (J / 4.0) * embed_two_site(d, N, i, j, np.kron(lam, lam))
    return H


def embed_two_site(d, N, i, j, op2):
    """Embed a (d^2 x d^2) two-site operator op2 (ordered site-i (x) site-j) into the full d^N space."""
    D = d ** N
    out = np.zeros((D, D), dtype=complex)
    sites = list(range(N))
    for a in itertools.product(range(d), repeat=N):
        ai, aj = a[i], a[j]
        rowin = ai * d + aj
        for bi in range(d):
            for bj in range(d):
                amp = op2[bi * d + bj, rowin]
                if amp != 0:
                    b = list(a); b[i] = bi; b[j] = bj
                    out[idx(b, d), idx(a, d)] += amp
    return out


def idx(digits, d):
    x = 0
    for v in digits:
        x = x * d + v
    return x


def L_H(rho, H):
    return -1j * (H @ rho - rho @ H)


def hamming_diag(d, N):
    """The vector of -2*Hamming(a,b) over the flattened (a,b) coherence basis, for the equidistant dephasing."""
    states = list(itertools.product(range(d), repeat=N))
    D = len(states)
    M = np.zeros((D, D))
    for ia, a in enumerate(states):
        for ib, b in enumerate(states):
            M[ia, ib] = -2.0 * sum(1 for x, y in zip(a, b) if x != y)
    return M


def L_dis(rho, deph_mask):
    return deph_mask * rho            # element-wise; deph_mask[a,b] = -2*Hamming(a,b)


def reduced_pair(rho, d, N, keep):
    """Partial trace keeping the qudits in `keep` (sorted), tracing the rest."""
    t = rho.reshape([d] * N + [d] * N)
    trace = [q for q in range(N) if q not in keep]
    for q in sorted(trace, reverse=True):
        offset = N - sum(1 for r in trace if r > q)
        t = np.trace(t, axis1=q, axis2=q + offset)
    nk = len(keep)
    return t.reshape((d ** nk, d ** nk))


def f94_coefficient(d, N=4, J=1.0, use_swap=False):
    """c = <00|_pair Tr_{1,3}[sym3 rho_0] |00>_pair / (6 * P_u0)."""
    plus = np.ones(d, dtype=complex) / np.sqrt(d)
    zero = np.zeros(d, dtype=complex); zero[0] = 1.0
    psi = zero
    for site in range(1, N):
        psi = np.kron(psi, plus if site % 2 == 1 else zero)
    rho0 = np.outer(psi, psi.conj())

    H = heisenberg_ring(d, N, J=J, use_swap=use_swap)
    mask = hamming_diag(d, N)

    A1 = L_H(L_H(L_dis(rho0, mask), H), H)        # L_H^2 L'_dis
    A2 = L_H(L_dis(L_H(rho0, H), mask), H)        # L_H L'_dis L_H
    A3 = L_dis(L_H(L_H(rho0, H), H), mask)        # L'_dis L_H^2
    sym3 = A1 + A2 + A3

    pair = reduced_pair(sym3, d, N, keep=[0, 2])
    elem = float(np.real(pair[0, 0]))             # <00| ... |00>
    pu0 = float(np.real(reduced_pair(rho0, d, N, keep=[0, 2])[0, 0]))
    return elem, pu0, elem / (6.0 * pu0)


# ====================================================================================================
# STAGE 0 -- PORT FIDELITY: d=2 reproduces F94's c = 4/3
# ====================================================================================================
print("=" * 100)
print("STAGE 0 -- PORT FIDELITY: the generalized build at d=2 reproduces F94's c = 4/3")
print("=" * 100)
elem2, pu2, c2 = f94_coefficient(2)
elem2s, pu2s, c2s = f94_coefficient(2, use_swap=True)
print(f"  d=2 Gell-Mann (J/4)*sum(lambda lambda): sym3_elem = {elem2:.6f}, P_u0 = {pu2:.6f}, c = {c2:.9f}")
print(f"  d=2 SWAP      (J/2)*SWAP             : sym3_elem = {elem2s:.6f}, P_u0 = {pu2s:.6f}, c = {c2s:.9f}")
assert abs(elem2 - 8.0) < 1e-9, f"STAGE 0 GATE: d=2 sym3 element {elem2} != 8 (F94 bit-exact)"
assert abs(c2 - 4.0 / 3.0) < 1e-9, f"STAGE 0 GATE: d=2 coefficient {c2} != 4/3"
assert abs(c2 - c2s) < 1e-9, "STAGE 0 GATE: Gell-Mann and SWAP builds disagree (the constant should drop)"
print(f"\nSTAGE 0 PASS: d=2 gives sym3 element = 8 and c = 4/3 EXACTLY (Gell-Mann == SWAP, the identity part "
      f"drops).\n  The faithful port reproduces F94.")

# ====================================================================================================
# STAGE 1 -- THE PROBE: d=3 -> c(3). Does F94's 4/3 become d^2/3 = 3 (family A)?
# ====================================================================================================
print("\n" + "=" * 100)
print("STAGE 1 -- THE QUTRIT: does F94's 4/3 -> d^2/3 = 3 (family A), stay 4/3, or something else?")
print("=" * 100)
elem3, pu3, c3 = f94_coefficient(3)
elem3s, pu3s, c3s = f94_coefficient(3, use_swap=True)
print(f"  d=3 Gell-Mann: sym3_elem = {elem3:.6f}, P_u0 = {pu3:.6f}, c = {c3:.9f}")
print(f"  d=3 SWAP     : sym3_elem = {elem3s:.6f}, P_u0 = {pu3s:.6f}, c = {c3s:.9f}")
assert abs(c3 - c3s) < 1e-9, "STAGE 1: Gell-Mann and SWAP disagree at d=3 (constant should still drop)"

cands = {"d^2/3 = 3 (family A)": 3.0, "4/3 (stays, family C/count)": 4.0 / 3.0,
         "2d/3 = 2 (family B)": 2.0, "d/3 = 1": 1.0}
print(f"\n  c(d=3) = {c3:.6f}.  Candidate readings:")
for name, val in cands.items():
    print(f"    {name:30} = {val:.6f}   {'<-- MATCH' if abs(val - c3) < 1e-6 else ''}")

# PRIMARY GATE (the family-A hypothesis, Tom's question): the gate FIRES -> family A is REFUTED for F94.
family_A = abs(c3 - 3.0) < 1e-6
print(f"\n  FAMILY-A HYPOTHESIS (c -> d^2/3 = 3): {'CONFIRMED' if family_A else 'REFUTED'} "
      f"(c(d=3) = {c3:.6f}, d^2/3 = 3.000000).")
assert not family_A, "unexpected: family A confirmed -- revisit the analysis"
print("  => F94's 4/3 does NOT generalize as d^2/3. The dynamics is (J/2)*SWAP (d-INDEPENDENT coupling, "
      "confirmed\n     Gell-Mann==SWAP), so the '4/3 = a_-1/3' reading is the d=2 coincidence (J/4)=(1/2)^2=1/d^2.")

print("  => F94's 4 is NOT the squared-dimension discriminant. Diagnosing the ACTUAL c(d) ...")

# ====================================================================================================
# STAGE 2 -- THE ACTUAL CLOSED FORM (re-pointing the gate at the validated truth, NOT loosening it)
# ====================================================================================================
print("\n" + "=" * 100)
print("STAGE 2 -- the actual c(d): closed form c(d) = 4(d+2)(d-1) / (3 d^2)")
print("=" * 100)
from fractions import Fraction


def c_closed(d):
    return Fraction(4 * (d + 2) * (d - 1), 3 * d * d)


print(f"  {'d':>2} {'computed c':>14} {'4(d+2)(d-1)/(3 d^2)':>22} {'match?':>7}")
worst = 0.0
for d in (2, 3, 4, 5, 6, 7):
    _, _, c = f94_coefficient(d)
    cf = c_closed(d)
    worst = max(worst, abs(c - float(cf)))
    print(f"  {d:>2} {c:>14.9f} {str(cf) + ' = ' + format(float(cf), '.6f'):>22} "
          f"{('YES' if abs(c - float(cf)) < 1e-9 else 'NO'):>7}")
assert worst < 1e-9, f"STAGE 2 GATE FIRED: c(d) != 4(d+2)(d-1)/(3 d^2) (worst {worst:.2e})"
print(f"\nSTAGE 2 PASS: c(d) = 4(d+2)(d-1)/(3 d^2) EXACTLY for d=2..7 (worst {worst:.1e}).")
print("  Reading: c = (4/3)*(1 + 1/d - 2/d^2). The '4/3' is the BASE (the d=2 value AND the d->oo limit);"
      "\n  the finite-d correction (d+2)(d-1)/d^2 PEAKS at d=4 (= 2^2; c(4)=3/2) and decays back to 4/3.")
print("\n  VERDICT: F94's 4/3 is NOT family A (d^2 -> 9/3=3). It is a SETUP-SPECIFIC coefficient with its own"
      "\n  qudit closed form 4(d+2)(d-1)/(3 d^2); the d=2 value 4/3 is the qubit point of that curve, and the"
      "\n  '4/3 = a_-1/3' discriminant reading is a d=2 coincidence ((J/4)=(1/2)^2=1/d^2 only at d=2; the"
      "\n  faithful dynamics is the d-INDEPENDENT (J/2)*SWAP). The proof's own universality remark is correct:"
      "\n  the coefficient counts setup-specific surviving diagrams, it is not the squared-dimension discriminant.")
print("\nDONE.")
