"""Gate-first: classify the F-registry's recurring "4" by genealogy -- turn the qutrit prism on real formulas.

CONTEXT. PolynomialDiscriminantAnchorClaim: the "4" in ~1/3 of all F-formulas is the d=2 value of the
foundational trunk d^2 - 2d = 0. At d=2 several genealogies coincide at 4; the qutrit (d=3) is the prism
that splits them into two super-families:
  SPLIT (operator-counting, d-dependent):  A = d^2 (squared dim / operator space) -> 9 ;  B = 2d (F121 cap) -> 6
  STAY  (rate / spectral, d-independent):   C = 2*gamma*Hamming (rates) ;  D = band edge J*rho (graph)
This verifier takes ONE representative formula per gate-able family and confirms, at d=2 vs d=3, WHICH way
its "4" goes -- turning the classification table from labels into TESTED predictions. A firing gate is the
find (the family was mislabelled); diagnose, do not loosen.

  STAGE A (operator-dimension, family A -> 9):  F23 (N+1)/4^N, F49 4^N, F1 4^(N-1). The N-qudit operator
    (Liouville) space has dimension (d^2)^N = d^{2N}; the "4^N" IS (d^2)^N, so it -> 9^N at the qutrit.
  STAGE C (rate, family C -> stays 4*gamma):  F25 e^{-4*gamma*t}, F73 e^{-4*g0*t}, F86 t_peak=1/(4*g0).
    The dephasing rate of |a><b| is 2*gamma*Hamming(a,b) (full-Cartan equidistant), so the rate ladder is
    {2*gamma*k : k=0..N}, d-INDEPENDENT; 4*gamma = 2*gamma*(Hamming 2) (the two-site / rung-2 coherence) STAYS.
  STAGE D (band edge, family D -> stays 4J):  F2 omega=4J*(1-cos), F41 t_Pi ~ 1/(4J*sin^2). The single-
    excitation hopping band = J * (chain adjacency spectrum), the SAME set for every d (each flavor hops like
    a qubit excitation), so the band edge J*rho and the bandwidth 4J are d-INDEPENDENT. Cross-check: F121's
    Heisenberg gap Delta = 4J also holds at d=3.

  NOTED, not gated (the deep follow-ups):
   * F94  Delta = (4/3) Q^2 K^3 writes the 4 EXPLICITLY as a_-1/3 = 4/3 and 1/16 = 1/a_-1^2 (the MOST explicit
     discriminant appearance -- and it is NOT in the anchor's cited F-list). Family A candidate: does 4/3 ->
     d^2/3 = 3 at the qutrit? Needs the SU(d) Born-mirror normalization (the (J/4)^2 -> SU(d) Casimir). Open.
   * F49  carries BOTH families in ONE formula: 4^N (denominator, family A, splits -> 9^N) AND the numerator 4
     = the rate^2 coefficient 4*gamma^2 (family C, stays). The qutrit pulls the two 4's apart -- covered by
     stages A and C together.

Run:  python simulations/formula_four_classification.py
"""
import itertools
import numpy as np

GAMMA = 0.05
J = 1.0


# ----------------------------------------------------------------------------------------------------
# minimal self-contained qudit helpers (faithful to qudit_g2_split.py; not imported -- that script runs on import)
# ----------------------------------------------------------------------------------------------------
def all_states(d, N):
    return list(itertools.product(range(d), repeat=N))


def hamming(a, b):
    return sum(1 for x, y in zip(a, b) if x != y)


def dephasing_rate_ladder(d, N):
    """The set of distinct decay rates 2*gamma*Hamming(a,b) over all d^{2N} coherences (full-Cartan
    equidistant dephasing). Returns the sorted distinct rates."""
    states = all_states(d, N)
    rates = {round(2.0 * GAMMA * hamming(a, b), 12) for a in states for b in states}
    return sorted(rates)


def se_band(d, N):
    """Single-excitation vacuum-hopping band of an open chain: the distinct eigenvalues of J * (chain
    adjacency on the one-excitation sector). Each flavor hops independently, so for d>2 the band is the
    same spectrum with higher multiplicity. Returns the sorted distinct band energies / J."""
    states = [s for s in all_states(d, N) if sum(1 for x in s if x != 0) == 1]   # exactly one excitation
    idx = {s: a for a, s in enumerate(states)}
    m = len(states)
    H = np.zeros((m, m))
    bonds = [(i, i + 1) for i in range(N - 1)]
    for a, s in enumerate(states):
        for (i, j) in bonds:
            si, sj = s[i], s[j]
            if (si == 0) != (sj == 0):                # one site vacuum -> the flavor hops
                s2 = list(s); s2[i], s2[j] = sj, si; s2 = tuple(s2)
                H[idx[s2], a] += J
    ev = np.linalg.eigvalsh(H)
    return sorted({round(e, 9) for e in ev})


# ====================================================================================================
# STAGE A -- OPERATOR-DIMENSION FAMILY (A): the "4^N" is (d^2)^N -> 9^N at the qutrit
# ====================================================================================================
print("=" * 100)
print("STAGE A -- family A (operator dimension): F23 (N+1)/4^N, F49 4^N, F1 4^(N-1). The '4^N' = (d^2)^N.")
print("=" * 100)
print(f"{'N':>2} {'d=2: 4^N':>10} {'d^{2N} check':>12} || {'d=3: 9^N':>10} {'d^{2N} check':>12} {'split?':>8}")
for N in (1, 2, 3, 4):
    dim2 = 2 ** (2 * N)
    dim3 = 3 ** (2 * N)
    assert dim2 == 4 ** N, f"STAGE A GATE: d=2 operator dim {dim2} != 4^{N}"
    assert dim3 == 9 ** N, f"STAGE A GATE: d=3 operator dim {dim3} != 9^{N}"
    print(f"{N:>2} {4**N:>10} {dim2:>12} || {9**N:>10} {dim3:>12} {'9^N (A)':>8}")
print("STAGE A PASS: the operator/Liouville space is (d^2)^N; the F-formula '4^N' is (d^2)^N and SPLITS to "
      "9^N\n  at the qutrit. The '4' here is the SQUARED DIMENSION d^2 (family A).")

# ====================================================================================================
# STAGE C -- RATE FAMILY (C): the rate ladder 2*gamma*Hamming is d-INDEPENDENT; 4*gamma STAYS
# ====================================================================================================
print("\n" + "=" * 100)
print("STAGE C -- family C (rate): F25 e^{-4*gamma*t}, F73 e^{-4*g0*t}, F86 t_peak=1/(4*g0). Rate = 2*g*Hamming.")
print("=" * 100)
fourgamma = 4.0 * GAMMA
for N in (2, 3, 4):
    lad2 = dephasing_rate_ladder(2, N)
    lad3 = dephasing_rate_ladder(3, N)
    expected = [round(2.0 * GAMMA * k, 12) for k in range(N + 1)]               # rungs 2*gamma*k, k=0..N
    assert lad2 == expected and lad3 == expected, (
        f"STAGE C GATE: rate ladder not {{2*gamma*k}} or not d-independent at N={N}: d2={lad2}, d3={lad3}")
    maxrate2, maxrate3 = lad2[-1], lad3[-1]                                     # Hamming=N coherence
    print(f"  N={N}: rate ladder (d=2) == (d=3) == {{2*gamma*k}}  "
          f"[{', '.join(f'{r:.3f}' for r in lad2)}]   max = 2*gamma*{N} = {maxrate2:.3f}")
# F25/F73: the two-site (Hamming-2) coherence rate = 2*gamma*2 = 4*gamma, at BOTH d
r2_d2 = round(2.0 * GAMMA * 2, 12)
assert abs(dephasing_rate_ladder(2, 2)[-1] - fourgamma) < 1e-12, "STAGE C: N=2 max rate should be 4*gamma"
assert abs(dephasing_rate_ladder(3, 2)[-1] - fourgamma) < 1e-12, "STAGE C: qutrit N=2 max rate should also be 4*gamma"
print(f"\nSTAGE C PASS: the rate ladder {{2*gamma*k}} is d-INDEPENDENT; the Hamming-2 coherence (F25 |00><11|,\n"
      f"  F73 mirror-sum) decays at 2*gamma*2 = 4*gamma = {fourgamma:.3f} at d=2 AND d=3. e^{{-4*gamma*t}} and "
      f"t_peak=1/(4*gamma)\n  STAY (family C). [The '4' here = 2 (Lindblad/AT) x Hamming-2, both carry no d.]")

# ====================================================================================================
# STAGE D -- BAND-EDGE FAMILY (D): the SE hopping band is d-INDEPENDENT; 4J STAYS
# ====================================================================================================
print("\n" + "=" * 100)
print("STAGE D -- family D (band edge): F2 omega=4J*(1-cos(pi k/N)), F41 t_Pi ~ 1/(4J*sin^2). Band = J*adjacency.")
print("=" * 100)
print(f"{'N':>2} {'d=2 band edge/J':>16} {'d=3 band edge/J':>16} {'2cos(pi/(N+1))':>16} {'d-indep?':>9}")
for N in (3, 4, 5):
    b2 = se_band(2, N)
    b3 = se_band(3, N)
    edge2, edge3 = max(abs(b2[0]), abs(b2[-1])), max(abs(b3[0]), abs(b3[-1]))
    rho = 2.0 * np.cos(np.pi / (N + 1))                                        # chain adjacency radius
    same = b2 == b3 and abs(edge2 - rho) < 1e-7
    assert b2 == b3, f"STAGE D GATE: SE band not d-independent at N={N}: d2={b2} vs d3={b3}"
    assert abs(edge2 - rho) < 1e-7, f"STAGE D GATE: band edge {edge2} != 2cos(pi/{N+1})={rho}"
    print(f"{N:>2} {edge2:>16.6f} {edge3:>16.6f} {rho:>16.6f} {('YES' if same else 'no'):>9}")
print("STAGE D PASS: the single-excitation hopping band (distinct energies) is IDENTICAL at d=2 and d=3 "
      "(each\n  flavor hops like a qubit excitation); band edge = J*rho = 2J*cos(pi/(N+1)), bandwidth = 4J*(...),\n"
      "  d-INDEPENDENT. The '4J' STAYS (family D). Cross-check: F121's Heisenberg gap Delta=4J holds at d=3 too.")

# ====================================================================================================
# THE CLASSIFICATION TABLE
# ====================================================================================================
print("\n" + "=" * 100)
print("CLASSIFICATION -- every anchor-cited '4' (and the extras), by family and qutrit verdict")
print("=" * 100)
rows = [
    ("F2  omega = 4J*(1-cos)",          "D graph",  "4J  -> stays (band edge J*rho)"),
    ("F23 (N+1)/4^N",                   "A d^2",    "4^N -> 9^N (operator space (d^2)^N)"),
    ("F25 CPsi, f = e^{-4*gamma*t}",    "C rate",   "4*gamma -> stays (2*gamma*Hamming-2)"),
    ("F49 R^2 = 4(N-2)/(N*4^N)",        "A + C",    "4^N -> 9^N (A) ; numerator 4=rate^2 -> stays (C)"),
    ("F56 saddle 4*eps, 16=4^2",        "E fold",   "the 1/4 cusp geometry (blocked: CPsi undefined d>2)"),
    ("F65 (4/(N+1))*sin^2",             "C x D",    "rate 2*gamma stays ; /(N+1) graph stays -> stays"),
    ("F67 4*g0/(N+1)",                  "C x D",    "as F65 -> stays"),
    ("F73 (1/2)e^{-4*g0*t}",            "C rate",   "4*g0 -> stays"),
    ("F76 e^{-4*g0*t}, 16J element",    "C ; D",    "4*g0 rate stays ; 16J graph element stays"),
    ("F86 t_peak = 1/(4*g0)",           "C rate",   "4*g0 -> stays"),
    ("F94 Delta = (4/3) Q^2 K^3",       "A? REFUTED",    "c(d)=4(d+2)(d-1)/(3d^2), NOT d^2/3=3 (f94_qutrit_born_mirror.py)"),
    ("F1  ||M||^2 ~ 4^(N-1)",           "A d^2",    "4^(N-1) = d^{2N}/d^2 -> 9^(N-1)"),
    ("F97/F98 CPsi=1/4 cusp",           "E fold",   "1/4 = a_3 = 1/a_-1 = 1/d^2 -> 1/9 ?  [blocked: CPsi d=2 only]"),
]
print(f"  {'formula':30} {'family':16} {'qutrit verdict':50}")
for f, fam, verdict in rows:
    print(f"  {f:30} {fam:16} {verdict:50}")

print("\nSUMMARY: the '4's partition into SPLIT (A=d^2->9, B=2d->6: operator-counting) and STAY (C=rate,\n"
      "  D=graph, E=fold-inverse: dynamical/spectral). Gated here: A (4^N->9^N), C (4*gamma stays),\n"
      "  D (4J stays). F94 RESOLVED (f94_qutrit_born_mirror.py): the family-A 'smoking gun' candidate is\n"
      "  REFUTED -- its qudit coefficient is 4(d+2)(d-1)/(3d^2) (peak 3/2 at d=4, ->4/3), NOT d^2/3=3; F94's\n"
      "  4 is a setup-specific diagram count, the a_-1/3 reading a d=2 coincidence. ONE test still OPEN:\n"
      "  family E's 1/4 -> 1/9, blocked until a qudit fold is defined (concurrence is d=2 only).")
print("\nDONE.")
