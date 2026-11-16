"""A post-flight reading of where arm B's joint failures can come from: the island, and what reaches it.

This is a reading without a verdict, written after the flight. None of its numbers is a pre-registered quantity of
experiments/X3_WINDOW_HARDWARE_PREDICTION.md or experiments/DRIVE_FRAME_HARDWARE_PREDICTION.md. It follows
drive_frame_colouring_reading.py, which found on both chips an excess of shots in which all three X outcomes of B's
pooled XX circuits fail together (signs normalised to the ideal |+,-,+>: U0 = X0, U1 = -X1, U2 = X2).

1. The island, exactly. Arm B's Hamiltonian is YZ + ZY on each bond. Written in the X basis (a Hadamard on every site
   carries Y to -Y and Z to X), each bond term flips an aligned pair (++ or --) and annihilates an anti-aligned one.
   The ideal configuration (+,-,+) and the all-fail configuration (-,+,-) are anti-aligned on both bonds, so the
   Hamiltonian has no matrix element into or out of either. The check is on Gaussian-integer matrices and is exact;
   the connected components of the eight configurations under H are printed beside it.
2. The island in the compiled circuits. Every one-site Pauli X, Y or Z is inserted after every gate of every B
   variant's XX circuit on both chips, in the ideal state-vector circuit, and the probability of the all-fail outcome
   is read. The same is done with letters on two and on three sites at one slot, with two one-site letters after two
   different gates (random draws on B50_XX), and, on B50_XX, with two-site letters sorted by where they act: inside a
   bond factor (from its first cz to its second) on that bond, inside it on a pair with the third qubit off the bond,
   or between factors.
3. A common Gaussian quasi-static Z field, through the colouring script's replay (T1 and T2 divided by 7) at standard
   deviations 0.2 to 0.8 rad/us: the sign of the third cumulant and the rate of all three failing it gives.
4. Arm B's pooled XX distribution against the local replay mixed with a fraction w of shots of one of five kinds:
   uniformly scrambled; all three losing their phase, independently (T2 = 0.05 us in those shots, T1 kept); all three
   turned by one common random detuning (uniform on [-6, 6] rad/us, 25 points); a common static detuning of 5 rad/us;
   a static detuning of 3 rad/us on the outer pair. The divisor of T1 and T2 runs over a grid and w is fitted. As in
   the colouring script, every B variant is read against the replay of B50_XX.
5. Arm Q's idle neighbours. In each Q circuit one qubit is active and the other two idle in |0> for up to about
   4 us; the rate at which both idle qubits read 1, against the product of their rates, raw counts.
6. Arm A, the bases with an even number of Y (XX, XZ, ZX, ZZ, YY; the odd ones carry the drive-frame page's odd-Y
   excess), 30 circuits: the local replay alone, each kind of item 4 at B's fitted w, and at a w fitted to A.

chi2 everywhere uses the multinomial covariance of the counts carried through the joint 8x8 inverse of the RO arm
(shot noise only; no RO drift factor, no uncertainty of the replays). A mixture with fixed component shapes and two
parameters is a reading, not a test (the colouring script's second caution); items 5 and 6 are the second routes.

    python drive_frame_island_reading.py

It reads the two counts files of data/ibm_drive_frame_november2026 (Marrakesh, Kingston), prints to stdout and writes
nothing. The output is kept in simulations/results/drive_frame/island_reading.txt.
"""
from __future__ import annotations

import json
import math
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import minimize_scalar

HERE = Path(__file__).resolve().parent
REPO = HERE.parent.parent
sys.path.insert(0, str(HERE))
import run_drive_frame as R  # noqa: E402  read-only use: Chip, compact, readout_matrix, replay_fn, ROLE, rl
import drive_frame_colouring_reading as CR  # noqa: E402  read-only use: U, readings, replay, quasi_static

rl = R.rl
DATA = REPO / "data" / "ibm_drive_frame_november2026"
FILES = {"ibm_marrakesh": "drive_frame_ibm_marrakesh_20261113_203627.json",
         "ibm_kingston": "drive_frame_ibm_kingston_20261113_203653.json"}
SEED = 20261116
N_PAIRS = 5000
DIVS = (3, 4, 5, 5.5, 6, 6.5, 7, 7.5, 8, 9, 10)
EVEN_Y = ("XX", "XZ", "ZX", "ZZ", "YY")
KEYS = ("U0", "U1", "U2", "cov01", "cov12", "cov02", "k3", "allfail")
IDEAL, ALL_FAIL = 1, 6   # runner index k = 4 q0 + 2 q2 + q1 (CR.U): all three read +1 at 1 and -1 at 6
assert all(CR.U[IDEAL] == 1) and all(CR.U[ALL_FAIL] == -1)

PAULI = {"I": np.eye(2, dtype=complex), "X": np.array([[0, 1], [1, 0]], dtype=complex),
         "Y": np.array([[0, -1j], [1j, 0]]), "Z": np.diag([1, -1]).astype(complex)}


def kron3(a, b, c):
    return np.kron(np.kron(PAULI[a], PAULI[b]), PAULI[c])


# ---------------------------------------------------------------- 1. the island, exactly
def island_exact(say):
    """Bond term Y_k Z_l + Z_k Y_l in the X basis is -(Y_k X_l + X_k Y_l); site order (0, 1, 2), bit 1 = '-'."""
    def bond(k, l):
        s = [["I"] * 3, ["I"] * 3]
        s[0][k], s[0][l] = "Y", "X"
        s[1][k], s[1][l] = "X", "Y"
        return -(kron3(*s[0]) + kron3(*s[1]))
    terms = {"(0,1)": bond(0, 1), "(1,2)": bond(1, 2)}
    for M in terms.values():                               # entries are exact Gaussian integers
        assert np.all(M.real == np.round(M.real)) and np.all(M.imag == np.round(M.imag))
    idx = lambda c: int("".join("1" if x == "-" else "0" for x in c), 2)
    x, xbar = idx("+-+"), idx("-+-")
    for name, M in terms.items():
        dead = [c for c in ("+-+", "-+-", "+--", "--+", "-++", "++-") if
                (c[int(name[1])] != c[int(name[3])]) and np.all(M[:, idx(c)] == 0)]
        anti = [c for c in ("+-+", "-+-", "+--", "--+", "-++", "++-") if c[int(name[1])] != c[int(name[3])]]
        say(f"  bond {name}: annihilates every configuration anti-aligned on it: {dead == anti} ({len(anti)} of 8)")
    Hx = sum(terms.values())
    ok = all(Hx[y, c] == 0 and Hx[c, y] == 0 for c in (x, xbar) for y in range(8))
    say(f"  H = sum of the two bond terms: every matrix element into and out of (+,-,+) and (-,+,-) is exactly 0: {ok}")
    assert ok
    names = ["".join("-" if (k >> (2 - i)) & 1 else "+" for i in range(3)) for k in range(8)]
    seen, comps = set(), []
    for k in range(8):
        if k in seen:
            continue
        stack, comp = [k], set()
        while stack:
            u = stack.pop()
            if u in comp:
                continue
            comp.add(u)
            stack += [w for w in range(8) if Hx[w, u] != 0 and w not in comp]
        seen |= comp
        comps.append(sorted(names[u] for u in comp))
    say(f"  connected components of the X-basis configurations under H: {comps}")
    return ok


# ---------------------------------------------------------------- 2. the island in the compiled circuits
def op_matrix(o):
    name, qs, prm, _ = o
    if name == "rz":
        return rl.on(qs[0], rl.rz(prm[0]))
    if name == "sx":
        return rl.on(qs[0], rl.rx(math.pi / 2))
    if name == "x":
        return rl.on(qs[0], rl.rx(math.pi))
    if name == "cz":
        return rl.czm(*sorted(qs))
    return None                                            # delay


def ideal_probs(ops, inject=()):
    """Ideal state-vector circuit with Pauli letters inserted after gate i: inject = [(i, qubit, letter), ...].
    Returns outcome probabilities in the runner's index."""
    psi = np.zeros(8, dtype=complex)
    psi[0] = 1
    for i, o in enumerate(ops):
        M = op_matrix(o)
        if M is not None:
            psi = M @ psi
        for (j, q, L) in inject:
            if j == i:
                psi = rl.on(q, PAULI[L]) @ psi
    return rl.probs_rho(np.outer(psi, psi.conj()))           # the runner's index, as run_drive_frame.logical_probs


def reach(ops, rng=None):
    slots = [i for i, o in enumerate(ops) if o[0] != "delay"]
    p0 = ideal_probs(ops)
    one = max(ideal_probs(ops, [(i, q, L)])[ALL_FAIL] for i in slots for q in range(3) for L in "XYZ")
    two = max(ideal_probs(ops, [(i, a, La), (i, b, Lb)])[ALL_FAIL]
              for i in slots for a, b in ((0, 1), (1, 2), (0, 2)) for La in "XYZ" for Lb in "XYZ")
    three = max(ideal_probs(ops, [(i, 0, L[0]), (i, 1, L[1]), (i, 2, L[2])])[ALL_FAIL]
                for i in slots for L in ("ZZZ", "YZY", "ZYZ", "YYY", "XXX"))
    return p0, one, two, three


def island_compiled(say):
    rng = np.random.default_rng(SEED)
    for b in FILES:
        ch = R.Chip(b)
        for v in ch.variants("B"):
            ops = ch.ops[f"{v}_XX"]
            p0, one, two, three = reach(ops)
            assert p0[IDEAL] > 1 - 1e-12 and np.array_equal(p0, R.logical_probs(ops))
            say(f"  {b:13s} {v}_XX: ideal all-fail {p0[ALL_FAIL]:.1e}; max over one-site letters at every slot "
                f"{one:.1e}; two-site letters at one slot {two:.3f}; three-site {three:.3f}")
        ops = ch.ops["B50_XX"]
        slots = [i for i, o in enumerate(ops) if o[0] != "delay"]
        hits, mx = 0, 0.0
        for _ in range(N_PAIRS):
            i, j = sorted(rng.choice(slots, 2, replace=False))
            (qa, qb), (La, Lb) = rng.integers(0, 3, 2), rng.choice(list("XYZ"), 2)
            a = ideal_probs(ops, [(i, int(qa), str(La)), (j, int(qb), str(Lb))])[ALL_FAIL]
            mx, hits = max(mx, a), hits + (a > 0.5)
        say(f"  {b:13s} B50_XX, two one-site letters after two different gates ({N_PAIRS} random draws): max "
            f"{mx:.3f}, draws above 0.5: {hits}")
        say(f"  {b:13s} B50_XX, two-site letters by where they act (max over letters and slots):")
        for line in slot_map(ops):
            say("      " + line)


def factors(ops):
    """Bond factors of the compiled circuit: (first cz index, second cz index, bond); consecutive cz on one bond."""
    cz = [(i, tuple(sorted(o[1]))) for i, o in enumerate(ops) if o[0] == "cz"]
    out = []
    for (i, b1), (j, b2) in zip(cz[0::2], cz[1::2]):
        assert b1 == b2
        out.append((i, j, b1))
    return out


def slot_map(ops):
    fac = factors(ops)
    inside = {}
    for i, j, bond in fac:
        for k in range(i, j):
            inside[k] = bond
    best = {}
    for k, o in enumerate(ops):
        if o[0] == "delay":
            continue
        for pair in ((0, 1), (1, 2), (0, 2)):
            if k in inside:
                cls = "inside a factor, on its own bond" if pair == inside[k] else "inside a factor, with the third qubit"
            else:
                cls = "between factors"
            a = max(ideal_probs(ops, [(k, pair[0], La), (k, pair[1], Lb)])[ALL_FAIL] for La in "XYZ" for Lb in "XYZ")
            best[cls] = max(best.get(cls, 0.0), a)
    return [f"{len(fac)} factors; {c}: {v:.1e}" for c, v in sorted(best.items())]


# ---------------------------------------------------------------- shared readers
def load_chip(b):
    d = json.load(open(DATA / FILES[b], encoding="utf-8"))
    ch = R.Chip(d["backend"])
    comp = R.compact(d["counts"], ch.nf)
    A_ro, _ = R.readout_matrix(comp)
    return ch, comp, np.linalg.inv(A_ro)


def pooled(comp, Ai, names):
    pc = sum(comp[n][0].sum(0) for n in names)
    N = pc.sum()
    pr = pc / N
    C = Ai @ ((np.diag(pr) - np.outer(pr, pr)) / N) @ Ai.T
    return Ai @ pr, np.linalg.pinv(C[:7, :7])


def clean(p):
    p = np.clip(p, 0, None)
    return p / p.sum()


def replay(ch, ops, f, L=None):
    T1 = [t / f for t in ch.T1]
    T2 = [t / f for t in ch.T2]
    return clean(rl.probs_rho(R.replay_fn(L)(ops, R.ROLE, T1, T2, L)))


KINDS = ("scrambled", "phase lost, independent", "common random detuning", "common detuning 5", "outer pair 3")


def component(ch, ops, f, kind):
    T1 = [t / f for t in ch.T1]
    T2 = [t / f for t in ch.T2]
    if kind == "scrambled":
        return np.ones(8) / 8
    if kind == "phase lost, independent":
        return clean(rl.probs_rho(R.replay_fn(None)(ops, R.ROLE, T1, [0.05] * 3, None)))
    if kind == "common random detuning":
        grid = np.linspace(-6, 6, 25)
        return sum(replay(ch, ops, f, {("det", R.ROLE[l]): x for l in range(3)}) for x in grid) / len(grid)
    pat, dd = {"common detuning 5": ((1, 1, 1), 5.0), "outer pair 3": ((1, 0, 1), 3.0)}[kind]
    return replay(ch, ops, f, {("det", R.ROLE[l]): c * dd for l, c in enumerate(pat)})


def fit_w(chi):
    best = None
    for f in DIVS:
        res = minimize_scalar(lambda w: chi(f, w), bounds=(0, 0.3), method="bounded", options={"xatol": 1e-5})
        if best is None or res.fun < best[0]:
            best = (res.fun, f, res.x)
    return best


# ---------------------------------------------------------------- 3. common Gaussian field
def gaussian_field(say):
    for b in FILES:
        ch = R.Chip(b)
        T1r, T2r = [t / 7 for t in ch.T1], [t / 7 for t in ch.T2]
        ref = CR.readings(CR.replay(ch, T1r, T2r))
        for s in (0.2, 0.4, 0.6, 0.8):
            r = CR.readings(CR.quasi_static(ch, T1r, T2r, (1, 1, 1), s))
            say(f"  {b:13s} s {s} rad/us: " + " ".join(f"{k} {r[k] - ref[k]:+.4f}" for k in KEYS))


# ---------------------------------------------------------------- 4 and 6. mixtures on arms B and A
def mixtures(say):
    out = {}
    for b in FILES:
        ch, comp, Ai = load_chip(b)
        pmB, WB = pooled(comp, Ai, [f"{v}_XX" for v in ch.variants("B")])
        Acirc = [f"{v}_{bs}" for v in ch.variants("A") for bs in EVEN_Y]
        A = {n: pooled(comp, Ai, [n]) for n in Acirc}
        cache = {}

        def parts(n, f, kind=None):
            key = (n, f, kind)
            if key not in cache:
                cache[key] = replay(ch, ch.ops[n], f) if kind is None else component(ch, ch.ops[n], f, kind)
            return cache[key]

        def chiB(f, w, kind):
            p = (1 - w) * parts("B50_XX", f) + w * parts("B50_XX", f, kind) if w else parts("B50_XX", f)
            r = (pmB - p)[:7]
            return float(r @ WB @ r)

        def chiA(f, w, kind):
            s = 0.0
            for n in Acirc:
                p = (1 - w) * parts(n, f) + w * parts(n, f, kind) if w else parts(n, f)
                r = (A[n][0] - p)[:7]
                s += float(r @ A[n][1] @ r)
            return s

        cB0, fB0 = min((chiB(f, 0, None), f) for f in DIVS)
        cA0, fA0 = min((chiA(f, 0, None), f) for f in DIVS)
        meas = CR.readings(pmB)
        say(f"\n  {b}")
        say(f"    measured B   " + " ".join(f"{k} {meas[k]:+.4f}" for k in KEYS) +
            f"   cov02|U1- {meas['cov02|U1-']:+.4f}")
        say(f"    B, local replay alone: chi2 {cB0:.1f} on 6 (divisor {fB0});   "
            f"A even-Y, local replay alone: chi2 {cA0:.1f} on {7 * len(Acirc)} (divisor {fA0}), "
            f"per dof {cA0 / (7 * len(Acirc)):.2f}")
        for kind in KINDS:
            cB, fB, wB = fit_w(lambda f, w: chiB(f, w, kind))
            m = CR.readings((1 - wB) * parts("B50_XX", fB) + wB * parts("B50_XX", fB, kind))
            cAfix = min(chiA(f, wB, kind) for f in DIVS)
            cA, fA, wA = fit_w(lambda f, w: chiA(f, w, kind))
            say(f"    + {kind:24s} B: chi2 {cB:7.1f} on 5, w {wB:.4f} (divisor {fB}) | A at B's w: chi2 {cAfix:7.1f} "
                f"({cAfix - cA0:+.1f}) | A's own w {wA:.4f}: chi2 {cA:7.1f} ({cA - cA0:+.1f})")
            say(f"      model B   " + " ".join(f"{k} {m[k]:+.4f}" for k in KEYS) + f"   cov02|U1- {m['cov02|U1-']:+.4f}")
            out[(b, kind)] = (cB, wB, cAfix - cA0, wA, cA - cA0)
    return out


# ---------------------------------------------------------------- 5. arm Q's idle neighbours
def bit(k, q):
    return (k >> {0: 2, 2: 1, 1: 0}[q]) & 1


def q_idle(say):
    for b in FILES:
        ch, comp, _ = load_chip(b)
        both = prod = n_all = 0.0
        worst = None
        for act in range(3):
            idle = [q for q in range(3) if q != act]
            for dly in range(4):
                pc = comp[f"Q_{act}X{dly}"][0].sum(0) + comp[f"Q_{act}Y{dly}"][0].sum(0)
                N = pc.sum()
                ma = sum(pc[k] for k in range(8) if bit(k, idle[0])) / N
                mc = sum(pc[k] for k in range(8) if bit(k, idle[1])) / N
                jj = sum(pc[k] for k in range(8) if bit(k, idle[0]) and bit(k, idle[1])) / N
                z = (jj - ma * mc) / math.sqrt(ma * mc / N) if ma * mc > 0 else 0.0
                both += jj * N
                prod += ma * mc * N
                n_all += N
                if worst is None or abs(z) > abs(worst[0]):
                    worst = (z, act, dly)
        say(f"  {b:13s} {int(n_all)} shots over 12 Q circuit pairs: both idle qubits read 1 in {int(round(both))} shots "
            f"({both / n_all:.6f}), the product of their rates gives {prod:.1f} ({prod / n_all:.6f}); largest single excess "
            f"{worst[0]:+.1f} sigma (active q{worst[1]}, delay {worst[2]})")


def main():
    lines = []

    def say(s):
        print(s, flush=True)
        lines.append(s)

    say("-- 1. the island in the X basis (exact)")
    island_exact(say)
    say("\n-- 2. the island in the compiled B circuits (ideal state vector; probability of all three failing)")
    island_compiled(say)
    say("\n-- 3. a common Gaussian quasi-static Z field, excess over the replay at T1 and T2 / 7 (B50_XX)")
    gaussian_field(say)
    say("\n-- 5. arm Q: idle neighbours, raw counts")
    q_idle(say)
    say("\n-- 4 and 6. arm B pooled XX and arm A even-Y bases: local replay mixed with a fraction w of shots of one kind")
    mixtures(say)


if __name__ == "__main__":
    main()
