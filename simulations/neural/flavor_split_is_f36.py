"""The neural flavor split is F36's mirror pairing, and it has a closed form.

The generator of FLAVOR_RESOLVED_T2_INHERITANCE's neural reading (neural_flavor_rule.build_jacobian, imported here,
not copied) is J = [[-gamma_E I + c A, -h I], [h I, -gamma_I I - c A]] on a graph adjacency A scaled to spectral
radius one. With Q the per-node E <-> I swap and s = (gamma_E + gamma_I)/2 it satisfies F36 identically,
Q J Q = -J - 2 s I, at every h, c, A and every pair of leaks: the leak condition d_i + d_Q(i) + 2s = 0 holds because
every E seat pairs with its own I seat, and the off-diagonal condition W[Q(i), Q(j)] = -W[i, j] holds because the
coupling enters as +cA on E and -cA on I while the h block is antisymmetric. Entry for entry the identity is a sign
flip plus one sum, gamma_E + gamma_I against 2 s, the same rounded sum on both sides, so the float residual is
exactly 0.0 at generic couplings too. Control: +cA on both blocks breaks the off-diagonal condition; a wrong s
breaks the leak condition.

CLOSED FORM. Per eigenvalue a of A (eigenvector w), the vectors (w, 0) and (0, w) span an invariant 2 x 2 block
-s I + [[Delta + c a, -h], [h, -Delta - c a]] with Delta = (gamma_I - gamma_E)/2, eigenvalues
lambda = -s +- sqrt((Delta + c a)^2 - h^2). A block is OVERDAMPED (two real rates, mirrored about s) iff
|Delta + c a| > h; its slow mode is E-dominant iff Delta + c a > 0. An underdamped block sits on Re lambda = -s with
E weight exactly one half ("mixed"). So the flavor split is the mirror pairing read on the overdamped blocks:
every E-dominant mode at rate r has an I-dominant partner at 2s - r. The experiment's ratio takes the slowest
E-dominant rate (the block a_max) against the slowest I-dominant rate, which is the FASTER rate s + sqrt(.) of the
weakest overdamped block, not the partner of the slowest mode; every row of neural_flavor_rule.txt and of the h = 0
control file is reproduced here from the closed form, with the classifier's own condition (the weakest block's fast
mode carries I weight above 0.65) read beside it. At h = 0 the reported ratio is (gamma_I + c a_min)/(gamma_E - c a_max),
2.333 on chain, ring and star (a_min = -1, a_max = 1), 2.6 on the complete graph (a_min = -1/5), 2 on the empty graph
(a = 0). The graph enters through which a are
overdamped, not through the slowest E rate, which is the same 0.940983 on every connected graph (a_max = 1).

What this bounds: F36 says which rates pair and that the weights exchange; which flavor is slow is the sign of
Delta + c a (control: gamma_E = 1, gamma_I = 2, h = 0.1, c = 1 on the chain has I-dominant modes slower than s);
and at equal leaks (Delta = 0) the "classifier-failed" reading of the experiment's control holds when h > c
(every block underdamped, every weight one half), while c > h gives overdamped blocks at equal leaks too; there the
blocks a and -a share their eigenvalues, so on a bipartite graph each overdamped rate is carried by one E- and one
I-dominant mode at once and the eigensolver's vectors inside that pair are basis-arbitrary (only a non-bipartite
graph splits genuinely). A block with Delta + c a = +-h is defective, a Jordan block at -s with eigenvector weight
one half; at the experiment's parameters that is every a = 0 block (the empty graph, the star, odd chains, rings
with 4 | N), read at the sqrt(eps) scale and not gated.
"""
import os
import sys
import numpy as np
from scipy.optimize import linear_sum_assignment
sys.stdout.reconfigure(encoding="utf-8")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import neural_flavor_rule as nf          # the producer; its __main__ guard keeps the import side-effect free

EPS = np.finfo(float).eps
FAILS = []
def check(name, ok, detail=""):
    print(f"  [{'ok' if ok else 'FAIL'}] {name}  {detail}")
    if not ok: FAILS.append(name)
def swap(N): return np.block([[np.zeros((N, N)), np.eye(N)], [np.eye(N), np.zeros((N, N))]])
def e_weight(v, N): return float(np.sum(np.abs(v[:N]) ** 2) / np.sum(np.abs(v) ** 2))
def bottleneck(ev, target):
    cost = np.abs(ev[:, None] - target[None, :]); r, c = linear_sum_assignment(cost); return float(cost[r, c].max())

print("=== F36 on the producer's generator, exactly, at generic couplings ===")
rng = np.random.default_rng(5)
for topo in ("chain", "ring", "star", "complete"):
    for trial in range(3):
        N = int(rng.integers(4, 9)); gE, gI, h, c = rng.uniform(0.1, 2.0, 4)
        J = nf.build_jacobian(N, topo, h, c, gE, gI); Q = swap(N); s = (gE + gI) / 2
        res = float(np.max(np.abs(Q @ J @ Q + J + 2 * s * np.eye(2 * N))))
        check(f"{topo:9s} N={N} generic (gamma_E, gamma_I, h, c): Q J Q + J + 2 s I == 0.0", res == 0.0, f"max {res}")
N = 6; gE, gI, h, c = 0.37, 0.91, 0.53, 0.29; Q = swap(N); s = (gE + gI) / 2
A = nf.adjacency(N, "ring"); I = np.eye(N)
Jb = np.block([[-gE * I + c * A, -h * I], [h * I, -gI * I + c * A]])
check("control, +cA on both blocks (off-diagonal condition broken): F36 residual != 0", float(np.max(np.abs(Q @ Jb @ Q + Jb + 2 * s * np.eye(2 * N)))) != 0.0)
J = nf.build_jacobian(N, "ring", h, c, gE, gI)
check("control, wrong centre s' = gamma_E (leak condition broken): residual != 0", float(np.max(np.abs(Q @ J @ Q + J + 2 * gE * np.eye(2 * N)))) != 0.0)

print("\n=== the closed form per adjacency eigenvalue, against the eigensolver (law: ratio to eps x spectral scale; exceptional rows read) ===")
def closed_form(N, topo, h, c, gE, gI):
    a = np.linalg.eigvalsh(nf.adjacency(N, topo)); s = (gE + gI) / 2; D = (gI - gE) / 2
    disc = (D + c * a) ** 2 - h * h
    return a, disc, np.concatenate([-s + np.sqrt(disc + 0j), -s - np.sqrt(disc + 0j)])
for topo, N in (("chain", 4), ("chain", 6), ("ring", 5), ("star", 6), ("complete", 6)):
    for (gE, gI, h, c) in ((1.0, 2.0, 0.5, 0.25), tuple(rng.uniform(0.1, 2.0, 4))):
        a, disc, lam = closed_form(N, topo, h, c, gE, gI); J = nf.build_jacobian(N, topo, h, c, gE, gI); ev = np.linalg.eigvals(J)
        bn = bottleneck(ev, lam); model = EPS * float(np.max(np.abs(ev))); near_ep = float(np.min(np.abs(disc)))
        if near_ep > 1e-6:
            check(f"{topo:8s} N={N} (gE, gI, h, c) = ({gE:.3f}, {gI:.3f}, {h:.3f}, {c:.3f}): spectrum == closed form within the model", bn < 1e3 * model, f"{bn:.1e} = {bn / model:.1f} x model")
        else:
            print(f"       {topo} N={N}: a block sits at an exceptional point (min |disc| = {near_ep:.1e}); closed form vs eig {bn:.1e}, the sqrt(eps) scale of a defective pair (read)")

print("\n=== every row of the experiment's data files from the closed form (the files are parsed, not retyped) ===")
import re
def parse_rows(path):
    rows = []
    for line in open(path, encoding="utf-8"):
        m = re.match(r"\s*(\d+)\s+(\w+)\s+\d+\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+\S+\s+(\S+)", line)
        if m: rows.append((int(m[1]), m[2], m[3], m[4], m[5], m[6], m[7], m[8]))
    return rows
def fast_mode_i_weight(x, h):
    """I weight of the FAST mode of the overdamped block [[x, -h], [h, -x]] (eigenvalue -sqrt(x^2 - h^2) relative to -s)."""
    mu = -np.sqrt(x * x - h * h); v = np.array([h, x - mu]); return float(v[1] ** 2 / (v @ v))   # (x - mu) alpha = h beta
for path, (gE, gI, h, c) in (("simulations/results/neural_flavor_rule.txt", (1.0, 2.0, 0.5, 0.25)),
                            ("simulations/results/neural_flavor_rule_h0_control.txt", (1.0, 2.0, 0.0, 0.25))):
    s = (gE + gI) / 2; D = (gI - gE) / 2
    for N, topo, kE, rE, kI, rI, ratio, status in parse_rows(path):
        A = nf.adjacency(N, topo) if topo != "empty" else np.zeros((N, N)); a = np.linalg.eigvalsh(A); x = D + c * a; disc = x * x - h * h
        over = disc > 1e-12
        if status != "ok":
            check(f"{path.split('/')[-1]} N={N} {topo:8s}: '{status}' and the closed form has no overdamped block with positive x", not np.any(over & (x > 0)) or np.all(np.abs(disc) < 1e-12), f"disc {np.round(disc, 4)}")
            continue
        roots = np.sqrt(disc[over]); xs = x[over]
        slow_E = s - roots[xs > 0].max(); weakest = np.argmin(roots); slow_I = s + roots[weakest]
        w_I = fast_mode_i_weight(xs[weakest], h)
        check(f"{path.split('/')[-1]} N={N} {topo:8s}: E rate {slow_E:.6f} (file {float(rE):.6f}), I rate {slow_I:.6f} (file {float(rI):.6f}), ratio {slow_I / slow_E:.6f} (file {float(ratio):.6f}); fast mode of the weakest block carries I weight {w_I:.3f} > 0.65",
              abs(slow_E - float(rE)) < 1e-6 and abs(slow_I - float(rI)) < 1e-6 and abs(slow_I / slow_E - float(ratio)) < 1e-6 and w_I > 0.65)
    if h == 0.0:
        for N, topo, kE, rE, kI, rI, ratio, status in parse_rows(path):
            A = nf.adjacency(N, topo) if topo != "empty" else np.zeros((N, N)); a_min = float(np.linalg.eigvalsh(A).min()); a_max = float(np.linalg.eigvalsh(A).max())
            check(f"h = 0 {topo:8s}: ratio == (gamma_I + c a_min)/(gamma_E - c a_max) = {(gI + c * a_min) / (gE - c * a_max):.6f} (file {float(ratio):.6f})", abs((gI + c * a_min) / (gE - c * a_max) - float(ratio)) < 1e-6)

print("\n=== the pairing carries the flavor: at a simple partner eigenvalue the eigensolver's vector has the block weights exchanged ===")
for topo in ("chain", "ring", "star", "complete"):
    N = 6; gE, gI, h, c = 0.37, 0.91, 0.53, 0.29; s = (gE + gI) / 2
    J = nf.build_jacobian(N, topo, h, c, gE, gI); w, V = np.linalg.eig(J); Q = swap(N)
    worst_transport = 0.0; worst_weight = 0.0; n_simple = 0
    for k in range(2 * N):
        v = V[:, k]; u = Q @ v; lam_p = -w[k] - 2 * s
        worst_transport = max(worst_transport, float(np.max(np.abs(J @ u - lam_p * u))))
        j = int(np.argmin(np.abs(w - lam_p)))
        if np.sum(np.abs(w - w[j]) < 1e-8) == 1:
            n_simple += 1; worst_weight = max(worst_weight, abs(e_weight(V[:, j], N) - (1 - e_weight(v, N))))
    model = EPS * float(np.max(np.abs(w)))
    check(f"{topo:9s}: J (Q v) = (-lambda - 2s) (Q v) for every eigenvector (the eigensolver's own residual, law)", worst_transport < 1e3 * model, f"{worst_transport:.1e} = {worst_transport / model:.1f} x model")
    check(f"{topo:9s}: at the {n_simple} simple partner eigenvalues the eigensolver's vector has E weight 1 - E weight of v (law)", n_simple > 0 and worst_weight < 1e3 * model, f"max {worst_weight:.1e}")

print("\n=== which flavor is slow is the sign of Delta + c a, not F36 ===")
N = 6; J = nf.build_jacobian(N, "chain", 0.1, 1.0, 1.0, 2.0); w, V = np.linalg.eig(J); s = 1.5
rates = -w.real; ew = np.array([e_weight(V[:, k], N) for k in range(2 * N)])
slow_I = sorted(r for r, e in zip(rates, ew) if e < 0.5 and r < s)
check("gamma_E = 1, gamma_I = 2, h = 0.1, c = 1, chain N=6: I-dominant modes slower than s exist (Delta + c a < 0 for a < -1/2)", len(slow_I) > 0, f"rates {np.round(slow_I, 4)}")

print("\n=== equal leaks: Delta = 0; h > c gives every block underdamped and every E weight one half; c > h gives flavored blocks ===")
for h, c, expect_half in ((0.5, 0.25, True), (0.5, 1.0, False)):
    J = nf.build_jacobian(6, "chain", h, c, 1.0, 1.0); w, V = np.linalg.eig(J); ew = np.array([e_weight(V[:, k], 6) for k in range(12)])
    a = np.linalg.eigvalsh(nf.adjacency(6, "chain")); all_under = bool(np.all((c * a) ** 2 < h * h))
    half = bool(np.max(np.abs(ew - 0.5)) < 1e3 * EPS)
    check(f"equal leaks h={h}, c={c}: every block underdamped is {all_under}; every E weight == 1/2 within the model is {half}", all_under == expect_half and half == expect_half,
          "" if expect_half else "(the eigensolver's weights inside the degenerate +-a pairs of this bipartite chain are basis-arbitrary and not printed)")
if True:
    a = np.linalg.eigvalsh(nf.adjacency(6, "chain")); x = 1.0 * a
    check("equal leaks, c = 1 > h = 0.5 on the bipartite chain: the overdamped blocks come in +-a pairs with equal eigenvalues (one E- and one I-dominant mode per rate)",
          bool(np.all(np.abs(np.sort(x[x > 0.5]) - np.sort(-x[x < -0.5])) < 1e-12)) and np.sum(x > 0.5) > 0)

print(f"\n{'ALL OK' if not FAILS else 'FAILED: ' + ', '.join(FAILS)}")
sys.exit(1 if FAILS else 0)
