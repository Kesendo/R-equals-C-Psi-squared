"""Gate for experiments/GAMMA_CONTROL.md (the 2026-08-23 two-lever repair).

Recomputes, through the committed N=5 primitives in mediator_bridge.py,
every number the repaired page's tables assert: the eight March profiles
(cross-checked against the March values), the nine matched-budget rows at
total 0.13, the one-site and multi-site DD rows, the one-factor mediator
controls, the lever crossing (all-of-0.25-on-centre beating uniform-0.13),
and the switch-control bound (a doubled centre rate imprints at most
0.013 in absolute MI at any lag out to ten units).

Exact-route policy: the propagator is an eigensolver-free expm route with
its own convergence, so values are gated at 1e-5 absolute against the
recorded six-decimal numbers (the two independent 2026-08-23 computations
agreed to all printed digits); ORDERINGS and the crossing are gated
strictly, with no tolerance.

Run: python simulations/gamma_control_two_lever_check.py  (~2 min)
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mediator_bridge as mb

TOL = 1e-5
checks = []

def mi(gammas, t=5.0, rho0=None):
    if rho0 is None:
        rho0 = mb.make_initial_state('bell_A_0M_pp_B')
    L, _, _ = mb.build_mediator_system(gammas=list(gammas))
    return mb.mutual_info(mb.evolve_rho(L, rho0, t), 5, [0, 1], [3, 4])

def gate(name, got, want):
    ok = abs(got - want) < TOL
    checks.append(ok)
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}: {got:.6f} vs {want:.6f}")

def gate_true(name, cond):
    checks.append(cond)
    print(f"  [{'PASS' if cond else 'FAIL'}] {name}")

print("March profiles (using_gamma.py set, values from using_gamma.txt)")
march = {
    "gate-closed":  ([.05, .05, .20, .05, .05], 0.192031),
    "uniform":      ([.05, .05, .05, .05, .05], 0.337552),
    "gate-open":    ([.05, .05, .01, .05, .05], 0.415907),
    "inv-V":        ([.05, .03, .01, .03, .05], 0.530100),
    "quiet-recv":   ([.05, .05, .05, .01, .01], 0.569031),
    "ramp-down":    ([.05, .04, .03, .02, .01], 0.628134),
    "ramp-up":      ([.01, .02, .03, .04, .05], 0.630704),
    "V-shape":      ([.01, .03, .05, .03, .01], 0.754963),
}
for name, (g, want) in march.items():
    gate(name, mi(g), want)

print("Matched budget 0.13 (the assignment axis)")
r = 13.0 / 17.0
matched = {
    "all on mediator": ([0, 0, .13, 0, 0], 1.039103),
    "sharp V":         ([0, .015, .10, .015, 0], 0.887763),
    "all on site 3":   ([0, 0, 0, .13, 0], 0.809086),
    "all on site 1":   ([0, .13, 0, 0, 0], 0.789445),
    "V-shape":         ([.01, .03, .05, .03, .01], 0.754963),
    "all on site 4":   ([0, 0, 0, 0, .13], 0.741572),
    "all on site 0":   ([.13, 0, 0, 0, 0], 0.740312),
    "uniform 0.026":   ([.026] * 5, 0.709348),
    "inv-V rescaled":  ([.05 * r, .03 * r, .01 * r, .03 * r, .05 * r], 0.689751),
}
vals = {}
for name, (g, want) in matched.items():
    vals[name] = mi(g)
    gate(name, vals[name], want)
gate_true("every single-site concentration beats uniform",
          all(vals[k] > vals["uniform 0.026"] for k in
              ("all on mediator", "all on site 0", "all on site 1",
               "all on site 3", "all on site 4")))

print("DD rows (x0.1 on treated sites, base 0.05)")
dd = {
    "site 0 alone": ([.005, .05, .05, .05, .05], 0.459097),
    "site 1 alone": ([.05, .005, .05, .05, .05], 0.448432),
    "mediator":     ([.05, .05, .005, .05, .05], 0.427889),
    "site 3 alone": ([.05, .05, .05, .005, .05], 0.441061),
    "site 4 alone": ([.05, .05, .05, .05, .005], 0.459422),
    "receiver 3+4": ([.05, .05, .05, .005, .005], 0.612053),
    "M+receiver":   ([.05, .05, .005, .005, .005], 0.784021),
    "everywhere":   ([.005] * 5, 1.534998),
}
dv = {}
for name, (g, want) in dd.items():
    dv[name] = mi(g)
    gate(name, dv[name], want)
gate_true("positional ranking: ends > sites 1/3 > mediator",
          dv["site 0 alone"] > dv["site 1 alone"] > dv["mediator"]
          and dv["site 4 alone"] > dv["site 3 alone"] > dv["mediator"])
gate_true("ends tie within 0.1%",
          abs(dv["site 0 alone"] - dv["site 4 alone"]) / dv["site 4 alone"] < 1e-3)

print("One-factor mediator controls")
gate("quiet edges, quiet mediator", mi([.01, .03, .01, .03, .01]), 0.935734)
gate("loud edges, loud mediator", mi([.05, .03, .05, .03, .05]), 0.428428)

print("The lever crossing")
gate("all-on-M at 0.25", mi([0, 0, .25, 0, 0]), 0.746888)
gate_true("0.25 concentrated beats 0.13 uniform",
          mi([0, 0, .25, 0, 0]) > vals["uniform 0.026"])
gate_true("within-shape monotone (all-on-M 0.13 > 0.25 > 0.40)",
          vals["all on mediator"] > mi([0, 0, .25, 0, 0]) > mi([0, 0, .40, 0, 0]))

print("Switch control (doubled centre rate at t=5, lag scan)")
rho0 = mb.make_initial_state('bell_A_0M_pp_B')
base, sw = [.05] * 5, [.05, .05, .10, .05, .05]
Lb, _, _ = mb.build_mediator_system(gammas=base)
Ls, _, _ = mb.build_mediator_system(gammas=sw)
r5 = mb.evolve_rho(Lb, rho0, 5.0)
diffs = []
for k in range(1, 21):
    dt = 0.5 * k
    a = mb.mutual_info(mb.evolve_rho(Ls, r5, dt), 5, [0, 1], [3, 4])
    b = mb.mutual_info(mb.evolve_rho(Lb, r5, dt), 5, [0, 1], [3, 4])
    diffs.append(abs(a - b))
gate_true(f"imprint <= 0.013 absolute MI at every lag (max {max(diffs):.4f})",
          max(diffs) <= 0.013)

n_ok = sum(checks)
print(f"\n{n_ok}/{len(checks)} gates pass")
sys.exit(0 if n_ok == len(checks) else 1)
