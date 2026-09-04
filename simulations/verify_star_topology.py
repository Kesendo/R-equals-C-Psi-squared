"""
STAR TOPOLOGY -- CLAIMS VERIFICATION
====================================
Verifies all key findings documented in STAR_TOPOLOGY_OBSERVERS.md
against star_topology_v2.py simulation.

Run: python verify_star_topology.py
Requires: numpy, star_topology_v2.py in same directory
"""
import sys
import numpy as np

# Windows' default console encoding cannot print this file's box-drawing and
# check marks, and the script died on its first print before any test ran. An
# outside reviewer found it in 2026-09; every session here had masked it by
# setting PYTHONIOENCODING=utf-8 out of habit, so the reproduction path this
# repository documents had been broken and invisible at once.
if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

# Import from star_topology_v2
from star_topology_v2 import (
    star_hamiltonian, dephasing_ops, make_state,
    rk4_step, ptrace, concurrence, psi_norm
)

def ab_cpsi_max(state="Bell_SA+B", J_SA=1.0, J_SB=1.0,
                gS=0.05, gA=0.05, gB=0.05, t_max=8.0, dt=0.001):
    """Run simulation, return maximum AB CΨ value."""
    H = star_hamiltonian(J_SA, J_SB)
    L = dephasing_ops(gS, gA, gB)
    rho = make_state(state)
    steps = int(t_max / dt)
    mx = 0.0
    for step in range(steps + 1):
        # Every step: the peak of CPsi is a supremum over continuous time, so
        # any skipped sample can only lower it and the thresholds read off it
        # come out too high. Sampling every 10th step at dt=0.005 gave 13 points
        # per oscillation period and was enough to move checks 3 and 4.
        if True:
            rp = ptrace(rho, 0)  # trace out S -> AB
            cpsi = concurrence(rp) * psi_norm(rp)
            if cpsi > mx:
                mx = cpsi
        if step < steps:
            rho = rk4_step(rho, H, L, dt)
    return mx

def check(name, condition):
    """Print check result, return success."""
    status = "✓" if condition else "✗"
    print(f"  {status} {name}")
    return condition

def main():
    print("╔══════════════════════════════════════════════╗")
    print("║   STAR TOPOLOGY v2 -- CLAIMS VERIFICATION    ║")
    print("╚══════════════════════════════════════════════╝\n")

    ok = 0
    total = 0

    def run(name, condition):
        nonlocal ok, total
        total += 1
        if check(name, condition):
            ok += 1

    # 1. Symmetric coupling: AB never crosses
    v = ab_cpsi_max(J_SB=1.0)
    run(f"Symmetric J=1,1: AB max={v:.4f} < 0.25", v < 0.25)

    # 2. Strong B: AB crosses
    v = ab_cpsi_max(J_SB=2.0)
    run(f"Strong B J=1,2: AB max={v:.4f} >= 0.25", v >= 0.25)

    # 3. Threshold at J_SB = 1.46295 (F29). The old bracket was 1.465/1.4655,
    #    which sits ABOVE the threshold on both sides and only read as a bracket
    #    because the coarse sampling pushed both peaks down.
    v1 = ab_cpsi_max(J_SB=1.4625)
    v2 = ab_cpsi_max(J_SB=1.4635)
    run(f"Threshold 1.46295: J=1.4625->{v1:.6f} NO, J=1.4635->{v2:.6f} YES",
        v1 < 0.25 and v2 >= 0.25)

    # 4. Receiver noise more destructive than sender noise. At a partner rate of
    #    0.05 the boundaries are gamma_A = 0.2699 and gamma_B = 0.4735, so 0.25
    #    is on the LIVING side for both and cannot separate them; 0.30 can.
    v_a = ab_cpsi_max(J_SB=2.0, gA=0.30, gB=0.05)
    v_b = ab_cpsi_max(J_SB=2.0, gA=0.05, gB=0.30)
    run(f"gamma_A=0.30 kills ({v_a:.4f}), gamma_B=0.30 doesn't ({v_b:.4f})",
        v_a < 0.25 and v_b >= 0.25)

    # 5. Only Bell state crosses, W and |0++⟩ don't
    vb = ab_cpsi_max("Bell_SA+B", J_SB=2.0, gA=0.001)
    vw = ab_cpsi_max("W", J_SB=2.0, gA=0.001)
    v0 = ab_cpsi_max("0++", J_SB=2.0, gA=0.001)
    run(f"Bell={vb:.3f} YES, W={vw:.3f} NO, |0++⟩={v0:.3f} NO",
        vb >= 0.25 and vw < 0.25 and v0 < 0.25)

    # 6. Weak A -> observers see each other
    v = ab_cpsi_max(J_SA=0.3, J_SB=1.0)
    run(f"Weak A (J=0.3,1): AB max={v:.4f} >= 0.25", v >= 0.25)

    # 7. Low gamma lowers threshold
    v = ab_cpsi_max(J_SB=1.2, gS=0.001, gA=0.001, gB=0.001)
    run(f"gamma=0.001, J_SB=1.2: AB max={v:.4f} >= 0.25", v >= 0.25)

    # 8. R_SA + R_SB not conserved.
    #    R = C * Psi^2, NOT (C * Psi)^2. The two agree at t=0 here because the
    #    Bell pair starts at C = 1, which is why the wrong exponent survived:
    #    the printed init was right and only the peak was wrong (0.1955 for
    #    1.8x, against 0.2705 for 2.43x, which is the figure Section 4.3 states).
    H = star_hamiltonian(1.0, 1.0)
    L = dephasing_ops(0.05, 0.05, 0.05)
    rho = make_state("Bell_SA+B")
    dt = 0.001
    max_sum = 0; init_sum = 0
    for step in range(int(5.0/dt) + 1):
        if True:
            sa = ptrace(rho, 2); sb = ptrace(rho, 1)
            r_sa = concurrence(sa) * psi_norm(sa) ** 2
            r_sb = concurrence(sb) * psi_norm(sb) ** 2
            s = r_sa + r_sb
            if step == 0: init_sum = s
            if s > max_sum: max_sum = s
        if step < int(5.0/dt):
            rho = rk4_step(rho, H, L, dt)
    # The fence is at 2.0 rather than 1.5 so that it can see the defect this
    # check was carrying: the wrong exponent (C*Psi)^2 gives 1.76x here and the
    # right one C*Psi^2 gives 2.43x, and a 1.5x fence passes both. It also
    # pins Section 4.3's stated 2.4x rather than merely "not conserved".
    run(f"R not conserved: init={init_sum:.4f}, peak={max_sum:.4f} ({max_sum/init_sum:.1f}x)",
        max_sum > init_sum * 2.0)

    # Summary
    print(f"\n  Result: {ok}/{total} claims verified")
    if ok == total:
        print("  ALL CLAIMS PASS ✓")
    else:
        print(f"  WARNING: {total - ok} claim(s) failed!")
    return 0 if ok == total else 1

if __name__ == "__main__":
    sys.exit(main())
