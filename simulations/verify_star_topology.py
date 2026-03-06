"""
STAR TOPOLOGY — CLAIMS VERIFICATION
====================================
Verifies all key findings documented in STAR_TOPOLOGY_OBSERVERS.md
against star_topology_v2.py simulation.

Run: python verify_star_topology.py
Requires: numpy, star_topology_v2.py in same directory
"""
import sys
import numpy as np

# Import from star_topology_v2
from star_topology_v2 import (
    star_hamiltonian, dephasing_ops, make_state,
    rk4_step, ptrace, concurrence, psi_norm
)

def ab_cpsi_max(state="Bell_SA+B", J_SA=1.0, J_SB=1.0,
                gS=0.05, gA=0.05, gB=0.05, t_max=8.0, dt=0.005):
    """Run simulation, return maximum AB CΨ value."""
    H = star_hamiltonian(J_SA, J_SB)
    L = dephasing_ops(gS, gA, gB)
    rho = make_state(state)
    steps = int(t_max / dt)
    mx = 0.0
    for step in range(steps + 1):
        if step % 10 == 0:
            rp = ptrace(rho, 0)  # trace out S → AB
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
    print("║   STAR TOPOLOGY v2 — CLAIMS VERIFICATION    ║")
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
    run(f"Strong B J=1,2: AB max={v:.4f} ≥ 0.25", v >= 0.25)

    # 3. Threshold at J_SB ≈ 1.466
    v1 = ab_cpsi_max(J_SB=1.465)
    v2 = ab_cpsi_max(J_SB=1.4655)
    run(f"Threshold: J=1.465→{v1:.6f} NO, J=1.4655→{v2:.6f} YES",
        v1 < 0.25 and v2 >= 0.25)

    # 4. Receiver noise more destructive than sender noise
    v_a = ab_cpsi_max(J_SB=2.0, gA=0.25, gB=0.05)
    v_b = ab_cpsi_max(J_SB=2.0, gA=0.05, gB=0.25)
    run(f"γ_A=0.25 kills ({v_a:.4f}), γ_B=0.25 doesn't ({v_b:.4f})",
        v_a < 0.25 and v_b >= 0.25)

    # 5. Only Bell state crosses, W and |0++⟩ don't
    vb = ab_cpsi_max("Bell_SA+B", J_SB=2.0, gA=0.001)
    vw = ab_cpsi_max("W", J_SB=2.0, gA=0.001)
    v0 = ab_cpsi_max("0++", J_SB=2.0, gA=0.001)
    run(f"Bell={vb:.3f} YES, W={vw:.3f} NO, |0++⟩={v0:.3f} NO",
        vb >= 0.25 and vw < 0.25 and v0 < 0.25)

    # 6. Weak A → observers see each other
    v = ab_cpsi_max(J_SA=0.3, J_SB=1.0)
    run(f"Weak A (J=0.3,1): AB max={v:.4f} ≥ 0.25", v >= 0.25)

    # 7. Low γ lowers threshold
    v = ab_cpsi_max(J_SB=1.2, gS=0.001, gA=0.001, gB=0.001)
    run(f"γ=0.001, J_SB=1.2: AB max={v:.4f} ≥ 0.25", v >= 0.25)

    # 8. R_SA + R_SB not conserved
    H = star_hamiltonian(1.0, 1.0)
    L = dephasing_ops(0.05, 0.05, 0.05)
    rho = make_state("Bell_SA+B")
    dt = 0.005
    max_sum = 0; init_sum = 0
    for step in range(int(5.0/dt) + 1):
        if step % 20 == 0:
            sa = ptrace(rho, 2); sb = ptrace(rho, 1)
            r_sa = (concurrence(sa) * psi_norm(sa)) ** 2
            r_sb = (concurrence(sb) * psi_norm(sb)) ** 2
            s = r_sa + r_sb
            if step == 0: init_sum = s
            if s > max_sum: max_sum = s
        if step < int(5.0/dt):
            rho = rk4_step(rho, H, L, dt)
    run(f"R not conserved: init={init_sum:.4f}, peak={max_sum:.4f} ({max_sum/init_sum:.1f}x)",
        max_sum > init_sum * 1.5)

    # Summary
    print(f"\n  Result: {ok}/{total} claims verified")
    if ok == total:
        print("  ALL CLAIMS PASS ✓")
    else:
        print(f"  WARNING: {total - ok} claim(s) failed!")
    return 0 if ok == total else 1

if __name__ == "__main__":
    sys.exit(main())
