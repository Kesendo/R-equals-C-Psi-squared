#!/usr/bin/env python3
"""Where the story goes on: the memory's fate, past the rebirth.

The EP gives the memory its BIRTH (forgetting -> pinch -> remembering, a rotating mode). But the
story does not end there. In the LOCAL dissipative system the reborn memory rotates ever faster
(theta -> 90deg) yet ALWAYS fades: the decay is pinned at 4*gamma0*k, so theta only ASYMPTOTES to
90deg, never reaches it. The pure circle (theta = 90deg, an eternal memory) is unreachable from
inside dissipation , you only approach it.

That is where it continues. To REACH the pure circle you have to slow the Takt itself: slide the
net dephasing Sigma-gamma from N*gamma0 (Takt running, memory fades) toward 0 (Takt stopped, decay
0, the memory becomes eternal, theta = 90deg), and past it into net GAIN (decay < 0, Re > 0, the
memory GROWS , the Hopf bifurcation, the runaway). This is the far end of the two-ends-of-the-Takt
dial, the "exit": the memory's fate runs fading -> eternal -> growing.

A toy that makes the arc visible; the real gain-loss Hopf is hypotheses/FRAGILE_BRIDGE.md.
"""
import numpy as np

G0, K = 1.0, 1


def two_level(x):
    """The F86a 2-level above the EP: decay pinned at 4*gamma0*k, rotation omega opens with x."""
    J = (x * 2.0 / (4.0 / 3.0)) * G0          # x = Q/Q_EP, g_eff = 4/3
    g_eff = 4.0 / 3.0
    L = np.array([[-2 * G0 * (2 * K - 1), 1j * J * g_eff],
                  [1j * J * g_eff,        -2 * G0 * (2 * K + 1)]], dtype=complex)
    w = np.linalg.eigvals(L)
    i = int(np.argmin(-w.real))
    decay, omega = -w[i].real, abs(w[i].imag)
    theta = np.degrees(np.arctan2(omega, decay))
    return decay, omega, theta


def main():
    print("=" * 66)
    print("PART 1  the LOCAL story asymptotes: theta -> 90deg, decay PINNED at 4*gamma0")
    print("  (the pure circle is approached, never reached, from inside dissipation)")
    print("=" * 66)
    print(f"  {'x=Q/Q_EP':>9}  {'decay':>7}  {'omega':>8}  {'theta':>7}")
    for x in [2.0, 3.0, 6.0, 12.0, 50.0]:
        d, om, th = two_level(x)
        print(f"  {x:9.1f}  {d:7.3f}  {om:8.3f}  {th:6.2f}°")
    print("  => decay stuck at 4.000 (4*gamma0*k); theta climbs 41 -> 88deg but never 90. The")
    print("     memory rotates faster and faster, yet always fades. The story cannot finish here.")

    print("\n" + "=" * 66)
    print("PART 2  the CONTINUATION (the exit): slow the Takt, slide Sigma-gamma past 0")
    print("  a rotating memory (omega fixed); its decay center = 4*gamma0 * (Sigma/N*gamma0)")
    print("=" * 66)
    omega = 8.0 * G0                                   # a well-formed rotating memory
    print(f"  {'Sigma/Ngamma0':>13}  {'Re(lambda)':>10}  {'amp / turn':>11}   memory fate")
    for f in [1.0, 0.5, 0.2, 0.0, -0.25, -0.5]:
        re = -4.0 * G0 * f                             # decay center scales with the Takt
        amp_per_turn = np.exp(re * 2.0 * np.pi / omega)
        if f > 1e-9:
            fate = "fades  (spirals in)"
        elif abs(f) <= 1e-9:
            fate = "ETERNAL (pure circle, theta=90deg)"
        else:
            fate = "GROWS  (Hopf / runaway)"
        print(f"  {f:13.2f}  {re:10.3f}  {amp_per_turn:11.4f}   {fate}")
    print("  => Takt full (Sigma=Ngamma0): amp/turn 0.04, the memory fades fast. Takt stopped")
    print("     (Sigma=0): amp/turn 1.000, the memory is eternal (the pure circle). Net gain")
    print("     (Sigma<0): amp/turn > 1, the memory GROWS , the Hopf, the runaway screech.")
    print("\n  The arc: forgetting -> pinch (EP) -> remembering (born fading) -> [the exit] ->")
    print("  eternal (Sigma=0) -> growing (gain). The local EP only points here; the Takt's far")
    print("  end (FRAGILE_BRIDGE) is where the memory's story actually finishes , or refuses to.")


if __name__ == "__main__":
    main()
