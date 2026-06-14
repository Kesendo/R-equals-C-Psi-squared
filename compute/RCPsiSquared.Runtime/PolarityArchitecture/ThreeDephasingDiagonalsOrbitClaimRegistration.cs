using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.PolarityArchitecture;

/// <summary>Schicht-1 wiring of <see cref="ThreeDephasingDiagonalsOrbitClaim"/> (2026-06-14): the
/// three dephasing diagonals {Q_X, Q_Y, Q_Z} as one orbit of the single-qubit Clifford basis-change
/// S₃ ⟨h_zx, h_yz⟩ (same spectrum), and the one-diagonal's three readings (rate = D-fix, mirror =
/// R·Q·R=−Q, judge = {D, 𝓕D} cell) as the mirror group D₄ acting within a diagonal; the structure is
/// S₃ ⋉ D₄.
///
/// <para>The two typed parents are the physics edge that welds the mirror-group and absorption
/// clusters (previously joined only at the d²−2d=0 foundation): <see cref="MirrorGroupD4Claim"/>
/// (the readings + the D₄ factor) and <see cref="AbsorptionTheoremClaim"/> (the dephasing diagonal
/// L_D = γ·(Q − N·I)). Anchor: <c>simulations/one_diagonal_mirror_group.py</c> (self-validating,
/// the physics-first gate that corrected the first hypothesis).</para></summary>
public static class ThreeDephasingDiagonalsOrbitClaimRegistration
{
    public static ClaimRegistryBuilder RegisterThreeDephasingDiagonalsOrbitClaim(
        this ClaimRegistryBuilder builder) =>
        builder.Register<ThreeDephasingDiagonalsOrbitClaim>(b =>
            new ThreeDephasingDiagonalsOrbitClaim(
                b.Get<MirrorGroupD4Claim>(),
                b.Get<AbsorptionTheoremClaim>()));
}
