using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.PolarityArchitecture;

/// <summary>Schicht-1 wiring of <see cref="MomentTowerPumpChannelClaim"/> (2026-06-11): the
/// F120 moment-tower pump channel. Amplitude damping is the unique non-unital piece of the
/// standard Lindbladian and pumps along pure local Z (D[σ⁻_l](I) = +Z_l,
/// D[σ⁺_l](I) = −Z_l), so d/dt ⟨A⟩|_{I/d} = (1/d)·Σ_l Δγ_l·Tr(A·Z_l) with
/// Δγ_l = γ↓_l − γ↑_l; with A = H^j the slope reads the girth-ladder tower
/// t_j(l) = Tr(Z_l·H^j) linearly. Dephasing-blind, evolution-blind, closed at detailed
/// balance; rung 1 is F113 (asymmetry = −4^N·slope⟨H⟩ exactly); the curvature is exactly
/// affine in the generator and fingerprints X/Y-flavored parasites against the commutator
/// probes [Z_l, H_p^j] while Z-flavored parasites stay exactly invisible; the girth
/// certificate is one-sided (a firing rung proves m* = 2ℓ+1; silence is not softness).
///
/// <para>Tier1Derived (one-line identities, exact; self-check battery at N = 2 and N = 3 in
/// the ctor). Typed parents: <see cref="LindbladBitBPiBreakMagnitude"/> (F113, whose closed
/// form is the channel's first rung) and
/// <see cref="F84ThermalAmplitudeDampingPi2Inheritance"/> (F84, whose vacuum rate Δγ is the
/// pump weight and whose detailed-balance regime is blindness #3). The girth-ladder
/// primitive (RCPsiSquared.Diagnostics/F87/GirthLadder.cs) is a compute primitive, not a
/// Claim, and is carried in prose. Anchor:
/// <c>docs/proofs/PROOF_MOMENT_TOWER_PUMP_CHANNEL.md</c> +
/// <c>simulations/moment_tower_pump_channel.py</c>.</para></summary>
public static class MomentTowerPumpChannelClaimRegistration
{
    public static ClaimRegistryBuilder RegisterMomentTowerPumpChannelClaim(
        this ClaimRegistryBuilder builder) =>
        builder.Register<MomentTowerPumpChannelClaim>(b =>
            new MomentTowerPumpChannelClaim(
                b.Get<LindbladBitBPiBreakMagnitude>(),
                b.Get<F84ThermalAmplitudeDampingPi2Inheritance>()));
}
