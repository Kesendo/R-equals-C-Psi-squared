using RCPsiSquared.Core.F86;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.PolarityArchitecture;

/// <summary>Schicht-1 wiring for <see cref="LEffMirrorAxisClaim"/>: the F86 2-level
/// reduction L_eff has real part −4γ₀ = −2γ₀·2, the mirror axis of the (−2γ₀, −6γ₀)
/// channel pair, with the Absorption Theorem holding exactly at the ⟨n_XY⟩ = 2 rung.
///
/// <para>Edge declared: <see cref="LEffMirrorAxisClaim"/> ← <see cref="AbsorptionTheoremClaim"/>.
/// The −4γ₀ real part is the Absorption Theorem's −2γ₀·⟨n_XY⟩ rate at ⟨n_XY⟩ = 2 (the
/// 50/50-hybridised pair on the integer rung at the mirror centre), so the mirror-axis
/// claim inherits from the absorption-quantum claim.</para>
///
/// <para>Tier consistency: both are Tier1Derived; the inheritance check (parent at least
/// as strong as child) passes (Tier1Derived ≥ Tier1Derived).</para></summary>
public static class F86LEffMirrorAxisRegistration
{
    public static ClaimRegistryBuilder RegisterF86LEffMirrorAxis(this ClaimRegistryBuilder builder) =>
        builder.Register<LEffMirrorAxisClaim>(b =>
        {
            _ = b.Get<AbsorptionTheoremClaim>();
            return LEffMirrorAxisClaim.Build();
        });
}
