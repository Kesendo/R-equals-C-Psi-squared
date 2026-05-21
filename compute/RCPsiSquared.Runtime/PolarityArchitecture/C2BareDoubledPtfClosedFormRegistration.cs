using RCPsiSquared.Core.F86.Item1Derivation;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.PolarityArchitecture;

/// <summary>Schicht-1 wiring of <see cref="C2BareDoubledPtfClosedForm"/>: the F86 c=2
/// bare-doubled-PTF closed form for the K_b observable (post-EP cos, pre-EP cosh, EP limit
/// −5·e⁻²/12), the analytical source of the universal constants x_peak = 2.196910 and
/// HWHM_left/x_peak = 0.671535.
///
/// <para>A standalone Tier1Derived closed-form claim: parameterless, block-independent, no
/// Claim parents (a registry root).</para></summary>
public static class C2BareDoubledPtfClosedFormRegistration
{
    public static ClaimRegistryBuilder RegisterC2BareDoubledPtfClosedForm(
        this ClaimRegistryBuilder builder) =>
        builder.Register<C2BareDoubledPtfClosedForm>(_ => new C2BareDoubledPtfClosedForm());
}
