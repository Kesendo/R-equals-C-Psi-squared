using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.PolarityArchitecture;

/// <summary>Schicht-1 wiring of <see cref="HeldLetterRoutingClaim"/> (the label layer, typed:
/// the dephasing routes by the letter it holds; L_P(S) = −2γ·n_anti(S, P)·S, the price
/// list letter-routed,
/// <c>docs/quantum/DEPHASING_TRANSLATED.md</c> §4). Two typed parent edges:
///
/// <list type="bullet">
///   <item><see cref="AbsorptionTheoremClaim"/>: the rate law Re λ = −2γ·⟨n_XY⟩ — the price
///         list itself, light as the letters anticommuting with the held letter.</item>
///   <item><see cref="Pi2KleinV4DephaseSwapGroup"/>: the {D, Q_zx, Q_yx} operator-space
///         involutions — the letter swap that relocates which cells pay.</item>
/// </list>
///
/// <para>Tier consistency: Tier 1 derived (both parents Tier 1 derived; the exact core is
/// machine-verified live at <c>inspect --root label</c>, <c>HeldLetterRoutingWitness</c>;
/// the Tier-4 canvas reading rides as prose children, never promoted).</para>
///
/// <para>Requires <see cref="AbsorptionTheoremClaimRegistration.RegisterAbsorptionTheoremClaim"/>
/// and the Pi2 Klein V₄ registration (the builder topo-resolves, so registration order is
/// free).</para></summary>
public static class HeldLetterRoutingClaimRegistration
{
    public static ClaimRegistryBuilder RegisterHeldLetterRoutingClaim(
        this ClaimRegistryBuilder builder) =>
        builder.Register<HeldLetterRoutingClaim>(b =>
            new HeldLetterRoutingClaim(
                b.Get<AbsorptionTheoremClaim>(),
                b.Get<Pi2KleinV4DephaseSwapGroup>()));
}
