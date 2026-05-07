using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Diagnostics.F87;

/// <summary>Registers the parameterless F87 trichotomy claims into the typed-knowledge
/// runtime. Three claims:
///
/// <list type="bullet">
///   <item><see cref="F87TrichotomyClassification"/>: the F87 truly/soft/hard split via
///         F1 residual (Tier1Derived).</item>
///   <item><see cref="DissipatorResonanceLaw"/>: SU(2)-symmetric Klein-cell alignment law
///         (Tier1Derived).</item>
///   <item><see cref="DissipatorAxisSelectsPolarityClaim"/>: the typed bridge claim that
///         declares F87's dissipator-axis choice corresponds to
///         <see cref="PolarityLayerOriginClaim"/>'s polarity axis (Tier1Derived).</item>
/// </list>
///
/// <para>The DissipatorAxisSelectsPolarityClaim claim declares an explicit cross-KB edge to
/// PolarityLayerOriginClaim, so the polarity family must be registered first.
/// F87CanonicalWitness is parameterised by chain + terms + expected class and is left
/// for a per-witness audit iteration.</para>
///
/// <para>Layer note: this extension lives in <c>RCPsiSquared.Diagnostics</c> because the
/// F87 typed claims live there. Diagnostics gains a project reference to Runtime to
/// reach <see cref="ClaimRegistryBuilder"/>; the production-side asymmetry remains
/// Core ← Runtime + Core ← Diagnostics → Runtime, no circular reference.</para></summary>
public static class F87FamilyRegistration
{
    public static ClaimRegistryBuilder RegisterF87Family(this ClaimRegistryBuilder builder) =>
        builder
            .Register<F87TrichotomyClassification>(_ => new F87TrichotomyClassification())
            .Register<DissipatorResonanceLaw>(_ => new DissipatorResonanceLaw())
            .Register<DissipatorAxisSelectsPolarityClaim>(b =>
            {
                _ = b.Get<PolarityLayerOriginClaim>();
                return new DissipatorAxisSelectsPolarityClaim();
            });
}
