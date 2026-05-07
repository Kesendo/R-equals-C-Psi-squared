using RCPsiSquared.Core.F86;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.PolarityArchitecture;

/// <summary>Registers <see cref="PolarityInheritanceLink"/> as a cross-KB consumer
/// of the Pi2 polarity-foundation family. Three edges declared:
/// PolynomialFoundationClaim (trunk), PolarityLayerOriginClaim (the +0/-0 source),
/// HalfAsStructuralFixedPointClaim (the three-faces synthesis).
///
/// <para>Tier consistency: PolarityInheritanceLink is Tier2Verified; all three
/// parents are Tier1Derived. The TierStrength inheritance check (parent at least as
/// strong as child) passes.</para>
///
/// <para>The F86KnowledgeBase aggregator's own PolarityInheritanceLink instance is
/// unrelated to the Runtime registration; both call PolarityInheritanceLink.Build()
/// and produce structurally-equivalent stateless claims.</para></summary>
public static class F86PolarityLinkRegistration
{
    public static ClaimRegistryBuilder RegisterF86PolarityLink(this ClaimRegistryBuilder builder) =>
        builder.Register<PolarityInheritanceLink>(b =>
        {
            _ = b.Get<PolynomialFoundationClaim>();
            _ = b.Get<PolarityLayerOriginClaim>();
            _ = b.Get<HalfAsStructuralFixedPointClaim>();
            return PolarityInheritanceLink.Build();
        });
}
