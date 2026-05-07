using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.PolarityArchitecture;

/// <summary>Registers the seven Pi2 polarity-foundation Claims with their internal
/// inheritance edges. Order is topology-driven (Builder resolves):
///
/// <code>
///   PolynomialFoundationClaim (trunk: d²−2d=0 ↔ R=CΨ²)
///     ├── QubitDimensionalAnchorClaim (1/2 = 1/d)
///     │   ├── PolarityLayerOriginClaim (+0/-0)
///     │   │   └── KleinFourCellClaim (F88 4-cell)
///     │   ├── BilinearApexClaim (apex 1/2)
///     │   └── HalfAsStructuralFixedPointClaim (three faces)
///     └── NinetyDegreeMirrorMemoryClaim (90° in F80's 2i)
/// </code>
///
/// <para>All seven claims are Tier1Derived with parameterless ctors. The Pi2KnowledgeBase
/// aggregator in Core remains untouched; the Runtime registry holds parallel instances
/// (stateless, so equivalence by content is sufficient).</para></summary>
public static class Pi2FamilyRegistration
{
    public static ClaimRegistryBuilder RegisterPi2Family(this ClaimRegistryBuilder builder) =>
        builder
            .Register<PolynomialFoundationClaim>(_ =>
                new PolynomialFoundationClaim())
            .Register<QubitDimensionalAnchorClaim>(b =>
            {
                _ = b.Get<PolynomialFoundationClaim>();
                return new QubitDimensionalAnchorClaim();
            })
            .Register<NinetyDegreeMirrorMemoryClaim>(b =>
            {
                _ = b.Get<PolynomialFoundationClaim>();
                return new NinetyDegreeMirrorMemoryClaim();
            })
            .Register<PolarityLayerOriginClaim>(b =>
            {
                _ = b.Get<QubitDimensionalAnchorClaim>();
                return new PolarityLayerOriginClaim();
            })
            .Register<BilinearApexClaim>(b =>
            {
                _ = b.Get<QubitDimensionalAnchorClaim>();
                return new BilinearApexClaim();
            })
            .Register<HalfAsStructuralFixedPointClaim>(b =>
            {
                _ = b.Get<QubitDimensionalAnchorClaim>();
                return new HalfAsStructuralFixedPointClaim();
            })
            .Register<KleinFourCellClaim>(b =>
            {
                _ = b.Get<PolarityLayerOriginClaim>();
                return new KleinFourCellClaim();
            });
}
