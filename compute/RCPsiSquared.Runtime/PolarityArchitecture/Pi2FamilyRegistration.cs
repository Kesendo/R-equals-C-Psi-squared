using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.PolarityArchitecture;

/// <summary>Registers the eight Pi2 polarity-foundation Claims with their internal
/// inheritance edges. Order is topology-driven (Builder resolves):
///
/// <code>
///   PolynomialFoundationClaim (trunk: d²−2d=0 ↔ R=CΨ²)
///     ├── QubitDimensionalAnchorClaim (1/2 = 1/d)
///     │   ├── PolarityLayerOriginClaim (+0/-0)
///     │   │   └── KleinFourCellClaim (F88 4-cell)
///     │   ├── BilinearApexClaim (apex 1/2)
///     │   ├── HalfAsStructuralFixedPointClaim (three faces)
///     │   └── QuarterAsBilinearMaxvalClaim (max p·(1−p) = 1/4)
///     └── NinetyDegreeMirrorMemoryClaim (90° in F80's 2i)
/// </code>
///
/// <para>All eight claims are Tier1Derived with parameterless ctors. The Pi2KnowledgeBase
/// aggregator in Core remains untouched; the Runtime registry holds parallel instances
/// (stateless, so equivalence by content is sufficient).</para>
///
/// <para>QuarterAsBilinearMaxvalClaim joined the eager registry on 2026-05-08 when
/// Pi2DyadicLadderClaim landed: the dyadic halving ladder a_n = 2^(1−n) has 1/4 as the
/// n=3 anchor, so the parent claim must be available in the registry for the ladder's
/// edges to resolve.</para>
///
/// <para>Two Pi2KnowledgeBase entries remain deliberately excluded: <c>Pi2InvolutionClaim</c>
/// (a consequence of F1, not a foundation) and <c>HalfIntegerMirrorClaim</c> (parameterised
/// by N, would need a ChainSystem injection that the polarity architecture does not require).
/// They can be added in a later iteration if the registry needs them.</para></summary>
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
            .Register<QuarterAsBilinearMaxvalClaim>(b =>
            {
                _ = b.Get<QubitDimensionalAnchorClaim>();
                return new QuarterAsBilinearMaxvalClaim();
            })
            .Register<KleinFourCellClaim>(b =>
            {
                _ = b.Get<PolarityLayerOriginClaim>();
                return new KleinFourCellClaim();
            });
}
