using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.PolarityArchitecture;

/// <summary>Schicht-1 registration for
/// <see cref="XGlobalEigenstateMirrorPi2Inheritance"/>: the F99 0°-anchor
/// Pi2 inheritance claim for X⊗N-eigenstates (α = 0 at γ = 1).
///
/// <para>Wired with two parents: the foundational
/// <see cref="HalfAsStructuralFixedPointClaim"/> (polarity 1/2 anchor)
/// and the optional BitA twin
/// <see cref="ZGlobalEigenstateMirrorBitAInheritance"/> via the new ctor
/// parameter introduced in Welle 7. With this registration in place, the
/// X-Mirror BitB ↔ Z-Mirror BitA twin edge is materialized in the registry-
/// built instance (BitATwinStatus = Filled), closing the orphan-partner
/// condition where Z-Mirror BitA was registered without its reciprocating
/// BitB sibling.</para></summary>
public static class XGlobalEigenstateMirrorPi2InheritanceRegistration
{
    public static ClaimRegistryBuilder RegisterXGlobalEigenstateMirrorPi2Inheritance(
        this ClaimRegistryBuilder builder) =>
        builder.Register<XGlobalEigenstateMirrorPi2Inheritance>(b =>
            new XGlobalEigenstateMirrorPi2Inheritance(
                b.Get<HalfAsStructuralFixedPointClaim>(),
                b.Get<ZGlobalEigenstateMirrorBitAInheritance>()));
}
