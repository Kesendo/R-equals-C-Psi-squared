using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.PolarityArchitecture;

/// <summary>Schicht-1 registration for
/// <see cref="ZGlobalEigenstateMirrorBitAInheritance"/>: the BitA twin of
/// <see cref="XGlobalEigenstateMirrorPi2Inheritance"/> at the Z⊗N-eigenstate
/// Mirror anchor.
///
/// <para>One typed parent: <see cref="HalfAsStructuralFixedPointClaim"/>. Wires
/// the polarity 1/2 anchor that grounds the γ_X = ±1 endpoint (same parent
/// pattern as the X version). The BitA twin edge from the BitB
/// X-version Claim is wired via that Claim's optional ctor parameter.</para></summary>
public static class ZGlobalEigenstateMirrorBitAInheritanceRegistration
{
    public static ClaimRegistryBuilder RegisterZGlobalEigenstateMirrorBitAInheritance(
        this ClaimRegistryBuilder builder) =>
        builder.Register<ZGlobalEigenstateMirrorBitAInheritance>(b =>
            new ZGlobalEigenstateMirrorBitAInheritance(b.Get<HalfAsStructuralFixedPointClaim>()));
}
