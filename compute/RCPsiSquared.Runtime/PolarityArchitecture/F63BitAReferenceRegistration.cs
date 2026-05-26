using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.PolarityArchitecture;

/// <summary>Schicht-1 registration for <see cref="F63BitAReference"/>: the
/// lightweight BitA-axis sibling Claim for <see cref="F63LCommutesPi2Pi2Inheritance"/>.
///
/// <para>Standalone Claim (no ctor parents). The F61 connection lives in the
/// docstring + inspector summary only, not as a typed ctor edge, to keep the
/// inheritance graph acyclic (F61 → F63 → F61 would otherwise close). The
/// canonical bit_a Π² conservation derivation lives in
/// <see cref="F61BitAParityPi2Inheritance"/>.</para></summary>
public static class F63BitAReferenceRegistration
{
    public static ClaimRegistryBuilder RegisterF63BitAReference(
        this ClaimRegistryBuilder builder) =>
        builder.Register<F63BitAReference>(_ => new F63BitAReference());
}
