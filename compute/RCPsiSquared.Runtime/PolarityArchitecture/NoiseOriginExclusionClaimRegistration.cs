using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.PolarityArchitecture;

/// <summary>Wiring of <see cref="NoiseOriginExclusionClaim"/>. Two parent edges:
/// <see cref="PolynomialFoundationClaim"/> (d²−2d=0, Candidate 5's algebra) and
/// <see cref="QubitDimensionalAnchorClaim"/> (the d=2 anchor), both registered earlier in
/// <c>BuildDefault</c> via the Pi2 family.
///
/// <para>Requires: RegisterPolynomialFoundationClaim, RegisterQubitDimensionalAnchorClaim (Pi2 family).</para></summary>
public static class NoiseOriginExclusionClaimRegistration
{
    public static ClaimRegistryBuilder RegisterNoiseOriginExclusionClaim(
        this ClaimRegistryBuilder builder) =>
        builder.Register<NoiseOriginExclusionClaim>(b =>
            new NoiseOriginExclusionClaim(b.Get<PolynomialFoundationClaim>(), b.Get<QubitDimensionalAnchorClaim>()));
}
