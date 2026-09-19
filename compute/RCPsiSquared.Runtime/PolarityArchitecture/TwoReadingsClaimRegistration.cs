using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.PolarityArchitecture;

/// <summary>
/// Registers the stable <see cref="TwoReadingsClaim"/> type as a parentless open
/// synthesis. Individual examples retain their own proofs; operator-space
/// vectorisation does not create a polynomial ancestry edge or a universal
/// exactly-two law.
/// </summary>
public static class TwoReadingsClaimRegistration
{
    public static ClaimRegistryBuilder RegisterTwoReadingsClaim(
        this ClaimRegistryBuilder builder) =>
        builder.Register<TwoReadingsClaim>(_ => new TwoReadingsClaim());
}
