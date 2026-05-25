using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.PolarityArchitecture;

/// <summary>Schicht-1 registration for <see cref="F108Part2Pi2XEvenAlwaysPalindromic"/>.
/// Standalone Claim (BitA-axis siblings do not point at BitB Claims; the twin edge
/// lives on F108 Part 1). Must be registered BEFORE F108 Part 1 so Part 1's ctor
/// dependency resolves.</summary>
public static class F108Part2Pi2XEvenAlwaysPalindromicRegistration
{
    public static ClaimRegistryBuilder RegisterF108Part2Pi2XEvenAlwaysPalindromic(
        this ClaimRegistryBuilder builder) =>
        builder.Register<F108Part2Pi2XEvenAlwaysPalindromic>(_ => new F108Part2Pi2XEvenAlwaysPalindromic());
}
