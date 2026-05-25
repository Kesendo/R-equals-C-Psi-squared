using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.PolarityArchitecture;

/// <summary>Schicht-1 registration for <see cref="F108Part2Pi2XEvenAlwaysPalindromic"/>
/// (F108 Part 2, BitA twin of F108 Part 1).
///
/// <para>Standalone Claim: no constructor parents (BitA-axis siblings do not point
/// at BitB-axis Claims; the twin relationship is encoded by F108 Part 1's BitATwin
/// slot pointing at this Claim, not the reverse). Must be registered BEFORE F108
/// Part 1 in <c>KnowledgeRegistryFactory.BuildDefault</c> so that Part 1's ctor
/// dependency on this Claim resolves correctly.</para></summary>
public static class F108Part2Pi2XEvenAlwaysPalindromicRegistration
{
    public static ClaimRegistryBuilder RegisterF108Part2Pi2XEvenAlwaysPalindromic(
        this ClaimRegistryBuilder builder) =>
        builder.Register<F108Part2Pi2XEvenAlwaysPalindromic>(_ => new F108Part2Pi2XEvenAlwaysPalindromic());
}
