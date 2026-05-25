using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.PolarityArchitecture;

/// <summary>Schicht-1 registration for <see cref="F108Part3Pi2YEvenAlwaysPalindromic"/>
/// (F108 Part 3, Y-dephasing sibling of F108 Part 1; same BitB axis).
///
/// <para>Standalone Claim: no ctor parents. Closes the Y-dephasing branch of
/// F108's Π²-even palindrome family and retroactively closes F109's Step 5
/// Y-dephasing branch closed-form (F109 now fully unconditional across all
/// three dephase letters).</para></summary>
public static class F108Part3Pi2YEvenAlwaysPalindromicRegistration
{
    public static ClaimRegistryBuilder RegisterF108Part3Pi2YEvenAlwaysPalindromic(
        this ClaimRegistryBuilder builder) =>
        builder.Register<F108Part3Pi2YEvenAlwaysPalindromic>(_ => new F108Part3Pi2YEvenAlwaysPalindromic());
}
