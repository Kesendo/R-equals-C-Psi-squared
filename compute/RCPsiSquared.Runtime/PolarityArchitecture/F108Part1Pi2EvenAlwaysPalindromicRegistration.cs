using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.PolarityArchitecture;

/// <summary>Schicht-1 registration for <see cref="F108Part1Pi2EvenAlwaysPalindromic"/>.
/// Typed ctor parent: <see cref="F108Part2Pi2XEvenAlwaysPalindromic"/> (BitA twin);
/// Part 2 must be registered BEFORE Part 1 so the <c>b.Get&lt;...&gt;()</c> resolves.</summary>
public static class F108Part1Pi2EvenAlwaysPalindromicRegistration
{
    public static ClaimRegistryBuilder RegisterF108Part1Pi2EvenAlwaysPalindromic(
        this ClaimRegistryBuilder builder) =>
        builder.Register<F108Part1Pi2EvenAlwaysPalindromic>(b =>
            new F108Part1Pi2EvenAlwaysPalindromic(
                b.Get<F108Part2Pi2XEvenAlwaysPalindromic>()));
}
