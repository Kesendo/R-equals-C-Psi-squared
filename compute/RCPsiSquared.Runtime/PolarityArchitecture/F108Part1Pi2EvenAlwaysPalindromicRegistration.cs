using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.PolarityArchitecture;

/// <summary>Schicht-1 registration for <see cref="F108Part1Pi2EvenAlwaysPalindromic"/>
/// (F108 Part 1).
///
/// <para>Standalone Claim: no constructor parents. First BitB-axis Claim added in the
/// post-F107 wave of the cubic-unpacking arc. Closes the previously open F108 Part 1
/// (Π²-even pairs are never F87-hard) via the Π_5bilinear phase-variant Π operator;
/// also retroactively promotes F109 (MotherSoftYParityOnePurity) from "Tier1Derived
/// modulo F108 Part 1" to fully unconditional Tier1Derived.</para></summary>
public static class F108Part1Pi2EvenAlwaysPalindromicRegistration
{
    public static ClaimRegistryBuilder RegisterF108Part1Pi2EvenAlwaysPalindromic(
        this ClaimRegistryBuilder builder) =>
        builder.Register<F108Part1Pi2EvenAlwaysPalindromic>(_ => new F108Part1Pi2EvenAlwaysPalindromic());
}
