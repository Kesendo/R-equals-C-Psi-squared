using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.PolarityArchitecture;

/// <summary>Schicht-1 registration for <see cref="F108Part1Pi2EvenAlwaysPalindromic"/>
/// (F108 Part 1).
///
/// <para>Typed ctor parent: <see cref="F108Part2Pi2XEvenAlwaysPalindromic"/> (BitA
/// twin, X-dephasing branch). Part 2 must be registered BEFORE Part 1 in the
/// factory chain so the <c>b.Get&lt;...&gt;()</c> dependency resolves.</para>
///
/// <para>First BitB-axis Claim added in the post-F107 wave of the cubic-unpacking
/// arc. Closes "Π²-even pairs are never F87-hard under Z-dephasing" via the
/// Π_5bilinear phase-variant Π operator; together with F108 Part 2 (which closes
/// the X-dephasing branch) retroactively closes F109's Step 5 for Z- and
/// X-dephasing. The Y-dephasing branch (F108 Part 3) remains open.</para></summary>
public static class F108Part1Pi2EvenAlwaysPalindromicRegistration
{
    public static ClaimRegistryBuilder RegisterF108Part1Pi2EvenAlwaysPalindromic(
        this ClaimRegistryBuilder builder) =>
        builder.Register<F108Part1Pi2EvenAlwaysPalindromic>(b =>
            new F108Part1Pi2EvenAlwaysPalindromic(
                b.Get<F108Part2Pi2XEvenAlwaysPalindromic>()));
}
