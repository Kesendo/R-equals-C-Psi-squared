using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.PolarityArchitecture;

/// <summary>Schicht-1 registration for <see cref="LindbladBitBPiBalance"/>
/// (F112, Tier1Derived for Hermitian H; Tier1Candidate non-Hermitian extension
/// documented in inspectables, not typed as a separate Claim).
///
/// <para>Typed ctor parent: <see cref="F108Part1Pi2EvenAlwaysPalindromic"/>
/// (shared bit_b axis foundation via F38 / F63 Π² eigenvalue formula on Pauli
/// strings). F108 Part 1 must be registered BEFORE F112 so the
/// <c>b.Get&lt;...&gt;()</c> resolves; the registration ordering in
/// <c>KnowledgeRegistryFactory.BuildDefault</c> places this call after F108 Part 1
/// (and Part 3) to satisfy that constraint while keeping the F108 family
/// contiguous in the factory.</para></summary>
public static class LindbladBitBPiBalanceRegistration
{
    public static ClaimRegistryBuilder RegisterLindbladBitBPiBalance(
        this ClaimRegistryBuilder builder) =>
        builder.Register<LindbladBitBPiBalance>(b =>
            new LindbladBitBPiBalance(
                b.Get<F108Part1Pi2EvenAlwaysPalindromic>()));
}
