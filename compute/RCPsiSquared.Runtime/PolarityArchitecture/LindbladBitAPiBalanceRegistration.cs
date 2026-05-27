using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.PolarityArchitecture;

/// <summary>Schicht-1 registration for <see cref="LindbladBitAPiBalance"/>
/// (F112-X, Tier1Derived universal N for both Hermitian and non-Hermitian H;
/// Welle 13, 2026-05-27; structural proof in
/// <c>docs/proofs/PROOF_F112_CROSS_DEPHASE_VIA_KLEIN_V4.md</c>).
///
/// <para>Typed ctor parent: <see cref="F108Part2Pi2XEvenAlwaysPalindromic"/>
/// (shared bit_a + X-dephase foundation via F38 Π_X² eigenvalue formula on Pauli
/// strings). F108 Part 2 must be registered BEFORE F112-X so the
/// <c>b.Get&lt;...&gt;()</c> resolves; the registration ordering in
/// <c>KnowledgeRegistryFactory.BuildDefault</c> places this call after F108 Part 2.</para></summary>
public static class LindbladBitAPiBalanceRegistration
{
    public static ClaimRegistryBuilder RegisterLindbladBitAPiBalance(
        this ClaimRegistryBuilder builder) =>
        builder.Register<LindbladBitAPiBalance>(b =>
            new LindbladBitAPiBalance(
                b.Get<F108Part2Pi2XEvenAlwaysPalindromic>()));
}
