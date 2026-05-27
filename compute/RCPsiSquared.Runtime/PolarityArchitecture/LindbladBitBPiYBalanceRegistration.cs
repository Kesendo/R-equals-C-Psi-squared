using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.PolarityArchitecture;

/// <summary>Schicht-1 registration for <see cref="LindbladBitBPiYBalance"/>
/// (F112-Y, Tier1Derived universal N for both Hermitian and non-Hermitian H;
/// Welle 13, 2026-05-27; structural proof in
/// <c>docs/proofs/PROOF_F112_CROSS_DEPHASE_VIA_KLEIN_V4.md</c>).
///
/// <para>Typed ctor parent: <see cref="F108Part3Pi2YEvenAlwaysPalindromic"/>
/// (shared bit_b + Y-dephase foundation via F38 Π_Y² = (−1)^bit_b eigenvalue formula
/// on Pauli strings). F108 Part 3 must be registered BEFORE F112-Y so the
/// <c>b.Get&lt;...&gt;()</c> resolves; the registration ordering in
/// <c>KnowledgeRegistryFactory.BuildDefault</c> places this call after F108 Part 3.</para></summary>
public static class LindbladBitBPiYBalanceRegistration
{
    public static ClaimRegistryBuilder RegisterLindbladBitBPiYBalance(
        this ClaimRegistryBuilder builder) =>
        builder.Register<LindbladBitBPiYBalance>(b =>
            new LindbladBitBPiYBalance(
                b.Get<F108Part3Pi2YEvenAlwaysPalindromic>()));
}
