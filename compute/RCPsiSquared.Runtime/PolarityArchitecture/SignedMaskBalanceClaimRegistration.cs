using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.PolarityArchitecture;

/// <summary>Wiring of <see cref="SignedMaskBalanceClaim"/> (F173). One parent edge, resolved topologically from
/// <c>BuildDefault</c>: <see cref="PalindromeComplementConnectionClaim"/> (Theorem 2 of the complement connection),
/// which the claim reads where every edge's ratio is a sign.
///
/// <para>Requires: RegisterPalindromeComplementConnectionClaim.</para></summary>
public static class SignedMaskBalanceClaimRegistration
{
    public static ClaimRegistryBuilder RegisterSignedMaskBalanceClaim(this ClaimRegistryBuilder builder) =>
        builder.Register<SignedMaskBalanceClaim>(b => new SignedMaskBalanceClaim(b.Get<PalindromeComplementConnectionClaim>()));
}
