using RCPsiSquared.Core.F1;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.PolarityArchitecture;

/// <summary>Wiring of <see cref="DirectSumDecompositionClaim"/>. Two parent edges,
/// both resolved earlier in <c>BuildDefault</c>:
///
/// <list type="bullet">
///   <item><see cref="F1PalindromeIdentity"/>: the global palindrome
///         Π·L·Π⁻¹ = −L − 2Σγ·I that the claim restricts to the parity sectors
///         (odd N: mirror-image sector exchange; even N: per-sector self-palindromy).</item>
///   <item><see cref="F61BitAParityPi2Inheritance"/>: the bit_a (n_XY) Z₂ grading
///         that defines V_even/V_odd and the superselection charge P_XY = Π²_X.</item>
/// </list>
///
/// <para>Requires: RegisterF1Family, RegisterF61BitAParityPi2Inheritance.</para></summary>
public static class DirectSumDecompositionClaimRegistration
{
    public static ClaimRegistryBuilder RegisterDirectSumDecompositionClaim(
        this ClaimRegistryBuilder builder) =>
        builder.Register<DirectSumDecompositionClaim>(b =>
        {
            var f1 = b.Get<F1PalindromeIdentity>();
            var f61 = b.Get<F61BitAParityPi2Inheritance>();
            return new DirectSumDecompositionClaim(f1, f61);
        });
}
