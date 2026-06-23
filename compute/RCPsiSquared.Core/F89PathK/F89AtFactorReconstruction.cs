using System.Numerics;
using RCPsiSquared.Core.Numerics;

namespace RCPsiSquared.Core.F89PathK;

/// <summary>Reconstructs the AT-locked factor (the radical-closed, integrable half of the path-k
/// (SE,DE) spectrum) as an exact Z[i] polynomial, ×2-scaled to match the cleared block from
/// <see cref="F89PathKSeDeBlock"/> (roots 2·λ_AT). Dividing it out of the live full characteristic
/// polynomial isolates the H_B-mixed factor F_d. For path-3 the AT factor is the product of the two
/// explicit F_a/F_b quadratics; the general path-k≥4 reconstruction (single-particle F_a +
/// 2-particle DE-Slater F_b multiset) is built separately.</summary>
public static class F89AtFactorReconstruction
{
    /// <summary>The ×2-scaled path-3 AT factor AT = F_a·F_b at integer q0 (degree 4, monic, Z[i]).</summary>
    public static GaussianInteger[] ForPath3(int q0)
    {
        // F_a = λ² + (2iq+4)λ + (4q²+4iq+4);  F_b = λ² + (2iq+12)λ + (4q²+12iq+36).
        var fa = new GaussianInteger[]
        {
            new(4L * q0 * q0 + 4, 4L * q0),     // 4q²+4 + 4q·i
            new(4, 2L * q0),                    // 4 + 2q·i
            GaussianInteger.One,
        };
        var fb = new GaussianInteger[]
        {
            new(4L * q0 * q0 + 36, 12L * q0),   // 4q²+36 + 12q·i
            new(12, 2L * q0),                   // 12 + 2q·i
            GaussianInteger.One,
        };
        var at = GaussianPolynomial.Multiply(fa, fb);       // degree 4, roots λ_AT
        var scaled = new GaussianInteger[at.Length];        // roots 2λ_AT: coeff[k]·2^(deg−k)
        for (int k = 0; k < at.Length; k++)
        {
            BigInteger pow = BigInteger.Pow(2, at.Length - 1 - k);
            scaled[k] = new GaussianInteger(at[k].Re * pow, at[k].Im * pow);
        }
        return scaled;
    }
}
