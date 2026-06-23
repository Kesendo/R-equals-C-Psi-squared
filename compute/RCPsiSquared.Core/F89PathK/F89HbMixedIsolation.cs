using System;
using RCPsiSquared.Core.Numerics;

namespace RCPsiSquared.Core.F89PathK;

/// <summary>Isolates the H_B-mixed factor F_d = charpoly / AT by exact division, gated by the
/// validation TRIPLE that makes the isolation trustworthy (the divide-out route is otherwise
/// vulnerable to shared roots / multiplicity errors): (i) the division is exact (remainder 0),
/// (ii) the quotient has the expected degree d, (iii) gcd(AT, F_d) = 1 (no AT∩F_d crossing at q0).
/// All three together are provably sufficient that the quotient is the true F_d.</summary>
public static class F89HbMixedIsolation
{
    /// <summary>F_d = charpoly ÷ atFactor, validated by the triple. Throws if any leg fails (a
    /// failure means the AT reconstruction was wrong or q0 is degenerate, not a normal branch).</summary>
    public static GaussianInteger[] Isolate(GaussianInteger[] charpoly, GaussianInteger[] atFactor, int expectedDegree)
    {
        var (quotient, remainder) = GaussianPolynomial.DivMod(charpoly, atFactor);
        if (remainder.Length != 0)
            throw new InvalidOperationException(
                "AT factor does not divide the characteristic polynomial exactly (remainder ≠ 0).");
        int degree = GaussianPolynomial.Degree(quotient);
        if (degree != expectedDegree)
            throw new InvalidOperationException(
                $"isolated F_d has degree {degree}, expected {expectedDegree}.");
        if (!GaussianPolynomial.AreCoprime(atFactor, quotient))
            throw new InvalidOperationException(
                "AT factor and F_d share a root (gcd ≠ 1) — q0 lies on an AT∩F_d crossing.");
        return quotient;
    }
}
