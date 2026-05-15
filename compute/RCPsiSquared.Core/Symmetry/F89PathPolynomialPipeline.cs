using System.Numerics;
using RCPsiSquared.Core.Numerics;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>F89 path-k closed form (P_k, D_k) native C# Chebyshev pipeline
/// (Tier-1-Derived 2026-05-15). Computes the integer-coefficient path polynomial
/// P_k(y) and the denominator D_k via the F_a eigenvector ansatz + Chebyshev
/// expansion + orbit-polynomial reduction. Native C# replacement for
/// simulations/f89_pathk_symbolic_derivation.py; bit-exact match across k=3..46
/// against the tabulation in <see cref="F89UnifiedFaClosedFormClaim"/> and extends
/// to arbitrary k beyond the int.MaxValue tabulation boundary (k ≥ 47).
///
/// <para>Pipeline:
/// (1) <c>p_n(c) = (2/(m²·k²))·(1−c²)²·A(c)²·B(c)</c> with
///     <c>A(c) = Σ_{j=0..k} U_j(c)·(k−2j)</c>,
///     <c>B(c) = Σ_{j=0..k} U_j(c)²·(k−2j)²</c>
///     using the Chebyshev <c>U_j</c> recursion <c>U_0=1, U_1=2c, U_{j+1}=2c·U_j−U_{j−1}</c>.
/// (2) Substitute <c>c = y/4</c> to get <c>p_n(y)</c> as a degree-(2k+4) rational polynomial in y.
/// (3) Compute the orbit polynomial in y from <c>4^m·(T_m(y/4)−1)</c> via the recurrence
///     <c>S_0=1, S_1=y, S_{j+1}=2y·S_j−16·S_{j−1}</c>; factor out (y−4) once, factor out
///     (y+4) once when m is even, halve once when m is even, then take the integer
///     polynomial square root and divide by 2^FA to get the monic orbit polynomial.
/// (4) Reduce <c>p_n(y)</c> mod the orbit polynomial; the result has degree FA−1 and
///     rational coefficients. D_k is the LCM of the coefficient denominators;
///     P_k(y) = D_k · (reduced polynomial) has integer coefficients.</para>
///
/// <para>All arithmetic is exact (BigInteger / BigRational); no floating-point
/// approximation. Performance: ~milliseconds per path for k ≤ 50, scaling as O(k³).</para></summary>
public static class F89PathPolynomialPipeline
{
    /// <summary>Compute (P_k coefficients low-to-high, D_k) for arbitrary path k ≥ 3.
    /// Bit-exact across k=3..46 against <see cref="F89UnifiedFaClosedFormClaim.PathPolynomial"/>;
    /// extends past the int.MaxValue tabulation boundary at k=47.</summary>
    public static (BigInteger[] CoefficientsLowToHigh, BigInteger Denominator) Compute(int k)
    {
        if (k < 3) throw new ArgumentOutOfRangeException(nameof(k), k, "k must be ≥ 3.");

        var pInY = DerivePnInY(k);
        var orbitPoly = BuildOrbitPolynomialInY(k);
        var reduced = pInY.Mod(orbitPoly);

        if (reduced.IsZero)
            return (new[] { BigInteger.Zero }, BigInteger.One);

        var denomLcm = reduced.DenominatorLcm();
        var coefs = new BigInteger[reduced.Degree + 1];
        for (int i = 0; i <= reduced.Degree; i++)
        {
            var c = reduced[i];
            coefs[i] = c.Numerator * (denomLcm / c.Denominator);
        }
        return (coefs, denomLcm);
    }

    /// <summary>Derive p_n(y) = (2/(m²k²))·(1−c²)²·A(c)²·B(c) substituted c = y/4.
    /// Result is a degree-(2k+4) polynomial in y with rational coefficients.</summary>
    private static RationalPolynomial DerivePnInY(int k)
    {
        int m = k + 2;
        var chebU = ChebyshevUInC(k);
        var A = RationalPolynomial.Zero;
        var B = RationalPolynomial.Zero;
        for (int j = 0; j <= k; j++)
        {
            var weight = new BigRational(k - 2 * j);
            A = A + (RationalPolynomial)weight * chebU[j];
            B = B + (RationalPolynomial)(weight * weight) * (chebU[j] * chebU[j]);
        }

        var oneMinusC2 = new RationalPolynomial(BigRational.One, BigRational.Zero, -BigRational.One);
        var oneMinusC2Squared = oneMinusC2 * oneMinusC2;
        var prefactor = new BigRational(
            new BigInteger(2),
            new BigInteger(m) * new BigInteger(m) * new BigInteger(k) * new BigInteger(k));
        var pInC = prefactor * oneMinusC2Squared * A * A * B;

        return SubstituteCEqualsYOver4(pInC);
    }

    /// <summary>U_j(c) Chebyshev polynomials of the second kind for j = 0..maxDegree.
    /// Recursion: U_0 = 1, U_1 = 2c, U_{j+1} = 2c·U_j − U_{j-1}.</summary>
    private static List<RationalPolynomial> ChebyshevUInC(int maxDegree)
    {
        var result = new List<RationalPolynomial>(maxDegree + 1) { RationalPolynomial.One };
        if (maxDegree == 0) return result;
        var twoC = new RationalPolynomial(BigRational.Zero, new BigRational(2));
        result.Add(twoC);
        for (int j = 2; j <= maxDegree; j++)
            result.Add(twoC * result[j - 1] - result[j - 2]);
        return result;
    }

    /// <summary>Substitute c = y/4 in a polynomial expressed in c: coef of c^j
    /// becomes coef · (1/4)^j as the coef of y^j.</summary>
    private static RationalPolynomial SubstituteCEqualsYOver4(RationalPolynomial pInC)
    {
        if (pInC.IsZero) return RationalPolynomial.Zero;
        var coefs = new BigRational[pInC.Degree + 1];
        var inv4 = new BigRational(BigInteger.One, new BigInteger(4));
        var power = BigRational.One;
        for (int i = 0; i <= pInC.Degree; i++)
        {
            coefs[i] = pInC[i] * power;
            power = power * inv4;
        }
        return new RationalPolynomial(coefs);
    }

    /// <summary>Build the monic orbit polynomial in y annihilating the S_2-anti
    /// orbit {y_n = 4·cos(πn/m) : n = 2, 4, ..., 2·FA}. Uses exact BigInteger
    /// arithmetic on the Chebyshev T_m identity: y_n is a root of
    /// T_m(y_n/4) − 1 = 0 (since m·(πn/m) = πn and cos(πn) = 1 for even n).</summary>
    internal static RationalPolynomial BuildOrbitPolynomialInY(int k)
    {
        int m = k + 2;
        int FA = (k + 1) / 2;

        var Q = BuildScaledChebyshevTMinusOne(m);

        var Q1 = SyntheticDivideByLinearFactor(Q, 4);
        var Q2 = (m % 2 == 0) ? SyntheticDivideByLinearFactor(Q1, -4) : Q1;

        if (m % 2 == 0)
        {
            for (int i = 0; i < Q2.Length; i++)
            {
                if (Q2[i] % 2 != 0)
                    throw new InvalidOperationException(
                        $"Orbit polynomial halving failed: coef at i={i} is {Q2[i]} (odd) for m={m}.");
                Q2[i] /= 2;
            }
        }

        var R = IntegerPolynomialSquareRoot(Q2);

        var twoPowFA = BigInteger.Pow(2, FA);
        var coefs = new BigRational[R.Length];
        for (int i = 0; i < R.Length; i++)
        {
            if (R[i] % twoPowFA != 0)
                throw new InvalidOperationException(
                    $"Orbit square root coef not divisible by 2^FA at i={i}: R[i]={R[i]}, 2^FA={twoPowFA}.");
            coefs[i] = new BigRational(R[i] / twoPowFA);
        }
        return new RationalPolynomial(coefs);
    }

    /// <summary>Compute Q(y) = 4^m · (T_m(y/4) − 1) exactly via the recurrence
    /// S_0 = 1, S_1 = y, S_{j+1} = 2y·S_j − 16·S_{j−1}, then subtract 4^m from the
    /// constant term. Result is a degree-m integer polynomial in y with integer
    /// coefficients and leading coefficient 2^(m−1).</summary>
    private static BigInteger[] BuildScaledChebyshevTMinusOne(int m)
    {
        var s0 = new BigInteger[] { BigInteger.One };
        var s1 = new BigInteger[] { BigInteger.Zero, BigInteger.One };
        if (m == 0) { var copy = new BigInteger[1]; copy[0] = BigInteger.Zero; return copy; }
        if (m == 1)
        {
            var q = new BigInteger[s1.Length];
            Array.Copy(s1, q, s1.Length);
            q[0] -= BigInteger.Pow(4, 1);
            return q;
        }

        BigInteger[] sPrev = s0;
        BigInteger[] sCur = s1;
        for (int j = 1; j < m; j++)
        {
            var sNext = new BigInteger[sCur.Length + 1];
            for (int i = 0; i < sCur.Length; i++)
                sNext[i + 1] = 2 * sCur[i];
            for (int i = 0; i < sPrev.Length; i++)
                sNext[i] -= 16 * sPrev[i];
            sPrev = sCur;
            sCur = sNext;
        }
        var qOut = new BigInteger[sCur.Length];
        Array.Copy(sCur, qOut, sCur.Length);
        qOut[0] -= BigInteger.Pow(4, m);
        return qOut;
    }

    /// <summary>Synthetic-divide a polynomial by (y − root). Throws if division is
    /// not clean.</summary>
    private static BigInteger[] SyntheticDivideByLinearFactor(BigInteger[] coefsLowToHigh, BigInteger root)
    {
        int n = coefsLowToHigh.Length - 1;
        if (n < 1)
            throw new InvalidOperationException("Cannot divide a constant polynomial by a linear factor.");
        var quot = new BigInteger[n];
        var carry = coefsLowToHigh[n];
        quot[n - 1] = carry;
        for (int i = n - 1; i >= 1; i--)
        {
            carry = coefsLowToHigh[i] + root * carry;
            quot[i - 1] = carry;
        }
        var remainder = coefsLowToHigh[0] + root * carry;
        if (!remainder.IsZero)
            throw new InvalidOperationException(
                $"Linear factor (y − {root}) does not divide cleanly; remainder = {remainder}.");
        return quot;
    }

    /// <summary>Compute the integer polynomial R(y) such that R(y)² equals the input,
    /// assuming the input is a perfect square of an integer polynomial. The leading
    /// coefficient of R is the integer square root of the input's leading coefficient.</summary>
    private static BigInteger[] IntegerPolynomialSquareRoot(BigInteger[] coefsLowToHigh)
    {
        int n = coefsLowToHigh.Length - 1;
        if (n % 2 != 0)
            throw new InvalidOperationException(
                $"Cannot take integer polynomial square root of odd-degree polynomial (degree {n}).");
        int d = n / 2;
        var r = new BigInteger[d + 1];

        var leadingSquared = coefsLowToHigh[n];
        if (leadingSquared.Sign <= 0)
            throw new InvalidOperationException($"Leading coef {leadingSquared} must be positive for square root.");
        var leading = IntegerSqrt(leadingSquared);
        if (leading * leading != leadingSquared)
            throw new InvalidOperationException($"Leading coef {leadingSquared} is not a perfect square.");
        r[d] = leading;

        for (int k = n - 1; k >= d; k--)
        {
            int targetIdx = k - d;
            var knownSum = BigInteger.Zero;
            int minI = Math.Max(0, k - d);
            int maxI = Math.Min(d, k);
            for (int i = minI; i <= maxI; i++)
            {
                int j = k - i;
                if (i == targetIdx || j == targetIdx) continue;
                knownSum += r[i] * r[j];
            }
            var numerator = coefsLowToHigh[k] - knownSum;
            var denom = 2 * r[d];
            var (quot, rem) = BigInteger.DivRem(numerator, denom);
            if (!rem.IsZero)
                throw new InvalidOperationException(
                    $"Polynomial not a perfect integer square at coef k={k}: remainder = {rem}.");
            r[targetIdx] = quot;
        }

        for (int k = 0; k < d; k++)
        {
            var sum = BigInteger.Zero;
            int minI = Math.Max(0, k - d);
            int maxI = Math.Min(d, k);
            for (int i = minI; i <= maxI; i++)
                sum += r[i] * r[k - i];
            if (sum != coefsLowToHigh[k])
                throw new InvalidOperationException(
                    $"Polynomial square-root verification failed at coef k={k}: expected {coefsLowToHigh[k]}, got {sum}.");
        }

        return r;
    }

    /// <summary>Newton-iteration integer square root. Used because <c>BigInteger.Sqrt</c>
    /// is documented for .NET 9+ but is not surfaced in the SDK this project resolves
    /// against; fall back to manual Newton until that's resolved (2026-05-15).</summary>
    private static BigInteger IntegerSqrt(BigInteger value)
    {
        if (value.Sign < 0)
            throw new ArgumentOutOfRangeException(nameof(value), "Cannot take square root of negative value.");
        if (value.IsZero) return BigInteger.Zero;
        if (value.IsOne) return BigInteger.One;
        int bitLength = (int)value.GetBitLength();
        var x = BigInteger.One << ((bitLength + 1) / 2);
        while (true)
        {
            var next = (x + value / x) / 2;
            if (next >= x) return x;
            x = next;
        }
    }
}
