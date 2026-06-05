namespace RCPsiSquared.Core.Symmetry;

/// <summary>The graded bit_a: a Pauli term's X/Y mask read as a tower of dyadic moments whose depth is
/// the (1 + x)-adic valuation of the mask polynomial. The Klein <c>bit_a</c> (#(X+Y) mod 2) is the
/// bottom rung, <c>bit_a = [Grade == 0]</c>. Inside the bit_a = 0 (even-popcount) cell the grade keeps
/// counting v = 1, 2, …, k−1, partitioning the masks into the F87 / F115 valuation classes
/// c_v = 2^(k-1-v); the F87 hard/soft "Z bit" is the pair-difference of grades (PROOF_F103 §7.7, F115).
///
/// <para>The reading. moment_j of an X/Y mask is the parity of the X/Y count on the sites whose index
/// has j as a submask (j AND i == j). moment_0 is the total X/Y parity, i.e. the Klein bit_a; moment_1
/// is the parity on odd sites; moment_2 the parity on sites with bit 1 set; and so on. The grade is the
/// number of LEADING vanishing moments (the smallest j with moment_j ≠ 0), and it equals the
/// (1 + x)-adic valuation: it measures how deeply the X/Y pattern is dyadically balanced. (1 + x) is an
/// adjacent pair (grade 1), (1 + x)^2 a distance-2 pair (grade 2), (1 + x)^4 a distance-4 pair (grade 4).
/// So the existing bit_a axis is not a bare Z₂: it is the v = 0 face of this graded valuation, and the
/// combinations the grade produces (the §7.8 classes) are the dyadic-balance strata of the X/Y pattern.
/// Verified bit-exact against the GF(2)[x] valuation (<c>WindowedObstructionScan.ValuationAtOnePlusX</c>)
/// and against the class sizes c_v.</para></summary>
public static class BitADyadicGrade
{
    /// <summary>moment_j: parity of the X/Y count on sites i (0 ≤ i &lt; k) whose index has j as a
    /// submask (j AND i == j). moment_0 is the total X/Y parity, i.e. the Klein bit_a.</summary>
    public static int Moment(ulong mask, int j, int k)
    {
        int parity = 0;
        for (int i = 0; i < k; i++)
            if (((mask >> i) & 1UL) != 0 && (i & j) == j)
                parity ^= 1;
        return parity;
    }

    /// <summary>The graded bit_a: the number of leading vanishing dyadic moments, equal to the
    /// (1 + x)-adic valuation of the k-bit X/Y mask. 0 for an odd-popcount mask (Klein bit_a = 1);
    /// for an even-popcount Mixed term it is the dyadic-balance depth v ∈ {1, …, k−1}; k for the
    /// empty mask (no X/Y, valuation conventionally capped at k).</summary>
    public static int Grade(ulong mask, int k)
    {
        for (int j = 0; j < k; j++)
            if (Moment(mask, j, k) != 0) return j;
        return k;
    }

    /// <summary>The Klein bit_a (#(X+Y) mod 2), the bottom rung of the grade: 1 iff the grade is 0.</summary>
    public static int BitA(ulong mask, int k) => Grade(mask, k) == 0 ? 1 : 0;
}
