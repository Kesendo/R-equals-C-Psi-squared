using RCPsiSquared.Core.Numerics;
using GM = RCPsiSquared.Core.Numerics.GaussianRationalMatrix;

namespace RCPsiSquared.Core.Quadrupolar;

/// <summary>The numerical-range bound of PROOF_QUDIT_CLOCK_MIRROR §4 for the quadrupolar dissipator.
/// The Hamiltonian part is skew-adjoint, so every eigenvalue of L has real part in the numerical range
/// of the dissipator D, which here is self-adjoint with eigenvalues in [−r_max, 0]; a palindrome's
/// centre −s is fixed by the trace, which H does not change, and the partner of the eigenvalue 0 sits
/// at −2s. So r_max &lt; 2s, that is D + 2s·1 positive definite, excludes a centre for every Hermitian H
/// (every splitting, offset, shift and lock at once). Decided exactly, block by coherence order, by
/// elimination without pivoting: the gauge is a diagonal similarity, which keeps every leading principal
/// minor of the symmetric physical block, so all pivots are positive exactly when the block is definite
/// (Sylvester). A false is silence, not a centre.</summary>
public static class NumericalRangeBound
{
    public static bool ExcludesEveryCentre(QuadrupolarGenerator g)
    {
        var twoS = 2 * g.MeanRate();
        int d = g.Space.Levels;
        for (int p = -(d - 1); p <= d - 1; p++)
        {
            var b = g.DissipatorBlock(p);
            for (int i = 0; i < b.GetLength(0); i++) b[i, i] += twoS;
            if (!AllPivotsPositive(b)) return false;
        }
        return true;
    }

    /// <summary>Gaussian elimination without pivoting; true when every pivot is strictly positive.</summary>
    public static bool AllPivotsPositive(BigRational[,] a)
    {
        var m = (BigRational[,])a.Clone();
        int n = m.GetLength(0);
        for (int k = 0; k < n; k++)
        {
            if (m[k, k].Sign <= 0) return false;
            for (int r = k + 1; r < n; r++)
            {
                if (m[r, k].IsZero) continue;
                var f = m[r, k] / m[k, k];
                for (int c = k; c < n; c++) m[r, c] -= f * m[k, c];
            }
        }
        return true;
    }

    /// <summary>At extreme narrowing (J₀ = J₁ = J₂ = J) the rate of the rank-k multipoles,
    /// R_k = (J/4)·k(k + 1)·(4I(I + 1) − k(k + 1) − 1), multiplicity 2k + 1 (at 3/2: 6J, 12J, 6J): the
    /// value of <see cref="ExtremeNarrowingFromCasimir"/> on rank k, where the rotation Casimir is k(k + 1);
    /// derived for every I on the sodium page. Abragam and Pound's attenuation of rank-k statistical
    /// tensors (Phys. Rev. 92, 943 (1953), read in a restatement), whose extra 1/(I²(2I − 1)²) is the
    /// normalization of the coupling, absent from these jumps.</summary>
    public static BigRational ExtremeNarrowingRankRate(SpinSpace s, int k, BigRational j)
    {
        BigRational kk = k * (k + 1);
        return j / 4 * kk * (4 * s.SpinSquared - kk - 1);
    }

    /// <summary>The dissipator at extreme narrowing as a polynomial in the rotation Casimir C of
    /// <see cref="SpinSpace.GaugedRotationCasimir"/>: D = −(J/4)·C·((4I(I + 1) − 1)·1 − C), in the gauge.
    /// With all five rates equal and the jumps paired, m with −m, D = −(J/2)·Σ_ab ad(T_ab)²; the trace
    /// part of T_ab drops out of ad; with A_a, B_a left and right multiplication by I_a, X = I(I + 1) =
    /// Σ_a A_a² = Σ_a B_a² and A·B = Σ_a A_a B_a, the contraction is 2X² − X − 2(A·B)² + A·B, while
    /// C = 2X − 2A·B.</summary>
    public static GaussianRational[,] ExtremeNarrowingFromCasimir(SpinSpace s, BigRational j)
    {
        var c = s.GaugedRotationCasimir();
        var shifted = GM.Subtract(GM.Scale(GM.Identity(c.GetLength(0)), 4 * s.SpinSquared - 1), c);
        return GM.Scale(GM.Multiply(c, shifted), -j / 4);
    }

    /// <summary>r_max/2s at extreme narrowing, independent of J: the largest R_k over twice the mean
    /// Σ_k (2k + 1)·R_k / d².</summary>
    public static BigRational ExtremeNarrowingRatio(SpinSpace s)
    {
        var max = BigRational.Zero;
        var sum = BigRational.Zero;
        for (int k = 1; k <= s.TwoI; k++)
        {
            var r = ExtremeNarrowingRankRate(s, k, 1);
            if ((r - max).Sign > 0) max = r;
            sum += (2 * k + 1) * r;
        }
        var mean = sum / (s.Levels * s.Levels);
        return max / (2 * mean);
    }
}
