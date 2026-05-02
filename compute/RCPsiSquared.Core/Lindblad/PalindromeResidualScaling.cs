namespace RCPsiSquared.Core.Lindblad;

public enum HamiltonianClass
{
    /// <summary>Bilinear nearest-neighbour bond Hamiltonians (XY, Heisenberg, etc.).</summary>
    Main,
    /// <summary>Single-body Hamiltonians (e.g. uniform field).</summary>
    SingleBody,
}

/// <summary>Closed-form scaling of the F1-palindrome residual ‖M‖²_F = ‖Π·L·Π⁻¹ + L + 2σ·I‖²_F.
///
/// Verified on chain, ring, star, K_N at N = 4, 5 to machine precision (see
/// experiments/OPERATOR_RIGIDITY_ACROSS_CUSP.md):
///
///   main class          ‖M(N, G)‖² = c_H · B(G) · 4^(N − 2)
///   single-body class   ‖M(N, G)‖² = c_H · (D2(G) / 2) · 4^(N − 2)
///
/// where B(G) = bond count, D2(G) = Σ_i deg_G(i)². For chains: B = N − 1 and D2 = 4N − 6.
/// </summary>
public static class PalindromeResidualScaling
{
    /// <summary>Graph-aware factor: F(N, G) such that ‖M(N, G)‖²_F = c_H · F(N, G).</summary>
    public static double FactorFromGraph(int N, int bondCount, int degreeSquaredSum, HamiltonianClass cls)
    {
        if (N < 2) throw new ArgumentOutOfRangeException(nameof(N), $"N must be ≥ 2; got {N}");
        if (bondCount < 1) throw new ArgumentOutOfRangeException(nameof(bondCount), $"B must be ≥ 1; got {bondCount}");
        if (degreeSquaredSum < 2 * bondCount)
            throw new ArgumentOutOfRangeException(nameof(degreeSquaredSum),
                $"D2 ({degreeSquaredSum}) inconsistent with B ({bondCount}); D2 ≥ 2B = {2 * bondCount}");
        double pow = Math.Pow(4, N - 2);
        return cls switch
        {
            HamiltonianClass.Main => bondCount * pow,
            HamiltonianClass.SingleBody => (degreeSquaredSum / 2.0) * pow,
            _ => throw new ArgumentOutOfRangeException(nameof(cls)),
        };
    }

    /// <summary>Chain-specific F(N): main class (N − 1)·4^(N−2); single-body (2N − 3)·4^(N−2).</summary>
    public static double FactorChain(int N, HamiltonianClass cls)
    {
        if (N < 2) throw new ArgumentOutOfRangeException(nameof(N), $"N must be ≥ 2; got {N}");
        double pow = Math.Pow(4, N - 2);
        return cls switch
        {
            HamiltonianClass.Main => (N - 1) * pow,
            HamiltonianClass.SingleBody => (2 * N - 3) * pow,
            _ => throw new ArgumentOutOfRangeException(nameof(cls)),
        };
    }

    /// <summary>Adjacent-N ratio ‖M(N+1)‖²/‖M(N)‖²: main 4·N/(N−1); single-body 4·(2N−1)/(2N−3).</summary>
    public static double AdjacentRatio(int N, HamiltonianClass cls)
    {
        if (N < 2) throw new ArgumentOutOfRangeException(nameof(N), $"N must be ≥ 2; got {N}");
        return cls switch
        {
            HamiltonianClass.Main => 4.0 * N / (N - 1),
            HamiltonianClass.SingleBody when 2 * N - 3 > 0 => 4.0 * (2 * N - 1) / (2 * N - 3),
            HamiltonianClass.SingleBody => throw new ArgumentOutOfRangeException(nameof(N), $"single-body formula requires N ≥ 2; got {N}"),
            _ => throw new ArgumentOutOfRangeException(nameof(cls)),
        };
    }
}
