namespace RCPsiSquared.Diagnostics.Foundation;

/// <summary>The 2D coordinate map of the birth-canal-surface witness: a symmetric per-site
/// γ-profile [w_edge, w_bulk, …, w_center, …, w_bulk, w_edge] summing to N, parameterized by the
/// two free coordinates (w_edge, w_center). The bulk sites absorb the normalization slack, so two
/// independent coordinates require a bulk reservoir (k_b ≥ 1), i.e. N ≥ 5. N=4 ([e,c,c,e], no bulk)
/// collapses to one free coordinate and is rejected. Pure; no eigendecomposition.</summary>
public static class SymmetricGammaSlice
{
    /// <summary>The number of center sites: 1 for odd N (the single middle site), 2 for even N
    /// (the two middle sites share w_center).</summary>
    public static int CenterCount(int n) => (n % 2 == 0) ? 2 : 1;

    /// <summary>The number of bulk sites: N − 2 edges − k_c centers. ≥ 1 exactly when N ≥ 5.</summary>
    public static int BulkCount(int n) => n - 2 - CenterCount(n);

    /// <summary>The bulk weight that makes the profile sum to N: (N − 2·w_edge − k_c·w_center)/k_b.</summary>
    private static double BulkWeight(int n, double wEdge, double wCenter) =>
        (n - 2.0 * wEdge - CenterCount(n) * wCenter) / BulkCount(n);

    /// <summary>True when the bulk weight is strictly positive (the admissible region). Requires
    /// N ≥ 5; returns false for N &lt; 5 (no bulk reservoir).</summary>
    public static bool IsAdmissible(int n, double wEdge, double wCenter) =>
        n >= 5 && wEdge > 0 && wCenter > 0 && BulkWeight(n, wEdge, wCenter) > 0;

    /// <summary>The length-N symmetric profile for (w_edge, w_center), bulk solved by Σ=N. Throws
    /// for N &lt; 5 (no bulk reservoir for two independent coords) or when the bulk weight is ≤ 0
    /// (outside the admissible region).</summary>
    public static double[] Profile(int n, double wEdge, double wCenter)
    {
        if (n < 5)
            throw new ArgumentOutOfRangeException(nameof(n), n,
                "the symmetric (w_edge, w_center) slice needs a bulk reservoir, i.e. N >= 5 " +
                "(N=4 has no bulk site and collapses to one free coordinate)");
        double wBulk = BulkWeight(n, wEdge, wCenter);
        if (wBulk <= 0)
            throw new ArgumentException(
                $"w_bulk = {wBulk} <= 0 for (w_edge={wEdge}, w_center={wCenter}, N={n}): outside the " +
                "admissible region (2*w_edge + k_c*w_center must be < N).");

        int kc = CenterCount(n);
        int centerLo = (n - kc) / 2;            // first center index
        int centerHi = centerLo + kc - 1;       // last center index
        var p = new double[n];
        for (int i = 0; i < n; i++)
        {
            if (i == 0 || i == n - 1) p[i] = wEdge;
            else if (i >= centerLo && i <= centerHi) p[i] = wCenter;
            else p[i] = wBulk;
        }
        return p;
    }
}
