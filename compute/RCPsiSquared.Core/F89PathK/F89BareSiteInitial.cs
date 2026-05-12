namespace RCPsiSquared.Core.F89PathK;

/// <summary>Closed-form (ρ_l)_{0,1}(0) for any single bare (out-of-block) site l in the F89
/// path-k initial condition: <c>(N - 1) / (2 · √(N²(N-1)/2))</c>.
///
/// <para>Same scalar for every bare site by site-uniformity; used as the per-bare-site amplitude
/// anchor in F89 multi-exponential and path-k vacuum/SE Parseval derivations.</para>
/// </summary>
public static class F89BareSiteInitial
{
    public static double BareSiteInitial01(int N)
    {
        if (N < 2)
            throw new ArgumentOutOfRangeException(nameof(N), N, "N must be ≥ 2.");
        return (N - 1) / (2.0 * Math.Sqrt((double)N * N * (N - 1) / 2.0));
    }
}
