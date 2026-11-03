namespace RCPsiSquared.Core.Resonance;

/// <summary>The F86 envelope max_t |K_b(Q, t)| read off a seed t-grid and refined in continuous t.
///
/// <para>F86's Q_peak is the peak of the ENVELOPE, argmax_Q max_t |K_b(Q, t)| (the origin,
/// <c>experiments/Q_SCALE_THREE_BANDS.md</c> Result 4, takes the maximum over t; a fixed-t slice
/// is a different observable). A grid maximum alone is not that envelope: it reads the fixed-t
/// slice K_b(Q, t_j) at the winning grid point t_j, and its Q_peak sits off the envelope's by
/// s·(t_j − t*), linear in the distance of the true maximiser t* from the grid point, where
/// s = dQ_slice/dt (about −0.49, −0.27, −0.23 on the three c = 2, N = 6 orbits). A finer grid
/// therefore does not converge as dt²; the reading jumps with where t* falls between grid points.</para>
///
/// <para>The refinement: take the grid winner t_j as a seed, then golden-section search for the
/// maximum of |K| on the bracket [t_(j−1), t_(j+1)] (clipped to the grid ends). Because
/// |K(t_j)| ≥ |K(t_(j±1))| the bracket contains a local maximum; the grid seed is what keeps the
/// search off a far local basin (the repo's lesson from the bounded-Brent trap, OpenArcs
/// <c>ptf_painter_pipeline</c>). The returned value is never below the grid maximum. A winner on
/// the first or last grid point means the maximiser may lie outside the window; that is reported,
/// not hidden.</para>
/// </summary>
public static class EnvelopeRefinement
{
    /// <summary>Bracket width at which the golden-section search stops. Near the maximum
    /// |K(t)| is flat to second order, so the envelope VALUE is accurate to O(tol²) while t*
    /// itself is resolved to tol.</summary>
    public const double TTolerance = 1e-9;

    private static readonly double InvPhi = (Math.Sqrt(5.0) - 1.0) / 2.0;

    /// <summary>Refine a grid maximum. <paramref name="tGrid"/> is ascending, <paramref name="iWin"/>
    /// the index of the grid point with the largest |K| (value <paramref name="kWin"/>),
    /// <paramref name="kAbs"/> evaluates |K| at any t.</summary>
    public static (double KMax, double TStar, bool AtWindowEdge) Refine(
        IReadOnlyList<double> tGrid, int iWin, double kWin, Func<double, double> kAbs)
    {
        if (tGrid is null || tGrid.Count == 0)
            throw new ArgumentException("tGrid must be non-empty.", nameof(tGrid));
        if (iWin < 0 || iWin >= tGrid.Count)
            throw new ArgumentOutOfRangeException(nameof(iWin));

        bool atEdge = iWin == 0 || iWin == tGrid.Count - 1;
        if (tGrid.Count == 1) return (kWin, tGrid[0], true);

        double a = tGrid[Math.Max(iWin - 1, 0)];
        double b = tGrid[Math.Min(iWin + 1, tGrid.Count - 1)];
        double c = b - InvPhi * (b - a);
        double d = a + InvPhi * (b - a);
        double fc = kAbs(c), fd = kAbs(d);
        while (b - a > TTolerance)
        {
            if (fc >= fd)
            {
                b = d; d = c; fd = fc;
                c = b - InvPhi * (b - a); fc = kAbs(c);
            }
            else
            {
                a = c; c = d; fc = fd;
                d = a + InvPhi * (b - a); fd = kAbs(d);
            }
        }

        double tStar = fc >= fd ? c : d;
        double kStar = Math.Max(fc, fd);
        if (kStar < kWin) return (kWin, tGrid[iWin], atEdge);
        return (kStar, tStar, atEdge);
    }
}
