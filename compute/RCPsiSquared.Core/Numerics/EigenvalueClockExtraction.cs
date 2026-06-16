using System;
using System.Collections.Generic;
using System.Numerics;

namespace RCPsiSquared.Core.Numerics;

/// <summary>Extract the Takt clock (the slowest decay rate + the coherence frequency at it) from a
/// Liouvillian eigenvalue spectrum. The Takt gap is the slowest nonzero decay rate (−Re λ, the inverse of
/// the longest coherence lifetime); ω_mem is the max |Im λ| among the eigenvalues whose decay rate is within
/// <see cref="GapMatch"/> of the strict gap.
///
/// <para>Shared so the reading is defined once: <c>Symphony.Clock</c> reads it off the one shared evolution
/// (capped at MaxN=5), while <c>TopologyBandEdgeWitness</c> builds the dense Liouvillian directly to reach
/// N=6 (above that cap) and reads the IDENTICAL clock here. Pure spectrum post-processing — no Liouvillian
/// build — so it lives in Core/Numerics next to the other spectral helpers.</para></summary>
public static class EigenvalueClockExtraction
{
    /// <summary>A decay rate −Re λ above this counts as a genuine (nonzero) decaying mode; at or below it the
    /// mode is steady (a conserved operator, not part of the gap).</summary>
    public const double GapFloor = 1e-9;

    /// <summary>ω_mem is read off the eigenvalues whose decay rate lies within this of the strict gap (the
    /// modes "at" the gap).</summary>
    public const double GapMatch = 1e-6;

    /// <summary>(Takt gap, ω_mem) from a Liouvillian spectrum: gap = the slowest nonzero decay rate
    /// (−Re λ, minimized over modes with rate &gt; <see cref="GapFloor"/>); ω_mem = max |Im λ| among the modes
    /// whose rate is within <see cref="GapMatch"/> of that gap. Returns (0, 0) when no mode decays above the
    /// floor (a stopped clock, e.g. γ = 0). The sequence is enumerated twice, so pass a materialized
    /// collection (an array or a MathNet vector), not a lazy query.</summary>
    public static (double Gap, double Omega) ExtractClockFromSpectrum(IEnumerable<Complex> eigenvalues)
    {
        if (eigenvalues is null) throw new ArgumentNullException(nameof(eigenvalues));
        double gap = double.PositiveInfinity;
        foreach (var e in eigenvalues) { double rate = -e.Real; if (rate > GapFloor && rate < gap) gap = rate; }
        if (double.IsInfinity(gap)) return (0.0, 0.0);
        double omega = 0.0;
        foreach (var e in eigenvalues)
            if (Math.Abs(-e.Real - gap) <= GapMatch) omega = Math.Max(omega, Math.Abs(e.Imaginary));
        return (gap, omega);
    }
}
