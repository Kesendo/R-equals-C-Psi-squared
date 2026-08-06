using System.Numerics;
using RCPsiSquared.Core.F1;

namespace RCPsiSquared.Compute;

/// <summary>
/// Mirror symmetry analysis of decay rate spectra: the F1 palindrome read on the
/// oscillatory rate subset, through the canonical check rather than a local one.
/// </summary>
public static class MirrorAnalysis
{
    /// <summary>Machine epsilon for double, the unit the eigensolver's backward error is measured in.</summary>
    private const double Eps = 2.220446049250313e-16;

    public record MirrorResult(
        double MaxPairingDistance,
        double DistanceInEps,
        double RateRadius,
        int Count,
        double Center
    );

    /// <summary>
    /// The F1 pairing distance of a rate multiset around a center, computed by the canonical
    /// check <see cref="F1SpectrumStatistics.MaxF1PairingDistance"/>: greedy nearest-neighbour
    /// matching WITH REMOVAL between the multiset and its own F1 reflection. Returns a DISTANCE;
    /// ~0 means the multiset is closed under the reflection, and a large value means it is not.
    ///
    /// <para>The embedding, since the canonical check reflects λ ↦ −2σ − λ on complex
    /// eigenvalues while this method receives decay RATES d = −Re λ: pass λ = −d, for which
    /// −2σ − λ = d − 2σ = −(2σ − d) is exactly the embedding of the mirrored rate. So the
    /// canon computes the rate palindrome with no local re-derivation. The equality of the two
    /// routes is EXACT, not approximate, because d ↦ −d is an isometry and preserves the index
    /// order the greedy matcher walks; <c>simulations/f1_rate_embedding_check.py</c> is the
    /// producer that reproduces it on chain and star at N=2..5.</para>
    ///
    /// <para>This is a PROJECTION and stays weaker than the full-spectrum check: the imaginary
    /// parts are discarded upstream by <c>Liouvillian.ExtractRates</c>, so two eigenvalues at
    /// the same rate and different frequency are indistinguishable here. The full-spectrum
    /// verification lives in the block-spectrum runs
    /// (<c>simulations/results/f1_n8_n9_metrics/</c>);
    /// <c>docs/proofs/MIRROR_SYMMETRY_PROOF.md</c> scopes this suite as a sanity check.</para>
    ///
    /// <para><see cref="MirrorResult.DistanceInEps"/> is a SCALE-FREE CONVENIENCE, not the
    /// eigensolver's floor, and must not be read as one. The backward-error scale of
    /// <c>zgeev</c> on a strongly non-normal L is eps·‖L‖ amplified by the eigenvalue condition
    /// number, and ‖L‖ tracks the coupling J, not γ; max|rate| is neither. Measured on the star
    /// at γ=0.05 the ratio grows by more than a decade across N=2..6, so it is not constant and
    /// there is no constant-ratio law here to gate against.
    /// <see cref="MirrorResult.MaxPairingDistance"/> in absolute units is the number to read;
    /// the ratio is for putting two runs at the same N side by side.</para>
    ///
    /// <para><paramref name="center"/> is σ, the sum of the per-site dephasing rates, NOT half
    /// of it: the rate reflection is d ↦ 2σ − d.</para>
    /// </summary>
    public static MirrorResult CheckSymmetry(List<double> rates, double center)
    {
        if (rates is null) throw new ArgumentNullException(nameof(rates));

        // An empty multiset has no distance. Returning 0.0 would print the BEST possible
        // reading for a spectrum that was never scored, which is how a metric fails quiet.
        // NOTE a divergence, not a shared contract: the Python sibling
        // fw.max_f1_pairing_distance returns 0.0 on an empty input rather than NaN, because
        // np.all(np.isfinite([])) is True and its loop never runs. This side is the stricter
        // of the two; closing the gap on the Python side is open.
        if (rates.Count == 0)
            return new MirrorResult(double.NaN, double.NaN, 0.0, 0, center);

        var embedded = new Complex[rates.Count];
        double radius = 0.0;
        for (int i = 0; i < rates.Count; i++)
        {
            embedded[i] = new Complex(-rates[i], 0.0);
            double m = Math.Abs(rates[i]);
            if (m > radius) radius = m;
        }

        double distance = F1SpectrumStatistics.MaxF1PairingDistance(embedded, center);

        // A zero radius means every rate is 0, so there is no eps·ρ to divide by. Report the
        // raw distance and let the caller see NaN rather than an invented unit.
        double inEps = radius > 0.0 ? distance / (Eps * radius) : double.NaN;

        return new MirrorResult(distance, inEps, radius, rates.Count, center);
    }
}
