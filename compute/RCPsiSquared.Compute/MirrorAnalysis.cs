using System.Numerics;
using RCPsiSquared.Core.F1;

namespace RCPsiSquared.Compute;

/// <summary>
/// Mirror symmetry analysis of decay rate spectra: the F1 palindrome read on the
/// oscillatory rate subset, through the canonical check rather than a local one.
/// </summary>
public static class MirrorAnalysis
{
    /// <summary>Machine epsilon for double. It is the unit <see cref="MirrorResult.DistanceInEps"/>
    /// is expressed in, and NOT a claim that the eigensolver's backward error is what that ratio
    /// measures; see the note on <see cref="CheckSymmetry"/> for why it is not.</summary>
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
    /// check <see cref="F1SpectrumStatistics.MaxF1RatePairingDistance"/>: greedy nearest-neighbour
    /// matching WITH REMOVAL between the multiset and its own F1 reflection. Returns a DISTANCE;
    /// ~0 means the multiset is closed under the reflection, and a large value means it is not.
    ///
    /// <para>The rate embedding is NOT this class's: it lives in Core beside the canon, as
    /// <see cref="F1SpectrumStatistics.MaxF1RatePairingDistance"/>, and is pinned there. What
    /// this method adds over that call is the eps ratio, the rate radius and count it is read
    /// beside, and the empty-input policy below.</para>
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

        // An empty multiset has no distance. Three answers are live and they are NOT one
        // contract, so state which is which. Core THROWS on empty, in both the rate and the
        // complex entry point. This method returns NaN instead, because it is a reporting
        // wrapper whose caller prints a table row and should see a hole rather than a crash.
        // The Python sibling fw.max_f1_pairing_distance still returns 0.0, the best possible
        // reading, because np.all(np.isfinite([])) is True and its loop never runs; that one
        // is a genuine gap and closing it is open. What no route may do is return 0.0 here.
        if (rates.Count == 0)
            return new MirrorResult(double.NaN, double.NaN, 0.0, 0, center);

        double radius = 0.0;
        for (int i = 0; i < rates.Count; i++)
        {
            double m = Math.Abs(rates[i]);
            if (m > radius) radius = m;
        }

        double distance = F1SpectrumStatistics.MaxF1RatePairingDistance(rates, center);

        // A zero radius means every rate is 0, so there is no eps·ρ to divide by. Report the
        // raw distance and let the caller see NaN rather than an invented unit.
        double inEps = radius > 0.0 ? distance / (Eps * radius) : double.NaN;

        return new MirrorResult(distance, inEps, radius, rates.Count, center);
    }
}
