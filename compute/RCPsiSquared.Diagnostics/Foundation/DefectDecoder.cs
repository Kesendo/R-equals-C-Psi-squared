using System;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;
using RCPsiSquared.Core.ChainSystems;
using RCPsiSquared.Core.Inspection;

namespace RCPsiSquared.Diagnostics.Foundation;

/// <summary>The defect decoder (M2b): inverting the painters' α-profile. The painters READ a defect
/// (its location in where the per-site α-profile breaks the mirror, its strength in the f-scale); this
/// decoder CLOSES the loop. Given an observed per-site α-profile, it returns the bond the defect sits
/// on and how strong it is, with the least-squares residual as a confidence.
///
/// <para>The mechanism is a linear dictionary. In the perturbative window α ≈ 1 + δJ·f⁽ᵇ⁾, so the
/// per-bond per-site f-profiles — calibrated once at δJ_cal through the SAME PaintersMovement pipeline
/// the Symphony uses (one source of α semantics, never a reimplemented fit) — ARE the dictionary. To
/// decode an observed α: for each candidate bond, project (α−1) onto that bond's f-profile by least
/// squares (δĴ = (α−1)·f / (f·f)), measure the residual ‖(α−1) − δĴ·f‖², and return the bond with the
/// smallest residual. Location = the best-fitting f-profile; strength = its projection coefficient; the
/// residual carries the decoder's own error bar (it grows orders of magnitude from the linear window's
/// center to its edge).</para>
///
/// <para>Honest scope (inherited from the painters): the canonical PTF protocol — XY chain, bonding-mode
/// carrier, Z-dephasing, |δJ| ≲ 0.1·J — and N in 3..5 (location needs an interior bond; the dense
/// d²×d² eigendecomposition caps N). Design:
/// docs/superpowers/specs/2026-06-12-defect-decoder-m2b-design.md.</para></summary>
public sealed class DefectDecoder
{
    private static readonly CultureInfo Inv = CultureInfo.InvariantCulture;

    /// <summary>The canonical calibration defect strength. The decoder calibrates its dictionary here
    /// (and a self-test/read-back must pass THIS, not the movement's own δJ — calibrating at the
    /// movement's δJ is circular).</summary>
    public const double DefaultDeltaJCal = 0.02;

    /// <summary>A reading is AMBIGUOUS when the runner-up bond's residual is within this factor of the
    /// winner's: the linear dictionary cannot cleanly separate the two letters. Calibrated from the two
    /// verified regimes at N=5: the near-anti-collinear mismatch case (edge bond 3 weakened reads almost
    /// like complementary interior bond 1 strengthened) sits at residual ratio ≈ 1.5 (ambiguous); the
    /// clean positive-defect cases sit at ratio ≫ 10 (unambiguous). The 3.0 threshold flags the former
    /// and clears the latter.</summary>
    public const double AmbiguityFactor = 3.0;

    public int N { get; }
    public double J { get; }
    public double Gamma { get; }
    public double DeltaJCal { get; }

    /// <summary>The dictionary: one f-profile (per-site f_i = (α_i − 1)/δJ_cal) per bond b in 0..N−2,
    /// read off a calibration painters run with a δJ_cal defect on bond b.</summary>
    public IReadOnlyList<double[]> Dictionary => _dictionary;
    private readonly double[][] _dictionary;

    /// <summary>How many calibration performances were built (N−1 painters movements, one per bond).
    /// The decoder owns this counter; it does NOT add to any host movement's BuildCount. A host that
    /// calibrates a decoder for its own read-back pays this many extra movements — counted here, openly,
    /// so the existing movement-level BuildCount accounting (4 trajectories per movement) stays exact.</summary>
    public int CalibrationBuildCount { get; }

    private DefectDecoder(int n, double j, double gamma, double deltaJCal, double[][] dictionary, int calibrationBuildCount)
    {
        N = n; J = j; Gamma = gamma; DeltaJCal = deltaJCal;
        _dictionary = dictionary;
        CalibrationBuildCount = calibrationBuildCount;
    }

    /// <summary>Calibrate a decoder for an (N, J, γ) system: for each bond b in 0..N−2, build a
    /// canonical-protocol painters run (XY chain, bonding-mode carrier, δJ = δJ_cal on bond b) and read
    /// its per-site f-profile through PaintersMovement. The resulting f-dictionary is the decoder's
    /// codebook. N is guarded to 3..5 (the painters' canonical range: location needs an interior bond,
    /// and the dense Liouvillian caps the top).</summary>
    public static DefectDecoder Calibrate(int n, double j, double gamma, double deltaJCal = DefaultDeltaJCal)
    {
        if (n < 3 || n > 5)
            throw new ArgumentOutOfRangeException(nameof(n),
                $"the defect decoder is specified for N in 3..5 (the painters' canonical range; location " +
                $"needs an interior bond and the dense Liouvillian caps the top); got {n}");

        int bonds = n - 1;
        var dictionary = new double[bonds][];
        int buildCount = 0;
        for (int b = 0; b < bonds; b++)
        {
            // The SAME construction the canonical painters tests use — one source of α semantics.
            var s = new Symphony(n: n, j: j, gamma: gamma, hType: HamiltonianType.XY,
                                 initialState: InitialStateKind.BondingMode, defectBond: b, deltaJ: deltaJCal);
            var pm = ((IInspectable)s).Children.OfType<PaintersMovement>().Single();
            if (!pm.HasLenses)
                throw new InvalidOperationException(
                    $"calibration bond {b} declined: {pm.DeclineReason} (the canonical PTF protocol must apply)");
            dictionary[b] = pm.F.ToArray();
            buildCount += pm.BuildCount;
        }

        return new DefectDecoder(n, j, gamma, deltaJCal, dictionary, buildCount);
    }

    /// <summary>The decode result: the identified (winning) bond, the recovered strength δĴ, and the
    /// least-squares residual ‖(α−1) − δĴ·f‖² of the winning fit (the confidence — smaller is more
    /// certain); plus the RUNNER-UP bond, its δĴ and residual, and the ambiguity flag. A linear
    /// dictionary cannot always separate sign+location: at N=5 an edge bond and the complementary
    /// interior bond have nearly anti-collinear f-profiles (weakening the edge reads almost like
    /// strengthening the interior). When the runner-up residual is within
    /// <see cref="AmbiguityFactor"/>× the winner's, <see cref="IsAmbiguous"/> is set and the honest
    /// reading reports BOTH candidate letters, never the winner alone.</summary>
    public readonly record struct DecodeResult(
        int Bond, double DeltaJ, double Residual,
        int RunnerUpBond, double RunnerUpDeltaJ, double RunnerUpResidual,
        bool IsAmbiguous);

    /// <summary>Decode an observed per-site α-profile: for each dictionary bond, least-squares project
    /// (α−1) onto that bond's f-profile (δĴ = (α−1)·f / (f·f)) and score the residual
    /// ‖(α−1) − δĴ·f‖²; return the minimum-residual bond with its δĴ and residual, the runner-up bond,
    /// and the ambiguity flag (runner-up residual within <see cref="AmbiguityFactor"/>× the winner's).</summary>
    public DecodeResult Decode(IReadOnlyList<double> alphaObserved)
    {
        if (alphaObserved is null) throw new ArgumentNullException(nameof(alphaObserved));
        if (alphaObserved.Count != N)
            throw new ArgumentException(
                $"the observed α-profile has {alphaObserved.Count} sites; the decoder expects {N}", nameof(alphaObserved));

        var delta = new double[N];
        for (int i = 0; i < N; i++) delta[i] = alphaObserved[i] - 1.0;
        return ProjectAndScore(delta, _dictionary);
    }

    /// <summary>The shared least-squares scoring core. Project the observed per-site vector onto each
    /// dictionary column (δĴ = obs·d / d·d), score the SQUARED residual ‖obs − δĴ·d‖², and return the
    /// minimum-residual bond with its δĴ, the runner-up, and the ambiguity flag (runner-up residual within
    /// <see cref="AmbiguityFactor"/>× the winner's). Both decode paths call this; the ONLY path-specific
    /// step is forming `obs` (the α path subtracts 1; the deviation path does not). The residual is squared,
    /// so the reported ratio is the SQUARE of the residual-norm ratio (review §11).</summary>
    private DecodeResult ProjectAndScore(double[] obs, double[][] dict)
    {
        int bestBond = -1, secondBond = -1;
        double bestDeltaJ = 0.0, secondDeltaJ = 0.0;
        double bestResidual = double.PositiveInfinity, secondResidual = double.PositiveInfinity;

        for (int b = 0; b < dict.Length; b++)
        {
            var f = dict[b];
            double dotFD = 0.0, dotFF = 0.0;
            for (int i = 0; i < N; i++) { dotFD += f[i] * obs[i]; dotFF += f[i] * f[i]; }
            double deltaJ = dotFF > 0.0 ? dotFD / dotFF : 0.0;

            double residual = 0.0;
            for (int i = 0; i < N; i++)
            {
                double e = obs[i] - deltaJ * f[i];
                residual += e * e;
            }

            if (residual < bestResidual)
            {
                secondResidual = bestResidual; secondBond = bestBond; secondDeltaJ = bestDeltaJ;
                bestResidual = residual; bestBond = b; bestDeltaJ = deltaJ;
            }
            else if (residual < secondResidual)
            {
                secondResidual = residual; secondBond = b; secondDeltaJ = deltaJ;
            }
        }

        bool ambiguous = secondBond >= 0 && secondResidual < AmbiguityFactor * bestResidual;
        return new DecodeResult(bestBond, bestDeltaJ, bestResidual,
                                secondBond, secondDeltaJ, secondResidual, ambiguous);
    }
}
