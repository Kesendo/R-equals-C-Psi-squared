using System.Numerics;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.BlockSpectrum.JordanWigner;

/// <summary>Two-shift F1 palindrome probe on the sparse JW Slater-pair Liouvillian block.
/// Runs <see cref="JwSlaterPairShiftInvertArnoldi"/> twice — once at σ_slow ≈ 0 to pull
/// the steady-state-adjacent (slow-decay) end, once at σ_fast = −2·Σγ − σ_slow to pull the
/// F1 mirror (fast-decay) end — and reports per-pair residuals on the F1 mirror map
/// <c>λ → −λ − 2·Σγ</c>.
///
/// <para><b>This is a probe, not a strict F1 prover.</b> Shift-invert Arnoldi returns the
/// largest-magnitude Ritz values of <c>(L − σI)^(−1)</c>, which approximates "K closest to
/// σ in L" but does not equal it at finite Krylov dimension — deeper Ritz values can capture
/// genuine eigenvalues that are NOT among the strict K-closest-to-σ. When that happens for
/// the slow run but not the fast run (or vice versa), the per-pair residual on the
/// non-K-closest stragglers is large, even though F1 is trivially satisfied at the
/// theorem-level for the full spectrum. Empirically at N=5 (2,2) dim=100, both runs
/// individually return genuine eigenvalues (checked against dense Evd) but their two top-K
/// sets are NOT F1-paired with each other. The pair-wise residual is therefore an Arnoldi-
/// convergence diagnostic, not a structural-F1 witness.</para>
///
/// <para><b>F1 (mirror-symmetry theorem, <c>docs/proofs/MIRROR_SYMMETRY_PROOF.md</c>):</b>
/// the spectrum of <c>L = −i[H, ·] + Σ_l γ_l·(Z_l ρ Z_l − ρ)</c> is closed under
/// <c>λ → −λ − 2·Σγ</c>. The chain XY Hamiltonian sits in the F87 "truly" class — no F1
/// Brecher — so the structural F1 prediction is independent of this probe. The probe's
/// value: when its per-pair residuals ARE small, the slow- and fast-end Arnoldi runs sampled
/// genuinely mirror-paired regions of the spectrum, providing a tight numerical sanity check
/// on the sparse pipeline. When residuals are large, the sample-pair sets diverged
/// (Arnoldi-depth issue), not the underlying F1.</para>
///
/// <para><b>Why this primitive exists.</b> At N=10 the half-filling sector (p_c=5, p_r=5)
/// has dimension 63 504; dense Evd needs ~64 GB and is infeasible. The two-shift sparse
/// path recovers the top-k slow modes AND their F1 mirror partners (subject to the Arnoldi
/// depth caveat above) — both ends inaccessible to the dense
/// <see cref="LiouvillianSectorSweep"/> which drops above-cap sectors entirely. For a
/// strict full-spectrum F1 verification at N=10 the structurally cleaner route is Prosen
/// third quantization (a 2N×2N Nambu spectral problem — Medvedyeva, Essler, Prosen 2016 —
/// not yet implemented in the repo; "Phase 3" of the N=10 plan).</para>
///
/// <para><b>L-conjugation in the residual.</b> The Liouvillian for Hermitian H + Hermitian-
/// preserving dissipator (Z-dephasing qualifies) has <c>λ ∈ spec ⟺ λ* ∈ spec</c>. When a
/// slow eigenvalue <c>a + bi</c> from a conjugate pair is captured by slow Arnoldi, the fast
/// Arnoldi (biased by its own σ) typically captures <c>−a − 2·Σγ − bi</c> (mirror of <c>λ*</c>),
/// not <c>−a − 2·Σγ + bi</c> (mirror of <c>λ</c>). Both are valid F1 partners. The residual
/// computation reports <c>min(|F1(λ) − μ|, |F1(λ*) − μ|)</c> across the fast set so the
/// conjugation choice does not pollute the result. Pure-real slow modes reduce trivially.</para>
///
/// <para><b>Tier outcome: <see cref="Tier.Tier1Derived"/>.</b> Composition of two Tier1
/// <see cref="JwSlaterPairShiftInvertArnoldi"/> extractions plus an arithmetic
/// nearest-neighbour distance over their union with conjugation flexibility. Each
/// extraction's eigenvalues are checked against dense Evd at small N in the test suite.</para>
///
/// <para>Anchor: <see cref="JwSlaterPairShiftInvertArnoldi"/> (per-end extraction) +
/// <see cref="F1.F1PalindromeIdentity"/> (the structural palindrome) +
/// <see cref="F1.F1OpenQuestions"/> (enumerated F1-Brecher mechanisms);
/// <c>docs/proofs/MIRROR_SYMMETRY_PROOF.md</c>.</para>
/// </summary>
public sealed class JwSlaterPairF1PalindromeProbe : Claim
{
    public JwSlaterPairSparseLBuilder Source { get; }
    public Complex SigmaSlow { get; }
    public Complex SigmaFast { get; }
    public int NumEigenvaluesRequested { get; }
    public int NumIterations { get; }
    public double InnerTolerance { get; }
    public int InnerMaxIterations { get; }
    public double SumGamma { get; }

    public JwSlaterPairShiftInvertArnoldi SlowExtraction { get; }
    public JwSlaterPairShiftInvertArnoldi FastExtraction { get; }

    /// <summary>Per-slow-λ residual on the F1 mirror map: for slow eigenvalue λ_i, the
    /// minimum of <c>|F1(λ_i) − μ|</c> and <c>|F1(λ_i*) − μ|</c> over the recovered fast set
    /// <c>{μ}</c>. Length equals <c>SlowExtraction.Eigenvalues.Length</c>.</summary>
    public double[] PerSlowResidual { get; }

    public double MaxF1Residual => PerSlowResidual.Length == 0 ? 0.0 : PerSlowResidual.Max();
    public double MeanF1Residual => PerSlowResidual.Length == 0 ? 0.0 : PerSlowResidual.Average();

    /// <summary>Number of slow eigenvalues whose per-pair residual is below
    /// <paramref name="tolerance"/>. A coarse "how many of the K pairs landed mirror-paired"
    /// counter — independent of which K stragglers diverged.</summary>
    public int CountMatchedPairs(double tolerance)
    {
        int n = 0;
        foreach (var r in PerSlowResidual) if (r < tolerance) n++;
        return n;
    }

    public static JwSlaterPairF1PalindromeProbe Build(JwSlaterPairSparseLBuilder source,
        int numEig, int numIter, int randomSeed,
        double innerTolerance, int innerMaxIter,
        Complex? sigmaSlow = null)
    {
        if (source is null) throw new ArgumentNullException(nameof(source));

        var sigmaS = sigmaSlow ?? new Complex(0.0, 1e-3);
        double sumGamma = source.GammaPerSite.Sum();
        // F1 mirror of the chosen slow shift: σ_fast = −σ_slow − 2·Σγ. This places σ_fast at
        // the same nearest-neighbour distance from its target eigenvalues as σ_slow is from
        // the slow modes, so the two shift-invert runs operate at matched conditioning.
        var sigmaF = new Complex(-2.0 * sumGamma - sigmaS.Real, -sigmaS.Imaginary);

        var slow = JwSlaterPairShiftInvertArnoldi.Build(source, sigmaS, numEig, numIter,
            randomSeed, innerTolerance, innerMaxIter);
        var fast = JwSlaterPairShiftInvertArnoldi.Build(source, sigmaF, numEig, numIter,
            randomSeed + 1, innerTolerance, innerMaxIter);

        var residuals = ComputeF1Residuals(slow.Eigenvalues, fast.Eigenvalues, sumGamma);

        return new JwSlaterPairF1PalindromeProbe(source, sigmaS, sigmaF, numEig, numIter,
            innerTolerance, innerMaxIter, sumGamma, slow, fast, residuals);
    }

    /// <summary>For each slow eigenvalue λ_i, the minimum of <c>|F1(λ_i) − μ|</c> and
    /// <c>|F1(λ_i*) − μ|</c> over the fast set. The dual-F1 form absorbs the L-conjugation
    /// degree of freedom (see class docstring).</summary>
    private static double[] ComputeF1Residuals(Complex[] slow, Complex[] fast, double sumGamma)
    {
        if (fast.Length == 0) return slow.Select(_ => double.PositiveInfinity).ToArray();
        var residuals = new double[slow.Length];
        for (int i = 0; i < slow.Length; i++)
        {
            var predDirect = new Complex(-slow[i].Real - 2.0 * sumGamma, -slow[i].Imaginary);
            var predConj = new Complex(-slow[i].Real - 2.0 * sumGamma, +slow[i].Imaginary);
            double best = double.MaxValue;
            for (int j = 0; j < fast.Length; j++)
            {
                double dDirect = (fast[j] - predDirect).Magnitude;
                double dConj = (fast[j] - predConj).Magnitude;
                double d = Math.Min(dDirect, dConj);
                if (d < best) best = d;
            }
            residuals[i] = best;
        }
        return residuals;
    }

    private JwSlaterPairF1PalindromeProbe(JwSlaterPairSparseLBuilder source,
        Complex sigmaSlow, Complex sigmaFast, int numEig, int numIter,
        double innerTol, int innerMaxIter, double sumGamma,
        JwSlaterPairShiftInvertArnoldi slow, JwSlaterPairShiftInvertArnoldi fast,
        double[] residuals)
        : base($"F1 two-shift palindrome probe on sparse L_JW (p_c={source.PCol}, p_r={source.PRow}, " +
               $"N={source.N}, dim={source.SectorDim}); top-{numEig} per end at σ_slow=({sigmaSlow.Real:G3}, {sigmaSlow.Imaginary:G3}), " +
               $"σ_fast=({sigmaFast.Real:G3}, {sigmaFast.Imaginary:G3}); 2·Σγ={2.0 * sumGamma:G4}; " +
               $"max F1 residual = {(residuals.Length == 0 ? 0.0 : residuals.Max()):G3}.",
               Tier.Tier1Derived,
               "compute/RCPsiSquared.Core/BlockSpectrum/JordanWigner/JwSlaterPairShiftInvertArnoldi.cs (per-end extraction) + " +
               "compute/RCPsiSquared.Core/F1/F1PalindromeIdentity.cs (structural palindrome) + " +
               "compute/RCPsiSquared.Core/F1/F1OpenQuestions.cs (enumerated F1-Brecher mechanisms); " +
               "docs/proofs/MIRROR_SYMMETRY_PROOF.md.")
    {
        Source = source;
        SigmaSlow = sigmaSlow;
        SigmaFast = sigmaFast;
        NumEigenvaluesRequested = numEig;
        NumIterations = numIter;
        InnerTolerance = innerTol;
        InnerMaxIterations = innerMaxIter;
        SumGamma = sumGamma;
        SlowExtraction = slow;
        FastExtraction = fast;
        PerSlowResidual = residuals;
    }

    public override string DisplayName =>
        $"F1 two-shift probe on L_JW (p_c={Source.PCol}, p_r={Source.PRow}, N={Source.N})";

    public override string Summary =>
        $"max F1 residual = {MaxF1Residual:G3}, mean = {MeanF1Residual:G3}, " +
        $"top-{NumEigenvaluesRequested}/end ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return SlowExtraction;
            yield return FastExtraction;
            yield return new InspectableNode("σ_slow", summary: $"({SigmaSlow.Real:G4}, {SigmaSlow.Imaginary:G4})");
            yield return new InspectableNode("σ_fast", summary: $"({SigmaFast.Real:G4}, {SigmaFast.Imaginary:G4})");
            yield return InspectableNode.RealScalar("2·Σγ", 2.0 * SumGamma, "G4");
            yield return InspectableNode.RealScalar("max F1 residual", MaxF1Residual, "G3");
            yield return InspectableNode.RealScalar("mean F1 residual", MeanF1Residual, "G3");
            for (int i = 0; i < PerSlowResidual.Length; i++)
            {
                var slow = SlowExtraction.Eigenvalues[i];
                var predDirect = new Complex(-slow.Real - 2.0 * SumGamma, -slow.Imaginary);
                var predConj = new Complex(-slow.Real - 2.0 * SumGamma, +slow.Imaginary);
                yield return new InspectableNode(
                    $"pair_{i}",
                    summary: $"slow λ=({slow.Real:F5}, {slow.Imaginary:F5}); F1 mirror candidates " +
                             $"({predDirect.Real:F5}, {predDirect.Imaginary:F5}) or " +
                             $"({predConj.Real:F5}, {predConj.Imaginary:F5}); min residual = {PerSlowResidual[i]:G3}");
            }
        }
    }
}
