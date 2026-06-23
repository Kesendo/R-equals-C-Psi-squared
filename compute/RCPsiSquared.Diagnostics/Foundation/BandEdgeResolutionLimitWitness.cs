using System;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;
using RCPsiSquared.Core.Inspection;
using F124 = RCPsiSquared.Diagnostics.Foundation.BandEdgeTransitionInvariantWitness;

namespace RCPsiSquared.Diagnostics.Foundation;

/// <summary>The defect-localization RESOLUTION LIMIT of F124's bond-to-mode map, computed live (typed home:
/// <c>BandEdgeResolutionLimitClaim</c>). F124's transition matrix M[b,k]=⟨ψ_k|V_b|ψ_1⟩ is the forward map of
/// a bond-recovery inverse problem (a signal/control engineer's ill-conditioned inverse; an optician's imaging
/// system). Read that way, F124's λ_min = E = σ_min²(M) is the optimal lower frame bound = the worst-case
/// reconstruction floor, and the condition number sets the resolution:
///
/// <para>  κ = λ_max/λ_min ~ N²  (the noise amplification),  contrast σ_max/σ_min = √κ ~ N,</para>
///
/// so a staggered (−1)^b zone-boundary (q=π) bond defect produces a mode-response √κ ~ N times weaker than a
/// band-edge defect of the same magnitude, i.e. it is √κ ~ N times harder to localize (matched-filter detection
/// SNR, confirmed in simulations/f124_inverse_problem_gate.py). The worst-conditioned bond direction is exactly
/// the staggered mode (F124's λ_min eigenvector), the q=π detail at the resolution cutoff, the optician's
/// diffraction limit. The reconstruction floor itself vanishes as σ_min = √E ~ (N+1)^(−3/2) (E·(N+1)³ → 4π²),
/// so the long chain goes singular along the K-direction.
///
/// <para>Three trades, one object (borrowing-a-discipline): signal/control = the ill-conditioned inverse / the
/// observability Gramian's least-observable direction (the K-partner is the unobservable channel); optics = the
/// modulation transfer function whose q=π detail transfers with vanishing contrast at the cutoff; the lower
/// frame bound = the distance to rank-collapse. F124 gives the floor (λ_min=E); this is its resolution reading.
/// The 1.5 DefectDecoder ambiguity is a SEPARATE matter (an α-rescaling parametrization artifact, NOT this
/// finite conditioning; see the f124_inverse_problem_resolution_seam arc).</para></summary>
public sealed class BandEdgeResolutionLimitWitness : IInspectable
{
    private static readonly CultureInfo Inv = CultureInfo.InvariantCulture;

    private static readonly int[] Ladder = { 3, 4, 5, 6, 7, 8, 9, 10, 11, 12 };

    /// <summary>σ_min(M) = √E, the worst-case reconstruction floor (the lower frame bound).</summary>
    public static double SigmaMin(int n) => Math.Sqrt(F124.Analyse(F124.Topo.Chain, n).LamMin);

    /// <summary>σ_max(M) = √λ_max, the best-conditioned (band-edge) direction's response magnitude.</summary>
    public static double SigmaMax(int n) => Math.Sqrt(F124.Analyse(F124.Topo.Chain, n).LamMax);

    /// <summary>κ = λ_max/λ_min (the Gram condition number, ~N²): the noise amplification of the inverse.</summary>
    public static double Kappa(int n) => F124.Analyse(F124.Topo.Chain, n).Kappa;

    /// <summary>The contrast / resolution ratio σ_max/σ_min = √κ ~ N: how many times weaker the staggered q=π
    /// defect's response is than a band-edge defect's, i.e. how many times harder it is to localize.</summary>
    public static double ContrastRatio(int n) => SigmaMax(n) / SigmaMin(n);

    /// <summary>The asymptotic floor constant E·(N+1)³ → 4π² = 39.478: σ_min = √E ~ (N+1)^(−3/2).</summary>
    public static double FloorConstant(int n) => F124.EndpointClosedForm(n) * Math.Pow(n + 1, 3);

    /// <summary>The worst-conditioned bond direction is the staggered (−1)^b mode (F124's λ_min eigenvector).</summary>
    public static bool WorstDirectionIsStaggered(int n) => F124.Analyse(F124.Topo.Chain, n).StaggeredIsGenuineMinimum;

    /// <summary>Log-log slope of f over the given N (least squares); the empirical scaling exponent.</summary>
    public static double LogLogSlope(IEnumerable<int> ns, Func<int, double> f, bool versusNPlus1 = false)
    {
        var xs = ns.Select(n => Math.Log(versusNPlus1 ? n + 1 : n)).ToArray();
        var ys = ns.Select(n => Math.Log(f(n))).ToArray();
        double mx = xs.Average(), my = ys.Average();
        double num = 0, den = 0;
        for (int i = 0; i < xs.Length; i++) { num += (xs[i] - mx) * (ys[i] - my); den += (xs[i] - mx) * (xs[i] - mx); }
        return num / den;
    }

    // ---- IInspectable ----

    public string DisplayName =>
        "BandEdgeResolutionLimitWitness (F124's conditioning = the defect-localization resolution limit √κ ~ N)";

    public string Summary
    {
        get
        {
            double k5 = Kappa(5), c5 = ContrastRatio(5);
            return
                "F124's bond-to-mode map M is a bond-recovery inverse problem; its conditioning is a defect-" +
                "LOCALIZATION RESOLUTION LIMIT (typed home: BandEdgeResolutionLimitClaim). σ_min=√E is the worst-" +
                "case reconstruction floor (the lower frame bound); κ=λ_max/λ_min ~ N² is the noise amplification; " +
                $"the contrast σ_max/σ_min=√κ ~ N (N=5: κ={k5.ToString("0.###", Inv)}, √κ={c5.ToString("0.###", Inv)}) " +
                "is how many times harder a staggered q=π zone-boundary defect is to localize than a band-edge one " +
                "(matched-filter SNR). The worst direction is the staggered mode (F124's λ_min eigenvector), the " +
                "q=π detail at the resolution cutoff (the optician's diffraction limit); the floor vanishes as " +
                "σ_min ~ (N+1)^(−3/2) (E·(N+1)³ → 4π²), the long chain going singular. One object in three trades: " +
                "the ill-conditioned inverse / the observability Gramian (signal/control) and the MTF cutoff (optics).";
        }
    }

    public IEnumerable<IInspectable> Children
    {
        get
        {
            yield return TheResolutionLadder();
            yield return TheAsymptoticFloor();
            yield return TheBorrowedReading();
        }
    }

    /// <summary>The resolution ladder: per N, the reconstruction floor σ_min=√E, the condition number κ~N², the
    /// contrast/resolution ratio √κ~N, and the worst direction (the staggered q=π mode).</summary>
    private InspectableNode TheResolutionLadder()
    {
        var rows = new List<IInspectable>();
        foreach (int n in Ladder)
        {
            double smin = SigmaMin(n), k = Kappa(n), c = ContrastRatio(n);
            bool worst = WorstDirectionIsStaggered(n);
            rows.Add(new InspectableNode($"chain N={n}",
                summary: $"σ_min=√E={smin.ToString("0.#####", Inv)} (the floor); κ=λ_max/λ_min={k.ToString("0.###", Inv)}; " +
                         $"contrast σ_max/σ_min=√κ={c.ToString("0.###", Inv)} (staggered defect is {c.ToString("0.##", Inv)}× harder " +
                         $"to localize); worst direction = staggered q=π: {worst}"));
        }
        double kSlope = LogLogSlope(new[] { 4, 5, 6, 7, 8, 9, 10, 11, 12 }, Kappa);
        double cSlope = LogLogSlope(new[] { 4, 5, 6, 7, 8, 9, 10, 11, 12 }, ContrastRatio);
        return new InspectableNode("the resolution ladder (κ ~ N², the contrast √κ ~ N, the staggered worst direction)",
            summary: $"the bond-localization resolution degrades with N: κ=λ_max/λ_min ~ N² (fitted exponent " +
                     $"{kSlope.ToString("0.##", Inv)}), the contrast σ_max/σ_min=√κ ~ N (fitted {cSlope.ToString("0.##", Inv)}). " +
                     $"The worst-conditioned bond direction is the staggered (−1)^b q=π zone-boundary mode (F124's λ_min " +
                     $"eigenvector) at every N: the longer the chain, the blinder it is to a fine staggered defect.");
    }

    /// <summary>The reconstruction floor's asymptotic form: σ_min=√E ~ (N+1)^(−3/2), E·(N+1)³ → 4π²=39.478.</summary>
    private InspectableNode TheAsymptoticFloor()
    {
        var rows = new List<IInspectable>();
        foreach (int n in new[] { 4, 8, 16, 32, 60 })
            rows.Add(new InspectableNode($"N={n}",
                summary: $"E·(N+1)³ = {FloorConstant(n).ToString("0.###", Inv)} (→ 4π² = {(4 * Math.PI * Math.PI).ToString("0.###", Inv)})"));
        double eSlope = LogLogSlope(Enumerable.Range(4, 57), F124.EndpointClosedForm, versusNPlus1: true);
        return new InspectableNode("the reconstruction floor σ_min = √E ~ (N+1)^(−3/2) (E·(N+1)³ → 4π²)",
            summary: $"the worst-case reconstruction floor E = (4/(N+1))sin²(π/(N+1)) vanishes as (N+1)⁻³ (fitted " +
                     $"exponent over N=4..60: {eSlope.ToString("0.###", Inv)}); E·(N+1)³ → 4π² = 39.478. So σ_min ~ " +
                     $"(N+1)^(−3/2) and the long-chain / ring limit goes singular along the staggered K-direction. " +
                     $"NOTE the small-N effective exponent in N (not N+1) is ≈ −2.5, a pre-asymptotic / wrong-variable " +
                     $"artifact (the robust law is the RATIO κ ~ N²); see the f124_inverse_problem_resolution_seam arc.");
    }

    private static InspectableNode TheBorrowedReading() =>
        new InspectableNode("one object, three trades (the borrowed reading)",
            summary: "the same conditioning seen by three disciplines that own this shape (borrowing-a-discipline). " +
                     "SIGNAL/inverse problems: an ill-conditioned inverse, the lower frame bound = 1/‖S⁻¹‖ = the " +
                     "worst-case reconstruction floor; κ the noise amplification; the discrete Picard condition. " +
                     "CONTROL: the observability Gramian S=MᵀM, λ_min its least-observable eigenvalue, the K-partner " +
                     "ψ_N the unobservable channel (the null column). OPTICS (Zeiss): the singular spectrum is a " +
                     "modulation transfer function whose q=π detail transfers with vanishing contrast at the cutoff " +
                     "(the diffraction limit). F124 gives the floor λ_min=E; this is its resolution-limit reading.");

    public InspectablePayload Payload => InspectablePayload.Empty;
}
