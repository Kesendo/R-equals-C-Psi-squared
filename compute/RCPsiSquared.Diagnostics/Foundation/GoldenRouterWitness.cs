using System.Globalization;
using System.Numerics;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Pauli;
using RCPsiSquared.Diagnostics.F87;

namespace RCPsiSquared.Diagnostics.Foundation;

/// <summary>The live lab for F116, the golden ceiling router (2026-06-11). Where the typed claim
/// <see cref="PalindromeSoftCertifierClaim"/> asserts the §7.12 soft-certifier ceiling is CLOSED at zero
/// (the 6 → 4 → 2 → 0 arc), this witness re-runs the existing router machinery at inspect time and watches
/// it close: it feeds the canonical Z-middle ceiling cases (XZX+XZY+YZX and its X↔Y sibling YZY+XZY+YZX)
/// through <see cref="PalindromeSoftCertifier"/> and shows the per-term lens DECLINING them while the
/// window-summed golden router (Stufe B′) CERTIFIES them, the 2 → 0 step made visible.
///
/// <para>It then evaluates the METALLIC FAMILY live (PROOF_CEILING_GOLDEN_ROUTER.md §8): for each weight
/// c the metallic period-4 router q_l[a, a, b, b] with frame ratio r is built and the window-summed
/// anticommutator residual ‖{W, S}‖_F is root-found (<see cref="KBodyPalindromeRouting.LiveMetallicRatio"/>)
/// to recover r(c) from the construction ALONE; that live ratio is then compared against the closed-form
/// metallic mean r(c) = (c + √(c² + 4))/2 (golden φ at c = 1, silver 1 + √2 at c = 2, bronze (3+√13)/2 at
/// c = 3, r = 1 at c = 0). The agreement is a genuine check: the live ratio comes from where the router's
/// anticommutator vanishes, the closed form from the metallic-mean quadratic.</para>
///
/// <para>Children: a ceiling-arc node (the 6 → 4 → 2 → 0 step ladder, each canonical case mapped to its
/// certifying strategy live), a per-c node for c ∈ {0, 1, 2, 3} with the live r(c) vs the closed form and a
/// MATCH/MISMATCH verdict, and a residual-curve node (‖{W, S}‖_F over a grid of r at c = 1, dipping to zero
/// at φ) for <c>--draw</c>. The c-grid is settable via <c>--router-c</c> (comma-separated). Summary, e.g.
/// "golden router closes the ceiling (Z-middle pair → RoutingWindowSummed); metallic family r(c) live =
/// closed form at c = 0, 1, 2, 3".</para>
///
/// <para>Anchors: <c>docs/proofs/PROOF_CEILING_GOLDEN_ROUTER.md</c> (incl. §8 metallic family) +
/// <c>simulations/metallic_router_family.py</c> (8 blocks) +
/// <see cref="KBodyPalindromeRouting"/> (the reused router) + <see cref="PalindromeSoftCertifier"/>.</para></summary>
public sealed class GoldenRouterWitness : IInspectable
{
    private static readonly CultureInfo Inv = CultureInfo.InvariantCulture;

    /// <summary>The tolerance for the live-vs-closed-form metallic ratio MATCH verdict. The live ratio is a
    /// bisection root (80 steps over a width-12 bracket ⟹ ≈ 1e-23 in principle), but the residual floor and
    /// float Kronecker products limit the practical agreement; 1e-9 is comfortably met.</summary>
    public const double MatchTolerance = 1e-9;

    /// <summary>The chain length used for the certifier ceiling re-run. Soft is established at N ≥ 4 for the
    /// Z-middle cases, and the golden router is N-free, so 4 suffices and is cheap.</summary>
    public const int CeilingN = 4;

    /// <summary>The weights c at which the metallic family is evaluated live (default {0, 1, 2, 3}).</summary>
    public IReadOnlyList<double> Weights { get; }

    public GoldenRouterWitness(IReadOnlyList<double>? weights = null)
    {
        Weights = weights is { Count: > 0 } ? weights : new[] { 0.0, 1.0, 2.0, 3.0 };
    }

    private static PauliTerm T(string label) => new(PauliLabel.Parse(label), Complex.One);
    private static List<PauliTerm> H(params string[] labels) => labels.Select(T).ToList();

    /// <summary>The two canonical Z-middle ceiling cases (the 2 → 0 survivors): the golden case and its
    /// X↔Y sibling, both soft and certified only by the window-summed router.</summary>
    private static readonly (string Name, string[] Labels)[] ZMiddleCases =
    {
        ("XZX+XZY+YZX (golden)", new[] { "XZX", "XZY", "YZX" }),
        ("YZY+XZY+YZX (golden sibling)", new[] { "YZY", "XZY", "YZX" }),
    };

    /// <summary>One metallic-family reading: the weight c, the live router-derived ratio, the closed-form
    /// metallic mean, and their absolute difference.</summary>
    public readonly record struct MetallicReading(double C, double LiveRatio, double ClosedForm)
    {
        public double Difference => Math.Abs(LiveRatio - ClosedForm);
        public bool Matches => Difference < MatchTolerance;
        public string Mean => C switch
        {
            1.0 => "golden φ",
            2.0 => "silver 1+√2",
            3.0 => "bronze (3+√13)/2",
            0.0 => "45° frame r=1",
            _ => "metallic",
        };
    }

    /// <summary>Evaluate the metallic family live across <see cref="Weights"/>: for each c, root-find the
    /// frame ratio from the router's window-summed anticommutator and compare to the closed-form mean.</summary>
    public IReadOnlyList<MetallicReading> Readings()
    {
        var readings = new List<MetallicReading>(Weights.Count);
        foreach (double c in Weights)
        {
            double live = KBodyPalindromeRouting.LiveMetallicRatio(c);
            double closed = KBodyPalindromeRouting.MetallicMean(c);
            readings.Add(new MetallicReading(c, live, closed));
        }
        return readings;
    }

    /// <summary>For the two canonical Z-middle cases: the per-term router verdict (declines) and the
    /// window-summed certificate (certifies). Returns (name, perTermRoutes, certified, strategy).</summary>
    public IReadOnlyList<(string Name, bool PerTermRoutes, bool Certified, PalindromeSoftCertifier.SoftStrategy Strategy)> CeilingClose()
    {
        var rows = new List<(string, bool, bool, PalindromeSoftCertifier.SoftStrategy)>();
        foreach (var (name, labels) in ZMiddleCases)
        {
            var terms = H(labels);
            bool perTerm = KBodyPalindromeRouting.Routes(terms, CeilingN);          // the per-term lens (declines)
            var cert = PalindromeSoftCertifier.Certify(terms, CeilingN);            // the full certifier (window-summed certifies)
            rows.Add((name, perTerm, cert.Certified, cert.Strategy));
        }
        return rows;
    }

    public string DisplayName => $"GoldenRouterWitness (F116 live lab, metallic family at c ∈ {{{string.Join(", ", Weights.Select(c => c.ToString("0.###", Inv)))}}})";

    public string Summary
    {
        get
        {
            var close = CeilingClose();
            bool allCertified = close.All(r => r.Certified) && close.All(r => !r.PerTermRoutes);
            var readings = Readings();
            int matched = readings.Count(r => r.Matches);
            string cList = string.Join(", ", readings.Select(r => r.C.ToString("0.###", Inv)));
            return $"golden router closes the ceiling ({(allCertified ? "Z-middle pair → RoutingWindowSummed, per-term lens declines" : "REGRESSION")}); " +
                   $"metallic family r(c) live = closed form at c = {cList} ({matched}/{readings.Count} MATCH). " +
                   $"Built the period-4 router live and root-found r(c) from ‖{{W, S}}‖_F = 0.";
        }
    }

    public IEnumerable<IInspectable> Children
    {
        get
        {
            // 1. The ceiling arc 6 → 4 → 2 → 0, with the 2 → 0 step re-run live.
            var close = CeilingClose();
            var arcChildren = new List<IInspectable>
            {
                new InspectableNode("6 → 4 (continuous-uniform router)",
                    summary: "XIX+XIY+YIX, YIY+XIY+YIX route via a continuous-uniform per-site Q (the once-counted pair)"),
                new InspectableNode("4 → 2 (single-site-field router)",
                    summary: "IXI+IIY+YII, IYI+IIX+XII are sums of single-site transverse fields (SingleSiteField strategy)"),
            };
            foreach (var (name, perTerm, certified, strategy) in close)
                arcChildren.Add(new InspectableNode(
                    displayName: $"2 → 0 LIVE: {name}",
                    summary: $"per-term router Routes = {perTerm} (declines, cross-template cancellation invisible to it); " +
                             $"full certifier: certified = {certified}, strategy = {strategy}; " +
                             (perTerm == false && certified && strategy == PalindromeSoftCertifier.SoftStrategy.RoutingWindowSummed
                                 ? "CLOSES (window-summed golden router certifies what the per-term lens cannot)."
                                 : "REGRESSION.")));
            yield return new InspectableNode(
                displayName: "the ceiling arc 6 → 4 → 2 → 0 (the 2 → 0 step re-run live)",
                summary: "the §7.12 soft-certifier structural ceiling, closed at zero (F116); the last 2 Z-middle cases " +
                         "certify via the window-summed golden router, re-run here at inspect time.",
                children: arcChildren);

            // 2. The metallic family, evaluated live per c.
            var readings = Readings();
            foreach (var r in readings)
                yield return new InspectableNode(
                    displayName: $"c = {r.C.ToString("0.###", Inv)} ({r.Mean})",
                    summary: $"live r(c) from the router (root of ‖{{W, S}}‖_F) = {r.LiveRatio.ToString("0.############", Inv)}; " +
                             $"closed-form metallic mean (c+√(c²+4))/2 = {r.ClosedForm.ToString("0.############", Inv)}; " +
                             $"|diff| = {r.Difference.ToString("0.##e+00", Inv)} ⟹ {(r.Matches ? "MATCH" : "MISMATCH")}.");

            // 3. The residual curve ‖{W, S}‖_F over a grid of r at c = 1, dipping to zero at φ (for --draw).
            const double cCurve = 1.0;
            const int gridPoints = 61;
            double rLo = 0.4, rHi = 3.0;
            var rAxis = new double[gridPoints];
            var residual = new double[gridPoints];
            for (int i = 0; i < gridPoints; i++)
            {
                double rr = rLo + (rHi - rLo) * i / (gridPoints - 1);
                rAxis[i] = rr;
                residual[i] = KBodyPalindromeRouting.WindowSummedAnticommutatorNorm(cCurve, rr);
            }
            double phi = KBodyPalindromeRouting.MetallicMean(cCurve);
            yield return new InspectableNode(
                displayName: "the residual ‖{W, S}‖_F over r (c = 1)",
                summary: $"the window-summed anticommutator norm as a function of the frame ratio r; it dips to zero " +
                         $"exactly at r = φ = {phi.ToString("0.######", Inv)} (the golden mean), the live root.",
                payload: new InspectablePayload.Curve("‖{W,S}‖_F vs r", rAxis, residual, "frame ratio r", "‖{W,S}‖_F"));
        }
    }

    public InspectablePayload Payload => InspectablePayload.Empty;
}
