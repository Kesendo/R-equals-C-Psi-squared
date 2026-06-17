using System;
using System.Collections.Generic;
using System.Globalization;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Symmetry;

namespace RCPsiSquared.Diagnostics.Foundation;

/// <summary>The Niven root computed live (typed home: <see cref="NivenRationalityRootClaim"/>). Niven's theorem
/// on the single-excitation cyclotomic angle π/(N+1) is the number-theoretic ceiling on the SE spectrum's closed
/// forms; this witness recomputes its three faces at inspect time and gates them:
///
/// <list type="bullet">
/// <item><b>RE-face</b> — the dissipator rates α_k/γ₀ = (4/(N+1))·sin²(kπ/(N+1)) are all rational iff
/// N+1 ∈ {1,2,3,4,6} (<see cref="NivenRationalityRootClaim.RatesAllRational"/>).</item>
/// <item><b>IM-face</b> — the band edge 2cos(π/(N+1)) has algebraic degree φ_euler(2(N+1))/2
/// (<see cref="NivenRationalityRootClaim.BandEdgeDegree"/>): rational iff N≤2, quadratic surd iff N≤5.</item>
/// <item><b>V-face</b> — the V-Effect gain 1+cos(π/N), golden at N=5 (docs/carbon/OFF_NIVEN_AS_WAVE_BREAKING.md).</item>
/// </list>
///
/// <para>N=4 is the first golden on both SE faces (band edge IS φ = 2cos(π/5); the rates carry √5). The exact
/// minimal-polynomial proof (degree = φ_euler(2(N+1))/2) is the sympy gate in
/// simulations/niven_rationality_root.py; this witness is the live recomputation, citing it. J-independent
/// (pure arithmetic), so it takes no J/γ.</para></summary>
public sealed class NivenRationalityRootWitness : IInspectable
{
    private static readonly CultureInfo Inv = CultureInfo.InvariantCulture;
    private const double Tol = 1e-12;

    /// <summary>The band edge 2cos(π/(N+1)).</summary>
    private static double BandEdge(int n) => 2.0 * Math.Cos(Math.PI / (n + 1));

    /// <summary>A dissipator rate α_k/γ₀ = (4/(N+1))·sin²(kπ/(N+1)) (= F65's SingleExcitationRate at γ₀=1).</summary>
    private static double Rate(int n, int k)
    {
        double s = Math.Sin(Math.PI * k / (n + 1));
        return 4.0 / (n + 1) * s * s;
    }

    /// <summary>The V-Effect gain 1 + cos(π/N).</summary>
    private static double VGain(int n) => 1.0 + Math.Cos(Math.PI / n);

    /// <summary>Name the band edge: rational value, √2/φ/√3, or "degree d".</summary>
    private static string BandEdgeIdentity(int n)
    {
        double be = BandEdge(n);
        int deg = NivenRationalityRootClaim.BandEdgeDegree(n);
        if (deg == 1) return $"{Math.Round(be):0} (rational)";
        if (Math.Abs(be - Math.Sqrt(2.0)) < Tol) return "√2";
        if (Math.Abs(be - (1.0 + Math.Sqrt(5.0)) / 2.0) < Tol) return "φ = (1+√5)/2";
        if (Math.Abs(be - Math.Sqrt(3.0)) < Tol) return "√3";
        return $"degree {deg}";
    }

    public string DisplayName => "NivenRationalityRootWitness (the SE cyclotomic arithmetic, three faces, live)";

    public string Summary =>
        "the Niven root computed live (typed home: NivenRationalityRootClaim): Niven's theorem on π/(N+1) is the " +
        "number-theoretic ceiling on the SE spectrum's closed forms. RE rates rational iff N+1∈{1,2,3,4,6}; IM " +
        "band edge degree φ_euler(2(N+1))/2 (rational N≤2, quadratic surd √2/φ/√3 N≤5); V-Effect golden at N=5. " +
        "N=4 = first golden on both SE faces (band edge = φ).";

    public IEnumerable<IInspectable> Children
    {
        get
        {
            yield return ImFace();
            yield return ReFace();
            yield return VFace();
            yield return Hinge();
            yield return TwoRoots();
        }
    }

    /// <summary>IM-face: the band edge degree ladder, gated rational-iff-N≤2 and quadratic-surd-iff-N≤5.</summary>
    private InspectableNode ImFace()
    {
        var rows = new List<IInspectable>();
        var rationalNs = new List<int>();
        var quadraticNs = new List<int>();
        for (int n = 1; n <= 8; n++)
        {
            int deg = NivenRationalityRootClaim.BandEdgeDegree(n);
            if (deg == 1) rationalNs.Add(n);
            if (deg <= 2) quadraticNs.Add(n);
            rows.Add(new InspectableNode($"N={n}",
                summary: $"2cos(π/{n + 1}) = {BandEdge(n).ToString("0.######", Inv)} = {BandEdgeIdentity(n)} (degree {deg})"));
        }
        bool ok = string.Join(",", rationalNs) == "1,2"
                  && string.Join(",", quadraticNs) == "1,2,3,4,5"
                  && NivenRationalityRootClaim.BandEdgeDegree(6) == 3;
        return new InspectableNode($"IM-face: the band edge 2cos(π/(N+1)) [{(ok ? "gate ✓" : "GATE FIRED")}]",
            summary: $"algebraic degree φ_euler(2(N+1))/2: rational iff N≤2 (last N=2 = 1), a single quadratic surd " +
                     $"(√2, φ, √3 at N=3,4,5) iff N≤5, degree ≥3 (first cubic) from N=6. " +
                     $"{(ok ? "All gates pass." : "A gate FIRED — investigate.")}",
            children: rows);
    }

    /// <summary>RE-face: the dissipator-rate rationality, gated against N+1 ∈ {1,2,3,4,6}.</summary>
    private InspectableNode ReFace()
    {
        var rows = new List<IInspectable>();
        var rationalNs = new List<int>();
        for (int n = 1; n <= 8; n++)
        {
            bool rat = NivenRationalityRootClaim.RatesAllRational(n);
            if (rat) rationalNs.Add(n);
            string sample = n <= 5
                ? "{" + string.Join(", ", RatesList(n)) + "}"
                : "(cyclotomic)";
            rows.Add(new InspectableNode($"N={n}",
                summary: $"rates {(rat ? "rational" : "IRRATIONAL")} (N+1={n + 1}); α_k/γ₀ = {sample}"));
        }
        bool ok = string.Join(",", rationalNs) == "1,2,3,5";   // N+1 ∈ {2,3,4,6} of {1,2,3,4,6}
        return new InspectableNode($"RE-face: the dissipator rates α_k (F65/F99) [{(ok ? "gate ✓" : "GATE FIRED")}]",
            summary: $"α_k/γ₀ = (4/(N+1))·sin²(kπ/(N+1)) all rational iff N+1 ∈ {{1,2,3,4,6}} — N=3 the last rational " +
                     $"before the gap, N=5 a rational island. {(ok ? "Gate passes." : "Gate FIRED — investigate.")}",
            children: rows);
    }

    private static string[] RatesList(int n)
    {
        var r = new string[n];
        for (int k = 1; k <= n; k++) r[k - 1] = Rate(n, k).ToString("0.####", Inv);
        return r;
    }

    /// <summary>V-face: the V-Effect gain, golden at N=5 (not N=4) — the third face, angle π/N.</summary>
    private InspectableNode VFace()
    {
        // golden at N=5: 1+cos(π/5) = (5+√5)/4
        double v5 = VGain(5), goldenV = (5.0 + Math.Sqrt(5.0)) / 4.0;
        double v4 = VGain(4), silverV = 1.0 + Math.Sqrt(2.0) / 2.0;
        bool ok = Math.Abs(v5 - goldenV) < Tol && Math.Abs(v4 - silverV) < Tol;
        var rows = new List<IInspectable>
        {
            new InspectableNode("N=4 (square)", summary: $"V = 1+cos(π/4) = {v4.ToString("0.######", Inv)} = 1+√2/2 (silver)"),
            new InspectableNode("N=5 (pentagon)", summary: $"V = 1+cos(π/5) = {v5.ToString("0.######", Inv)} = (5+√5)/4 (GOLDEN)"),
            new InspectableNode("N=6 (hexagon)", summary: $"V = 1+cos(π/6) = {VGain(6).ToString("0.######", Inv)} = 1+√3/2 (√3)"),
        };
        return new InspectableNode($"V-face: the V-Effect gain 1+cos(π/N) [{(ok ? "gate ✓" : "GATE FIRED")}]",
            summary: "Niven-rational iff N ∈ {2,3}; for N ≥ 4 the named constant of the next ring polygon — silver " +
                     "(N=4), GOLDEN (N=5), √3 (N=6). Its golden is N=5 because the angle is π/N, not π/(N+1) " +
                     "(docs/carbon/OFF_NIVEN_AS_WAVE_BREAKING.md). 'First golden' is angle-convention-dependent.",
            children: rows);
    }

    /// <summary>The hinge: N=4 first golden on both SE faces (band edge = φ exactly; rates carry √5).</summary>
    private InspectableNode Hinge()
    {
        double be4 = BandEdge(4), phi = (1.0 + Math.Sqrt(5.0)) / 2.0;
        double rate41 = Rate(4, 1), goldenRate = (5.0 - Math.Sqrt(5.0)) / 10.0;   // (4/5)sin²(π/5) = (5−√5)/10
        bool ok = Math.Abs(be4 - phi) < Tol && Math.Abs(rate41 - goldenRate) < Tol && !NivenRationalityRootClaim.RatesAllRational(4);
        return new InspectableNode($"the hinge: N=4 = first golden [{(ok ? "gate ✓" : "GATE FIRED")}]",
            summary: $"N=4 is the FIRST GOLDEN on both SE faces: the band edge IS φ (2cos(π/5) = {be4.ToString("0.######", Inv)} " +
                     $"= (1+√5)/2 = {phi.ToString("0.######", Inv)}), and the rates carry √5 (α_1/γ₀ = (4/5)sin²(π/5) = " +
                     $"{rate41.ToString("0.######", Inv)} = (5−√5)/10). The golden ratio is forced by the cyclotomic " +
                     $"geometry, not chosen. {(ok ? "Gate passes." : "Gate FIRED — investigate.")}");
    }

    private static InspectableNode TwoRoots() =>
        new InspectableNode("two roots, not one (the n3_special_cases answer)",
            summary: "this number-theoretic root (Niven on π/(N+1)) governs the band-edge / dissipator-rate / V-Effect " +
                     "family. The (n,n)/{0,2} filling-maximality specials (the N=3 extra modes, the N=4 (2,2) anomalies) " +
                     "are a SEPARATE combinatorial root; star-3 = path P_3 a third, graph-theoretic coincidence. So the " +
                     "arc's 'do they share a root' = no: at least two real roots (arithmetic + combinatorial).");

    public InspectablePayload Payload => InspectablePayload.Empty;
}
