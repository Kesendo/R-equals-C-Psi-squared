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
/// <item><b>RE-face</b> — for the uniform open XX chain with one dephased endpoint, the exact arithmetic of
/// the first-order-in-γ₀/J H-eigenmode rate comb is a_k = (4/(N+1))·sin²(kπ/(N+1)), with
/// α_k^full = γ₀·a_k + O(γ₀³/J²) and α_k^full/γ₀ = a_k + O((γ₀/J)²); its coefficients are
/// all rational iff N+1 ∈ {1,2,3,4,6}
/// (<see cref="NivenRationalityRootClaim.FirstOrderRateCombAllRational"/>). Exact Niven rationality belongs to the first-order coefficient comb. At finite γ₀/J the relative full-L
/// rate shift is O((γ₀/J)²), equivalently the absolute shift δα_k = O(γ₀³/J²), so no exact finite-γ₀/J
/// full-L rationality is claimed.</item>
/// <item><b>IM-face</b> — the band edge 2cos(π/(N+1)) has algebraic degree φ_euler(2(N+1))/2
/// (<see cref="NivenRationalityRootClaim.BandEdgeDegree"/>): degree at most 2 iff N≤5, with a rational
/// value for N≤2 and a quadratic surd exactly for N∈{3,4,5}.</item>
/// <item><b>F6 face</b> — the F6 Q-edge gain 1+cos(π/N), historically called the V-Effect gain, golden at N=5
/// (docs/carbon/OFF_NIVEN_AS_WAVE_BREAKING.md).</item>
/// </list>
///
/// <para>N=4 is the first golden on both SE faces (band edge IS φ = 2cos(π/5); the first-order comb carries √5). The exact
/// minimal-polynomial proof (degree = φ_euler(2(N+1))/2) is the sympy gate in
/// simulations/niven_rationality_root.py; this witness is the live recomputation, citing it. J-independent
/// (pure arithmetic), so it takes no J/γ.</para></summary>
public sealed class NivenRationalityRootWitness : IInspectable
{
    private static readonly CultureInfo Inv = CultureInfo.InvariantCulture;
    private const double Tol = 1e-12;

    /// <summary>The band edge 2cos(π/(N+1)).</summary>
    private static double BandEdge(int n) => 2.0 * Math.Cos(Math.PI / (n + 1));

    /// <summary>An exact first-order F65 coefficient
    /// a_k = (4/(N+1))·sin²(kπ/(N+1)), not a finite-γ₀/J full-L rate.</summary>
    private static double FirstOrderRateCoefficient(int n, int k)
    {
        double s = Math.Sin(Math.PI * k / (n + 1));
        return 4.0 / (n + 1) * s * s;
    }

    /// <summary>The F6 Q-edge gain 1 + cos(π/N).</summary>
    private static double VGain(int n) => 1.0 + Math.Cos(Math.PI / n);

    /// <summary>The independent numeric route: recognize the band-edge value 2cos(π/(N+1)) and return its name
    /// plus its algebraic degree (0 = unrecognized, i.e. a higher-degree surd). This is the second route the
    /// IM-face cross-checks against the totient degree <see cref="NivenRationalityRootClaim.BandEdgeDegree"/> —
    /// value recognition vs φ_euler(2(N+1))/2 must agree wherever the value is a named surd (the sibling-witness
    /// two-route convention, not a literal-vs-static check).</summary>
    private static (string Name, int NumericDegree) RecognizeBandEdge(double be)
    {
        if (Math.Abs(be - Math.Round(be)) < Tol) return ($"{Math.Round(be):0} (rational)", 1);
        if (Math.Abs(be - Math.Sqrt(2.0)) < Tol) return ("√2", 2);
        if (Math.Abs(be - (1.0 + Math.Sqrt(5.0)) / 2.0) < Tol) return ("φ = (1+√5)/2", 2);
        if (Math.Abs(be - Math.Sqrt(3.0)) < Tol) return ("√3", 2);
        return ("", 0);   // a higher-degree surd, not a named quadratic
    }

    /// <summary>The shared gate label: a tag for the node header and a sentence for its summary.</summary>
    private static (string Tag, string Note) Gate(bool ok) =>
        ok ? ("gate ✓", "Gate passes.") : ("GATE FIRED", "Gate FIRED — investigate.");

    public string DisplayName => "NivenRationalityRootWitness (the SE band edge + first-order endpoint comb, three faces, live)";

    public string Summary =>
        "the Niven root computed live (typed home: NivenRationalityRootClaim): Niven's theorem on π/(N+1) is the " +
        "number-theoretic ceiling on the SE cyclotomic forms. For the uniform open XX chain with one dephased endpoint, " +
        "a_k = (4/(N+1))·sin²(kπ/(N+1)); α_k^full = γ₀·a_k + O(γ₀³/J²), equivalently " +
        "α_k^full/γ₀ = a_k + O((γ₀/J)²). Exact Niven rationality belongs to the first-order coefficient comb, " +
        "which is rational for every k iff N+1∈{1,2,3,4,6}; at finite γ₀/J the " +
        "relative full-L rate shift is O((γ₀/J)²), equivalently the absolute shift δα_k = O(γ₀³/J²), so no exact " +
        "finite-γ₀/J full-L rationality is claimed. IM " +
        "band edge degree φ_euler(2(N+1))/2 (degree at most 2 iff N≤5: rational for N≤2, " +
        "quadratic surd exactly for N∈{3,4,5}, namely √2/φ/√3); F6 Q-edge gain golden at N=5. " +
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

    /// <summary>IM-face: the band edge degree ladder. A genuine two-route cross-check (the sibling-witness
    /// convention): the totient degree φ_euler(2(N+1))/2 must agree with the independent numeric value
    /// recognition (√2/φ/√3/rational) wherever the value is named; the degree-at-most-two cutoff is N≤5,
    /// with rational values at N≤2 and quadratic surds exactly at N∈{3,4,5}.</summary>
    private InspectableNode ImFace()
    {
        var rows = new List<IInspectable>();
        var rationalNs = new List<int>();
        var degreeAtMostTwoNs = new List<int>();
        bool routesAgree = true;
        for (int n = 1; n <= 8; n++)
        {
            double be = BandEdge(n);
            int deg = NivenRationalityRootClaim.BandEdgeDegree(n);   // route A: the totient φ_euler(2(N+1))/2
            var (name, numericDeg) = RecognizeBandEdge(be);          // route B: independent value recognition
            if (deg == 1) rationalNs.Add(n);
            if (deg <= 2) degreeAtMostTwoNs.Add(n);
            if (numericDeg != 0 && numericDeg != deg) routesAgree = false;   // the two routes disagree
            rows.Add(new InspectableNode($"N={n}",
                summary: $"2cos(π/{n + 1}) = {be.ToString("0.######", Inv)} = {(name.Length > 0 ? name : $"degree {deg}")} (degree {deg})"));
        }
        bool ok = rationalNs.SequenceEqual(new[] { 1, 2 })
                  && degreeAtMostTwoNs.SequenceEqual(new[] { 1, 2, 3, 4, 5 })
                  && NivenRationalityRootClaim.BandEdgeDegree(6) == 3
                  && routesAgree;
        var (tag, note) = Gate(ok);
        return new InspectableNode($"IM-face: the band edge 2cos(π/(N+1)) [{tag}]",
            summary: $"algebraic degree φ_euler(2(N+1))/2 (agreeing with the independent numeric value recognition): " +
                     $"degree at most 2 iff N≤5; rational iff N≤2 (last N=2 = 1), quadratic surd exactly for " +
                     $"N∈{{3,4,5}} (√2, φ, √3), " +
                     $"degree ≥3 (first cubic) from N=6. {note}",
            children: rows);
    }

    /// <summary>RE-face: the first-order endpoint-comb coefficient rationality, gated against
    /// N+1 ∈ {1,2,3,4,6}; it does not classify finite-γ₀/J full-L rates.</summary>
    private InspectableNode ReFace()
    {
        var rows = new List<IInspectable>();
        var rationalNs = new List<int>();
        for (int n = 1; n <= 8; n++)
        {
            bool rat = NivenRationalityRootClaim.FirstOrderRateCombAllRational(n);
            if (rat) rationalNs.Add(n);
            string sample = n <= 5
                ? "{" + string.Join(", ", FirstOrderCombList(n)) + "}"
                : "(cyclotomic)";
            rows.Add(new InspectableNode($"N={n}",
                summary: $"first-order coefficient comb: {(rat ? "all rational" : "contains irrational coefficients")} (N+1={n + 1}); a_k = {sample}"));
        }
        bool ok = rationalNs.SequenceEqual(new[] { 1, 2, 3, 5 });   // N+1 ∈ {2,3,4,6} of {1,2,3,4,6}
        var (tag, note) = Gate(ok);
        return new InspectableNode($"RE-face: the first-order endpoint coefficient comb a_k (F65/F99) [{tag}]",
            summary: $"For the uniform open XX chain with one dephased endpoint, a_k = (4/(N+1))·sin²(kπ/(N+1)); " +
                     $"α_k^full = γ₀·a_k + O(γ₀³/J²), equivalently α_k^full/γ₀ = a_k + O((γ₀/J)²). Exact Niven rationality " +
                     $"belongs to the first-order coefficient comb, whose coefficients are all rational iff N+1 ∈ {{1,2,3,4,6}} — N=3 the last rational " +
                     $"before the gap, N=5 a rational island. At finite γ₀/J the relative full-L rate shift is O((γ₀/J)²), " +
                     $"equivalently the absolute shift δα_k = O(γ₀³/J²), so " +
                     $"no exact finite-γ₀/J full-L rationality is claimed (the coefficient proof is the sympy gate in " +
                     $"niven_rationality_root.py). {note}",
            children: rows);
    }

    private static string[] FirstOrderCombList(int n)
    {
        var r = new string[n];
        for (int k = 1; k <= n; k++) r[k - 1] = FirstOrderRateCoefficient(n, k).ToString("0.####", Inv);
        return r;
    }

    /// <summary>F6 face: the Q-edge gain, golden at N=5 (not N=4) — the third face, angle π/N.</summary>
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
        var (tag, note) = Gate(ok);
        return new InspectableNode($"F6 face: the Q-edge gain 1+cos(π/N) [{tag}]",
            summary: "Niven-rational iff N ∈ {2,3}; for N ≥ 4 the named constant of the next ring polygon — silver " +
                     "(N=4), GOLDEN (N=5), √3 (N=6). Its golden is N=5 because the angle is π/N, not π/(N+1) " +
                     $"(docs/carbon/OFF_NIVEN_AS_WAVE_BREAKING.md). 'First golden' is angle-convention-dependent. {note}",
            children: rows);
    }

    /// <summary>The hinge: N=4 first golden on both SE faces (band edge = φ exactly; first-order comb carries √5).</summary>
    private InspectableNode Hinge()
    {
        double be4 = BandEdge(4), phi = (1.0 + Math.Sqrt(5.0)) / 2.0;
        double rate41 = FirstOrderRateCoefficient(4, 1), goldenRate = (5.0 - Math.Sqrt(5.0)) / 10.0;   // (4/5)sin²(π/5) = (5−√5)/10
        bool ok = Math.Abs(be4 - phi) < Tol && Math.Abs(rate41 - goldenRate) < Tol && !NivenRationalityRootClaim.FirstOrderRateCombAllRational(4);
        var (tag, note) = Gate(ok);
        return new InspectableNode($"the hinge: N=4 = first golden [{tag}]",
            summary: $"N=4 is the FIRST GOLDEN on both SE faces: the band edge IS φ (2cos(π/5) = {be4.ToString("0.######", Inv)} " +
                     $"= (1+√5)/2 = {phi.ToString("0.######", Inv)}), and the first-order comb carries √5 (a_1 = (4/5)sin²(π/5) = " +
                     $"{rate41.ToString("0.######", Inv)} = (5−√5)/10). The golden ratio is forced by the cyclotomic " +
                     $"geometry, not chosen. {note}");
    }

    private static InspectableNode TwoRoots() =>
        new InspectableNode("two roots, not one (the n3_special_cases answer)",
            summary: "this number-theoretic root (Niven on π/(N+1)) governs the band-edge / first-order endpoint-rate-comb / F6 Q-edge-gain " +
                     "family. The (n,n)/{0,2} filling-maximality specials (the N=3 extra modes, the N=4 (2,2) anomalies) " +
                     "are a SEPARATE combinatorial root; star-3 = path P_3 a third, graph-theoretic coincidence. So the " +
                     "arc's 'do they share a root' = no: at least two real roots (arithmetic + combinatorial).");

    public InspectablePayload Payload => InspectablePayload.Empty;
}
