using System.Globalization;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>The qudit partial palindrome (Tier1Derived, 2026-06-11, F121): for d > 2 the
/// dephasing spectrum is not random but partial, and the partial pairing is the symmetric
/// overlap of the disagreement-count multiplicity. Under full-Cartan dephasing the d levels
/// are equidistant, so the decay rate of a coherence |i⟩⟨j| is −2γ·Hamming(i, j), the SAME
/// rate ladder as the qubit (the Absorption Theorem one dimension up). Only the multiplicity
/// per Hamming rung differs:
/// <code>
///   c_k = d^N · C(N, k) · (d−1)^k,   Σ_k c_k = d^{2N}.
/// </code>
/// The palindrome reflects rung k ↔ N−k, so the dissipator's paired ceiling is
/// <code>
///   paired(d, N) = Σ_k d^N · C(N, k) · (d−1)^{min(k, N−k)},
/// </code>
/// which equals d^{2N} (everything pairs) IFF d = 2: the (d−1)^k tilt base equals 1 only
/// there. This is the d² − 2d = 0 necessity of <see cref="QubitNecessityPi2Inheritance"/>
/// re-seen as the unique fully-paired column of an N-indexed family. For d = 3, N = 2 the
/// counts are c = [9, 36, 36], paired = 54/81, excess = 27. The tilt base d − 1 is exactly
/// the parent's per-site decaying : immune ratio (d² − d) : d = (d − 1) : 1, raised to the
/// number of disagreeing sites.
///
/// <para><b>Open (not claimed here):</b> the ceiling is the dissipator's partial palindrome.
/// The full L = L_H + L_D does not respect it (the SU(3) Heisenberg redistributes real parts;
/// N = 2 reaches 60/81 at center −3γ, exceeding the ceiling). That is the open follow-up,
/// carried in <c>docs/proofs/PROOF_QUDIT_PARTIAL_PALINDROME.md</c> §4, not in this claim.</para>
///
/// <para><b>Self-check battery (integer-exact, built in the ctor):</b> c_k vs brute
/// enumeration at d = 3, N = 2 and d = 2, N = 3; the closed-form ceiling vs an independent
/// combinatorial pairing on the (d, N) ∈ {2,3,4}×{1,2,3} grid; d = 2 full in every column
/// while d > 2 is strictly partial; the qutrit 54/81 with excess 27; the tilt base equal to
/// the parent's per-site decaying : immune ratio at d = 2. Mirrors
/// <c>simulations/qutrit_partial_palindrome.py</c>.</para>
///
/// <para>Anchors: <c>docs/proofs/PROOF_QUDIT_PARTIAL_PALINDROME.md</c> +
/// <c>simulations/qutrit_partial_palindrome.py</c>.</para></summary>
public sealed class QuditPartialPalindromeCeiling : Claim
{
    /// <summary>Typed parent: the d² − 2d = 0 necessity, whose per-site d : (d²−d) split is
    /// the k = 0 / k = 1 rung of this ladder and whose unique solution d = 2 is the unique
    /// fully-paired column here.</summary>
    public QubitNecessityPi2Inheritance QubitNecessity { get; }

    /// <summary>One integer-exact self-check.</summary>
    public readonly record struct BatteryCase(string Name, string Detail, string Expected, string Actual)
    {
        public bool Passes => string.Equals(Expected, Actual, StringComparison.Ordinal);
    }

    public IReadOnlyList<BatteryCase> Cases { get; }
    public int PassCount => Cases.Count(c => c.Passes);

    public QuditPartialPalindromeCeiling(QubitNecessityPi2Inheritance qubitNecessity)
        : base("The qudit partial palindrome: under full-Cartan dephasing the d levels are " +
               "equidistant so the rate is −2γ·Hamming(i,j) (the qubit ladder), and the partial " +
               "pairing for d > 2 is the symmetric overlap of the disagreement-count multiplicity " +
               "c_k = d^N·C(N,k)·(d−1)^k under k ↔ N−k: paired(d,N) = Σ_k d^N·C(N,k)·(d−1)^min(k,N−k), " +
               "= d^{2N} iff d = 2 (the unique fully-paired column; d² − 2d = 0 re-seen). " +
               "d = 3, N = 2: 54/81, excess 27. Tier1Derived (closed-form combinatorial identity, exact)",
               Tier.Tier1Derived,
               "docs/proofs/PROOF_QUDIT_PARTIAL_PALINDROME.md + " +
               "simulations/qutrit_partial_palindrome.py")
    {
        QubitNecessity = qubitNecessity ?? throw new ArgumentNullException(nameof(qubitNecessity));
        Cases = BuildBattery(qubitNecessity);
    }

    /// <summary>Binomial coefficient C(n, k) as an exact long.</summary>
    public static long Binom(int n, int k)
    {
        if (k < 0 || k > n) return 0;
        if (k > n - k) k = n - k;
        long r = 1;
        for (int i = 0; i < k; i++)
        {
            r = r * (n - i) / (i + 1);
        }
        return r;
    }

    /// <summary>Integer power b^e for non-negative e.</summary>
    public static long IntPow(long b, int e)
    {
        long r = 1;
        for (int i = 0; i < e; i++) r *= b;
        return r;
    }

    /// <summary>Coherence multiplicity at Hamming distance k: c_k = d^N·C(N,k)·(d−1)^k.</summary>
    public static long Multiplicity(int d, int N, int k) =>
        IntPow(d, N) * Binom(N, k) * IntPow(d - 1, k);

    /// <summary>Total Liouville dimension d^{2N}.</summary>
    public static long Total(int d, int N) => IntPow(d, 2 * N);

    /// <summary>The dissipator's paired ceiling: Σ_k d^N·C(N,k)·(d−1)^{min(k, N−k)}.</summary>
    public static long Ceiling(int d, int N)
    {
        long sum = 0;
        for (int k = 0; k <= N; k++)
            sum += IntPow(d, N) * Binom(N, k) * IntPow(d - 1, Math.Min(k, N - k));
        return sum;
    }

    /// <summary>The paired fraction Ceiling / Total.</summary>
    public double PairedFraction(int d, int N) => (double)Ceiling(d, N) / Total(d, N);

    /// <summary>Independent pairing count over a rung-multiplicity array (rung k ↔ N−k):
    /// each two-rung pair contributes 2·min(c_k, c_{N−k}); the middle rung (N even) self-mirrors.
    /// Used to cross-check the closed-form <see cref="Ceiling"/>.</summary>
    public static long CombinatorialPairing(IReadOnlyList<long> c)
    {
        int N = c.Count - 1;
        long paired = 0;
        for (int k = 0; k <= N; k++)
        {
            int m = N - k;
            if (k < m) paired += 2 * Math.Min(c[k], c[m]);
            else if (k == m) paired += c[k];
        }
        return paired;
    }

    /// <summary>Brute-force rung counts by enumerating all d^{2N} (i, j) basis pairs.</summary>
    private static long[] BruteCounts(int d, int N)
    {
        var c = new long[N + 1];
        int total = (int)IntPow(d, N);
        var digits = new int[N];
        for (int i = 0; i < total; i++)
        {
            DigitsOf(i, d, N, digits);
            var iDig = (int[])digits.Clone();
            for (int j = 0; j < total; j++)
            {
                DigitsOf(j, d, N, digits);
                int k = 0;
                for (int l = 0; l < N; l++) if (iDig[l] != digits[l]) k++;
                c[k]++;
            }
        }
        return c;
    }

    private static void DigitsOf(int value, int d, int N, int[] outDigits)
    {
        for (int l = 0; l < N; l++)
        {
            outDigits[l] = value % d;
            value /= d;
        }
    }

    private static IReadOnlyList<BatteryCase> BuildBattery(QubitNecessityPi2Inheritance parent)
    {
        var cases = new List<BatteryCase>();

        // (a) c_k vs brute enumeration, d = 3, N = 2.
        var brute32 = BruteCounts(3, 2);
        var form32 = new[] { Multiplicity(3, 2, 0), Multiplicity(3, 2, 1), Multiplicity(3, 2, 2) };
        cases.Add(new BatteryCase(
            Name: "c_k vs brute enumeration (d=3, N=2)",
            Detail: "c_k = d^N·C(N,k)·(d−1)^k against full enumeration of the 81 basis pairs",
            Expected: "9,36,36",
            Actual: brute32.SequenceEqual(form32)
                ? string.Join(",", form32)
                : "brute=" + string.Join(",", brute32) + " form=" + string.Join(",", form32)));

        // (b) c_k vs brute enumeration, d = 2, N = 3 (the qubit binomial × 8).
        var brute23 = BruteCounts(2, 3);
        var form23 = Enumerable.Range(0, 4).Select(k => Multiplicity(2, 3, k)).ToArray();
        cases.Add(new BatteryCase(
            Name: "c_k vs brute enumeration (d=2, N=3)",
            Detail: "qubit rungs c_k = 8·C(3,k); the (d−1)^k tilt is 1",
            Expected: "8,24,24,8",
            Actual: brute23.SequenceEqual(form23)
                ? string.Join(",", form23)
                : "brute=" + string.Join(",", brute23) + " form=" + string.Join(",", form23)));

        // (c) Ceiling closed form vs independent combinatorial pairing on the grid.
        int gridOk = 0, gridTot = 0;
        foreach (int d in new[] { 2, 3, 4 })
            foreach (int N in new[] { 1, 2, 3 })
            {
                gridTot++;
                var c = Enumerable.Range(0, N + 1).Select(k => Multiplicity(d, N, k)).ToList();
                if (Ceiling(d, N) == CombinatorialPairing(c)) gridOk++;
            }
        cases.Add(new BatteryCase(
            Name: "ceiling = combinatorial pairing on (d,N) grid {2,3,4}×{1,2,3}",
            Detail: "the closed-form Σ d^N·C(N,k)·(d−1)^min(k,N−k) vs the min-based rung pairing",
            Expected: "9/9",
            Actual: gridOk.ToString(CultureInfo.InvariantCulture) + "/" + gridTot.ToString(CultureInfo.InvariantCulture)));

        // (d) d = 2 is the unique fully-paired column.
        bool d2Full = Enumerable.Range(1, 3).All(N => Ceiling(2, N) == Total(2, N));
        bool dGt2Partial = new[] { 3, 4 }.All(d => Enumerable.Range(1, 3).All(N => Ceiling(d, N) < Total(d, N)));
        cases.Add(new BatteryCase(
            Name: "d = 2 unique fully-paired column (d² − 2d = 0 re-seen)",
            Detail: "Ceiling(2,N) = 4^N for N = 1..3; Ceiling(d,N) < d^{2N} strictly for d ∈ {3,4}, N = 1..3",
            Expected: "d=2 full, d>2 partial",
            Actual: (d2Full && dGt2Partial) ? "d=2 full, d>2 partial"
                : $"d2Full={d2Full}, dGt2Partial={dGt2Partial}"));

        // (e) The qutrit anchor: 54/81 with excess 27.
        cases.Add(new BatteryCase(
            Name: "qutrit N=2: paired = 54/81, excess = 27",
            Detail: "rung 0 (×9) pairs into rung 2 (×36) leaving 27; rung 1 (×36) self-mirrors",
            Expected: "54/81 excess 27",
            Actual: Ceiling(3, 2) == 54 && Total(3, 2) == 81 && (Total(3, 2) - Ceiling(3, 2)) == 27
                ? "54/81 excess 27"
                : $"{Ceiling(3, 2)}/{Total(3, 2)} excess {Total(3, 2) - Ceiling(3, 2)}"));

        // (f) The tilt base equals the parent's per-site decaying : immune ratio at d = 2.
        double parentRatio = parent.DecayingOpsPerSite / parent.ImmuneOpsPerSite;   // = 1 at d = 2
        bool tiltOk = Math.Abs(parentRatio - (2 - 1)) < 1e-12;
        cases.Add(new BatteryCase(
            Name: "tilt base (d−1) = parent's per-site decaying : immune ratio at d = 2",
            Detail: "the multiplicity tilt base d−1 is (d²−d)/d; at d = 2 it is 1, the parent's 2:2 balance",
            Expected: "ratio = 1",
            Actual: tiltOk ? "ratio = 1" : "ratio = " + parentRatio.ToString("0.###", CultureInfo.InvariantCulture)));

        return cases;
    }

    public override string DisplayName =>
        "The qudit partial palindrome: paired = Σ d^N·C(N,k)·(d−1)^min(k,N−k), full iff d=2";

    public override string Summary =>
        "under full-Cartan dephasing the d levels are equidistant (rate −2γ·Hamming, the qubit ladder); " +
        "the partial pairing for d > 2 is the symmetric overlap of c_k = d^N·C(N,k)·(d−1)^k under k ↔ N−k, " +
        "ceiling Σ_k d^N·C(N,k)·(d−1)^min(k,N−k) = d^{2N} iff d = 2 (the d²−2d=0 column of an N-family); " +
        $"d = 3, N = 2 gives {Ceiling(3, 2)}/{Total(3, 2)}; {PassCount}/{Cases.Count} battery PASS ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("rate ladder",
                summary: "full-Cartan dephasing makes the d levels equidistant: rate(|i⟩⟨j|) = −2γ·Hamming(i,j), identical to the qubit (Absorption Theorem)");
            yield return new InspectableNode("multiplicity tilt",
                summary: "c_k = d^N·C(N,k)·(d−1)^k coherences at Hamming distance k; tilt base d−1 = per-site (d²−d):d ratio");
            yield return new InspectableNode("ceiling",
                summary: "paired(d,N) = Σ_k d^N·C(N,k)·(d−1)^min(k,N−k); = d^{2N} iff d = 2");
            yield return InspectableNode.RealScalar("PairedFraction(d=3, N=2)", PairedFraction(3, 2));
            yield return InspectableNode.RealScalar("PairedFraction(d=4, N=2)", PairedFraction(4, 2));
            yield return new InspectableNode("qutrit anchor",
                summary: $"d = 3, N = 2: c = [9,36,36], paired = {Ceiling(3, 2)}/{Total(3, 2)}, excess = {Total(3, 2) - Ceiling(3, 2)}");
            yield return new InspectableNode("d = 2 uniqueness",
                summary: "the unique fully-paired column; (d−1)^k = 1 only at d = 2; this IS d² − 2d = 0 (parent QubitNecessityPi2Inheritance) as an N-family");
            yield return new InspectableNode("open: interacting spectrum",
                summary: "the full L = L_H + L_D exceeds the ceiling (H redistributes real parts; N = 2 reaches 60/81 at −3γ); its closed form is open (PROOF §4)");
            foreach (var c in Cases)
                yield return new InspectableNode(c.Name,
                    summary: $"{c.Detail}; expected {c.Expected}, got {c.Actual}, " + (c.Passes ? "PASS" : "FAIL"));
        }
    }
}
