using System.Numerics;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Numerics;
using RCPsiSquared.Core.Symmetry;
using T = RCPsiSquared.Core.Symmetry.ThermalThirdMomentClaim.Term;

namespace RCPsiSquared.Diagnostics.Foundation;

/// <summary>The live lab for <see cref="ConditionedFieldPairingClaim"/> (F137, the conditioned field): H = h(X₀ + X₀Z₁) +
/// w·Z₁, site 0 cooled at a and heated at b, site 1 cooled alone at g, the generator 2L + c in exact Gaussian integers
/// (<see cref="ExactThermalGenerator"/>).
///
/// <para>The proof, on a road of its own. L is linear in H and in the rates, so 2L + c = h·P_h + w·P_w + a·P_a + b·P_b +
/// g·P_g with five integer pieces (by construction; the generator is compared with it at two points as a sanity check); the triangularity is read on each
/// piece, entry by entry, so it holds at every parameter. The block identities are
/// polynomial identities in (h, w, a, b, g) of degree at most four in each variable (every entry is linear in the
/// parameters, the blocks are 4 × 4), so agreement on the grid {0, …, 4}⁵ of 3125 points proves them: the
/// field-free block's characteristic polynomial is (y − (a + b − g))·(y + g)²·(y + a + b + g), the one-sided blocks
/// satisfy p₀₁(y) = p₁₀(−y), and the field's block vanishes at y = a + b + g and y = g with trace 4g; the two roots
/// differ as polynomials (a + b is not the zero polynomial), so both linear factors divide it and the trace fixes the
/// remaining pair's sum (in 2L + c these
/// are the claim's −g, −s₀ − g twice, −2s₀ − g; 0, −s₀ and a pair of mean −3s₀/2). The registry's gate proves the same
/// with h, w and the rates as symbols. The converse is proven on the same grid: Tr((2L + c)³) has degree three, and at
/// every grid point it is the parent's law, so the coefficient of y¹³, −Tr((2L + c)³)/3 by Newton's identities, is
/// −64h²(a + b − 2g) at every parameter.</para>
///
/// <para>Then the second route: at two points with a + b = 2g, h ≠ 0 and w ≠ 0, for the mirrored setting (w ≠ 0) and
/// for the turn to Y (w = 0), every odd coefficient of the full characteristic polynomial and every odd power sum through the dimension 16
/// vanish; at points off the line the coefficient of y¹³ is the parent's law, −64h²(a + b − 2g), and not zero. The
/// controls cancel the third moment and break at the fifth: a ZZ term, heating on site 1, the field left on in the
/// other sector, a second neighbour conditioned on site 1.</para></summary>
public sealed class ConditionedFieldPairingWitness : IInspectable
{
    public sealed record Setting(string Label, IReadOnlyList<T> Terms, IReadOnlyList<long> Cool, IReadOnlyList<long> Heat,
        (long H, long A, long B, long G)? XForm);

    public sealed record PairingReading(Setting Setting, IReadOnlyList<int> NonzeroOddCoefficients, IReadOnlyList<int> NonzeroOddSums,
        BigInteger Thirteenth, BigInteger? Law);

    public sealed record ControlReading(Setting Setting, BigInteger Third, BigInteger Fifth, BigInteger Law)
    {
        public bool CancelsThenBreaks => Third.IsZero && Law.IsZero && !Fifth.IsZero;
    }

    private static Setting X(string label, long h, long w, long a, long b, long g) =>
        new(label, new[] { new T("XI", h), new T("XZ", h), new T("IZ", w) }, new[] { a, g }, new[] { b, 0L }, (h, a, b, g));

    public static IReadOnlyList<Setting> OnLine { get; } = new[]
    {
        X("h = 2, w = 3, site 0 at 3/1, site 1 cooled at 2", 2, 3, 3, 1, 2),
        X("h = 1, w = −1, site 0 cooled alone at 4, site 1 cooled at 2", 1, -1, 4, 0, 2),
        new Setting("mirrored: X₀ − X₀Z₁ + 2Z₁, site 0 at 1/3, site 1 heated alone at 2",
            new[] { new T("XI", 1), new T("XZ", -1), new T("IZ", 2) }, new long[] { 1, 0 }, new long[] { 3, 2 }, null),
        new Setting("turned: 3(Y₀ + Y₀Z₁), site 0 at 2/2, site 1 cooled at 2",
            new[] { new T("YI", 3), new T("YZ", 3) }, new long[] { 2, 2 }, new long[] { 2, 0 }, null),
    };

    public static IReadOnlyList<Setting> OffLine { get; } = new[]
    {
        X("h = 2, w = 3, site 0 at 3/1, site 1 cooled at 1", 2, 3, 3, 1, 1),
        X("h = 1, w = 0, site 0 at 1/1, site 1 cooled at 3", 1, 0, 1, 1, 3),
    };

    public static IReadOnlyList<Setting> Controls { get; } = new[]
    {
        new Setting("ZZ added: X₀ + X₀Z₁ + Z₀Z₁, site 0 at 1/1, site 1 cooled at 1",
            new[] { new T("XI", 1), new T("XZ", 1), new T("ZZ", 1) }, new long[] { 1, 1 }, new long[] { 1, 0 }, null),
        new Setting("site 1 cooled at 2 and heated at 1: X₀ + X₀Z₁, site 0 at 1/1",
            new[] { new T("XI", 1), new T("XZ", 1) }, new long[] { 1, 2 }, new long[] { 1, 1 }, null),
        new Setting("the field on in both sectors: 2X₀ + X₀Z₁, site 0 at 1/1, site 1 cooled at 2",
            new[] { new T("XI", 2), new T("XZ", 1) }, new long[] { 1, 2 }, new long[] { 1, 0 }, null),
        new Setting("two neighbours conditioned on site 1: X₀(1 + Z₁) + X₂(1 + Z₁), the ends at 1/1, site 1 cooled at 1",
            new[] { new T("XII", 1), new T("XZI", 1), new T("IIX", 1), new T("IZX", 1) }, new long[] { 1, 1, 1 }, new long[] { 1, 0, 1 }, null),
    };

    public static readonly long[] GridValues = { 0, 1, 2, 3, 4 };

    /// <summary>The decomposition 2L + c = Σ parameter·piece reproduces the generator at a test point.</summary>
    public bool PiecesExact { get; }
    /// <summary>Entries of any piece that leave site 1's sector, besides site 1's jump (1, 1) → (0, 0).</summary>
    public int Leaks { get; }
    /// <summary>Entries of P_g carrying site 1's jump from (1, 1) to (0, 0).</summary>
    public int ArrowEntries { get; }
    public int GridPoints { get; }
    public int FieldFreeFailures { get; }
    public int OneSidedFailures { get; }
    public int FieldBlockFailures { get; }
    /// <summary>Grid points where Tr((2L + c)³) is not the parent's law.</summary>
    public int ThirdFailures { get; }
    public IReadOnlyList<PairingReading> OnLineReadings { get; }
    public IReadOnlyList<PairingReading> OffLineReadings { get; }
    public IReadOnlyList<ControlReading> ControlReadings { get; }

    public ConditionedFieldPairingWitness()
    {
        var none = Array.Empty<T>();
        var zero = new long[] { 0, 0 };
        var pieces = new (string Name, GaussianInteger[,] M)[]
        {
            ("h", ExactThermalGenerator.Build(new[] { new T("XI", 1), new T("XZ", 1) }, zero, zero)),
            ("w", ExactThermalGenerator.Build(new[] { new T("IZ", 1) }, zero, zero)),
            ("a", ExactThermalGenerator.Build(none, new long[] { 1, 0 }, zero)),
            ("b", ExactThermalGenerator.Build(none, zero, new long[] { 1, 0 })),
            ("g", ExactThermalGenerator.Build(none, new long[] { 0, 1 }, zero)),
        };

        // a test point with every parameter nonzero and distinct from the others: h = 2, w = −3, a = 5, b = 1, g = 4
        var direct = ExactThermalGenerator.Build(new[] { new T("XI", 2), new T("XZ", 2), new T("IZ", -3) }, new long[] { 5, 4 }, new long[] { 1, 0 });
        var second = ExactThermalGenerator.Build(new[] { new T("XI", -1), new T("XZ", -1), new T("IZ", 2) }, new long[] { 0, 1 }, new long[] { 3, 0 });
        PiecesExact = Equal(direct, Combine(pieces, new long[] { 2, -3, 5, 1, 4 })) && Equal(second, Combine(pieces, new long[] { -1, 2, 0, 3, 1 }));

        int leaks = 0, arrows = 0;
        foreach (var (name, m) in pieces)
            for (int r = 0; r < 16; r++)
                for (int c = 0; c < 16; c++)
                {
                    if (m[r, c] == GaussianInteger.Zero) continue;
                    var sr = ExactThermalGenerator.SiteOneSector(r);
                    var sc = ExactThermalGenerator.SiteOneSector(c);
                    if (sr == sc) continue;
                    if (name == "g" && sc == (1, 1) && sr == (0, 0)) arrows++;
                    else leaks++;
                }
        Leaks = leaks;
        ArrowEntries = arrows;

        var idx = new Dictionary<(int, int), int[]>();
        foreach (var kb in new[] { (0, 0), (0, 1), (1, 0), (1, 1) })
            idx[kb] = Enumerable.Range(0, 16).Where(i => ExactThermalGenerator.SiteOneSector(i) == kb).ToArray();
        int points = 0, ff = 0, os = 0, fb = 0, tf = 0;
        foreach (var h in GridValues)
            foreach (var w in GridValues)
                foreach (var a in GridValues)
                    foreach (var b in GridValues)
                        foreach (var g in GridValues)
                        {
                            points++;
                            var m = Combine(pieces, new[] { h, w, a, b, g });
                            var p11 = ExactThermalGenerator.CharPoly(ExactThermalGenerator.Extract(m, idx[(1, 1)]));
                            var expected = Product(new long[] { -(a + b - g), g, g, a + b + g });
                            if (!p11.SequenceEqual(expected)) ff++;
                            var p01 = ExactThermalGenerator.CharPoly(ExactThermalGenerator.Extract(m, idx[(0, 1)]));
                            var p10 = ExactThermalGenerator.CharPoly(ExactThermalGenerator.Extract(m, idx[(1, 0)]));
                            if (Enumerable.Range(0, 5).Any(k => p01[k] != (k % 2 == 0 ? p10[k] : -p10[k]))) os++;
                            var block00 = ExactThermalGenerator.Extract(m, idx[(0, 0)]);
                            var p00 = ExactThermalGenerator.CharPoly(block00);
                            if (ExactThermalGenerator.Evaluate(p00, new GaussianInteger(a + b + g, 0)) != GaussianInteger.Zero
                                || ExactThermalGenerator.Evaluate(p00, new GaussianInteger(g, 0)) != GaussianInteger.Zero
                                || ExactThermalGenerator.Trace(block00) != new GaussianInteger(4 * g, 0)) fb++;
                            var third = ExactThermalGenerator.TraceOfProduct(ExactThermalGenerator.Multiply(m, m), m);
                            var law = ThermalThirdMomentClaim.CubedCentredTrace(new[] { new T("XI", h), new T("XZ", h), new T("IZ", w) },
                                new[] { a, g }, new[] { b, 0L });
                            if (third != new GaussianInteger(law, 0)) tf++;
                        }
        GridPoints = points;
        FieldFreeFailures = ff;
        OneSidedFailures = os;
        FieldBlockFailures = fb;
        ThirdFailures = tf;

        OnLineReadings = OnLine.Select(ReadPairing).ToArray();
        OffLineReadings = OffLine.Select(ReadPairing).ToArray();
        ControlReadings = Controls.Select(s =>
        {
            var sums = ExactThermalGenerator.OddPowerSums(ExactThermalGenerator.Build(s.Terms, s.Cool, s.Heat), 5);
            return new ControlReading(s, sums[3], sums[5], ThermalThirdMomentClaim.CubedCentredTrace(s.Terms, s.Cool, s.Heat));
        }).ToArray();
    }

    private static PairingReading ReadPairing(Setting s)
    {
        var m = ExactThermalGenerator.Build(s.Terms, s.Cool, s.Heat);
        var cp = ExactThermalGenerator.CharPoly(m);
        var oddCoefficients = Enumerable.Range(0, 17).Where(k => k % 2 == 1 && cp[k] != GaussianInteger.Zero).ToArray();
        var sums = ExactThermalGenerator.OddPowerSums(m, 16);
        var oddSums = sums.Where(kv => !kv.Value.IsZero).Select(kv => kv.Key).ToArray();
        if (cp[13].Im != 0) throw new InvalidOperationException("the y¹³ coefficient has an imaginary part");
        BigInteger? law = s.XForm is { } x ? ConditionedFieldPairingClaim.ThirteenthCoefficientOfTwoLPlusC(x.H, x.A, x.B, x.G) : null;
        return new PairingReading(s, oddCoefficients, oddSums, cp[13].Re, law);
    }

    private static GaussianInteger[,] Combine((string Name, GaussianInteger[,] M)[] pieces, long[] values)
    {
        var r = new GaussianInteger[16, 16];
        for (int p = 0; p < pieces.Length; p++)
        {
            if (values[p] == 0) continue;
            var v = new GaussianInteger(values[p], 0);
            for (int i = 0; i < 16; i++)
                for (int j = 0; j < 16; j++)
                    if (pieces[p].M[i, j] != GaussianInteger.Zero) r[i, j] += v * pieces[p].M[i, j];
        }
        return r;
    }

    private static bool Equal(GaussianInteger[,] x, GaussianInteger[,] y)
    {
        for (int i = 0; i < x.GetLength(0); i++)
            for (int j = 0; j < x.GetLength(1); j++)
                if (x[i, j] != y[i, j]) return false;
        return true;
    }

    /// <summary>The coefficients of Π (y + r) over the given r, index k for y^k.</summary>
    private static GaussianInteger[] Product(long[] roots)
    {
        var c = new GaussianInteger[] { GaussianInteger.One };
        foreach (var r in roots)
        {
            var next = new GaussianInteger[c.Length + 1];
            for (int k = 0; k < c.Length; k++)
            {
                next[k + 1] += c[k];
                next[k] += new GaussianInteger(r, 0) * c[k];
            }
            c = next;
        }
        return c;
    }

    public bool TriangularEverywhere => PiecesExact && Leaks == 0 && ArrowEntries > 0;
    public bool BlocksProven => GridPoints == 3125 && FieldFreeFailures == 0 && OneSidedFailures == 0 && FieldBlockFailures == 0;
    public bool ConverseProven => GridPoints == 3125 && ThirdFailures == 0;
    public bool PairsOnTheLine => OnLineReadings.Count == OnLine.Count && OnLineReadings.All(x => x.NonzeroOddCoefficients.Count == 0 && x.NonzeroOddSums.Count == 0 && x.Thirteenth.IsZero);
    public bool ThirteenthIsTheLaw => OnLineReadings.Concat(OffLineReadings).All(x => x.Law is not BigInteger l || x.Thirteenth == l);
    public bool BreaksOffTheLine => OffLineReadings.Count >= 2 && OffLineReadings.All(x => !x.Thirteenth.IsZero && x.NonzeroOddCoefficients.Contains(13));
    public bool ControlsCancelThenBreak => ControlReadings.Count >= 4 && ControlReadings.All(x => x.CancelsThenBreaks);

    public string DisplayName => "ConditionedFieldPairingWitness (h(X₀ + X₀Z₁) + w·Z₁, site 0 at a/b, site 1 cooled at g)";

    public string Summary =>
        $"triangular on every piece {(TriangularEverywhere ? "(no leak beside site 1's jump)" : "NOT")}; the block identities and the converse on {GridPoints} grid points " +
        $"{(BlocksProven && ConverseProven ? "hold, a proof at every parameter" : "FAIL somewhere")}; on a + b = 2g the odd coefficients and odd sums " +
        $"{(PairsOnTheLine ? "vanish" : "do NOT all vanish")}, off it the y¹³ coefficient {(ThirteenthIsTheLaw && BreaksOffTheLine ? "is the law and not zero" : "is NOT the law or is zero")}; " +
        $"the controls {(ControlsCancelThenBreak ? "cancel the third and break at the fifth" : "do NOT all cancel then break")}.";

    public IEnumerable<IInspectable> Children
    {
        get
        {
            yield return new InspectableNode(
                displayName: "the triangularity",
                summary: $"2L + c = h·P_h + w·P_w + a·P_a + b·P_b + g·P_g {(PiecesExact ? "exactly at both test points" : "NOT at a test point")} (L is linear in H and the rates); entries leaving site 1's ket-bra sector: " +
                         $"{Leaks} beside site 1's jump (1, 1) → (0, 0), which has {ArrowEntries}.",
                provenance: NodeProvenance.Live);
            yield return new InspectableNode(
                displayName: "the blocks, proven on the grid",
                summary: $"{GridPoints} points of {{0, …, 4}}⁵ in (h, w, a, b, g), one wider than every degree: field-free block " +
                         $"(y − (a + b − g))(y + g)²(y + a + b + g) failed at {FieldFreeFailures}; p₀₁(y) = p₁₀(−y) failed at {OneSidedFailures}; " +
                         $"the field's block zero at a + b + g and g with trace 4g failed at {FieldBlockFailures}; Tr((2L + c)³) against the parent's " +
                         $"law, the converse, failed at {ThirdFailures}.",
                provenance: NodeProvenance.Live);
            foreach (var x in OnLineReadings)
                yield return new InspectableNode(
                    displayName: $"on the line: {x.Setting.Label}",
                    summary: $"odd coefficients nonzero at {(x.NonzeroOddCoefficients.Count == 0 ? "none" : string.Join(", ", x.NonzeroOddCoefficients))}, " +
                             $"odd power sums through 16 nonzero at {(x.NonzeroOddSums.Count == 0 ? "none" : string.Join(", ", x.NonzeroOddSums))}" +
                             (x.Law is BigInteger l ? $"; y¹³ coefficient {x.Thirteenth}, the law {l}" : "") + ".",
                    provenance: NodeProvenance.Live);
            foreach (var x in OffLineReadings)
                yield return new InspectableNode(
                    displayName: $"off the line: {x.Setting.Label}",
                    summary: $"y¹³ coefficient {x.Thirteenth}, the law −64h²(a + b − 2g) = {x.Law}: " +
                             $"{(x.Law is BigInteger l && l == x.Thirteenth && !l.IsZero ? "equal, and not zero" : "DIFFERENT or zero")}.",
                    provenance: NodeProvenance.Live);
            foreach (var x in ControlReadings)
                yield return new InspectableNode(
                    displayName: $"control: {x.Setting.Label}",
                    summary: $"Tr((2L + c)³) = {x.Third} (law {x.Law}), Tr((2L + c)⁵) = {x.Fifth}: " +
                             $"{(x.CancelsThenBreaks ? "the third cancels, the fifth breaks" : "NOT the expected cancel-then-break")}.",
                    provenance: NodeProvenance.Live);
        }
    }

    public InspectablePayload Payload => InspectablePayload.Empty;
}
