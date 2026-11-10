using System.Numerics;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Numerics;
using RCPsiSquared.Core.Symmetry;
using B = RCPsiSquared.Core.Symmetry.FluxAreaFifthMomentClaim.Bond;
using T = RCPsiSquared.Core.Symmetry.ThermalThirdMomentClaim.Term;

namespace RCPsiSquared.Diagnostics.Foundation;

/// <summary>The live lab for <see cref="SeventhMomentDegreeFourClaim"/> (Theorems 1 and 2 of
/// <c>docs/proofs/PROOF_SEVENTH_TRACE_DEGREE_FOUR_TRIANGLE.md</c>; F137): the part of degree four in the rates of
/// Tr((2L + c)⁷), read from the generator in exact Gaussian integers (<see cref="ExactThermalGenerator"/>) against the
/// claim's sum over the triangles, and the whole trace against the two closed forms together.
///
/// <para>The road: 2L + c = A + B with A the generator at zero rates and B its part linear in the rates. The part of
/// degree four of Tr((A + B)⁷) is the sum of the 35 words with four letters B and three letters A; seven being prime,
/// each word has seven distinct rotations, and the classes are the five necklaces A·B^g₁·A·B^g₂·A·B^g₃ with the gaps
/// (4, 0, 0), (3, 1, 0), (3, 0, 1), (2, 2, 0), (2, 1, 1), so the part is 7 times the sum of their traces. The proof's
/// gate scales the rates and solves a small linear system instead, and its fit and walk sum are other roads again.</para>
///
/// <para>The rows, at N = 3 and 4: the two recorded values (F137's temperature gradient, 430080, and the transport
/// proof's triangle with unequal weights, 322560, both all of degree four), a triangle with DM and ZZ, a site hanging
/// on a triangle by a hopping bond with ZZ to another corner, two triangles sharing an edge, the complete graph on four
/// sites (four triangles); three settings where the part is zero (a ring of four with no triangle, one ratio on every
/// site, flux 0); and two controls that leave the class, whose whole trace must miss the closed forms together. Read
/// there: 3·X₀X₁ moves the part of degree four off the sum over the triangles, while a Z field 3 on one site leaves it
/// on it at that row and moves the trace through the other degrees (whether a Z field always leaves it, and why, is not
/// shown).</para></summary>
public sealed class SeventhMomentDegreeFourWitness : IInspectable
{
    public sealed record Row(string Label, int N, IReadOnlyList<B> Bonds, IReadOnlyList<long> Cool, IReadOnlyList<long> Heat,
        BigInteger? RecordedTrace = null, IReadOnlyList<T>? Extra = null)
    {
        /// <summary>A control: H carries a term outside the class, and the closed form of the bonds alone must miss.</summary>
        public bool IsControl => Extra is not null;
    }

    public sealed record Reading(Row Row, BigInteger DegreeFour, BigInteger Trace, BigInteger Law, BigInteger WholeLaw)
    {
        /// <summary>The generator's part of degree four equals the sum over the triangles, the whole trace equals the
        /// closed forms of degree two and four together, and the recorded value where there is one; a control's whole
        /// trace is not zero and misses the closed forms together (its part of degree four may or may not miss: the Z
        /// field below leaves it at the sum over the triangles and moves the trace through the other degrees).</summary>
        public bool Holds => Row.IsControl
            ? Trace != WholeLaw && !Trace.IsZero
            : DegreeFour == Law && Trace == WholeLaw && (Row.RecordedTrace is not BigInteger r || Trace == r);

        public bool DegreeFourMisses => DegreeFour != Law;
    }

    public static IReadOnlyList<Row> Rows { get; } = BuildRows();

    private static IReadOnlyList<Row> BuildRows()
    {
        var unit = new[] { new B(0, 1, 1, 0, 0), new B(1, 2, 1, 0, 0), new B(2, 0, 0, 1, 0) };
        var skew = new[] { new B(0, 1, 1, 0, 0), new B(1, 2, 1, 1, 0), new B(2, 0, 0, 2, 0) };
        var dmZZ = new[] { new B(0, 1, 2, 1, 3), new B(1, 2, 1, -2, -1), new B(2, 0, 3, 1, 2) };
        var hanging = new[] { new B(0, 1, 1, 1, 0), new B(1, 2, 2, 0, 1), new B(2, 0, 0, 1, 0), new B(3, 0, 1, -1, 0), new B(3, 2, 0, 0, 2) };
        var shared = new[] { new B(0, 1, 1, 0, 2), new B(1, 2, 2, 1, 0), new B(2, 0, 0, 1, -1), new B(1, 3, 1, 2, 1), new B(3, 2, 2, -1, 0) };
        var k4 = new[]
        {
            new B(0, 1, 1, 2, 1), new B(0, 2, 2, -1, -2), new B(0, 3, 1, 1, 3),
            new B(1, 2, 3, 1, 2), new B(3, 1, -1, 2, -1), new B(2, 3, 2, 3, 4),
        };
        var ring = new[]
        {
            new B(0, 1, 1, 0, 3), new B(1, 2, 2, 1, -1), new B(2, 3, 1, -1, 2), new B(3, 0, 0, 2, 5),
            new B(0, 2, 0, 0, 7), new B(1, 3, 0, 0, -4),
        };
        var real = new[] { new B(0, 1, 1, 0, 1), new B(1, 2, 2, 0, 0), new B(2, 0, 3, 0, 2) };
        long[] c4 = { 2, 1, 3, 1 }, h4 = { 1, 3, 1, 2 };
        return new[]
        {
            new Row("F137's temperature gradient on the unit triangle, γ↓ = 2(n̄ + 1), γ↑ = 2n̄, n̄ = (1, 2, 3): 430080, all of degree four", 3, unit,
                new long[] { 4, 6, 8 }, new long[] { 2, 4, 6 }, 430080 * (BigInteger)128),
            new Row("hoppings 1, 1 + i, 2i at γ↓ = (3, 5, 2), γ↑ = γ↓ − 2: 322560, all of degree four", 3, skew,
                new long[] { 3, 5, 2 }, new long[] { 1, 3, 0 }, 322560 * (BigInteger)128),
            new Row("a triangle with DM on every bond and ZZ 3, −1, 2", 3, dmZZ, new long[] { 3, 1, 4 }, new long[] { 1, 5, 2 }),
            new Row("a triangle with a site hanging on corner 0 by a hopping bond and ZZ to corner 2", 4, hanging, c4, h4),
            new Row("two triangles sharing the edge 1–2, ZZ on three bonds", 4, shared, new long[] { 3, 1, 2, 1 }, new long[] { 1, 2, 4, 3 }),
            new Row("the complete graph on four sites: four triangles, complex hoppings, six ZZ", 4, k4, c4, h4),
            new Row("a ring of four with flux and ZZ on its edges and diagonals, no triangle: degree four zero", 4, ring, c4, h4),
            new Row("one ratio on every site of the complete graph: every ω zero, the whole trace zero", 4, k4,
                new long[] { 2, 4, 6, 8 }, new long[] { 1, 2, 3, 4 }, BigInteger.Zero),
            new Row("flux 0 (real hoppings 1, 2, 3) with ZZ, three ratios: degree four zero", 3, real, new long[] { 2, 1, 3 }, new long[] { 1, 3, 1 }),
            new Row("control: the triangle with DM and ZZ and a Z field 3 on site 0", 3, dmZZ, new long[] { 3, 1, 4 }, new long[] { 1, 5, 2 },
                Extra: new[] { new T("ZII", 3) }),
            new Row("control: the triangle with DM and ZZ and 3·X₀X₁ (the excitation number not conserved)", 3, dmZZ, new long[] { 3, 1, 4 }, new long[] { 1, 5, 2 },
                Extra: new[] { new T("XXI", 3) }),
        };
    }

    /// <summary>The five necklaces of three A and four B, as the gaps of B between the three A.</summary>
    public static IReadOnlyList<(int, int, int)> Necklaces { get; } = new[] { (4, 0, 0), (3, 1, 0), (3, 0, 1), (2, 2, 0), (2, 1, 1) };

    public IReadOnlyList<Reading> Readings { get; }

    public SeventhMomentDegreeFourWitness()
    {
        Readings = Rows.Select(Read).ToArray();
    }

    public static Reading Read(Row row)
    {
        var terms = FluxAreaFifthMomentClaim.Terms(row.N, row.Bonds).Concat(row.Extra ?? Array.Empty<T>()).ToArray();
        var m = ExactThermalGenerator.Build(terms, row.Cool, row.Heat);
        var a = ExactThermalGenerator.Build(terms, new long[row.N], new long[row.N]);
        int dim = m.GetLength(0);
        var b = new GaussianInteger[dim, dim];
        for (int i = 0; i < dim; i++)
            for (int j = 0; j < dim; j++) b[i, j] = m[i, j] - a[i, j];

        var pow = new GaussianInteger[5][,];
        pow[1] = b;
        for (int g = 2; g <= 4; g++) pow[g] = ExactThermalGenerator.Multiply(pow[g - 1], b);
        GaussianInteger[,] Times(GaussianInteger[,] x, int g) => g == 0 ? x : ExactThermalGenerator.Multiply(x, pow[g]);

        var four = GaussianInteger.Zero;
        foreach (var (g1, g2, g3) in Necklaces)
        {
            var x = Times(a, g1);
            x = Times(ExactThermalGenerator.Multiply(x, a), g2);
            x = ExactThermalGenerator.Multiply(x, a);
            four += g3 == 0 ? ExactThermalGenerator.Trace(x) : ExactThermalGenerator.TraceOfProduct(x, pow[g3]);
        }
        if (!four.Im.IsZero) throw new InvalidOperationException($"{row.Label}: the part of degree four has an imaginary part");

        var m2 = ExactThermalGenerator.Multiply(m, m);
        var m3 = ExactThermalGenerator.Multiply(m2, m);
        var m4 = ExactThermalGenerator.Multiply(m3, m);
        var seventh = ExactThermalGenerator.TraceOfProduct(m3, m4);
        if (!seventh.Im.IsZero) throw new InvalidOperationException($"{row.Label}: Tr(M⁷) has an imaginary part");

        return new Reading(row, 7 * four.Re, seventh.Re,
            SeventhMomentDegreeFourClaim.SeventhDegreeFour(row.N, row.Bonds, row.Cool, row.Heat),
            SeventhMomentDegreeFourClaim.SeventhCentredTrace(row.N, row.Bonds, row.Cool, row.Heat));
    }

    public bool AllHold => Readings.All(x => x.Holds);
    public int Speaking => Readings.Count(x => !x.Row.IsControl && !x.Law.IsZero);
    public int Silent => Readings.Count(x => !x.Row.IsControl && x.Law.IsZero);
    public int Controls => Readings.Count(x => x.Row.IsControl);
    public int Recorded => Readings.Count(x => x.Row.RecordedTrace is not null);
    public bool NotVacuous => Speaking >= 6 && Silent >= 3 && Controls >= 2 && Recorded >= 2;

    public string DisplayName => "SeventhMomentDegreeFourWitness (the part of degree four of Tr((2L + c)⁷) against the sum over the triangles)";

    public string Summary =>
        $"{Readings.Count(x => x.Holds)} of {Readings.Count} rows: the generator's part of degree four of Tr((2L + c)⁷), as the words with four rate " +
        $"letters, {(AllHold ? "equals" : "does NOT always equal")} the sum over the triangles on {Speaking} speaking rows (the 430080 and the 322560 " +
        $"among them, graphs with one, two and four triangles) and on {Silent} silent ones, the whole trace equals the closed forms of degree two and " +
        $"four together, and the whole trace misses them on {Controls} controls outside the class {(NotVacuous ? "(not vacuous)" : "(VACUOUS)")}.";

    public IEnumerable<IInspectable> Children
    {
        get
        {
            foreach (var x in Readings)
                yield return new InspectableNode(
                    displayName: x.Row.Label,
                    summary: $"degree four {x.DegreeFour}, sum over the triangles {x.Law}, the whole trace {x.Trace}, the closed forms together {x.WholeLaw}: " +
                             (x.Row.IsControl
                                 ? (x.Holds ? $"the whole trace different, as a control must be; degree four {(x.DegreeFourMisses ? "off" : "on")} the sum over the triangles" : "EQUAL, the control does not discriminate")
                                 : (x.Holds ? "equal" : "DIFFERENT")) + ".",
                    provenance: NodeProvenance.Live);
        }
    }

    public InspectablePayload Payload => InspectablePayload.Empty;
}
