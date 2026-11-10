using System.Numerics;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Numerics;
using RCPsiSquared.Core.Symmetry;
using B = RCPsiSquared.Core.Symmetry.FluxAreaFifthMomentClaim.Bond;
using T = RCPsiSquared.Core.Symmetry.ThermalThirdMomentClaim.Term;

namespace RCPsiSquared.Diagnostics.Foundation;

/// <summary>The live lab for <see cref="TransportSeventhMomentClaim"/> (Theorem B of
/// <c>docs/proofs/PROOF_TRANSPORT_FORMULA_LOSS_AND_GAIN.md</c>; F137): the part of degree two in the rates of
/// Tr((2L + c)⁷), read from the generator in exact Gaussian integers (<see cref="ExactThermalGenerator"/>) against
/// the claim's closed form, read from the bonds and the rates alone.
///
/// <para>The road: 2L + c = A + B with A the generator at zero rates and B its part linear in the rates, and the part
/// of degree two of Tr((A + B)⁷) is the sum of the 21 words with two letters B, 7·(Tr(B·B·A⁵) + Tr(B·A·B·A⁴) +
/// Tr(B·A²·B·A³)). This is not the road of the proof's gate <c>simulations/seventh_moment_transport.py</c>,
/// where the rates are scaled and a small linear system is solved. The whole trace Tr((2L + c)⁷) is read beside it, so a
/// row can say that the trace equals its part of degree two, or that this part is zero while the trace is not; the
/// witness does not read the part of degree four by itself, and "all of degree two" and "of degree four" below lean on
/// the trace having no parts of odd degree or of degree zero or six (F137's maps and its argument at first order in H,
/// proof §6; degrees zero and six checked by the gate).</para>
///
/// <para>The rows, at N = 3 and 4: a triangle with unequal hoppings and ZZ and one with equal hoppings and one Δ, a site hanging on a triangle by a hopping
/// bond (the copy), ZZ bonds from a triangle to a fourth site, a ring of four with ZZ on its edges and diagonals, the
/// complete graph on four sites; the recorded values (gate M11's two breaks, −44800 and −20380416⁄25, the second
/// with H times 10 and the rates times 840; the ring of four of the fifth trace's gate, −68239360, and zero without
/// its ZZ; the temperature gradient, 430080 all of degree four; the triangle with unequal weights, 322560 all of
/// degree four); three settings where the part is zero (rate points on a line at equal weights, cooling alone, flux 0);
/// and three controls that leave the class and must miss the closed form. A ring of five and a ZZ bond off a
/// triangle need five sites: no row reads them against a generator here, the gate does.</para></summary>
public sealed class TransportSeventhMomentWitness : IInspectable
{
    public sealed record Row(string Label, int N, IReadOnlyList<B> Bonds, IReadOnlyList<long> Cool, IReadOnlyList<long> Heat,
        BigInteger? RecordedDegreeTwo = null, BigInteger? RecordedTrace = null, IReadOnlyList<T>? Extra = null)
    {
        /// <summary>A control: H carries a term outside the class, and the closed form of the bonds alone must miss.</summary>
        public bool IsControl => Extra is not null;
    }

    public sealed record Reading(Row Row, BigInteger DegreeTwo, BigInteger Trace, BigInteger Law)
    {
        /// <summary>The generator's part of degree two equals the closed form, and the recorded values where there are
        /// any; for a control, it is not zero and does not equal the closed form (a zero can come from a pairing,
        /// as with 2·X₀X₁ below, which is a miss for another reason).</summary>
        public bool Holds => Row.IsControl
            ? DegreeTwo != Law && !DegreeTwo.IsZero
            : DegreeTwo == Law
              && (Row.RecordedDegreeTwo is not BigInteger r || DegreeTwo == r)
              && (Row.RecordedTrace is not BigInteger f || Trace == f);

        public bool AllOfDegreeTwo => Trace == DegreeTwo;
    }

    public static IReadOnlyList<Row> Rows { get; } = BuildRows();

    private static IReadOnlyList<Row> BuildRows()
    {
        var unit = new[] { new B(0, 1, 1, 0, 0), new B(1, 2, 1, 0, 0), new B(2, 0, 0, 1, 0) };                  // loop i, flux π/2
        var zzTriangle = new[] { new B(0, 1, 1, 2, 5), new B(1, 2, 3, -1, -2), new B(0, 2, 2, 1, 1) };          // |t|² = 5, 10, 5
        var hanging = new[] { new B(0, 1, 1, 0, 0), new B(1, 2, 2, 1, 0), new B(2, 0, 0, 1, 0), new B(3, 1, 2, -3, 0) };
        var zzOut = new[] { new B(0, 1, 1, 0, 0), new B(1, 2, 2, 1, 0), new B(2, 0, 0, 1, 0), new B(0, 3, 0, 0, 2), new B(1, 3, 0, 0, -3), new B(2, 3, 0, 0, 1) };
        var ringDiagonals = new[]
        {
            new B(0, 1, 1, 0, 3), new B(1, 2, 2, 1, -1), new B(2, 3, 1, -1, 2), new B(3, 0, 0, 2, 5),
            new B(0, 2, 0, 0, 7), new B(1, 3, 0, 0, -4),
        };
        var k4 = new[]
        {
            new B(0, 1, 1, 2, 1), new B(0, 2, 2, -1, -2), new B(0, 3, 1, 1, 3),
            new B(1, 2, 3, 1, 2), new B(3, 1, -1, 2, -1), new B(2, 3, 2, 3, 4),
        };
        var uniformDelta = new[] { new B(0, 1, 1, 0, 3), new B(1, 2, 1, 0, 3), new B(2, 0, 0, 1, 3) };
        var skew = new[] { new B(0, 1, 1, 0, 0), new B(1, 2, 1, 1, 0), new B(2, 0, 0, 2, 0) };                   // weights 98, 91, 77
        var pendant = new[] { new B(0, 1, 1, 0, 0), new B(1, 2, 1, 0, 0), new B(2, 0, 0, 1, 0), new B(3, 0, 1, 0, 1) };
        var ringFifth = new[] { new B(0, 1, 1, 0, 7), new B(1, 2, 2, 0, 7), new B(2, 3, 3, 0, 7), new B(3, 0, 0, 17, 7) };
        var ringFifthBare = new[] { new B(0, 1, 1, 0, 0), new B(1, 2, 2, 0, 0), new B(2, 3, 3, 0, 0), new B(3, 0, 0, 17, 0) };
        var ringM11 = new[] { new B(0, 1, 10, 0, 7), new B(1, 2, 20, 0, 7), new B(2, 3, 30, 0, 7), new B(3, 0, 0, 17, 7) };   // H times 10
        long[] c3 = { 2, 1, 3 }, h3 = { 1, 3, 1 }, c4 = { 3, 1, 2, 5 }, h4 = { 1, 2, 2, 1 };
        BigInteger scale72 = 128 * 72 * 72;
        BigInteger scaleM11 = BigInteger.Pow(10, 5) * 840 * 840 * 128;
        return new[]
        {
            new Row("a triangle with |t|² = 5, 10, 5 and ZZ 5, −2, 1", 3, zzTriangle, c3, h3),
            new Row("a triangle with a site hanging on one corner by a hopping bond: the copy", 4, hanging, c4, h4),
            new Row("a triangle with ZZ bonds to a fourth site and no hopping there", 4, zzOut, c4, h4),
            new Row("a ring of four with flux, ZZ on its edges and both diagonals", 4, ringDiagonals, new long[] { 2, 1, 3, 1 }, new long[] { 1, 3, 1, 2 }),
            new Row("the complete graph on four sites, complex hoppings, six different ZZ", 4, k4, new long[] { 2, 1, 3, 1 }, new long[] { 1, 3, 1, 2 }),
            new Row("the unit triangle at flux π/2 with one Δ on every bond, three rate points off a line", 3, uniformDelta, c3, h3),
            new Row("gate M11's triangle at ratio ¼ with a fourth site at ratio 1⁄9, rates times 72 (recorded −44800, all of degree two)", 4, pendant,
                new long[] { 36, 72, 108, 72 }, new long[] { 9, 18, 27, 8 }, -44800 * scale72, -44800 * scale72),
            new Row("gate M11's ring of four at Δ = 7⁄10 and flux π/2, H times 10, rates times 840 (recorded −20380416⁄25, all of degree two)", 4, ringM11,
                new long[] { 420, 840, 1260, 1680 }, new long[] { 105, 280, 168, 120 }, -20380416 * scaleM11 / 25, -20380416 * scaleM11 / 25),
            new Row("the ring of four of the fifth trace's gate, hoppings 1, 2, 3, 17i with ZZ = 7, rates times 72 (recorded −68239360, all of degree two)", 4, ringFifth,
                new long[] { 36, 72, 108, 72 }, new long[] { 9, 24, 20, 8 }, -68239360 * scale72, -68239360 * scale72),
            new Row("the same ring without its ZZ: the whole seventh trace is zero", 4, ringFifthBare,
                new long[] { 36, 72, 108, 72 }, new long[] { 9, 24, 20, 8 }, BigInteger.Zero, BigInteger.Zero),
            new Row("a temperature gradient on the unit triangle, γ↓ = 2(n̄ + 1), γ↑ = 2n̄, n̄ = (1, 2, 3): degree two zero, the trace 430080 of degree four", 3, unit,
                new long[] { 4, 6, 8 }, new long[] { 2, 4, 6 }, BigInteger.Zero, 430080 * 128),
            new Row("hoppings 1, 1 + i, 2i (weights 98, 91, 77), γ↓ = (3, 5, 2), γ↑ = γ↓ − 2: degree two zero, the trace 322560 of degree four", 3, skew,
                new long[] { 3, 5, 2 }, new long[] { 1, 3, 0 }, BigInteger.Zero, 322560 * 128),
            new Row("the unit triangle with one Δ on every bond, rate points on the line γ↑ = γ↓ − 1: degree two zero", 3, uniformDelta,
                new long[] { 2, 5, 3 }, new long[] { 1, 4, 2 }, BigInteger.Zero),
            new Row("cooling alone on the triangle with ZZ: degree two zero", 3, zzTriangle, c3, new long[] { 0, 0, 0 }, BigInteger.Zero),
            new Row("flux 0 (real hoppings 1, 2, 3) with ZZ, three ratios: degree two zero", 3,
                new[] { new B(0, 1, 1, 0, 1), new B(1, 2, 2, 0, 0), new B(2, 0, 3, 0, 2) }, c3, h3, BigInteger.Zero),
            // 3·X₀X₁ and not 2: the part of 2·X₀X₁ that conserves the number moves t₀₁ to 2 + 2i, the loop becomes real (20), and
            // every odd trace is zero by a pairing, not by leaving the class
            new Row("control: the triangle with ZZ and 3·X₀X₁ (the excitation number not conserved)", 3, zzTriangle, c3, h3, Extra: new[] { new T("XXI", 3) }),
            new Row("control: the triangle with ZZ and 2·Y₁Z₂", 3, zzTriangle, c3, h3, Extra: new[] { new T("IYZ", 2) }),
            new Row("control: the triangle with ZZ and a Z field 3 on site 0", 3, zzTriangle, c3, h3, Extra: new[] { new T("ZII", 3) }),
        };
    }

    public IReadOnlyList<Reading> Readings { get; }

    public TransportSeventhMomentWitness()
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

        // x[p] = B·A^p; the words with two letters B: 7·(Tr(x0·x5) + Tr(x1·x4) + Tr(x2·x3))
        var x = new GaussianInteger[6][,];
        x[0] = b;
        for (int p = 1; p < 6; p++) x[p] = ExactThermalGenerator.Multiply(x[p - 1], a);
        var two = ExactThermalGenerator.TraceOfProduct(x[0], x[5]) + ExactThermalGenerator.TraceOfProduct(x[1], x[4])
                  + ExactThermalGenerator.TraceOfProduct(x[2], x[3]);
        if (!two.Im.IsZero) throw new InvalidOperationException($"{row.Label}: the part of degree two has an imaginary part");

        var m2 = ExactThermalGenerator.Multiply(m, m);
        var m3 = ExactThermalGenerator.Multiply(m2, m);
        var m4 = ExactThermalGenerator.Multiply(m3, m);
        var seventh = ExactThermalGenerator.TraceOfProduct(m3, m4);
        if (!seventh.Im.IsZero) throw new InvalidOperationException($"{row.Label}: Tr(M⁷) has an imaginary part");

        return new Reading(row, 7 * two.Re, seventh.Re,
            TransportSeventhMomentClaim.SeventhDegreeTwo(row.N, row.Bonds, row.Cool, row.Heat));
    }

    public bool AllHold => Readings.All(x => x.Holds);
    public int Speaking => Readings.Count(x => !x.Row.IsControl && !x.Law.IsZero);
    public int Silent => Readings.Count(x => !x.Row.IsControl && x.Law.IsZero);
    public int Controls => Readings.Count(x => x.Row.IsControl);
    public int Recorded => Readings.Count(x => x.Row.RecordedTrace is not null);
    public bool NotVacuous => Speaking >= 8 && Silent >= 5 && Controls >= 3 && Recorded >= 5;

    public string DisplayName => "TransportSeventhMomentWitness (the part of degree two of Tr((2L + c)⁷) against 3584·4^N·Σ s_l·a_k·3β_lk)";

    public string Summary =>
        $"{Readings.Count(x => x.Holds)} of {Readings.Count} rows: the generator's part of degree two of Tr((2L + c)⁷) " +
        $"{(AllHold ? "equals" : "does NOT always equal")} the closed form on {Speaking} speaking rows (the recorded −44800, −20380416⁄25 and " +
        $"−68239360 among them, each all of degree two) and on {Silent} silent ones (the 430080 and the 322560 of degree four among them), and misses " +
        $"it on {Controls} controls outside the class {(NotVacuous ? "(not vacuous)" : "(VACUOUS)")}.";

    public IEnumerable<IInspectable> Children
    {
        get
        {
            foreach (var x in Readings)
                yield return new InspectableNode(
                    displayName: x.Row.Label,
                    summary: $"degree two {x.DegreeTwo}, closed form {x.Law}, the whole trace {x.Trace}" +
                             $"{(x.AllOfDegreeTwo && !x.Trace.IsZero ? " (all of degree two)" : "")}: " +
                             (x.Row.IsControl
                                 ? (x.Holds ? "different, as a control must be" : "EQUAL, the control does not discriminate")
                                 : (x.Holds ? "equal" : "DIFFERENT")) + ".",
                    provenance: NodeProvenance.Live);
        }
    }

    public InspectablePayload Payload => InspectablePayload.Empty;
}
