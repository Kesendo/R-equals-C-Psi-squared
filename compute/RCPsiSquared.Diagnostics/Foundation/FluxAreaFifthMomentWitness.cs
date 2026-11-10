using System.Numerics;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Symmetry;
using B = RCPsiSquared.Core.Symmetry.FluxAreaFifthMomentClaim.Bond;
using T = RCPsiSquared.Core.Symmetry.ThermalThirdMomentClaim.Term;

namespace RCPsiSquared.Diagnostics.Foundation;

/// <summary>The live lab for <see cref="FluxAreaFifthMomentClaim"/> (F137, the fifth centred trace under cooling and
/// heating as flux times the area of the rate triangle): Tr((2L + c)⁵) from the generator in exact Gaussian integers
/// (<see cref="ExactThermalGenerator"/>, matrix units, column-major) against the claim's closed form
/// 1920·4^N·Σ_triangles Im(loop)·det, read from the bonds and the rates alone.
///
/// <para>The rows: the registry's two recorded values with the rates times 120 (F137's triangle at flux π/2, −656,
/// and with loss and gain exchanged, +656; gate M11's control with hoppings 1, 2, 3i, −3936, the trace being of
/// degree two in the rates); graphs with several triangles sharing edges (the complete graph on four sites, two
/// plaquettes) and a triangle with a pendant site carrying its own rate point, which does not enter; and the silent
/// settings (flux 0, one ratio, cooling alone, rate points on a line that misses the origin, a temperature gradient
/// at one coupling, a ring of four, a triangle closed by a ZZ bond without hopping, and M11's pendant setting, the
/// unit triangle at one ratio with a pendant site at another, silent because the triangle's own rate points lie on
/// one line). Three controls leave the class and must miss the closed form: a longitudinal field on one site, a
/// uniform one, and an XX term, which does not conserve the excitation number.</para></summary>
public sealed class FluxAreaFifthMomentWitness : IInspectable
{
    public sealed record Row(string Label, int N, IReadOnlyList<B> Bonds, IReadOnlyList<long> Cool, IReadOnlyList<long> Heat,
        BigInteger? Recorded = null, IReadOnlyList<T>? Extra = null)
    {
        /// <summary>A control: H carries a term outside the class, and the closed form of the bonds alone must miss.</summary>
        public bool IsControl => Extra is not null;
    }

    public sealed record Reading(Row Row, BigInteger Generator, BigInteger Law, IReadOnlyList<FluxAreaFifthMomentClaim.Triangle> Triangles)
    {
        /// <summary>The generator's fifth trace equals the closed form (and the recorded value, where there is one);
        /// for a control, it does not.</summary>
        public bool Holds => Row.IsControl
            ? Generator != Law
            : Generator == Law && (Row.Recorded is not BigInteger r || Generator == r);
    }

    // rates times 120: γ↓ = (½, 1, 3⁄2), γ↑ = (⅛, ⅓, ⅕); Tr((2L + c)⁵) = 32·120²·Tr((L + c/2)⁵) at the rates as recorded
    private static readonly long[] Down120 = { 60, 120, 180 };
    private static readonly long[] Up120 = { 15, 40, 24 };
    private static readonly BigInteger Scale120 = 32 * 120 * 120;

    public static IReadOnlyList<Row> Rows { get; } = BuildRows();

    private static IReadOnlyList<Row> BuildRows()
    {
        var unit = new[] { new B(0, 1, 1, 0, 0), new B(1, 2, 1, 0, 0), new B(2, 0, 0, 1, 0) };          // loop i, flux π/2
        var m11 = new[] { new B(0, 1, 1, 0, 7), new B(1, 2, 2, 0, 7), new B(2, 0, 0, 3, 7) };            // "the 1, 2, 3i triangle": loop 6i, ZZ on every bond
        var k4 = new[]
        {
            new B(0, 1, 1, 2, 1), new B(0, 2, 2, -1, -2), new B(0, 3, 1, 1, 0),
            new B(1, 2, 3, 1, 2), new B(3, 1, -1, 2, 1), new B(2, 3, 2, 3, -1),
        };
        var plaquettes = new[] { new B(0, 1, 1, 0, 1), new B(1, 2, 0, 1, 0), new B(0, 2, 1, 1, 2), new B(2, 3, 1, 0, -1), new B(0, 3, 2, 1, 0) };
        var pendant = new[] { new B(0, 1, 1, 0, 0), new B(1, 2, 1, 0, 0), new B(2, 0, 0, 1, 0), new B(3, 0, 1, 0, 1) };
        var ring = new[] { new B(0, 1, 1, 0, 1), new B(1, 2, 2, 0, 1), new B(2, 3, 3, 0, 1), new B(3, 0, 0, 2, 1) };
        return new[]
        {
            new Row("F137's triangle at flux π/2, rates times 120 (recorded −656)", 3, unit, Down120, Up120, -656 * Scale120),
            new Row("the same with loss and gain exchanged (recorded +656)", 3, unit, Up120, Down120, 656 * Scale120),
            new Row("gate M11's control, hoppings 1, 2, 3i with ZZ = 7 on every bond, rates times 120 (recorded −3936)", 3, m11, Down120, Up120, -3936 * Scale120),
            new Row("the complete graph on four sites, complex hoppings, unequal ZZ: four triangles", 4, k4,
                new long[] { 2, 1, 3, 1 }, new long[] { 1, 3, 1, 2 }),
            new Row("two plaquettes sharing the edge (0, 2)", 4, plaquettes, new long[] { 3, 1, 2, 5 }, new long[] { 1, 2, 2, 1 }),
            new Row("a triangle with a pendant site, four different ratios: the pendant does not enter", 4, pendant,
                new long[] { 2, 1, 3, 7 }, new long[] { 1, 3, 1, 2 }),
            new Row("flux 0 (real hoppings 1, 2, 3), three ratios: silent", 3,
                new[] { new B(0, 1, 1, 0, 1), new B(1, 2, 2, 0, 0), new B(2, 0, 3, 0, 2) }, new long[] { 2, 1, 3 }, new long[] { 1, 3, 1 }),
            new Row("the 1, 2, 3i triangle at one ratio γ↑ = γ↓/2 on every site: silent", 3, m11, new long[] { 2, 4, 6 }, new long[] { 1, 2, 3 }),
            new Row("cooling alone: silent", 3, m11, new long[] { 2, 1, 3 }, new long[] { 0, 0, 0 }),
            new Row("rate points on the line γ↑ = γ↓/3 + 20, three ratios: silent", 3, m11, Down120, new long[] { 40, 60, 80 }),
            new Row("a temperature gradient at one coupling, γ↓ = 2(n̄ + 1), γ↑ = 2n̄, n̄ = (1, 2, 3): silent", 3, unit,
                new long[] { 4, 6, 8 }, new long[] { 2, 4, 6 }),
            new Row("a ring of four with flux π/2 and ZZ, four ratios: silent", 4, ring, new long[] { 2, 1, 3, 1 }, new long[] { 1, 3, 1, 2 }),
            new Row("a triangle whose third edge carries ZZ alone: silent", 3,
                new[] { new B(0, 1, 1, 2, 1), new B(1, 2, 2, -1, 0), new B(2, 0, 0, 0, 5) }, new long[] { 2, 1, 3 }, new long[] { 1, 3, 1 }),
            new Row("M11's pendant setting, the unit triangle at ratio ¼ with a pendant site at ratio 1⁄9 joined by hopping and ZZ, rates times 72: silent", 4, pendant,
                new long[] { 36, 72, 108, 72 }, new long[] { 9, 18, 27, 8 }),
            new Row("control: the 1, 2, 3i triangle with a Z field 3 on site 0", 3, m11, Down120, Up120, Extra: new[] { new T("ZII", 3) }),
            new Row("control: the 1, 2, 3i triangle with a uniform Z field 1", 3, m11, Down120, Up120,
                Extra: new[] { new T("ZII", 1), new T("IZI", 1), new T("IIZ", 1) }),
            new Row("control: the 1, 2, 3i triangle with 2·X₀X₁ (the excitation number not conserved)", 3, m11, Down120, Up120, Extra: new[] { new T("XXI", 2) }),
        };
    }

    public IReadOnlyList<Reading> Readings { get; }

    public FluxAreaFifthMomentWitness()
    {
        Readings = Rows.Select(Read).ToArray();
    }

    public static Reading Read(Row row)
    {
        var terms = FluxAreaFifthMomentClaim.Terms(row.N, row.Bonds).Concat(row.Extra ?? Array.Empty<T>()).ToArray();
        var m = ExactThermalGenerator.Build(terms, row.Cool, row.Heat);
        var m2 = ExactThermalGenerator.Multiply(m, m);
        var m3 = ExactThermalGenerator.Multiply(m, m2);
        var fifth = ExactThermalGenerator.TraceOfProduct(m2, m3);
        if (!fifth.Im.IsZero) throw new InvalidOperationException($"{row.Label}: Tr(M⁵) has an imaginary part");
        return new Reading(row, fifth.Re,
            FluxAreaFifthMomentClaim.FifthCentredTrace(row.N, row.Bonds, row.Cool, row.Heat),
            FluxAreaFifthMomentClaim.Triangles(row.N, row.Bonds, row.Cool, row.Heat));
    }

    public bool AllHold => Readings.All(x => x.Holds);
    public int Speaking => Readings.Count(x => !x.Row.IsControl && !x.Law.IsZero);
    public int SeveralTriangles => Readings.Count(x => !x.Row.IsControl && !x.Law.IsZero && x.Triangles.Count(t => !t.Product.IsZero) >= 2);
    public int Silent => Readings.Count(x => !x.Row.IsControl && x.Law.IsZero);
    public int Controls => Readings.Count(x => x.Row.IsControl);
    public bool NotVacuous => Speaking >= 5 && SeveralTriangles >= 2 && Silent >= 5 && Controls >= 2;

    public string DisplayName => "FluxAreaFifthMomentWitness (Tr((2L + c)⁵) against 1920·4^N·Σ Im(loop)·det)";

    public string Summary =>
        $"{Readings.Count(x => x.Holds)} of {Readings.Count} rows: the generator's Tr((2L + c)⁵) {(AllHold ? "equals" : "does NOT always equal")} " +
        $"the closed form on {Speaking} speaking rows ({SeveralTriangles} with several triangles; the recorded −656, +656 and −3936 among them) and " +
        $"on {Silent} silent ones, and misses it on {Controls} controls outside the class " +
        $"{(NotVacuous ? "(not vacuous)" : "(VACUOUS)")}.";

    public IEnumerable<IInspectable> Children
    {
        get
        {
            foreach (var x in Readings)
            {
                string tris = x.Triangles.Count == 0
                    ? "no triangle"
                    : string.Join("; ", x.Triangles.Select(t => $"({t.U}, {t.V}, {t.W}): Im(loop) {t.ImLoop}, det {t.RateDeterminant}"));
                yield return new InspectableNode(
                    displayName: x.Row.Label,
                    summary: $"generator {x.Generator}, closed form {x.Law} [{tris}]: " +
                             (x.Row.IsControl
                                 ? (x.Holds ? "different, as a control must be" : "EQUAL, the control does not discriminate")
                                 : (x.Holds ? "equal" : "DIFFERENT")) + ".",
                    provenance: NodeProvenance.Live);
            }
        }
    }

    public InspectablePayload Payload => InspectablePayload.Empty;
}
