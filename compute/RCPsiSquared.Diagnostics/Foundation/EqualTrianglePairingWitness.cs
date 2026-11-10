using System.Numerics;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Numerics;
using RCPsiSquared.Core.Symmetry;
using B = RCPsiSquared.Core.Symmetry.FluxAreaFifthMomentClaim.Bond;
using V = RCPsiSquared.Core.Symmetry.EqualTrianglePairingClaim.Verdict;

namespace RCPsiSquared.Diagnostics.Foundation;

/// <summary>The live lab for <see cref="EqualTrianglePairingClaim"/> (Theorem 3 of
/// <c>docs/proofs/PROOF_SEVENTH_TRACE_DEGREE_FOUR_TRIANGLE.md</c>; F137's open question on flux reversal, on this
/// triangle): the verdict the claim reads from the rates alone, against the generator in exact Gaussian integers
/// (<see cref="ExactThermalGenerator"/>).
///
/// <para>Two roads per row. The odd centred traces of M = 2L + c up to the ninth, each against its closed form where
/// the verdict names it (the fifth, <see cref="FluxAreaFifthMomentClaim"/>; the seventh, the parent
/// <see cref="SeventhMomentDegreeFourClaim"/> with the part of degree two; the ninth on equal heights, the claim's own);
/// and the characteristic polynomial det(y − M) by Samuelson and Berkowitz (<see cref="GaussianMatrixCharpoly"/>,
/// without division), whose coefficients of odd degree (M has an even dimension, 64) are all zero exactly when the
/// spectrum of M is symmetric about zero, that is when L's is paired about −c/2. The second road passes through no
/// trace, no closed form and no carrier, and it reads both halves of the theorem: even on the carrier rows, not even on
/// the break rows.</para>
///
/// <para>The rows: on the triangle with every hop 2 + i and ZZ 1 (Im = 11), every case of the theorem, one ratio,
/// cooling alone, heating alone, an undamped site beside one ratio, two equal rate points, rate points off a line, on
/// a line of each kind and of negative slope, and on equal heights; on the unit triangle (t = 1, 1, i, no ZZ, Im = 1)
/// equal heights, two equal rate points and points off a line.</para></summary>
public sealed class EqualTrianglePairingWitness : IInspectable
{
    public sealed record Row(string Label, IReadOnlyList<B> Bonds, IReadOnlyList<long> Cool, IReadOnlyList<long> Heat);

    public sealed record Reading(Row Row, V Verdict, BigInteger Fifth, BigInteger Seventh, BigInteger Ninth, bool Paired, BigInteger? Law)
    {
        /// <summary>A carrier: the three odd traces zero and the characteristic polynomial even. A break: the odd traces
        /// below the named one zero, the named one not zero and equal to its closed form, the characteristic polynomial
        /// not even.</summary>
        public bool Holds => Verdict switch
        {
            V.PairedOneRatio or V.PairedSiteReflection => Fifth.IsZero && Seventh.IsZero && Ninth.IsZero && Paired,
            V.BrokenAtFifth => !Fifth.IsZero && Fifth == Law && !Paired,
            V.BrokenAtSeventh => Fifth.IsZero && !Seventh.IsZero && Seventh == Law && !Paired,
            V.BrokenAtNinth => Fifth.IsZero && Seventh.IsZero && !Ninth.IsZero && Ninth == Law && !Paired,
            _ => false,
        };
    }

    public static IReadOnlyList<Row> Rows { get; } = BuildRows();

    private static IReadOnlyList<Row> BuildRows()
    {
        var e = new[] { new B(0, 1, 2, 1, 1), new B(1, 2, 2, 1, 1), new B(2, 0, 2, 1, 1) };     // t = 2 + i, Δ = 1, loop (2 + i)³, Im 11
        var unit = new[] { new B(0, 1, 1, 0, 0), new B(1, 2, 1, 0, 0), new B(2, 0, 0, 1, 0) };   // t = 1, 1, i, Im 1
        Row R(string label, B[] bonds, long[] cool, long[] heat) => new(label, bonds, cool, heat);
        return new[]
        {
            R("t = 2 + i, Δ = 1: one ratio, γ↑ = γ↓/2 on every site", e, new long[] { 2, 4, 6 }, new long[] { 1, 2, 3 }),
            R("t = 2 + i, Δ = 1: cooling alone", e, new long[] { 1, 2, 4 }, new long[] { 0, 0, 0 }),
            R("t = 2 + i, Δ = 1: heating alone", e, new long[] { 0, 0, 0 }, new long[] { 3, 1, 2 }),
            R("t = 2 + i, Δ = 1: site 0 undamped, the others at one ratio", e, new long[] { 0, 2, 4 }, new long[] { 0, 1, 2 }),
            R("t = 2 + i, Δ = 1: sites 0 and 1 at one rate point", e, new long[] { 3, 3, 1 }, new long[] { 1, 1, 4 }),
            R("t = 2 + i, Δ = 1: three rate points off a line", e, new long[] { 2, 1, 3 }, new long[] { 1, 3, 2 }),
            R("t = 2 + i, Δ = 1: on the line γ↑ = 1 + γ↓", e, new long[] { 1, 2, 4 }, new long[] { 2, 3, 5 }),
            R("t = 2 + i, Δ = 1: on the line γ↓ = 2", e, new long[] { 2, 2, 2 }, new long[] { 1, 3, 4 }),
            R("t = 2 + i, Δ = 1: on the line γ↑ = 7 − 2γ↓", e, new long[] { 1, 2, 3 }, new long[] { 5, 3, 1 }),
            R("t = 2 + i, Δ = 1: equal heights γ↓ + γ↑ = 6, γ↓ − γ↑ = 4, −2, 2", e, new long[] { 5, 2, 4 }, new long[] { 1, 4, 2 }),
            R("unit triangle: equal heights γ↓ + γ↑ = 4, γ↓ − γ↑ = 4, −2, 2", unit, new long[] { 4, 1, 3 }, new long[] { 0, 3, 1 }),
            R("unit triangle: sites 0 and 2 at one rate point", unit, new long[] { 2, 1, 2 }, new long[] { 1, 1, 1 }),
            R("unit triangle: three rate points off a line", unit, new long[] { 2, 1, 3 }, new long[] { 1, 3, 2 }),
        };
    }

    public IReadOnlyList<Reading> Readings { get; }

    public EqualTrianglePairingWitness()
    {
        Readings = Rows.Select(Read).ToArray();
    }

    public static Reading Read(Row row)
    {
        var verdict = EqualTrianglePairingClaim.Decide(row.Bonds, row.Cool, row.Heat);
        var m = ExactThermalGenerator.Build(FluxAreaFifthMomentClaim.Terms(3, row.Bonds), row.Cool, row.Heat);
        var m2 = ExactThermalGenerator.Multiply(m, m);
        var m3 = ExactThermalGenerator.Multiply(m2, m);
        var m4 = ExactThermalGenerator.Multiply(m3, m);
        var m5 = ExactThermalGenerator.Multiply(m4, m);
        BigInteger Real(GaussianInteger z, int k) =>
            z.Im.IsZero ? z.Re : throw new InvalidOperationException($"{row.Label}: Tr(M^{k}) has an imaginary part");
        var fifth = Real(ExactThermalGenerator.Trace(m5), 5);
        var seventh = Real(ExactThermalGenerator.TraceOfProduct(m3, m4), 7);
        var ninth = Real(ExactThermalGenerator.TraceOfProduct(m4, m5), 9);

        var poly = GaussianMatrixCharpoly.Characteristic(m);
        bool paired = true;
        for (int k = 1; k < poly.Length; k += 2) paired &= poly[k] == GaussianInteger.Zero;

        BigInteger? law = verdict switch
        {
            V.BrokenAtFifth => FluxAreaFifthMomentClaim.FifthCentredTrace(3, row.Bonds, row.Cool, row.Heat),
            V.BrokenAtSeventh => SeventhMomentDegreeFourClaim.SeventhCentredTrace(3, row.Bonds, row.Cool, row.Heat),
            V.BrokenAtNinth => EqualTrianglePairingClaim.NinthOnEqualHeights(row.Bonds, row.Cool, row.Heat),
            _ => null,
        };
        return new Reading(row, verdict, fifth, seventh, ninth, paired, law);
    }

    public bool AllHold => Readings.All(x => x.Holds);
    public int Paired => Readings.Count(x => EqualTrianglePairingClaim.IsPaired(x.Verdict));
    public int CountOf(V v) => Readings.Count(x => x.Verdict == v);
    public bool NotVacuous => CountOf(V.PairedOneRatio) >= 3 && CountOf(V.PairedSiteReflection) >= 2 && CountOf(V.BrokenAtFifth) >= 2
                              && CountOf(V.BrokenAtSeventh) >= 3 && CountOf(V.BrokenAtNinth) >= 2;

    public string DisplayName => "EqualTrianglePairingWitness (the verdict from the rates against the odd traces and the characteristic polynomial)";

    public string Summary =>
        $"{Readings.Count(x => x.Holds)} of {Readings.Count} rows on two equal triangles: the verdict read from the rates alone " +
        $"{(AllHold ? "agrees" : "does NOT always agree")} with the generator, {Paired} rows paired by a carrier (characteristic polynomial even, " +
        $"fifth, seventh and ninth zero), {CountOf(V.BrokenAtFifth)} broken at the fifth, {CountOf(V.BrokenAtSeventh)} at the seventh and " +
        $"{CountOf(V.BrokenAtNinth)} at the ninth, each named trace equal to its closed form {(NotVacuous ? "(not vacuous)" : "(VACUOUS)")}.";

    public IEnumerable<IInspectable> Children
    {
        get
        {
            foreach (var x in Readings)
                yield return new InspectableNode(
                    displayName: x.Row.Label,
                    summary: $"verdict {x.Verdict}; Tr M⁵ = {x.Fifth}, Tr M⁷ = {x.Seventh}, Tr M⁹ = {x.Ninth}" +
                             (x.Law is BigInteger l ? $", closed form {l}" : "") +
                             $"; characteristic polynomial {(x.Paired ? "even (paired)" : "not even (unpaired)")}: " +
                             (x.Holds ? "agrees" : "DISAGREES") + ".",
                    provenance: NodeProvenance.Live);
        }
    }

    public InspectablePayload Payload => InspectablePayload.Empty;
}
