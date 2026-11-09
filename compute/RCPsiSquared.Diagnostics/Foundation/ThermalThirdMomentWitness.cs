using System.Numerics;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Symmetry;
using T = RCPsiSquared.Core.Symmetry.ThermalThirdMomentClaim.Term;

namespace RCPsiSquared.Diagnostics.Foundation;

/// <summary>The live lab for <see cref="ThermalThirdMomentClaim"/> (F137, the third moment under cooling and heating):
/// Tr((2L + c)³) from the generator in exact Gaussian integers (<see cref="ExactThermalGenerator"/>) against the
/// claim's closed form 12·4^N·(F − 2·R), read from the Pauli strings alone.
///
/// <para>The rows: the four values the registry cites at Heisenberg N = 3 under cooling 1, as Tr((L + c/2)³) (288 for a
/// uniform X field, a formal setting at zero splitting; −768 for a uniform Z field as written, the lab frame of the
/// field-free chain; 96, CAUGHT_ERRORS' swap counterexample, and −384/5 as written for two Z profiles, read with H scaled to integers, the moment being of degree two
/// in H); settings of its own carrying a field part and an arrow part
/// together, one of them through a string pair two sites apart and beside an anticommuting pair the law must ignore;
/// a constant term the law drops; a bond of all nine letters of WHICH_BONDS_PAIR_UNDER_COOLING, where both parts
/// vanish; equal rates on every site, where the arrow is gone; the conditioned field off its cancellation, −48 in the
/// registry's gate. Rows with a field part alone, an arrow part alone and both pin each factor and their sign.</para></summary>
public sealed class ThermalThirdMomentWitness : IInspectable
{
    public sealed record Row(string Label, IReadOnlyList<T> Terms, IReadOnlyList<long> Cool, IReadOnlyList<long> Heat, BigInteger? RecordedMoment)
    {
        public int N => Cool.Count;
    }

    public sealed record Reading(Row Row, BigInteger Generator, BigInteger Law, BigInteger Field, BigInteger Arrow)
    {
        /// <summary>The generator's cubed trace equals the closed form, and, where the registry records a value of
        /// Tr((L + c/2)³), eight times it.</summary>
        public bool Holds => Generator == Law && (Row.RecordedMoment is not BigInteger r || Generator == 8 * r);
    }

    private static IReadOnlyList<T> Chain(int n, params (string Template, long Coefficient)[] templates)
    {
        var list = new List<T>();
        foreach (var (tpl, c) in templates)
            for (int k = 0; k + tpl.Length <= n; k++)
                list.Add(new T(new string('I', k) + tpl + new string('I', n - k - tpl.Length), c));
        return list;
    }

    public static IReadOnlyList<Row> Rows { get; } = BuildRows();

    private static IReadOnlyList<Row> BuildRows()
    {
        var heis3 = Chain(3, ("XX", 1), ("YY", 1), ("ZZ", 1));
        return new[]
        {
            new Row("Heisenberg N = 3, uniform X field 1, cooling 1, formal at zero splitting (recorded 288)",
                heis3.Concat(Chain(3, ("X", 1))).ToArray(), new long[] { 1, 1, 1 }, new long[] { 0, 0, 0 }, 288),
            new Row("Heisenberg N = 3, uniform Z field 1 as written, cooling 1 (recorded −768)",
                heis3.Concat(Chain(3, ("Z", 1))).ToArray(), new long[] { 1, 1, 1 }, new long[] { 0, 0, 0 }, -768),
            new Row("Heisenberg N = 3, Z fields (1/2, −1/2, 0), cooling 1, H doubled (recorded 96, read as 4·96)",
                Chain(3, ("XX", 2), ("YY", 2), ("ZZ", 2)).Concat(new[] { new T("ZII", 1), new T("IZI", -1) }).ToArray(),
                new long[] { 1, 1, 1 }, new long[] { 0, 0, 0 }, 4 * 96),
            new Row("Heisenberg N = 3, Z fields (3/10, −1/5, 1/2), cooling 1, H times 10 (recorded −384/5, read as 100·(−384/5))",
                Chain(3, ("XX", 10), ("YY", 10), ("ZZ", 10)).Concat(new[] { new T("ZII", 3), new T("IZI", -2), new T("IIZ", 5) }).ToArray(),
                new long[] { 1, 1, 1 }, new long[] { 0, 0, 0 }, -7680),
            new Row("N = 2, XX + 3ZZ, X field 2 on site 0, Y field −1 and Z field 1 on site 1",
                new[] { new T("XX", 1), new T("ZZ", 3), new T("XI", 2), new T("IY", -1), new T("IZ", 1) },
                new long[] { 2, 1 }, new long[] { 1, 3 }, null),
            new Row("N = 3, X₀X₂ beside −2·X₀Z₁X₂, Z₀ beside Z₀Z₂, the anticommuting pair X₀Z₁ and Y₀Z₁, Y field 3 on site 1",
                new[] { new T("XIX", 1), new T("XZX", -2), new T("ZII", 2), new T("ZIZ", 1), new T("XZI", 1), new T("YZI", 2), new T("IYI", 3) },
                new long[] { 3, 1, 2 }, new long[] { 1, 2, 0 }, null),
            new Row("N = 3, X field 1 on site 0, Z field −2 on site 2 beside Z₁Z₂, XY − YX on (0, 1)",
                new[] { new T("XII", 1), new T("IIZ", -2), new T("IZZ", 3), new T("XYI", 1), new T("YXI", -1) },
                new long[] { 2, 3, 1 }, new long[] { 1, 0, 2 }, null),
            new Row("N = 2, a constant 5 beside Z₀ and 2·Z₀Z₁ (the law reads the traceless part)",
                new[] { new T("II", 5), new T("ZI", 1), new T("ZZ", 2) },
                new long[] { 1, 2 }, new long[] { 0, 1 }, null),
            new Row("N = 3, the nine letters on every bond, XX + 2XY − YX + 3YY − ZZ + XZ − 2YZ + ZX + 3ZY, cooling (2, 1, 3)",
                Chain(3, ("XX", 1), ("XY", 2), ("YX", -1), ("YY", 3), ("ZZ", -1), ("XZ", 1), ("YZ", -2), ("ZX", 1), ("ZY", 3)),
                new long[] { 2, 1, 3 }, new long[] { 0, 0, 0 }, null),
            new Row("N = 2, equal rates 2 on both sites, X₀ + 2·Y₁ + Z₀ + Z₀Z₁ (the arrow is gone)",
                new[] { new T("XI", 1), new T("IY", 2), new T("ZI", 1), new T("ZZ", 1) },
                new long[] { 2, 2 }, new long[] { 2, 2 }, null),
            new Row("the conditioned field X₀ + X₀Z₁ off its cancellation, site 0 at 1/1, site 1 cooled at 2 (the gate's −48)",
                new[] { new T("XI", 1), new T("XZ", 1) }, new long[] { 1, 2 }, new long[] { 1, 0 }, -48),
        };
    }

    public IReadOnlyList<Reading> Readings { get; }

    public ThermalThirdMomentWitness()
    {
        Readings = Rows.Select(Read).ToArray();
    }

    public static Reading Read(Row row)
    {
        var m = ExactThermalGenerator.Build(row.Terms, row.Cool, row.Heat);
        var cubed = ExactThermalGenerator.OddPowerSums(m, 3)[3];
        var (f, r) = ThermalThirdMomentClaim.Parts(row.Terms, row.Cool, row.Heat);
        return new Reading(row, cubed, ThermalThirdMomentClaim.CubedCentredTrace(row.Terms, row.Cool, row.Heat), f, r);
    }

    public bool AllHold => Readings.All(x => x.Holds);
    public int BothPartsNonzero => Readings.Count(x => !x.Field.IsZero && !x.Arrow.IsZero);
    public bool NotVacuous => BothPartsNonzero >= 3 && Readings.Any(x => !x.Generator.IsZero && x.Arrow.IsZero)
                              && Readings.Any(x => !x.Generator.IsZero && x.Field.IsZero);

    public string DisplayName => "ThermalThirdMomentWitness (Tr((2L + c)³) against 12·4^N·(F − 2R))";

    public string Summary =>
        $"{Readings.Count(x => x.Holds)} of {Readings.Count} settings: the generator's Tr((2L + c)³) {(AllHold ? "equals" : "does NOT always equal")} " +
        $"the closed form 12·4^N·(F − 2R), the registry's four recorded values and the gate's −48 included; {BothPartsNonzero} carry a field " +
        $"part and an arrow part together " +
        $"{(NotVacuous ? "(not vacuous)" : "(VACUOUS)")}.";

    public IEnumerable<IInspectable> Children
    {
        get
        {
            foreach (var x in Readings)
                yield return new InspectableNode(
                    displayName: x.Row.Label,
                    summary: $"generator {x.Generator}, law {x.Law} (F = {x.Field}, R = {x.Arrow})" +
                             (x.Row.RecordedMoment is BigInteger rec ? $", Tr((L + c/2)³) = {x.Generator}/8 against {rec}" : "") +
                             $": {(x.Holds ? "equal" : "DIFFERENT")}.",
                    provenance: NodeProvenance.Live);
        }
    }

    public InspectablePayload Payload => InspectablePayload.Empty;
}
