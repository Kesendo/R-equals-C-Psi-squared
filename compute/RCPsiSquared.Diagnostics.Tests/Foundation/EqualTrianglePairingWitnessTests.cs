using System.Linq;
using System.Numerics;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Diagnostics.Foundation;
using Xunit;
using B = RCPsiSquared.Core.Symmetry.FluxAreaFifthMomentClaim.Bond;
using V = RCPsiSquared.Core.Symmetry.EqualTrianglePairingClaim.Verdict;

namespace RCPsiSquared.Diagnostics.Tests.Foundation;

/// <summary>The live witness of <see cref="EqualTrianglePairingClaim"/> (Theorem 3 of
/// PROOF_SEVENTH_TRACE_DEGREE_FOUR_TRIANGLE, F137): the verdict from the rates alone against the generator's odd
/// traces and its characteristic polynomial. Every comparison is exact.</summary>
public class EqualTrianglePairingWitnessTests
{
    private static readonly EqualTrianglePairingWitness W = new();
    private static readonly B[] E = { new(0, 1, 2, 1, 1), new(1, 2, 2, 1, 1), new(2, 0, 2, 1, 1) };

    [Fact]
    public void EveryRow_Holds()
    {
        foreach (var x in W.Readings)
            Assert.True(x.Holds, $"{x.Row.Label}: {x.Verdict}, Tr5 {x.Fifth}, Tr7 {x.Seventh}, Tr9 {x.Ninth}, law {x.Law}, paired {x.Paired}");
    }

    [Fact]
    public void Rows_AreNotVacuous() => Assert.True(W.NotVacuous, W.Summary);

    /// <summary>The ninth on equal heights as a number written here: on the unit triangle at γ↓ + γ↑ = 4 (S = 2) with
    /// γ↓ − γ↑ = (4, −2, 2), the arrows δ = (2, −1, 1) and V(δ) = 6, 36864·Im·(r² + 7Δ²)·S·V(δ) = 36864·1·1·2·6, times 512.</summary>
    [Fact]
    public void Ninth_IsTheNumberOfTheProof()
    {
        var x = W.Readings.Single(r => r.Row.Label.StartsWith("unit triangle: equal heights"));
        Assert.Equal(V.BrokenAtNinth, x.Verdict);
        Assert.Equal((BigInteger)36864 * 2 * 6 * 512, x.Ninth);
    }

    /// <summary>Theorem 2 with b = 0 by hand: on γ↑ = 1 + γ↓ (e = 1, m = 1) at γ↓ = (1, 2, 4) the seventh trace is
    /// −3360·e·Im·(1 + m)²·V(γ↓) = −3360·11·4·6, times 128.</summary>
    [Fact]
    public void Seventh_IsTheLineForm()
    {
        var x = W.Readings.Single(r => r.Row.Label.EndsWith("on the line γ↑ = 1 + γ↓"));
        Assert.Equal(V.BrokenAtSeventh, x.Verdict);
        Assert.Equal((BigInteger)(-3360) * 11 * 4 * 6 * 128, x.Seventh);
    }

    /// <summary>Every carrier row and every break row agree between the two roads: paired by the characteristic
    /// polynomial exactly where the verdict says paired.</summary>
    [Fact]
    public void CharacteristicPolynomial_AgreesWithTheVerdict()
    {
        foreach (var x in W.Readings)
            Assert.Equal(EqualTrianglePairingClaim.IsPaired(x.Verdict), x.Paired);
        Assert.Contains(W.Readings, x => x.Paired);
        Assert.Contains(W.Readings, x => !x.Paired);
    }

    [Fact]
    public void Decide_RefusesWhatTheTheoremDoesNotCover()
    {
        long[] cool = { 2, 1, 3 }, heat = { 1, 3, 2 };
        // unequal |t|
        Assert.Throws<System.ArgumentException>(() => EqualTrianglePairingClaim.Decide(
            new[] { new B(0, 1, 1, 0, 0), new B(1, 2, 2, 0, 0), new B(2, 0, 0, 1, 0) }, cool, heat));
        // two different Δ
        Assert.Throws<System.ArgumentException>(() => EqualTrianglePairingClaim.Decide(
            new[] { new B(0, 1, 2, 1, 1), new B(1, 2, 2, 1, 1), new B(2, 0, 2, 1, 2) }, cool, heat));
        // flux 0
        Assert.Throws<System.ArgumentException>(() => EqualTrianglePairingClaim.Decide(
            new[] { new B(0, 1, 1, 0, 0), new B(1, 2, 1, 0, 0), new B(2, 0, 1, 0, 0) }, cool, heat));
        // a negative rate
        Assert.Throws<System.ArgumentException>(() => EqualTrianglePairingClaim.Decide(E, new long[] { 2, -1, 3 }, heat));
        // the ninth's closed form only on equal heights
        Assert.Throws<System.ArgumentException>(() => EqualTrianglePairingClaim.NinthOnEqualHeights(E, cool, heat));
    }

    /// <summary>The claim's verdict depends on the rate points only: permuting the sites with their rates, by a cycle
    /// and by a transposition, changes nothing in the verdict, on every row.</summary>
    [Fact]
    public void Verdict_IsBlindToTheOrderOfTheSites()
    {
        foreach (var row in EqualTrianglePairingWitness.Rows)
        {
            var v = EqualTrianglePairingClaim.Decide(row.Bonds, row.Cool, row.Heat);
            foreach (var p in new[] { new[] { 1, 2, 0 }, new[] { 1, 0, 2 } })
            {
                var cool = p.Select(i => row.Cool[i]).ToArray();
                var heat = p.Select(i => row.Heat[i]).ToArray();
                Assert.Equal(v, EqualTrianglePairingClaim.Decide(row.Bonds, cool, heat));
            }
        }
    }

    /// <summary>The ninth's constants can fail: 36864·Im·(r² + 7Δ²)·S·V(δ) with 36864 or the 7 moved by one misses a
    /// ninth row of the generator (the Δ = 1 row exercises the 7).</summary>
    [Fact]
    public void WrongNinth_MissesARow()
    {
        var ninthRows = W.Readings.Where(x => x.Verdict == V.BrokenAtNinth).ToArray();
        Assert.NotEmpty(ninthRows);
        BigInteger Form(EqualTrianglePairingWitness.Reading x, BigInteger a, BigInteger seven)
        {
            var t = FluxAreaFifthMomentClaim.HopAmplitudes(3, x.Row.Bonds);
            var im = (t[0, 1] * t[1, 2] * t[2, 0]).Im;
            var r2 = t[0, 1].Re * t[0, 1].Re + t[0, 1].Im * t[0, 1].Im;
            BigInteger delta = x.Row.Bonds[0].Delta, s = x.Row.Cool[0] + x.Row.Heat[0];
            BigInteger a0 = x.Row.Cool[0] - x.Row.Heat[0], a1 = x.Row.Cool[1] - x.Row.Heat[1], a2 = x.Row.Cool[2] - x.Row.Heat[2];
            return 32 * a * im * (r2 + seven * delta * delta) * s * (a0 - a1) * (a1 - a2) * (a2 - a0);
        }
        Assert.All(ninthRows, x => Assert.Equal(x.Ninth, Form(x, 36864, 7)));
        Assert.Contains(ninthRows, x => x.Ninth != Form(x, 36865, 7));
        Assert.Contains(ninthRows, x => x.Ninth != Form(x, 36864, 8));
    }
}
