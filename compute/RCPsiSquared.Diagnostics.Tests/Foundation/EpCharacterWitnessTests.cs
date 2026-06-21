using System;
using System.Linq;
using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Numerics;
using RCPsiSquared.Diagnostics.Foundation;
using Xunit;

namespace RCPsiSquared.Diagnostics.Tests.Foundation;

/// <summary>The live artifact-free EP-character witness. The gate-first tests assert the instrument
/// reads the KNOWN toys correctly (a defective Jordan → DEFECTIVE, a diabolic diag → DIABOLIC); the
/// horizon tests assert the single-excitation {0,2} freezing pair reads DEFECTIVE at every N=2..5 with
/// departure-from-normality ≈ 4 and geometric mult 1 &lt; algebraic mult 2 — reproducing
/// <c>simulations/_review_coherence_horizon_ep.py</c> in C#, without an eig eigenvector pairing.</summary>
public class EpCharacterWitnessTests
{
    private static System.Collections.Generic.List<IInspectable> Children(EpCharacterWitness w) =>
        ((IInspectable)w).Children.ToList();

    // ---- GATE: the instrument self-validates on known answers before any horizon claim ----

    [Fact]
    public void Gate_Passes_ToyJordanDefective_DiagDiabolic()
    {
        Assert.True(new EpCharacterWitness().GatePassed,
            "the gate must pass: the toy Jordan reads DEFECTIVE and diag(-2γ,-2γ) reads DIABOLIC");
    }

    [Fact]
    public void Gate_ToyDefective_ReadsDefective_WithDepartureFour()
    {
        var r = EpCharacter.Characterize(EpCharacterWitness.ToyDefective(), new Complex(-4.0, 0), radius: 0.5);
        Assert.Equal(EpCharacter.EpKind.Defective, r.Kind);
        Assert.Equal(1, r.Geometric);
        Assert.Equal(2, r.Algebraic);
        Assert.Equal(4.0, r.Departure, 3);   // the F86a toy_Leff EP: dep = 4 exactly
    }

    [Fact]
    public void Gate_ToyDiabolic_ReadsDiabolic_WithDepartureZero()
    {
        var r = EpCharacter.Characterize(EpCharacterWitness.ToyDiabolic(), new Complex(-2.0, 0), radius: 0.5);
        Assert.Equal(EpCharacter.EpKind.Diabolic, r.Kind);
        Assert.Equal(2, r.Geometric);
        Assert.Equal(2, r.Algebraic);
        Assert.True(r.Departure < 1e-6, $"diabolic departure {r.Departure} should be ≈ 0");
    }

    [Fact]
    public void Gate_Node_IsSurfaced_AndReportsPass()
    {
        var gate = Children(new EpCharacterWitness()).First(c => c.DisplayName.Contains("the gate"));
        Assert.Contains("GATE PASSED", gate.Summary);
        var leaves = gate.Children.ToList();
        Assert.Contains(leaves, l => l.DisplayName.Contains("known-DEFECTIVE") && l.Summary.Contains("PASS"));
        Assert.Contains(leaves, l => l.DisplayName.Contains("known-DIABOLIC") && l.Summary.Contains("PASS"));
    }

    // ---- the Q* ladder + SE-block correctness ----

    [Theory]
    [InlineData(2, 1.0)]
    [InlineData(3, 1.4142135623730951)] // sqrt(2)
    [InlineData(4, 1.87874)]
    [InlineData(5, 2.37367)]
    public void Qstar_ReproducesTheCoherenceHorizonLadder(int n, double expected)
    {
        double q = new EpCharacterWitness().Qstar(n);
        Assert.Equal(expected, q, 4);
    }

    [Fact]
    public void SeBlock_SpectrumIsASubSpectrumOfTheFullLiouvillian_N3()
    {
        // The single-excitation block is an invariant sub-block, so its spectrum is a subset of the
        // full 4^N Liouvillian's (the Python coherence_horizon_se_block._assert_full_reduction).
        int n = 3; double j = 1.0, g = 0.5;
        var seEvals = EpCharacterWitness.Lse(n, j, g).Evd().EigenValues.ToArray();
        var fullEvals = FullLiouvillianXy(n, j, g).Evd().EigenValues.ToArray();
        foreach (var s in seEvals)
        {
            double nearest = fullEvals.Min(f => (f - s).Magnitude);
            Assert.True(nearest < 1e-9, $"SE eigenvalue {s} is not in the full spectrum (nearest {nearest:E2})");
        }
    }

    // ---- the horizon verdict: artifact-free DEFECTIVE at every N=2..5 ----

    [Theory]
    [InlineData(2)]
    [InlineData(3)]
    [InlineData(4)]
    [InlineData(5)]
    public void HorizonVerdict_IsDefective_AtEveryN(int n)
    {
        var w = new EpCharacterWitness();
        var r = w.CharacterizeAtQ(n, w.Qstar(n) * EpCharacterWitness.ProbeDetuning);
        Assert.NotNull(r);
        var reading = r!.Value;
        Assert.Equal(EpCharacter.EpKind.Defective, reading.Kind);
        Assert.Equal(2, reading.Algebraic);          // a coalescing pair is enclosed
        Assert.Equal(1, reading.Geometric);          // geo < alg ⟹ a Jordan block (DEFECTIVE), not diabolic
        Assert.True(reading.Departure > 3.5,         // dep ≈ 4 (N=2,3) drifting to ≈3.75 (N=5)
            $"N={n}: departure-from-normality {reading.Departure} should be bounded away from 0 (≈4)");
        Assert.True(reading.EigenvectorMergeCos > 0.99,
            $"N={n}: the two compression eigenvectors should be merging (|cos| {reading.EigenvectorMergeCos} → 1)");
    }

    [Theory]
    [InlineData(2, 4.0)]
    [InlineData(3, 4.0)]
    [InlineData(4, 3.894)]
    [InlineData(5, 3.759)]
    public void HorizonVerdict_DepartureMatchesThePythonReference(int n, double expectedDep)
    {
        // The departure-from-normality at the 1.002·Q* probe, bit-for-bit against
        // simulations/_review_coherence_horizon_ep.py (4.0000 / 4.0000 / 3.8941 / 3.7589).
        var w = new EpCharacterWitness();
        var reading = w.CharacterizeAtQ(n, w.Qstar(n) * EpCharacterWitness.ProbeDetuning)!.Value;
        Assert.Equal(expectedDep, reading.Departure, 2);
    }

    [Fact]
    public void HorizonVerdict_FreezingPairSitsAtReMinusTwoGamma()
    {
        // The freezing {0,2}-coherence pair sits on the Re = −2γ band (γ=1 in the probe convention):
        // Re ≈ −2 at N=2,3 (the clean 2×2), drifting slightly more negative as the pair is dressed.
        var w = new EpCharacterWitness();
        foreach (int n in new[] { 2, 3 })
        {
            var reading = w.CharacterizeAtQ(n, w.Qstar(n) * EpCharacterWitness.ProbeDetuning)!.Value;
            double meanRe = reading.CompressionEigenvalues.Average(e => e.Real);
            Assert.Equal(-2.0, meanRe, 2);
        }
    }

    // ---- the limit Q → Q* (the defectiveness lives in the limit) ----

    [Fact]
    public void Limit_DepartureStaysBoundedAsSplitShrinks()
    {
        // As Q → Q*(4) the pair-split shrinks but the departure-from-normality stays bounded away from
        // 0 (the DEFECTIVE signature; a diabolic point would have dep → 0 with the split).
        int n = 4;
        var w = new EpCharacterWitness();
        double qstar = w.Qstar(n);
        var wide = w.CharacterizeAtQ(n, qstar * 1.02)!.Value;
        var tight = w.CharacterizeAtQ(n, qstar * 1.002)!.Value;
        Assert.True(tight.Departure > 3.5, $"the tighter probe departure {tight.Departure} should stay ≈ 4");
        // the eigenvectors are MORE merged at the tighter detuning (closer to coalescence)
        Assert.True(tight.EigenvectorMergeCos > wide.EigenvectorMergeCos,
            "the eigenvectors should merge further as Q → Q*");
    }

    // ---- children + JSON surfacing ----

    [Fact]
    public void Witness_SurfacesGateHorizonLimitAndEigContrast()
    {
        var labels = Children(new EpCharacterWitness()).Select(c => c.DisplayName).ToList();
        Assert.Contains(labels, l => l.Contains("the gate"));
        Assert.Contains(labels, l => l.Contains("horizon verdict"));
        Assert.Contains(labels, l => l.Contains("limit Q → Q*"));
        Assert.Contains(labels, l => l.Contains("eig contrast"));
    }

    [Fact]
    public void Witness_SummaryReports_GatePassedAndDefectiveHorizon()
    {
        var s = new EpCharacterWitness().Summary;
        Assert.Contains("GATE PASSED", s);
        Assert.Contains("N=2:Defective", s);
        Assert.Contains("N=5:Defective", s);
    }

    [Fact]
    public void Witness_RendersToJson()
    {
        var json = InspectionJsonExporter.ToJson(new EpCharacterWitness());
        Assert.Contains("departure-from-normality", json);
        Assert.Contains("Defective", json);
    }

    // ---- a full N=2 XY Liouvillian for the SE-block sub-spectrum check ----

    private static Matrix<Complex> FullLiouvillianXy(int n, double j, double g)
    {
        int d = 1 << n;
        // XY Hamiltonian H = (J/2) Σ (X_b X_{b+1} + Y_b Y_{b+1}); built as a real symmetric hopping.
        var h = new double[d, d];
        for (int b = 0; b < n - 1; b++)
            for (int s = 0; s < d; s++)
            {
                int bit0 = (s >> b) & 1, bit1 = (s >> (b + 1)) & 1;
                if (bit0 != bit1) { int s2 = s ^ (1 << b) ^ (1 << (b + 1)); h[s2, s] += j; } // (J/2)(XX+YY) flips a differing pair with weight J
            }
        int dim = d * d;
        var l = Matrix<Complex>.Build.Dense(dim, dim);
        var negI = new Complex(0, -1);
        for (int a = 0; a < d; a++)
            for (int c = 0; c < d; c++)
            {
                double hv = h[a, c];
                if (hv != 0.0) for (int b = 0; b < d; b++) l[a * d + b, c * d + b] += negI * hv;
            }
        for (int b = 0; b < d; b++)
            for (int e = 0; e < d; e++)
            {
                double hv = h[b, e];
                if (hv != 0.0) for (int a = 0; a < d; a++) l[a * d + b, a * d + e] += -negI * hv;
            }
        for (int a = 0; a < d; a++)
            for (int b = 0; b < d; b++)
            {
                int hamm = System.Numerics.BitOperations.PopCount((uint)(a ^ b));
                l[a * d + b, a * d + b] += new Complex(-2.0 * g * hamm, 0.0);
            }
        return l;
    }
}
