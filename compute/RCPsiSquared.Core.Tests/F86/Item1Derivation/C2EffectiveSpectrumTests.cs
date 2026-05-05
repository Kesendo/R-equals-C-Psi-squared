using System;
using System.Linq;
using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.CoherenceBlocks;
using RCPsiSquared.Core.Decomposition;
using RCPsiSquared.Core.F86.Item1Derivation;
using RCPsiSquared.Core.Knowledge;
using Xunit;

namespace RCPsiSquared.Core.Tests.F86.Item1Derivation;

public class C2EffectiveSpectrumTests
{
    [Theory]
    [InlineData(5, 0.5)]
    [InlineData(5, 1.0)]
    [InlineData(5, 1.5)]
    [InlineData(5, 2.0)]
    [InlineData(5, 2.5)]
    [InlineData(5, 3.0)]
    [InlineData(7, 0.5)]
    [InlineData(7, 1.5)]
    [InlineData(7, 3.0)]
    public void Eigenvalues_MatchNumericalEvd_AtMachinePrecision(int N, double Q)
    {
        var block = new CoherenceBlock(N, n: 1, gammaZero: 0.05);
        var spectrum = C2EffectiveSpectrum.Build(block);
        var effective = FourModeEffective.Build(block);

        for (int b = 0; b < block.NumBonds; b++)
        {
            var ours = spectrum.Eigenvalues(Q, b);
            var theirs = effective.LEffAtQ(Q).Evd().EigenValues.ToArray();

            // Sort both by Re desc, then Im asc — stable comparison
            var oursSorted = ours.OrderByDescending(z => z.Real).ThenBy(z => z.Imaginary).ToArray();
            var theirsSorted = theirs.OrderByDescending(z => z.Real).ThenBy(z => z.Imaginary).ToArray();

            for (int i = 0; i < 4; i++)
            {
                var diff = (oursSorted[i] - theirsSorted[i]).Magnitude;
                Assert.True(diff < 1e-10,
                    $"λ_{i}({Q}, b={b}) at N={N}: ours={oursSorted[i]}, theirs={theirsSorted[i]}, diff={diff:E2}");
            }
        }
    }

    [Fact]
    public void Tier_IsTier1DerivedOrTier2Verified()
    {
        var block = new CoherenceBlock(N: 5, n: 1, gammaZero: 0.05);
        var spectrum = C2EffectiveSpectrum.Build(block);
        Assert.True(spectrum.Tier == Tier.Tier1Derived || spectrum.Tier == Tier.Tier2Verified);
    }

    [Fact]
    public void IsAnalyticallyDerived_ConsistentWithTier()
    {
        var block = new CoherenceBlock(N: 5, n: 1, gammaZero: 0.05);
        var spectrum = C2EffectiveSpectrum.Build(block);
        if (spectrum.IsAnalyticallyDerived) Assert.Equal(Tier.Tier1Derived, spectrum.Tier);
        else Assert.Equal(Tier.Tier2Verified, spectrum.Tier);
    }

    [Fact]
    public void Build_ThrowsIfNotC2()
    {
        var block = new CoherenceBlock(N: 5, n: 2, gammaZero: 0.05);  // c=3
        Assert.Throws<ArgumentException>(() => C2EffectiveSpectrum.Build(block));
    }

    [Fact]
    public void Eigenvalues_OutOfRangeBond_Throws()
    {
        var block = new CoherenceBlock(N: 5, n: 1, gammaZero: 0.05);
        var spectrum = C2EffectiveSpectrum.Build(block);
        Assert.Throws<ArgumentOutOfRangeException>(() => spectrum.Eigenvalues(1.0, -1));
        Assert.Throws<ArgumentOutOfRangeException>(() => spectrum.Eigenvalues(1.0, block.NumBonds));
    }

    [Theory]
    [InlineData(5)]
    [InlineData(6)]
    [InlineData(7)]
    [InlineData(8)]
    public void Eigenvalues_BondIndependent_AtCEquals2(int N)
    {
        // At c=2 the spectrum is bond-independent by construction (uniform-J means
        // L_eff = D_eff + Q·γ₀·Σ_b V_b, and the Σ_b is the same for any bond index passed).
        // This property test pins the contract.
        var block = new CoherenceBlock(N, n: 1, gammaZero: 0.05);
        var spectrum = C2EffectiveSpectrum.Build(block);
        double Q = 1.5;
        var reference = spectrum.Eigenvalues(Q, bond: 0);
        for (int b = 1; b < block.NumBonds; b++)
        {
            var actual = spectrum.Eigenvalues(Q, bond: b);
            for (int i = 0; i < 4; i++)
            {
                var diff = (actual[i] - reference[i]).Magnitude;
                Assert.True(diff < 1e-12,
                    $"Eigenvalues should be bond-independent at c=2: bond=0 vs bond={b}, λ_{i} diff={diff:E2}");
            }
        }
    }

    [Theory]
    [InlineData(5)]
    [InlineData(6)]
    [InlineData(7)]
    [InlineData(8)]
    public void LEffAtQ_MatchesFourModeEffective_AtMachinePrecision(int N)
    {
        // Cross-check: our bond-summed L_eff(Q) assembly must match FourModeEffective's
        // direct projection L_eff(Q) = D_eff + (Q·γ₀)·M_h_total_eff entry-by-entry.
        var block = new CoherenceBlock(N, n: 1, gammaZero: 0.05);
        var spectrum = C2EffectiveSpectrum.Build(block);
        var effective = FourModeEffective.Build(block);

        foreach (double Q in new[] { 0.5, 1.0, 2.0 })
        {
            var ours = spectrum.LEffAtQ(Q);
            var theirs = effective.LEffAtQ(Q);
            var diff = (ours - theirs).FrobeniusNorm();
            Assert.True(diff < 1e-12, $"L_eff(Q={Q}) at N={N}: ‖ours − theirs‖_F = {diff:E2}");
        }
    }

    [Theory]
    [InlineData(5)]
    [InlineData(7)]
    public void Eigenvalues_TraceMatchesDEffPlusQGammaTraceMh(int N)
    {
        // Sanity: Σ λ_i = trace(L_eff(Q)) = trace(D_eff) + Q·γ₀·trace(M_h_total_eff).
        var block = new CoherenceBlock(N, n: 1, gammaZero: 0.05);
        var spectrum = C2EffectiveSpectrum.Build(block);
        double Q = 1.0;
        var eigs = spectrum.Eigenvalues(Q, bond: 0);
        var traceFromEigs = eigs.Aggregate(Complex.Zero, (acc, z) => acc + z);
        var traceFromMatrix = spectrum.LEffAtQ(Q).Diagonal().Sum();
        Assert.True((traceFromEigs - traceFromMatrix).Magnitude < 1e-10,
            $"trace mismatch at N={N}: from eigs = {traceFromEigs}, from matrix = {traceFromMatrix}");
    }

    [Fact]
    public void Tier_TracksBondCouplingTier()
    {
        // The class Tier is currently independent of BondCoupling's Tier (which is Tier2
        // due to A3's even-N σ_0 obstruction). C2's own Tier comes from whether the 4×4
        // factorisation lands; the two are independent obstructions.
        // This test just pins that BondCoupling is reachable and is its expected Tier.
        var block = new CoherenceBlock(N: 5, n: 1, gammaZero: 0.05);
        var spectrum = C2EffectiveSpectrum.Build(block);
        Assert.NotNull(spectrum.BondCoupling);
        Assert.Equal(Tier.Tier2Verified, spectrum.BondCoupling.Tier);
    }

    [Fact]
    public void PendingDerivationNote_IsNullOnTier1Derived_NonNullOnTier2()
    {
        var block = new CoherenceBlock(N: 5, n: 1, gammaZero: 0.05);
        var spectrum = C2EffectiveSpectrum.Build(block);
        if (spectrum.Tier == Tier.Tier1Derived)
            Assert.Null(spectrum.PendingDerivationNote);
        else
        {
            Assert.NotNull(spectrum.PendingDerivationNote);
            Assert.Contains("Tier2Verified", spectrum.PendingDerivationNote!);
        }
    }
}
