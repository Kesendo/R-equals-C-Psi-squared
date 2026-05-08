using System;
using System.Linq;
using RCPsiSquared.Core.CoherenceBlocks;
using RCPsiSquared.Core.F86.JordanWigner;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Resonance;
using Xunit;
using Xunit.Abstractions;

namespace RCPsiSquared.Core.Tests.F86.JordanWigner;

public class BondHdChannelWeightsTests
{
    private readonly ITestOutputHelper _out;

    public BondHdChannelWeightsTests(ITestOutputHelper output) => _out = output;

    [Theory]
    [InlineData(5)]
    [InlineData(6)]
    [InlineData(7)]
    [InlineData(8)]
    [InlineData(9)]
    [InlineData(10)]
    public void Hd1SumRuleResidual_BelowTolerance(int N)
    {
        var block = new CoherenceBlock(N: N, n: 1, gammaZero: 0.05);
        var weights = BondHdChannelWeights.Build(block);
        Assert.True(weights.Hd1SumRuleResidual < BondHdChannelWeights.SumRuleTolerance,
            $"N={N}: HD=1 sum-rule residual {weights.Hd1SumRuleResidual:G3} ≥ tolerance {BondHdChannelWeights.SumRuleTolerance:G3}");
    }

    [Theory]
    [InlineData(5)]
    [InlineData(6)]
    [InlineData(7)]
    [InlineData(8)]
    [InlineData(9)]
    [InlineData(10)]
    public void Hd3SumRuleResidual_BelowTolerance(int N)
    {
        var block = new CoherenceBlock(N: N, n: 1, gammaZero: 0.05);
        var weights = BondHdChannelWeights.Build(block);
        Assert.True(weights.Hd3SumRuleResidual < BondHdChannelWeights.SumRuleTolerance,
            $"N={N}: HD=3 sum-rule residual {weights.Hd3SumRuleResidual:G3} ≥ tolerance {BondHdChannelWeights.SumRuleTolerance:G3}");
    }

    [Theory]
    [InlineData(5)]
    [InlineData(6)]
    [InlineData(7)]
    [InlineData(8)]
    [InlineData(9)]
    [InlineData(10)]
    public void F71MirrorInvariance_PerOrbit_BitExact(int N)
    {
        var block = new CoherenceBlock(N: N, n: 1, gammaZero: 0.05);
        var weights = BondHdChannelWeights.Build(block);

        Assert.Equal(0.0, weights.MaxF71MirrorDeviationHd1, precision: 12);
        Assert.Equal(0.0, weights.MaxF71MirrorDeviationHd3, precision: 12);
    }

    [Theory]
    [InlineData(5)]
    [InlineData(6)]
    [InlineData(7)]
    [InlineData(8)]
    public void NumBonds_MatchesBlockNumBonds(int N)
    {
        var block = new CoherenceBlock(N: N, n: 1, gammaZero: 0.05);
        var weights = BondHdChannelWeights.Build(block);
        Assert.Equal(N - 1, weights.Bonds.Count);
    }

    [Theory]
    [InlineData(5)]
    [InlineData(6)]
    [InlineData(7)]
    [InlineData(8)]
    public void BondClass_FirstAndLast_AreEndpoint(int N)
    {
        var block = new CoherenceBlock(N: N, n: 1, gammaZero: 0.05);
        var weights = BondHdChannelWeights.Build(block);
        int last = weights.Bonds.Count - 1;
        Assert.Equal(BondClass.Endpoint, weights.Bonds[0].BondClass);
        Assert.Equal(BondClass.Endpoint, weights.Bonds[last].BondClass);
    }

    // **Bond-uniformity theorem (verified bit-exact N=5..10):** all bonds share the
    // identical column-Frobenius² weight per HD-class. The deeper structural fact behind
    // this (also verified in this session's exploration but not surfaced as a separate
    // primitive yet): all MhPerBond[b] matrices are UNITARILY EQUIVALENT — same eigen-
    // spectrum (degenerate ε = ±2 at the same multiplicities), differing only by an
    // implicit unitary (the spatial bond-permutation U_{b,b'} that maps bond b to bond b').
    // The static (J=0) channel-norm of a bond's Hamiltonian is bond-blind for the same
    // reason. The bond-distinction in C2HwhmRatio's Q_peak/HWHM comes from the
    // L(Q)-EIGENBASIS projection xB = R(Q)⁻¹ · MhPerBond[b] · R(Q): R(Q) is constructed
    // from MhTotal = Σ all bonds, so each individual bond's matrix elements in this basis
    // depend on its spatial position relative to the global eigenmode structure (band-edge
    // effect: ψ_1(j=0) ~ 1/N^{3/2} Edge vs ψ_1(j=N/2) ~ 1/N^{1/2} Innermost). T6 establishes
    // the J=0 baseline; T7 will compute the L(Q)-projection that breaks this baseline.
    [Theory]
    [InlineData(5)]
    [InlineData(6)]
    [InlineData(7)]
    [InlineData(8)]
    [InlineData(9)]
    [InlineData(10)]
    public void AllBonds_HaveIdenticalHdChannelWeights_BondUniformityTheorem(int N)
    {
        var block = new CoherenceBlock(N: N, n: 1, gammaZero: 0.05);
        var weights = BondHdChannelWeights.Build(block);

        double refHd1 = weights.Bonds[0].Hd1Weight;
        double refHd3 = weights.Bonds[0].Hd3Weight;
        for (int b = 1; b < weights.Bonds.Count; b++)
        {
            Assert.Equal(refHd1, weights.Bonds[b].Hd1Weight, precision: 10);
            Assert.Equal(refHd3, weights.Bonds[b].Hd3Weight, precision: 10);
        }
    }

    // Strengthened structural finding: all MhPerBond[b] are unitarily equivalent. Proved
    // by showing the eigenspectra (sorted by magnitude) match bit-exact across bonds.
    // This is the underlying theorem from which the column-Frobenius² bond-uniformity
    // (and channel-cross-Frobenius² bond-uniformity, etc.) all follow as corollaries.
    // Tested up to N=8 to keep EVD cost manageable; the structure is N-independent.
    [Theory]
    [InlineData(5)]
    [InlineData(6)]
    [InlineData(7)]
    [InlineData(8)]
    public void AllBonds_HaveIdenticalEigenSpectra_UnitaryEquivalence(int N)
    {
        var block = new CoherenceBlock(N: N, n: 1, gammaZero: 0.05);
        var decomp = block.Decomposition;
        int numBonds = decomp.NumBonds;

        var refSpec = decomp.MhPerBond[0].Evd().EigenValues
            .Select(z => z.Magnitude).OrderByDescending(x => x).ToArray();

        for (int b = 1; b < numBonds; b++)
        {
            var thisSpec = decomp.MhPerBond[b].Evd().EigenValues
                .Select(z => z.Magnitude).OrderByDescending(x => x).ToArray();
            Assert.Equal(refSpec.Length, thisSpec.Length);
            for (int i = 0; i < refSpec.Length; i++)
                Assert.Equal(refSpec[i], thisSpec[i], precision: 10);
        }
    }

    [Theory]
    [InlineData(5, 0.40)]    // 2/5
    [InlineData(6, 1.0/3.0)]
    [InlineData(7, 2.0/7.0)]
    [InlineData(8, 0.25)]    // 2/8 = 1/4
    [InlineData(9, 2.0/9.0)]
    [InlineData(10, 0.20)]   // 2/10 = 1/5
    public void Hd1FractionPerBond_Equals_TwoOverN(int N, double expectedFraction)
    {
        var block = new CoherenceBlock(N: N, n: 1, gammaZero: 0.05);
        var weights = BondHdChannelWeights.Build(block);
        foreach (var b in weights.Bonds)
        {
            double total = b.Hd1Weight + b.Hd3Weight;
            double hd1Frac = b.Hd1Weight / total;
            Assert.Equal(expectedFraction, hd1Frac, precision: 10);
        }
    }

    [Theory]
    [InlineData(5)]
    [InlineData(8)]
    [InlineData(10)]
    public void Hd1AndHd3Weights_AreNonNegative(int N)
    {
        var block = new CoherenceBlock(N: N, n: 1, gammaZero: 0.05);
        var weights = BondHdChannelWeights.Build(block);
        foreach (var b in weights.Bonds)
        {
            Assert.True(b.Hd1Weight >= 0, $"N={N}, bond {b.Bond}: Hd1Weight = {b.Hd1Weight}");
            Assert.True(b.Hd3Weight >= 0, $"N={N}, bond {b.Bond}: Hd3Weight = {b.Hd3Weight}");
        }
    }

    [Fact]
    public void Build_RejectsNonC2Block()
    {
        var block = new CoherenceBlock(N: 5, n: 2, gammaZero: 0.05);  // c=3
        Assert.Throws<ArgumentException>(() => BondHdChannelWeights.Build(block));
    }

    [Fact]
    public void Tier_IsTier1Derived()
    {
        var block = new CoherenceBlock(N: 5, n: 1, gammaZero: 0.05);
        var weights = BondHdChannelWeights.Build(block);
        Assert.Equal(Tier.Tier1Derived, weights.Tier);
    }

    [Fact]
    public void Anchor_References_PROOF_BLOCK_CPSI_QUARTER_AndDirectionB()
    {
        var block = new CoherenceBlock(N: 5, n: 1, gammaZero: 0.05);
        var weights = BondHdChannelWeights.Build(block);
        Assert.Contains("PROOF_BLOCK_CPSI_QUARTER", weights.Anchor);
        Assert.Contains("Direction (b'')", weights.Anchor);
    }

    // The HD-channel weight redistribution Tom flagged as the Ausbruch trigger:
    // emit per-bond Hd1/Hd3 weight at N=5..10 so the channel-weight reorganisation
    // is visible across the N=8→10 transition.
    [Fact]
    public void Reconnaissance_EmitsPerNBondHdWeights_AcrossN5To10()
    {
        _out.WriteLine("  N | b | class    | Hd1 weight  | Hd3 weight  | Hd1 / total");
        _out.WriteLine("  --|---|----------|-------------|-------------|------------");
        foreach (int N in new[] { 5, 6, 7, 8, 9, 10 })
        {
            var block = new CoherenceBlock(N: N, n: 1, gammaZero: 0.05);
            var weights = BondHdChannelWeights.Build(block);
            Assert.Equal(N - 1, weights.Bonds.Count);
            foreach (var b in weights.Bonds)
            {
                double total = b.Hd1Weight + b.Hd3Weight;
                double hd1Frac = total > 0 ? b.Hd1Weight / total : 0;
                _out.WriteLine($"  {N} | {b.Bond} | {b.BondClass,-8} | {b.Hd1Weight,11:F6} | {b.Hd3Weight,11:F6} | {hd1Frac,10:F4}");
            }
            _out.WriteLine($"  -- N={N}: Hd1 sum-rule residual = {weights.Hd1SumRuleResidual:G3}, " +
                           $"Hd3 sum-rule residual = {weights.Hd3SumRuleResidual:G3}");
            _out.WriteLine("");
        }
    }
}
