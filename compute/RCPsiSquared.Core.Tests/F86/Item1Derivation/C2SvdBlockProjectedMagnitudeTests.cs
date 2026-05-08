using System.Numerics;
using RCPsiSquared.Core.CoherenceBlocks;
using RCPsiSquared.Core.F86.Item1Derivation;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Resonance;

namespace RCPsiSquared.Core.Tests.F86.Item1Derivation;

public class C2SvdBlockProjectedMagnitudeTests
{
    private static CoherenceBlock C2Block(int N) => new(N: N, n: 1, gammaZero: 0.05);

    [Theory]
    [InlineData(5)]
    [InlineData(6)]
    [InlineData(7)]
    [InlineData(8)]
    public void Build_AcrossN5To8_LandsTier1Derived(int N)
    {
        var m = C2SvdBlockProjectedMagnitude.Build(C2Block(N));
        Assert.Equal(Tier.Tier1Derived, m.Tier);
    }

    [Theory]
    [InlineData(5)]
    [InlineData(6)]
    [InlineData(7)]
    [InlineData(8)]
    public void WitnessCount_EqualsNumBonds(int N)
    {
        var m = C2SvdBlockProjectedMagnitude.Build(C2Block(N));
        Assert.Equal(m.Block.NumBonds, m.Witnesses.Count);
    }

    [Theory]
    [InlineData(5)]
    [InlineData(6)]
    [InlineData(7)]
    [InlineData(8)]
    public void EveryMagnitude_IsRealNonNegative(int N)
    {
        var m = C2SvdBlockProjectedMagnitude.Build(C2Block(N));
        foreach (var w in m.Witnesses)
        {
            Assert.True(w.MagnitudeSquared >= 0,
                $"bond {w.Bond} magnitude must be non-negative; got {w.MagnitudeSquared:G6}");
            Assert.True(w.NormalisedRatio >= 0,
                $"bond {w.Bond} normalised ratio must be non-negative; got {w.NormalisedRatio:G6}");
        }
    }

    [Theory]
    [InlineData(5)]
    [InlineData(6)]
    [InlineData(7)]
    [InlineData(8)]
    public void SumRuleResidual_BelowTolerance(int N)
    {
        var m = C2SvdBlockProjectedMagnitude.Build(C2Block(N));
        Assert.True(m.SumRuleResidual < C2SvdBlockProjectedMagnitude.SumRuleTolerance,
            $"N={N}: F73 sum-rule residual {m.SumRuleResidual:G3} should be < tolerance {C2SvdBlockProjectedMagnitude.SumRuleTolerance:G3}");
    }

    [Theory]
    [InlineData(5)]
    [InlineData(7)]
    public void OddN_NonDegenerate_ProjectedMagnitude_AgreesWithSingleVectorReference(int N)
    {
        // At odd N (RankTop=1) the projector lift reduces to the single-vector outer products,
        // so |V_b|²_proj = |⟨u_0|M_b|v_0⟩|² which equals C2BondCoupling.SvdBlockEntry(b, 2, 3)
        // squared in magnitude.
        var block = C2Block(N);
        var m = C2SvdBlockProjectedMagnitude.Build(block);
        var bondCoupling = C2BondCoupling.Build(block);
        Assert.Equal(1, m.Projector.RankTop);

        for (int b = 0; b < block.NumBonds; b++)
        {
            Complex vEntry = bondCoupling.SvdBlockEntry(b, 2, 3);
            double singleVectorSquared = vEntry.Magnitude * vEntry.Magnitude;
            double projectedSquared = m.Witnesses[b].MagnitudeSquared;
            Assert.Equal(singleVectorSquared, projectedSquared, precision: 10);
        }
    }

    [Theory]
    [InlineData(6)]
    [InlineData(8)]
    public void EvenN_DegenerateTop_ProjectedMagnitude_DiffersFromSingleVectorAtSomeBond(int N)
    {
        // At even N (RankTop=2) the projector includes both directions of the degenerate
        // top eigenspace. The single-vector V_b[2,3]² captures only the library-tiebreaker
        // direction; the projected version is generally LARGER (or at minimum different)
        // because it sums the contributions across the 2D top eigenspace. At least one bond
        // must show a non-trivial difference.
        var block = C2Block(N);
        var m = C2SvdBlockProjectedMagnitude.Build(block);
        var bondCoupling = C2BondCoupling.Build(block);
        Assert.Equal(2, m.Projector.RankTop);

        bool foundDifference = false;
        for (int b = 0; b < block.NumBonds; b++)
        {
            Complex vEntry = bondCoupling.SvdBlockEntry(b, 2, 3);
            double singleVectorSquared = vEntry.Magnitude * vEntry.Magnitude;
            double projectedSquared = m.Witnesses[b].MagnitudeSquared;
            if (Math.Abs(projectedSquared - singleVectorSquared) > 1e-6)
            {
                foundDifference = true;
                break;
            }
        }
        Assert.True(foundDifference,
            $"N={N}: at degenerate even N the projector lift must differ from the single-vector form on at least one bond (otherwise the lift adds nothing)");
    }

    [Theory]
    [InlineData(5)]
    [InlineData(6)]
    [InlineData(7)]
    [InlineData(8)]
    public void EndpointMean_DiffersFromInteriorMean(int N)
    {
        // The bond-class signature is the empirical anchor of (a''): Endpoint and Interior
        // must produce distinguishable mean magnitudes. Existing C2BondCoupling
        // CrossBlockWitnesses already shows Endpoint < Interior in |V_b cross|_F at c=2
        // N=5..8 (PROOF_F86_QPEAK Statement B2). Our SVD-block magnitude is a different
        // observable but should still split per bond class (the question (a'') asks is
        // whether HWHM-lift maps to it via a single function).
        var m = C2SvdBlockProjectedMagnitude.Build(C2Block(N));
        double endpoint = m.NormalisedRatioMean(BondClass.Endpoint);
        double interior = m.NormalisedRatioMean(BondClass.Interior);
        Assert.False(double.IsNaN(endpoint), $"N={N}: Endpoint mean must be defined");
        Assert.False(double.IsNaN(interior), $"N={N}: Interior mean must be defined");
        Assert.NotEqual(endpoint, interior);
    }

    [Theory]
    [InlineData(5)]
    [InlineData(6)]
    [InlineData(7)]
    [InlineData(8)]
    public void EndpointBonds_AreFirstAndLast(int N)
    {
        var m = C2SvdBlockProjectedMagnitude.Build(C2Block(N));
        int last = m.Block.NumBonds - 1;
        Assert.Equal(BondClass.Endpoint, m.Witnesses[0].BondClass);
        Assert.Equal(BondClass.Endpoint, m.Witnesses[last].BondClass);
        for (int b = 1; b < last; b++)
            Assert.Equal(BondClass.Interior, m.Witnesses[b].BondClass);
    }

    [Fact]
    public void Build_RejectsNonC2Block()
    {
        var block = new CoherenceBlock(N: 5, n: 2, gammaZero: 0.05);
        Assert.Equal(3, block.C);
        Assert.Throws<ArgumentException>(() => C2SvdBlockProjectedMagnitude.Build(block));
    }

    [Fact]
    public void Anchor_References_F86Item1()
    {
        var m = C2SvdBlockProjectedMagnitude.Build(C2Block(5));
        Assert.Contains("PROOF_F86_QPEAK", m.Anchor);
    }

    [Theory]
    [InlineData(5)]
    [InlineData(8)]
    public void Children_IncludeWitnessGroupAndSumRuleNode(int N)
    {
        IInspectable claim = C2SvdBlockProjectedMagnitude.Build(C2Block(N));
        var labels = claim.Children.Select(c => c.DisplayName).ToList();
        Assert.Contains("SumRuleResidual", labels);
        Assert.Contains("Witnesses (per bond)", labels);
    }
}
