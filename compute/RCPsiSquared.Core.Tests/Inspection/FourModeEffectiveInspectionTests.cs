using RCPsiSquared.Core.CoherenceBlocks;
using RCPsiSquared.Core.Decomposition;
using RCPsiSquared.Core.Decomposition.Views;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Resonance;
using Xunit.Abstractions;

namespace RCPsiSquared.Core.Tests.Inspection;

public class FourModeEffectiveInspectionTests(ITestOutputHelper output)
{
    [Fact]
    public void IsInspectable_ExposesFiveStructuralChildren()
    {
        var block = new CoherenceBlock(N: 5, n: 1, gammaZero: 0.05);
        var eff = FourModeEffective.Build(block);

        IInspectable root = eff;
        var children = root.Children.ToList();
        // Block + Basis + DEffView + MhTotalView + per-bond group + ProbeView + SKernelView
        Assert.Equal(7, children.Count);
    }

    [Fact]
    public void DEffView_PerModeRates_AreMinusTwoGammaAndMinusSixGamma()
    {
        var block = new CoherenceBlock(N: 5, n: 1, gammaZero: 0.05);
        var eff = FourModeEffective.Build(block);

        Assert.Equal(-0.1, eff.DEffView.RateOnC1.Real, 10);
        Assert.Equal(-0.3, eff.DEffView.RateOnC3.Real, 10);
        Assert.Equal(-0.1, eff.DEffView.RateOnU0.Real, 10);
        Assert.Equal(-0.3, eff.DEffView.RateOnV0.Real, 10);
        Assert.True(eff.DEffView.OffDiagonalResidual < 1e-10);
    }

    [Fact]
    public void ProbeView_SvdTopFraction_IsZero_ChannelUniformFraction_IsOne()
    {
        var block = new CoherenceBlock(N: 5, n: 1, gammaZero: 0.05);
        var eff = FourModeEffective.Build(block);

        Assert.True(eff.ProbeView.SvdTopFraction < 1e-20,
            $"probe SVD-top fraction should be 0; got {eff.ProbeView.SvdTopFraction}");
        Assert.Equal(1.0, eff.ProbeView.ChannelUniformFraction, 10);
    }

    [Fact]
    public void MhPerBondViews_BondClasses_MatchEndpointInteriorPattern()
    {
        var block = new CoherenceBlock(N: 6, n: 1, gammaZero: 0.05);
        var eff = FourModeEffective.Build(block);
        var bonds = eff.MhPerBondViews;

        Assert.Equal(5, bonds.Count);
        Assert.Equal(BondClass.Endpoint, bonds[0].BondClass);
        Assert.Equal(BondClass.Interior, bonds[1].BondClass);
        Assert.Equal(BondClass.Interior, bonds[2].BondClass);
        Assert.Equal(BondClass.Interior, bonds[3].BondClass);
        Assert.Equal(BondClass.Endpoint, bonds[4].BondClass);
    }

    [Fact]
    public void MhPerBondViews_CrossBlockFrobenius_IsNonZeroForAllBonds()
    {
        // The bond-position-dependent cross-block coupling drives the F86 universal-shape
        // split. Every bond carries non-trivial cross-block content.
        var block = new CoherenceBlock(N: 5, n: 1, gammaZero: 0.05);
        var eff = FourModeEffective.Build(block);

        foreach (var bond in eff.MhPerBondViews)
            Assert.True(bond.CrossBlockFrobenius > 1e-3,
                $"bond {bond.BondIndex} ({bond.BondClass}) cross-block Frobenius unexpectedly small: {bond.CrossBlockFrobenius}");
    }

    [Fact]
    public void DiagonalRatesIn4Mode_RejectsNon4x4()
    {
        var bad = MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex>.Build.Dense(3, 3);
        Assert.Throws<ArgumentException>(() => new DiagonalRatesIn4Mode(bad));
    }

    [Fact]
    public void Block2x2_RejectsNon2x2()
    {
        var bad = MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex>.Build.Dense(3, 3);
        Assert.Throws<ArgumentException>(() => new Block2x2("bad", bad));
    }

    [Fact]
    public void BondCouplingIn4Mode_RejectsOutOfRangeIndex()
    {
        var matrix = MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex>.Build.Dense(4, 4);
        Assert.Throws<ArgumentOutOfRangeException>(() => new BondCouplingIn4Mode(5, 4, matrix));
        Assert.Throws<ArgumentOutOfRangeException>(() => new BondCouplingIn4Mode(-1, 4, matrix));
    }

    [Fact]
    public void TreeWalk_DiscoversExpectedStructuralLeaves()
    {
        // End-to-end check that an external walker can enumerate every named structural
        // sub-property of the F86 4-mode object via IInspectable alone. We assert presence
        // of a few specific labels — if the wrappers stop self-naming or their children
        // change, this test pins the regression.
        var block = new CoherenceBlock(N: 5, n: 1, gammaZero: 0.05);
        var eff = FourModeEffective.Build(block);
        var lines = EnumerateLines(eff, depth: 0, maxDepth: 4).ToList();
        foreach (var line in lines) output.WriteLine(line);

        Assert.Contains(lines, l => l.Contains("FourModeEffective (c=2, N=5"));
        Assert.Contains(lines, l => l.Contains("D_eff (diagonal rates)"));
        Assert.Contains(lines, l => l.Contains("M_h_per_bond_eff"));
        Assert.Contains(lines, l => l.Contains("bond 0 (Endpoint)"));
        Assert.Contains(lines, l => l.Contains("bond 1 (Interior)"));
        Assert.Contains(lines, l => l.Contains("bond 3 (Endpoint)"));    // last bond at N=5 → b=3
        Assert.Contains(lines, l => l.Contains("channel-uniform 2×2"));
        Assert.Contains(lines, l => l.Contains("SVD-top 2×2"));
        Assert.Contains(lines, l => l.Contains("cross (CU → SVD-top) 2×2"));
        Assert.Contains(lines, l => l.Contains("probe_eff"));
        Assert.Contains(lines, l => l.Contains("S_kernel_eff"));
    }

    private static IEnumerable<string> EnumerateLines(IInspectable node, int depth, int maxDepth)
    {
        string indent = new(' ', depth * 4);
        string summary = string.IsNullOrEmpty(node.Summary) ? "" : $"  —  {node.Summary}";
        yield return indent + node.DisplayName + summary;
        if (depth >= maxDepth) yield break;
        foreach (var c in node.Children)
            foreach (var l in EnumerateLines(c, depth + 1, maxDepth)) yield return l;
    }
}
