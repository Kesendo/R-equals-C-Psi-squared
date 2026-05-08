using RCPsiSquared.Core.CoherenceBlocks;
using RCPsiSquared.Core.F86.Item1Derivation;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;
using Xunit.Abstractions;

namespace RCPsiSquared.Core.Tests.F86.Item1Derivation;

public class C2FullBlockEigenAnatomyTests
{
    private readonly ITestOutputHelper _out;

    public C2FullBlockEigenAnatomyTests(ITestOutputHelper output) => _out = output;

    private static CoherenceBlock C2Block(int N) => new(N: N, n: 1, gammaZero: 0.05);

    [Theory]
    [InlineData(5)]
    [InlineData(6)]
    [InlineData(7)]
    [InlineData(8)]
    public void Build_AcrossN5To8_LandsTier2Verified(int N)
    {
        var anatomy = C2FullBlockEigenAnatomy.Build(C2Block(N));
        Assert.Equal(Tier.Tier2Verified, anatomy.Tier);
    }

    [Theory]
    [InlineData(5, 50)]   // dim = 5 · 10 = 50
    [InlineData(6, 90)]   // dim = 6 · 15 = 90
    [InlineData(7, 147)]  // dim = 7 · 21 = 147
    [InlineData(8, 224)]  // dim = 8 · 28 = 224
    public void EigenSpectrum_HasFullBlockLDimension(int N, int expectedDim)
    {
        var anatomy = C2FullBlockEigenAnatomy.Build(C2Block(N));
        Assert.Equal(expectedDim, anatomy.EigenSpectrum.Count);
        Assert.Equal(expectedDim, anatomy.Block.Basis.MTotal);
    }

    [Theory]
    [InlineData(5)]
    [InlineData(6)]
    [InlineData(7)]
    [InlineData(8)]
    public void K90_IsFinite_AndAtLeastOne(int N)
    {
        var anatomy = C2FullBlockEigenAnatomy.Build(C2Block(N));
        Assert.True(anatomy.K90 >= 1, $"N={N}: K_90 must be ≥ 1 (sanity)");
        Assert.True(anatomy.K90 <= anatomy.EigenSpectrum.Count,
            $"N={N}: K_90 ≤ block dim");
        Assert.True(anatomy.K99 >= anatomy.K90, $"N={N}: K_99 ≥ K_90 (cumulative)");
    }

    [Fact]
    public void DirectionBPrimePrimeReconnaissance_AcrossN5To8_EmitTable()
    {
        // The decisive question for Direction (b''): does K_90 stay roughly constant
        // (4-8) as N grows, or does it scale with the block dim? Emit the per-N table
        // for human inspection in CI logs, plus the top-3 eigenmodes per N to reveal
        // structure.
        _out.WriteLine("Direction (b'') reconnaissance: full block-L eigenanatomy at Q = Q_EP");
        _out.WriteLine("");
        _out.WriteLine("  N | dim | K_90 | K_99 | top-mode Re(λ) | top-mode |c|² | top-mode weight");
        _out.WriteLine("  --|-----|------|------|----------------|--------------|----------------");
        foreach (int N in new[] { 5, 6, 7, 8 })
        {
            var anatomy = C2FullBlockEigenAnatomy.Build(C2Block(N));
            var top = anatomy.EigenSpectrum[0];
            _out.WriteLine(
                $"  {N} | {anatomy.Block.Basis.MTotal,3} | {anatomy.K90,4} | {anatomy.K99,4} | " +
                $"{top.EigenvalueReal,+13:F4}  | {top.ProbeOverlapSquared,11:G4}   | {top.DiagonalWeight,11:G4}");
        }
        _out.WriteLine("");
        _out.WriteLine("Reading: K_90 small + N-stable ⟹ HWHM lives in a low-rank truncation,");
        _out.WriteLine("Tier1Candidate path with explicit error term is realistic.");
        _out.WriteLine("K_90 grows with N ⟹ no clean closed form via this anatomy; alternative");
        _out.WriteLine("structural reduction (XY-mode decomposition, JW free-fermion picture) needed.");
    }

    [Theory]
    [InlineData(5)]
    [InlineData(7)]
    public void TopMode_HasNegativeRealPart_AtQEp(int N)
    {
        // At Q = Q_EP the slowest-pair eigenvalue Re ≈ −4γ₀ (PROOF_F86_QPEAK Statement 1).
        // The top-weight eigenmode at this Q should sit in the slow-decay regime
        // (Re(λ) close to or above −10γ₀).
        var anatomy = C2FullBlockEigenAnatomy.Build(C2Block(N));
        var top = anatomy.EigenSpectrum[0];
        Assert.True(top.EigenvalueReal < 0, $"N={N}: top eigenmode must have Re(λ) < 0 (decay)");
        Assert.True(top.EigenvalueReal > -10.0 * anatomy.Block.GammaZero,
            $"N={N}: top eigenmode Re(λ) = {top.EigenvalueReal:F4} should be above −10γ₀ = {-10.0 * anatomy.Block.GammaZero:F4}");
    }

    [Fact]
    public void Build_RejectsNonC2Block()
    {
        var block = new CoherenceBlock(N: 5, n: 2, gammaZero: 0.05);
        Assert.Equal(3, block.C);
        Assert.Throws<ArgumentException>(() => C2FullBlockEigenAnatomy.Build(block));
    }

    [Fact]
    public void Build_AcceptsCustomQ_DifferentFromQEp()
    {
        var anatomy = C2FullBlockEigenAnatomy.Build(C2Block(5), q: 0.5);
        Assert.Equal(0.5, anatomy.Q);
    }

    [Theory]
    [InlineData(5)]
    [InlineData(6)]
    public void Children_IncludeK90AndTopModesGroup(int N)
    {
        IInspectable claim = C2FullBlockEigenAnatomy.Build(C2Block(N));
        var labels = claim.Children.Select(c => c.DisplayName).ToList();
        Assert.Contains("K_90", labels);
        Assert.Contains("K_99", labels);
        Assert.Contains(labels, l => l.StartsWith("top ") && l.Contains("eigenmodes"));
    }
}
