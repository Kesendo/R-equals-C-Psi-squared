using RCPsiSquared.Core.ChainSystems;
using RCPsiSquared.Core.Pauli;
using RCPsiSquared.Diagnostics.Foundation;

namespace RCPsiSquared.Diagnostics.Tests.Foundation;

public class MemoryAxisChainTests
{
    [Fact]
    public void Compute_WithoutHamiltonian_ReportsBalancedPi2Partition()
    {
        var chain = new ChainSystem(N: 3, J: 1.0, GammaZero: 0.05);
        var result = MemoryAxisChain.Compute(chain);

        // F88 bilinear apex 1/2: 4^3 / 2 = 32 strings on each side
        Assert.Equal(32L, result.Pi2EvenStringCount);
        Assert.Equal(32L, result.Pi2OddStringCount);
        Assert.True(result.BilinearApexHolds);
        Assert.Null(result.F80ImaginarySpectrum);
        Assert.True(result.F80SpectrumIsMirrorSymmetric); // vacuously true with no H
    }

    [Theory]
    [InlineData(2, 8L, 8L)]      // 4^2 = 16, half = 8
    [InlineData(3, 32L, 32L)]    // 4^3 = 64, half = 32
    [InlineData(4, 128L, 128L)]  // 4^4 = 256, half = 128
    public void BilinearApex_HoldsAtAllTestedN(int N, long expectedEven, long expectedOdd)
    {
        var chain = new ChainSystem(N, J: 1.0, GammaZero: 0.05);
        var result = MemoryAxisChain.Compute(chain);

        Assert.Equal(expectedEven, result.Pi2EvenStringCount);
        Assert.Equal(expectedOdd, result.Pi2OddStringCount);
        Assert.True(result.BilinearApexHolds);
    }

    [Fact]
    public void Compute_WithPi2OddH_ReportsMirrorSymmetricF80Spectrum()
    {
        // XY is chain Π²-odd: F80 lift gives a non-trivial ±-symmetric imaginary spectrum
        var chain = new ChainSystem(N: 3, J: 1.0, GammaZero: 0.05);
        var terms = new[] { new PauliPairBondTerm(PauliLetter.X, PauliLetter.Y) };
        var result = MemoryAxisChain.Compute(chain, terms);

        Assert.NotNull(result.F80ImaginarySpectrum);
        Assert.True(result.F80SpectrumIsMirrorSymmetric);

        // Total multiplicity sums to 4^N = 64
        Assert.Equal(64, result.F80ImaginarySpectrum!.Values.Sum());

        // Non-trivial: more than one cluster, OR a single non-zero cluster
        bool nontrivial = result.F80ImaginarySpectrum.Count > 1
                       || (result.F80ImaginarySpectrum.Count == 1 && !result.F80ImaginarySpectrum.ContainsKey(0.0));
        Assert.True(nontrivial);
    }

    [Fact]
    public void Compute_WithTrulyOnlyH_ReportsTrivialF80Spectrum()
    {
        // All-truly H has no F80-non-truly content: spectrum is {0 → 4^N}
        var chain = new ChainSystem(N: 3, J: 1.0, GammaZero: 0.05);
        var terms = new[]
        {
            new PauliPairBondTerm(PauliLetter.X, PauliLetter.X),
            new PauliPairBondTerm(PauliLetter.Y, PauliLetter.Y),
            new PauliPairBondTerm(PauliLetter.Z, PauliLetter.Z),
        };
        var result = MemoryAxisChain.Compute(chain, terms);

        Assert.NotNull(result.F80ImaginarySpectrum);
        Assert.Single(result.F80ImaginarySpectrum!);
        Assert.Equal(64, result.F80ImaginarySpectrum[0.0]);
        Assert.True(result.F80SpectrumIsMirrorSymmetric); // 0 = -0
    }

    [Fact]
    public void StructuralPair_IsPlusMinusHalf()
    {
        // The static memory-axis anchor: same pair at every trio layer
        Assert.Equal(2, MemoryAxisChainResult.StructuralPair.Count);
        Assert.Equal(+0.5, MemoryAxisChainResult.StructuralPair[0]);
        Assert.Equal(-0.5, MemoryAxisChainResult.StructuralPair[1]);
        Assert.Equal(0.0, MemoryAxisChainResult.StructuralPair.Sum());
    }

    [Fact]
    public void Compute_DephaseLetterX_StillBalances()
    {
        // F88 holds for any dephase letter (the 4^N/2 split is a combinatorial property
        // of bit_a or bit_b parity, both balance to 4^N/2)
        var chain = new ChainSystem(N: 3, J: 1.0, GammaZero: 0.05);
        var result = MemoryAxisChain.Compute(chain, hamiltonianTerms: null, dephaseLetter: PauliLetter.X);

        Assert.Equal(32L, result.Pi2EvenStringCount);
        Assert.Equal(32L, result.Pi2OddStringCount);
        Assert.True(result.BilinearApexHolds);
    }
}
