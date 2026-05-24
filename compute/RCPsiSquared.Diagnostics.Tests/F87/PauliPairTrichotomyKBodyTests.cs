using System.Numerics;
using RCPsiSquared.Core.ChainSystems;
using RCPsiSquared.Core.Pauli;
using RCPsiSquared.Diagnostics.F87;
using Xunit;

namespace RCPsiSquared.Diagnostics.Tests.F87;

public class PauliPairTrichotomyKBodyTests
{
    private static ChainSystem MakeChainN4() => new(N: 4, J: 1.0, GammaZero: 0.05);

    [Fact]
    public void Classify_K3_TrulyWitness_UnderXDephase_ReturnsTruly()
    {
        // Pair (XXZ, ZZZ) at Klein (0, 1), y_par=0. Per F103 PROOF table, X-dephase
        // in Klein (0,1) is the TRULY cell (count 55 y_par=0, 0 y_par=1).
        var templates = new[]
        {
            new PauliTerm(new[] { PauliLetter.X, PauliLetter.X, PauliLetter.Z }, Complex.One),
            new PauliTerm(new[] { PauliLetter.Z, PauliLetter.Z, PauliLetter.Z }, Complex.One),
        };
        var result = PauliPairTrichotomy.Classify(MakeChainN4(), templates,
            dephaseLetter: PauliLetter.X);
        Assert.Equal(TrichotomyClass.Truly, result);
    }

    [Fact]
    public void Classify_K3_HardWitness_UnderZDephase_ReturnsHard()
    {
        // Same pair (XXZ, ZZZ) at Klein (0, 1), y_par=0. Under Z-dephase Klein (0,1)
        // is the HARD diagonal cell (F103 frozen: Z-deph hard = (42, 8); this pair
        // is one of the 42 y_par=0 hard pairs).
        var templates = new[]
        {
            new PauliTerm(new[] { PauliLetter.X, PauliLetter.X, PauliLetter.Z }, Complex.One),
            new PauliTerm(new[] { PauliLetter.Z, PauliLetter.Z, PauliLetter.Z }, Complex.One),
        };
        var result = PauliPairTrichotomy.Classify(MakeChainN4(), templates,
            dephaseLetter: PauliLetter.Z);
        Assert.Equal(TrichotomyClass.Hard, result);
    }

    [Fact]
    public void Classify_K3_SoftWitness_MotherSector_UnderZDephase_ReturnsSoft()
    {
        // Pair (XYZ, YXZ) at Klein (0, 0), y_par=1. Per F103 PROOF table, mother
        // sector (0,0) soft under any dephase is y_par=1-pure (21 pairs per letter).
        // This pair is one of the 21 Z-dephase mother soft pairs.
        var templates = new[]
        {
            new PauliTerm(new[] { PauliLetter.X, PauliLetter.Y, PauliLetter.Z }, Complex.One),
            new PauliTerm(new[] { PauliLetter.Y, PauliLetter.X, PauliLetter.Z }, Complex.One),
        };
        var result = PauliPairTrichotomy.Classify(MakeChainN4(), templates,
            dephaseLetter: PauliLetter.Z);
        Assert.Equal(TrichotomyClass.Soft, result);
    }
}
