using RCPsiSquared.Core.ChainSystems;
using RCPsiSquared.Core.Pauli;
using RCPsiSquared.Diagnostics.F49;
using RCPsiSquared.Diagnostics.F80;

namespace RCPsiSquared.Diagnostics.Tests.F80;

public class BlochSignWalkTests
{
    [Fact]
    public void TrulyOnly_ReturnsTrivialSpectrum()
    {
        // All terms truly → spectrum is {0 → 4^N}.
        var chain = new ChainSystem(N: 3, J: 1.0, GammaZero: 0.05);
        var terms = new[]
        {
            new PauliPairBondTerm(PauliLetter.X, PauliLetter.X),
            new PauliPairBondTerm(PauliLetter.Y, PauliLetter.Y),
        };
        var spec = BlochSignWalk.PredictMSpectrumImaginaryParts(chain, terms);
        Assert.Single(spec);
        Assert.Equal(64, spec[0.0]); // 4^3 = 64
    }

    [Fact]
    public void Pi2EvenNonTruly_IsRejected()
    {
        // YZ+ZY is Π²-even non-truly — rejected (richer cluster, not in F80 scope).
        var chain = new ChainSystem(N: 3, J: 1.0, GammaZero: 0.05);
        var terms = new[]
        {
            new PauliPairBondTerm(PauliLetter.Y, PauliLetter.Z),
            new PauliPairBondTerm(PauliLetter.Z, PauliLetter.Y),
        };
        Assert.Throws<ArgumentException>(() => BlochSignWalk.PredictMSpectrumImaginaryParts(chain, terms));
    }

    [Fact]
    public void XYAlone_GivesOpenChainSineSpectrumScaled()
    {
        // For pure XY (Π²-odd), F80 says imaginary parts of M = 2·λ_H. The XY chain has
        // OBC sine dispersion 2J·cos(πk/(N+1)). With J=1 here. So M imaginary parts come
        // from H_xy at the input bilinear coupling c.
        // Just check that the spectrum is non-trivial and total multiplicity = 4^N.
        var chain = new ChainSystem(N: 3, J: 1.0, GammaZero: 0.05);
        var terms = new[] { new PauliPairBondTerm(PauliLetter.X, PauliLetter.Y) };
        var spec = BlochSignWalk.PredictMSpectrumImaginaryParts(chain, terms);
        int totalMult = spec.Values.Sum();
        Assert.Equal(64, totalMult); // 4^3
        Assert.True(spec.Count > 1, "non-trivial Π²-odd term should produce multiple eigenvalue clusters");
    }

    [Fact]
    public void NonTrulyIdentityLetter_IsRejected()
    {
        // (I, Y) has nY=1 odd → not truly → falls through to the I-check, which raises.
        // (I, X) is truly (nY=nZ=0 even) and is filtered earlier.
        var chain = new ChainSystem(N: 3, J: 1.0, GammaZero: 0.05);
        var terms = new[] { new PauliPairBondTerm(PauliLetter.I, PauliLetter.Y) };
        Assert.Throws<ArgumentException>(() => BlochSignWalk.PredictMSpectrumImaginaryParts(chain, terms));
    }
}
