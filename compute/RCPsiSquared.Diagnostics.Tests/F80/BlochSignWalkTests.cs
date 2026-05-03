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
        // YZ+ZY is Π²-even non-truly: rejected (richer cluster, not in F80 scope).
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
        // For pure XY (Π²-odd), F80 says imaginary parts of M = ±2·λ_H, with multiplicity ×2^N
        // per H eigenvalue. Diagnostic: print the cluster contents on failure so the test
        // documents what the implementation actually produces.
        var chain = new ChainSystem(N: 3, J: 1.0, GammaZero: 0.05);
        var terms = new[] { new PauliPairBondTerm(PauliLetter.X, PauliLetter.Y) };
        var spec = BlochSignWalk.PredictMSpectrumImaginaryParts(chain, terms);

        string dump = string.Join(", ", spec.OrderBy(kv => kv.Key).Select(kv => $"{kv.Key:G6}→{kv.Value}"));

        int totalMult = spec.Values.Sum();
        Assert.True(totalMult == 64, $"expected total mult 4^N=64, got {totalMult}; spec = [{dump}]");

        // Spectrum must be ±-symmetric: every key k has its mirror −k with equal multiplicity.
        foreach (var (key, mult) in spec)
        {
            double mirror = Math.Round(-key, 10);
            Assert.True(spec.TryGetValue(mirror, out var mirrorMult) && mirrorMult == mult,
                $"asymmetric cluster: {key:G6}→{mult} but {mirror:G6}→{(spec.TryGetValue(mirror, out var m) ? m.ToString() : "missing")}; spec = [{dump}]");
        }

        // Multiple non-zero clusters expected (Π²-odd → non-trivial spectrum).
        Assert.True(spec.Count > 1 || (spec.Count == 1 && !spec.ContainsKey(0.0)),
            $"trivial spectrum for non-trivial Π²-odd term; spec = [{dump}]");
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
