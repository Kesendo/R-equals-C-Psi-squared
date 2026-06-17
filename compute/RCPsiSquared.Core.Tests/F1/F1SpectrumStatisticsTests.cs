using System.Numerics;
using RCPsiSquared.Core.F1;
using Xunit;

namespace RCPsiSquared.Core.Tests.F1;

/// <summary>The F1 symmetry distance (<see cref="F1SpectrumStatistics.MaxF1PairingDistance"/>): it
/// must test MULTISET-closure under λ ↦ −2σ − λ, not set-closure. A set/Hausdorff distance is blind
/// to a dropped/duplicated eigenvalue that still has a same-valued neighbour — the exact failure a
/// spectrum-reconstruction bug produces. These pin the multiplicity-aware behaviour at its Core home.</summary>
public class F1SpectrumStatisticsTests
{
    [Fact]
    public void MaxF1PairingDistance_IsZero_ForAnF1SymmetricMultiset()
    {
        // {0, -2, -1+i, -1-i} is closed under λ ↦ -2σ - λ with σ=1 (center -1).
        var spectrum = new[]
        {
            new Complex(0, 0), new Complex(-2, 0),
            new Complex(-1, 1), new Complex(-1, -1),
        };
        Assert.True(F1SpectrumStatistics.MaxF1PairingDistance(spectrum, sigma: 1.0) < 1e-12);
    }

    [Fact]
    public void MaxF1PairingDistance_DetectsAMultiplicityDefect()
    {
        // {+1 x3, -1 x1} is closed under λ ↦ -λ (σ=0) as a SET but NOT as a multiset
        // (reflected = {-1 x3, +1 x1}). A set distance reports ~0; the greedy NN with removal must
        // see the surplus +1 forced onto a -1 at distance 2.
        var defective = new[]
        {
            new Complex(1, 0), new Complex(1, 0), new Complex(1, 0), new Complex(-1, 0),
        };
        Assert.Equal(2.0, F1SpectrumStatistics.MaxF1PairingDistance(defective, sigma: 0.0), 9);
    }
}
