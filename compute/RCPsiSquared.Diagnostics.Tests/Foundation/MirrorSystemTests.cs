using System;
using System.Linq;
using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Diagnostics.Foundation;

namespace RCPsiSquared.Diagnostics.Tests.Foundation;

public class MirrorSystemTests
{
    private static Matrix<Complex> Heisenberg2(double J)
    {
        var X = Matrix<Complex>.Build.DenseOfArray(new Complex[,] { { 0, 1 }, { 1, 0 } });
        var Y = Matrix<Complex>.Build.DenseOfArray(new Complex[,] { { 0, -Complex.ImaginaryOne }, { Complex.ImaginaryOne, 0 } });
        var Z = Matrix<Complex>.Build.DenseOfArray(new Complex[,] { { 1, 0 }, { 0, -1 } });
        return (X.KroneckerProduct(X) + Y.KroneckerProduct(Y) + Z.KroneckerProduct(Z)).Multiply((Complex)J);
    }

    [Fact]
    public void MirrorSystem_HoldsSpectrumAndPalindromeAsLiveProperties()
    {
        var channels = new[] { new ChannelRate("a", 0.05), new ChannelRate("b", 0.20) };
        var sys = new MirrorSystem(2, Heisenberg2(1.0), channels);

        // σ = Σγ = 0.25, the palindrome centre.
        Assert.Equal(0.25, sys.TotalDephasing, precision: 9);

        // Voice 1 (inner law): all 4^2 modes present as portfolios.
        Assert.Equal(16, sys.Spectrum.Modes.Count);

        // Voice 2 (F1, live): every decay rate r has its mirror partner 2σ − r in the spectrum.
        Assert.True(sys.PalindromeHolds);
        // and the pairing is the one the palindrome dictates: r + partner = 2σ for each mode.
        foreach (var p in sys.PalindromePartners)
            Assert.Equal(2.0 * sys.TotalDephasing, p.Rate + p.PartnerRate, precision: 9);
    }

    [Fact]
    public void MovingTheCarrier_MovesTheWholeBox()
    {
        var sys = new MirrorSystem(2, Heisenberg2(1.0),
            new[] { new ChannelRate("a", 0.05), new ChannelRate("b", 0.20) });
        var moved = sys.WithChannels(new[] { new ChannelRate("a", 0.50), new ChannelRate("b", 0.50) });

        // The dephasing centre moved (0.25 -> 1.0)...
        Assert.Equal(0.25, sys.TotalDephasing, precision: 9);
        Assert.Equal(1.00, moved.TotalDephasing, precision: 9);

        // ...and so did the spectrum: the readings recompute from the moved input.
        var before = sys.Spectrum.Modes.Select(m => m.ActualDecayRate).OrderBy(x => x).ToArray();
        var after = moved.Spectrum.Modes.Select(m => m.ActualDecayRate).OrderBy(x => x).ToArray();
        double maxShift = before.Zip(after, (x, y) => Math.Abs(x - y)).Max();
        Assert.True(maxShift > 1e-6, $"the box did not move: max rate shift {maxShift:E3}");

        // F1 still holds for the moved system: the palindrome is a property of the box, not the input.
        Assert.True(moved.PalindromeHolds);
    }
}
