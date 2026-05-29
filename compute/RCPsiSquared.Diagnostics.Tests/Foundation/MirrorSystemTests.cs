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

    private static readonly Matrix<Complex> X2 = Matrix<Complex>.Build.DenseOfArray(new Complex[,] { { 0, 1 }, { 1, 0 } });
    private static readonly Matrix<Complex> Y2 = Matrix<Complex>.Build.DenseOfArray(new Complex[,] { { 0, -Complex.ImaginaryOne }, { Complex.ImaginaryOne, 0 } });
    private static readonly Matrix<Complex> Z2 = Matrix<Complex>.Build.DenseOfArray(new Complex[,] { { 1, 0 }, { 0, -1 } });
    private static readonly Matrix<Complex> I2 = Matrix<Complex>.Build.DenseIdentity(2);

    // Place A on site l and B on site l+1, identity elsewhere; one chain bond.
    private static Matrix<Complex> Bond(int n, int l, Matrix<Complex> A, Matrix<Complex> B)
    {
        Matrix<Complex>? r = null;
        for (int s = 0; s < n; s++)
        {
            var op = s == l ? A : (s == l + 1 ? B : I2);
            r = r is null ? op : r.KroneckerProduct(op);
        }
        return r!;
    }

    // Truly: H = J * sum_l (X_lX_{l+1} + Y_lY_{l+1} + Z_lZ_{l+1}) on an open chain.
    private static Matrix<Complex> HeisenbergChain(int n, double J)
    {
        var H = Matrix<Complex>.Build.Dense(1 << n, 1 << n);
        for (int l = 0; l < n - 1; l++)
            H += (Bond(n, l, X2, X2) + Bond(n, l, Y2, Y2) + Bond(n, l, Z2, Z2)).Multiply((Complex)J);
        return H;
    }

    // Non-truly chain Pi2-odd: H = sum_l X_l Y_{l+1} (the F80 case).
    private static Matrix<Complex> ChainXYNonTruly(int n)
    {
        var H = Matrix<Complex>.Build.Dense(1 << n, 1 << n);
        for (int l = 0; l < n - 1; l++)
            H += Bond(n, l, X2, Y2);
        return H;
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

    [Fact]
    public void Evolve_UnrollsTheStaticSpectrumIntoTime()
    {
        var sys = new MirrorSystem(2, Heisenberg2(1.0),
            new[] { new ChannelRate("a", 0.05), new ChannelRate("b", 0.20) });

        // K = 0: nothing has decayed yet, every mode is fully present (survival = 1).
        Assert.All(sys.Evolve(0.0), s => Assert.Equal(1.0, s.Survival, precision: 9));

        // Late time: the slow (memory) mode outlives the fast one, and the fastest is essentially gone.
        var late = sys.Evolve(5.0);
        var slowestNonzero = late.Where(s => s.Rate > 1e-9).OrderBy(s => s.Rate).First();
        var fastest = late.OrderByDescending(s => s.Rate).First();
        Assert.True(slowestNonzero.Survival > fastest.Survival,
            "the slow memory mode should outlive the fast mode at late time");
        Assert.True(fastest.Survival < 0.01,
            $"the fastest mode should be essentially gone by K=5, got {fastest.Survival:E2}");

        // survival is exactly e^(-(rate/sigma)*K): the spectrum unrolled, no hidden state.
        double sigma = sys.TotalDephasing;
        foreach (var s in late)
            Assert.Equal(Math.Exp(-(s.Rate / sigma) * 5.0), s.Survival, precision: 9);
    }

    [Fact]
    public void MemoryRotation_TrulyHamiltonian_IsPerfectMirror()
    {
        // non-uniform carrier, to exercise the non-uniform F1 palindrome too.
        var sys = new MirrorSystem(3, HeisenbergChain(3, 1.0),
            new[] { new ChannelRate("a", 0.05), new ChannelRate("b", 0.10), new ChannelRate("c", 0.20) });
        var mr = sys.MemoryRotation;

        // Truly H + Z-dephasing: the F1 palindrome holds exactly, M = 0, a perfect mirror.
        Assert.True(mr.PerfectMirror, $"truly H should give M=0, got defect {mr.DefectNorm:E3}");
        Assert.True(mr.DefectNorm < 1e-9);
        // M=0 carries no spectrum: a mirror-symmetric wave needs no memory-defect.
        Assert.False(mr.MemoryCarriesEnergy);
    }

    [Fact]
    public void MemoryRotation_NonTrulyChainBilinear_CarriesEnergyRotated90()
    {
        var H = ChainXYNonTruly(3);  // sum_l X_l Y_{l+1}, the F80 non-truly chain Pi2-odd case
        var sys = new MirrorSystem(3, H,
            new[] { new ChannelRate("a", 0.05), new ChannelRate("b", 0.10), new ChannelRate("c", 0.20) });
        var mr = sys.MemoryRotation;

        // The palindrome breaks: a real mirror-defect M, and it is not noise , it carries H.
        Assert.False(mr.PerfectMirror);
        Assert.True(mr.MemoryCarriesEnergy, "F80: non-truly chain Pi2-odd H gives Spec(M) = +-2i*Spec(H)");

        // F80 Frobenius relation: ||M||^2_F = 4 * ||H||^2_F * 2^N  (M = -2i * H (x) I_bra).
        double hNormSq = H.FrobeniusNorm() * H.FrobeniusNorm();
        Assert.Equal(4.0 * hNormSq * (1 << 3), mr.DefectNorm * mr.DefectNorm, precision: 4);

        // The quarter turn itself: every energy lambda maps to 2*lambda on the imaginary memory axis.
        Assert.NotEmpty(mr.Rotation);
        Assert.All(mr.Rotation, p => Assert.Equal(2.0 * p.Energy, p.MemoryAxisValue, precision: 9));
    }

    [Fact]
    public void Takt_AnchorsToTauMax_GapIsTwoGamma_TauIsTen()
    {
        // N=3 Heisenberg chain, uniform gamma = 0.05, J = 1: TAU_MAX's verified clock.
        var sys = new MirrorSystem(3, HeisenbergChain(3, 1.0),
            new[] { new ChannelRate("a", 0.05), new ChannelRate("b", 0.05), new ChannelRate("c", 0.05) });

        // The slowest felt motion is the 2*gamma floor; the longest breath is 1/(2*gamma) = 10.
        Assert.Equal(0.1, sys.Takt.Gap, precision: 6);
        Assert.Equal(10.0, sys.Takt.Tau, precision: 6);
        Assert.False(sys.Takt.Stopped);
    }

    [Fact]
    public void Takt_ClockStops_WhenAllChannelsAreZero()
    {
        // gamma = 0 everywhere: no decay clock, only frozen Hamiltonian oscillation.
        var sys = new MirrorSystem(3, HeisenbergChain(3, 1.0),
            new[] { new ChannelRate("a", 0.0), new ChannelRate("b", 0.0), new ChannelRate("c", 0.0) });

        Assert.True(sys.Takt.Stopped);
        Assert.Equal(0.0, sys.Takt.Gap);
        Assert.True(double.IsPositiveInfinity(sys.Takt.Tau));
    }

    [Fact]
    public void Takt_Gap_DelegatesToSpectrumSlowestRate()
    {
        // The voice does not re-derive the floor: Takt.Gap IS the spectrum's SlowestRate,
        // the framework's shared memory floor (single source of truth, non-uniform carrier, N=3).
        var sys = new MirrorSystem(3, HeisenbergChain(3, 1.0),
            new[] { new ChannelRate("a", 0.05), new ChannelRate("b", 0.10), new ChannelRate("c", 0.20) });

        Assert.Equal(sys.Spectrum.SlowestRate, sys.Takt.Gap, precision: 12);
    }

    [Fact]
    public void Takt_GammaScalesTheClock_HalvingGammaDoublesTau()
    {
        // Halving the uniform carrier doubles the slowest felt duration (Tau scales as 1/(2*gamma)).
        var fast = new MirrorSystem(3, HeisenbergChain(3, 1.0),
            new[] { new ChannelRate("a", 0.10), new ChannelRate("b", 0.10), new ChannelRate("c", 0.10) });
        var slow = fast.WithChannels(
            new[] { new ChannelRate("a", 0.05), new ChannelRate("b", 0.05), new ChannelRate("c", 0.05) });

        Assert.Equal(2.0 * fast.Takt.Tau, slow.Takt.Tau, precision: 6);
    }

    [Fact]
    public void Rotation_CapturesImLambda_TheChainOscillates()
    {
        // N=3 Heisenberg, uniform gamma = 0.05, J = 1: the chain rotates, so at least one mode
        // carries a nonzero oscillation frequency omega = Im(lambda) in its captured data.
        var sys = new MirrorSystem(3, HeisenbergChain(3, 1.0),
            new[] { new ChannelRate("a", 0.05), new ChannelRate("b", 0.05), new ChannelRate("c", 0.05) });

        Assert.Contains(sys.Spectrum.Modes, m => Math.Abs(m.OscillationFrequency) > 1e-9);
    }

    [Fact]
    public void Rotation_NoTurning_WhenJIsZero()
    {
        // J = 0: no Hamiltonian rotation, only pure radial decay. The angular hand does not turn.
        var sys = new MirrorSystem(3, HeisenbergChain(3, 0.0),
            new[] { new ChannelRate("a", 0.05), new ChannelRate("b", 0.05), new ChannelRate("c", 0.05) });

        Assert.False(sys.Rotation.Turning);
        Assert.Equal(0.0, sys.Rotation.Angle);
    }

    [Fact]
    public void Rotation_Angle_IsArctanOfOmegaOverGap()
    {
        // N=3 Heisenberg, gamma = 0.05, J = 1: the memory mode's angle theta = arctan(omega / Gap),
        // the F95 angle (= arctan Q for the 2-level case), wired off Takt.Gap and the captured omega.
        var sys = new MirrorSystem(3, HeisenbergChain(3, 1.0),
            new[] { new ChannelRate("a", 0.05), new ChannelRate("b", 0.05), new ChannelRate("c", 0.05) });

        Assert.Equal(Math.Atan2(sys.Rotation.Frequency, sys.Takt.Gap), sys.Rotation.Angle, precision: 9);
    }

    [Fact]
    public void Rotation_AngleGrowsWithJ_MoreRotationPerUnitDecay()
    {
        // At fixed gamma = 0.05, raising J turns the angular hand further: more rotation per unit decay.
        var weak = new MirrorSystem(3, HeisenbergChain(3, 0.5),
            new[] { new ChannelRate("a", 0.05), new ChannelRate("b", 0.05), new ChannelRate("c", 0.05) });
        var strong = new MirrorSystem(3, HeisenbergChain(3, 2.0),
            new[] { new ChannelRate("a", 0.05), new ChannelRate("b", 0.05), new ChannelRate("c", 0.05) });

        Assert.True(strong.Rotation.Angle > weak.Rotation.Angle,
            $"theta should grow with J: J=2 gave {strong.Rotation.Angle:F6}, J=0.5 gave {weak.Rotation.Angle:F6}");
    }

    [Fact]
    public void Rotation_PureCircleLimit_GammaZero_AngleIsNinetyDegrees()
    {
        // gamma = 0 everywhere, J = 1, N = 3: the radial hand stops (Takt.Stopped) and the angular
        // hand turns forever at the pure circle, theta = pi/2 (no inward pull, only rotation).
        var sys = new MirrorSystem(3, HeisenbergChain(3, 1.0),
            new[] { new ChannelRate("a", 0.0), new ChannelRate("b", 0.0), new ChannelRate("c", 0.0) });

        Assert.True(sys.Takt.Stopped);
        Assert.True(sys.Rotation.Turning);
        Assert.Equal(Math.PI / 2.0, sys.Rotation.Angle, tolerance: 1e-9);
    }
}
