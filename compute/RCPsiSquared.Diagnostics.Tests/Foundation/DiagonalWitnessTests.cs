using RCPsiSquared.Diagnostics.Foundation;

namespace RCPsiSquared.Diagnostics.Tests.Foundation;

public class DiagonalWitnessTests
{
    [Theory]
    [InlineData(2)]
    [InlineData(3)]
    [InlineData(4)]
    public void LH_is_an_even_step_ladder_parity_conserved(int n)
    {
        var w = new DiagonalWitness(n);
        // L_H connects rungs only by even steps {0,2} -- never k+-1
        Assert.All(w.RungStepsOfLH(), s => Assert.True(s == 0 || s == 2, $"unexpected rung step {s}"));
        Assert.DoesNotContain(1, w.RungStepsOfLH());
        // hence the disagreement-count parity is conserved by the full Lindbladian
        Assert.True(w.DisagreementParityConserved());
    }

    [Theory]
    [InlineData(2)]
    [InlineData(3)]
    public void three_diagonals_are_one_orbit_same_spectrum(int n)
    {
        var w = new DiagonalWitness(n);
        Assert.Equal(3, w.OrbitSizeOfQZ());          // {Q_X, Q_Y, Q_Z} one orbit of the basis-S3
        Assert.True(w.ThreeDiagonalsSameSpectrum());  // conjugate => co-spectral
    }

    [Theory]
    [InlineData(2)]
    [InlineData(3)]
    [InlineData(4)]
    public void mirror_group_acts_on_LH_per_PROOF_PI_FACTORS_section3(int n)
    {
        // the §3 Hamiltonian column of the palindrome split: D flips L_H, R fixes L_H
        var w = new DiagonalWitness(n);
        var (dFlip, rFix) = w.MirrorActionOnLH();
        Assert.True(dFlip < 1e-9, $"§3 row 1 D·L_H·D = −L_H failed: dev {dFlip}");
        Assert.True(rFix < 1e-9, $"§3 row 3 R·L_H·R = +L_H failed: dev {rFix}");
    }

    [Fact]
    public void guards_against_too_large_N()
    {
        // 4^6 = 4096 > MaxDim (1024)
        Assert.Throws<ArgumentOutOfRangeException>(() => new DiagonalWitness(6));
    }
}
