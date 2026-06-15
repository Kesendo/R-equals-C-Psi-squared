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
    [InlineData(3)]
    [InlineData(4)]
    public void U1_parity_is_a_tested_fact_a_transverse_field_breaks_it(int n)
    {
        // The U(1) boundary as a TESTED fact (gate-first): k mod 2 = (m_i − m_j) mod 2 is a Z₂ shadow of
        // the U(1) magnetization. A transverse field h·Σ X_l flips ONE bit (Δk = ±1, an ODD step) and MUST
        // break the parity. This test fires if it does not.
        var w = new DiagonalWitness(n);
        // h = 0: even steps only, parity conserved
        Assert.DoesNotContain(1, w.RungStepsWithField(0.0));
        Assert.True(w.EvenOddBlockMaxWithField(0.0) < 1e-12);
        // h ≠ 0: the odd step appears and the even↔odd block coupling becomes non-zero — the parity breaks
        Assert.Contains(1, w.RungStepsWithField(0.5));
        Assert.True(w.EvenOddBlockMaxWithField(0.5) > 1e-3,
            $"the transverse field failed to break the U(1) parity: even↔odd max {w.EvenOddBlockMaxWithField(0.5)}");
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
