using MirrorWorld;

namespace MirrorWorldTests;

// From-below guard for the within-block self-dual seed (adopted 2026-07-07 from
// experiments/F89_SEED_EXISTENCE_REDUCTION.md + docs/proofs/PROOF_CODIM1_BY_ADDITIVITY.md). The seed
// count is the nullity surplus r(0+) - r(inf) of the (1,2) pencil A + q*C, taken over GF(p) with no
// eigensolver. Pinned against F89's proven value: N-1 at odd N, and 0 at even N (the unmirrorable-seat
// asymmetry the census reads and this holds algebraically).
public class SeedTests
{
    const double G = 0.5;
    static readonly World W = new();

    // the proven count: a real defective seed forced on the (1,2) block at every odd N, N-1 of them.
    [Theory]
    [InlineData(3, 2)]
    [InlineData(5, 4)]
    [InlineData(7, 6)]
    [InlineData(9, 8)]
    public void Surplus_Is_N_Minus_1_At_Odd_N(int n, int expected)
    {
        var (_, _, surplus, _) = new Seed(W, n, G).Count();
        Assert.Equal(expected, surplus);
        Assert.Equal(n - 1, surplus);
    }

    // the mirror-image half: at even N there is no unmirrorable seat, every seat mirrors to another, so
    // nothing is forced real and the surplus is 0 (the retraction's "no real-axis seed" is right HERE).
    [Theory]
    [InlineData(4)]
    [InlineData(6)]
    [InlineData(8)]
    public void Surplus_Is_Zero_At_Even_N(int n)
    {
        var (_, _, surplus, _) = new Seed(W, n, G).Count();
        Assert.Equal(0, surplus);
    }

    // a sharper pin than the surplus alone (which a compensating error in both halves could fake): the
    // two halves and the per-rung nullities at N=5, matching F89's recorded reduction.
    [Fact]
    public void Halves_And_Rungs_Match_F89_At_N5()
    {
        var (rInf, r0, surplus, parts) = new Seed(W, 5, G).Count();
        Assert.Equal(6, rInf);
        Assert.Equal(10, r0);
        Assert.Equal(4, surplus);
        Assert.Contains(parts, p => p.NDiff == 3 && p.Dim == 30 && p.Nullity == 6);   // rung A = -6
        Assert.Contains(parts, p => p.NDiff == 1 && p.Dim == 20 && p.Nullity == 4);   // rung A = -2
    }

    // the fusion-resonance count, closed (adopted 2026-07-12): r(inf) = 3*Z3 at EVERY N (the cyclotomic
    // Step-4 theorem; even N included), and at odd N Z3 matches the Conway-Jones closed form. N = 8 is
    // the smallest even N with resonances (Z3 = 2, the 2cos(pi/9) family), yet seedless.
    [Theory]
    [InlineData(3, 1)]
    [InlineData(4, 0)]
    [InlineData(5, 2)]
    [InlineData(6, 0)]
    [InlineData(7, 3)]
    [InlineData(8, 2)]
    [InlineData(9, 4)]
    public void RInf_Is_3_Z3_And_Odd_N_Matches_ConwayJones(int n, int z3)
    {
        var seed = new Seed(W, n, G);
        var (rInf, _, _, _) = seed.Count();
        Assert.Equal(0, rInf % 3);
        Assert.Equal(z3, Seed.Z3FromRInf(rInf));
        if (n % 2 == 1) Assert.Equal(z3, Seed.Z3ClosedFormOdd(n));
    }

    // the criterion and the closed form past the rank range (pure arithmetic, F89's 10c pins):
    // resonant (odd N) iff 3 | N+1 and N >= 11; the next resonant N after 17 is 23, NOT 29; N = 29 is the
    // first with the PENT family (+2).
    [Fact]
    public void Resonance_Criterion_And_ClosedForm_Pins()
    {
        Assert.Equal(7, Seed.Z3ClosedFormOdd(11));    // (11-1)/2 + (12/3 - 2) = 5 + 2
        Assert.Equal(6, Seed.Z3ClosedFormOdd(13));    // (13-1)/2, 14 has no family
        Assert.Equal(24, Seed.Z3ClosedFormOdd(29));   // 14 + (10 - 2) + 2: the first PENT N
        Assert.False(Seed.IsResonant(5));             // 3 | 6 but N < 11 (there the CJ multiplicity n/3 - 2 = 0)
        Assert.False(Seed.IsResonant(8));             // even N: resonances exist, no seeds, criterion odd-N only
        Assert.False(Seed.IsResonant(14));            // even N with 3 | N+1: the parity guard, not the mod-3 test
        Assert.True(Seed.IsResonant(11));
        Assert.False(Seed.IsResonant(13));
        Assert.True(Seed.IsResonant(17));
        Assert.False(Seed.IsResonant(19));
        Assert.True(Seed.IsResonant(23));
        Assert.True(Seed.IsResonant(29));
    }

    // the count is gamma-independent (gamma scales A, not its kernel structure; the pencil's real-strand
    // count is a property of the hop C and the rungs, not the rate magnitude).
    [Fact]
    public void Count_Is_Independent_Of_Gamma()
    {
        var a = new Seed(W, 7, 0.05).Count();
        var b = new Seed(W, 7, 1.7).Count();
        Assert.Equal(a.Surplus, b.Surplus);
        Assert.Equal(6, a.Surplus);
    }
}
