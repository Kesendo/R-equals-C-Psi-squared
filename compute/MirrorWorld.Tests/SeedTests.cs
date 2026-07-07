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
