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

    // ---- the coupled-level law (adopted 2026-07-12 from the 10d section of the same doc) ----

    // the soundness tie: the cyclotomic triple enumeration can only OVERcount (a true zero is zero mod
    // every prime), and the proved r(inf) = 3*Z3 makes the independent GF(p) rank count the exact
    // ceiling -- equality pins the enumeration at every N in the rank range, even N included.
    [Theory]
    [InlineData(3)]
    [InlineData(4)]
    [InlineData(5)]
    [InlineData(6)]
    [InlineData(7)]
    [InlineData(8)]
    [InlineData(9)]
    [InlineData(10)]
    [InlineData(11)]
    public void Triple_Count_Ties_To_RInf_Over_3(int n)
    {
        var seed = new Seed(W, n, G);
        Assert.Equal(Seed.Z3FromRInf(seed.Count().RInf), seed.VanishingTriples().Length);
    }

    // the odd-N family decomposition behind the Conway-Jones closed form: (N-1)/2 TRIV + [3|n](n/3 - 2)
    // ROT3 + 2[15|n] PENT (the even-n counting; odd n counts differently -- see the N=8 pin below).
    [Theory]
    [InlineData(5)]
    [InlineData(7)]
    [InlineData(9)]
    [InlineData(11)]
    [InlineData(13)]
    [InlineData(17)]
    [InlineData(23)]
    [InlineData(29)]
    public void Odd_N_Family_Counts_Match_ConwayJones_Decomposition(int nSites)
    {
        int n = nSites + 1;
        var triples = new Seed(W, nSites, G).VanishingTriples();
        Assert.Equal((nSites - 1) / 2, triples.Count(t => t.Family == TripleFamily.Triv));
        Assert.Equal(n % 3 == 0 ? n / 3 - 2 : 0, triples.Count(t => t.Family == TripleFamily.Rot3));
        Assert.Equal(n % 15 == 0 ? 2 : 0, triples.Count(t => t.Family == TripleFamily.Pent));
        Assert.Equal(Seed.Z3ClosedFormOdd(nSites), triples.Length);
    }

    // the smallest even laboratory (run 2026-07-12): at N = 8 (n = 9) the vanishing set is exactly the
    // one mode-disjoint mirror pair {2,4,8} / {1,5,7}, both ROT3 (cos(2pi/9) + cos(4pi/9) = cos(pi/9)),
    // no TRIV floor (n odd), pair levels (18/9, 6/9) = (2, 2/3) -- a*n = 6 extends to even N.
    [Fact]
    public void Even_N_Inventory_At_N8()
    {
        var triples = new Seed(W, 8, G).VanishingTriples();
        Assert.Equal(2, triples.Length);
        Assert.Contains((1, 5, 7, TripleFamily.Rot3), triples);
        Assert.Contains((2, 4, 8, TripleFamily.Rot3), triples);
    }

    // mirror closure: lam_{n-k} = -lam_k, so the mirror k -> n-k of a vanishing triple vanishes too.
    // TRIV is exactly the self-mirror class (the class-split theorem); every non-TRIV triple pairs
    // with a DISTINCT mirror partner in the inventory.
    [Theory]
    [InlineData(8)]
    [InlineData(11)]
    [InlineData(14)]
    [InlineData(17)]
    public void NonTriv_Triples_Pair_Under_Mirror_Triv_Is_SelfMirror(int nSites)
    {
        int n = nSites + 1;
        var triples = new Seed(W, nSites, G).VanishingTriples();
        var set = triples.Select(t => (t.K1, t.K2, t.K3)).ToHashSet();
        foreach (var t in triples)
        {
            var mirror = (n - t.K3, n - t.K2, n - t.K1);
            Assert.Contains(mirror, set);
            bool selfMirror = mirror == (t.K1, t.K2, t.K3);
            Assert.Equal(t.Family == TripleFamily.Triv, selfMirror);
        }
    }

    // the law's two exact factors: the structural coupling-shape spectrum (3,1,0) = roots of
    // x^3 - 4x^2 + 3x (integer char poly, no eigensolver), and the integer pair levels times n --
    // a*n = 12 - sum lam^2 = 6 (ROT3) and 8 (PENT), spec(X) = (3a, a, 0). TRIV has no pair to twin
    // with (self-mirror) and asking for its pair levels throws.
    [Fact]
    public void Coupling_Shape_And_Pair_Levels_Are_Exact()
    {
        Assert.Equal((3, 1, 0), Seed.CouplingShapeSpectrum());
        Assert.Equal((18, 6, 0), Seed.PairLevelsTimesN(TripleFamily.Rot3));
        Assert.Equal((24, 8, 0), Seed.PairLevelsTimesN(TripleFamily.Pent));
        Assert.Throws<InvalidOperationException>(() => Seed.PairLevelsTimesN(TripleFamily.Triv));
    }
}
