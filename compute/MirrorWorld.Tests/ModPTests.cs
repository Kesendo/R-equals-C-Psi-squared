using System.Numerics;
using MirrorWorld;

namespace MirrorWorldTests;

// From-below guard for ModP, the world's one exact-arithmetic primitive, which Seed, LevelCollision,
// BlindSeat, Divisor and Crack are all written in.
//
// Every gate judges ModP against a route independent OF IT: System.Numerics.BigInteger for the modular
// arithmetic, trial division for primality, a membership question over every prime up to m for
// factorisation (deliberately not the peel-off loop the implementation uses, which would catch only a
// typo), all proper divisors for the multiplicative order. No gate here is ModP compared with itself.
//
// The inputs are chosen to BREAK the arithmetic rather than to pass it: the top residues p-1, p-2,
// which are enough to separate anything evaluated in int (the trap that cost four hours on
// 2026-07-03, (a % p + p) widened too late) from a long path, though they are not the only inputs
// that would; a modulus near 2^61, which is what separates the long path from the UInt128 one; an
// entry equal to p, which is zero in the field and is what the pivot search must not take; and, on
// each end of the rank, a matrix whose deficiency is known by construction.
//
// Three of these gates exist because the first version of this file did not have them, and three
// mutations (UInt128 to long, the 1 mod 4 guard deleted, the reduction in Rank deleted) all survived
// it. They fail it now.
public class ModPTests
{
    static long Big => ModP.Primes.Max();

    // ---- the prime list itself: the reason it moved off the old 3 mod 4 pair ----

    // Both primes are prime, by trial division, not by ModP.IsPrime (which is under test here).
    [Fact]
    public void ThePrimes_AreActuallyPrime_ByTrialDivision()
    {
        Assert.Equal(2, ModP.Primes.Length);
        Assert.Equal(2, ModP.Primes.Distinct().Count());
        foreach (long p in ModP.Primes)
        {
            Assert.True(p > 1);
            for (long d = 2; d * d <= p; d++)
                Assert.True(p % d != 0, $"{p} is divisible by {d}");
        }
    }

    // The whole point of the list: p = 1 mod 4, so -1 is a square and the Gaussian ranks (Divisor)
    // and the plain ranks (Seed, BlindSeat) can share one list. 2^31 - 1 is 3 mod 4 and cannot serve.
    [Fact]
    public void ThePrimes_AreOneModFour_SoMinusOneIsASquare()
    {
        // (2147483647 is 3 mod 4, which is why the list the other objects carried could not be the
        // shared one; that is arithmetic, stated here rather than asserted against itself.)
        foreach (long p in ModP.Primes) Assert.Equal(1, p % 4);
    }

    // ---- the modular atoms, against BigInteger, at the residues that break int arithmetic ----

    [Fact]
    public void MulMod_MatchesBigInteger_AtTheTopResidues()
    {
        long p = Big;
        foreach (var (a, b) in new[] { (p - 1, p - 1), (p - 1, p - 2), (p - 2, p - 2), (p / 2, p - 1), (1L, p - 1) })
        {
            long expected = (long)(new BigInteger(a) * b % p);
            Assert.Equal(expected, ModP.MulMod(a, b, p));
        }
    }

    // The list's primes are below 2^30, so a reduced product fits a long and the UInt128 route is never
    // needed by anything the world runs today. It is a contract of the primitive rather than a property
    // of the current callers, and F161/F162 are meant to be written in it, so it is read where it bites:
    // at the Mersenne prime 2^61 - 1 a plain long product of two top residues wraps and BigInteger does not.
    [Fact]
    public void MulMod_IsExact_AtAModulusWhereALongProductWouldWrap()
    {
        const long p = 2305843009213693951L;   // 2^61 - 1, prime
        Assert.True(ModP.IsPrime(p));
        foreach (var (a, b) in new[] { (p - 1, p - 1), (p - 1, p - 2), (p / 2, p - 1) })
        {
            long expected = (long)(new BigInteger(a) * b % p);
            Assert.Equal(expected, ModP.MulMod(a, b, p));
        }
    }

    [Fact]
    public void ModPow_MatchesBigInteger_AtTheTopResidues()
    {
        long p = Big;
        foreach (long b in new[] { p - 1, p - 2, p / 2, 2L })
        {
            long expected = (long)BigInteger.ModPow(b, p - 2, p);
            Assert.Equal(expected, ModP.ModPow(b, p - 2, p));
        }
    }

    // The inverse is a claim about a PRODUCT, so the gate reads the product, not the inverse.
    [Fact]
    public void ModInverse_TimesItsInput_IsOne_IncludingNegativeInputs()
    {
        long p = Big;
        foreach (long x in new[] { 1L, 2L, p / 2, p - 2, p - 1, -1L, -(p - 1), -(p * 3 + 7) })
        {
            BigInteger reduced = ((new BigInteger(x) % p) + p) % p;
            Assert.Equal(BigInteger.One, reduced * ModP.ModInverse(x, p) % p);
        }
    }

    // Mod is the trap's own line: it must reduce into [0, p) for inputs on both sides and far outside.
    [Fact]
    public void Mod_LandsInTheResidueRange_OnBothSigns()
    {
        long p = Big;
        foreach (long x in new[] { 0L, 1L, -1L, p, -p, p - 1, -(p - 1), 5 * p + 17, -(5 * p + 17) })
        {
            long r = ModP.Mod(x, p);
            Assert.InRange(r, 0, p - 1);
            Assert.Equal((long)(((new BigInteger(x) % p) + p) % p), r);
        }
    }

    // ---- the square root of -1: what the 1 mod 4 list buys ----

    [Fact]
    public void SqrtMinusOne_Squares_To_MinusOne_AtBothPrimes()
    {
        foreach (long p in ModP.Primes)
        {
            long r = ModP.SqrtMinusOne(p);
            Assert.Equal(new BigInteger(p - 1), new BigInteger(r) * r % p);
        }
    }

    // and refuses where it does not exist, rather than returning something that squares to anything.
    [Fact]
    public void SqrtMinusOne_Throws_AtAThreeModFourPrime()
    {
        var ex = Assert.Throws<InvalidOperationException>(() => ModP.SqrtMinusOne(2147483647L));
        // SqrtMinusOne throws on two paths; without this the search-exhausted throw would pass the test
        // just as well, and the mod-4 guard, the one line this gate exists for, could be deleted unseen.
        Assert.Contains("1 mod 4", ex.Message);
    }

    // ---- primality and factorisation, against trial division ----

    [Fact]
    public void IsPrime_AgreesWithTrialDivision_OnASmallRange()
    {
        for (long m = 2; m <= 2000; m++)
        {
            bool trial = true;
            for (long d = 2; d * d <= m; d++) if (m % d == 0) { trial = false; break; }
            Assert.Equal(trial, ModP.IsPrime(m));
        }
        Assert.False(ModP.IsPrime(1));
        Assert.False(ModP.IsPrime(0));
        Assert.False(ModP.IsPrime(-7));
    }

    // the large end, where trial division is still cheap but the Miller-Rabin path is the one taken.
    [Fact]
    public void IsPrime_AgreesWithTrialDivision_OnTheLargeCandidatesInPlay()
    {
        foreach (long m in new[] { 998244353L, 1004535809L, 2147483647L, 999999937L, 998244353L * 2, 1004535809L + 1 })
        {
            bool trial = m >= 2;
            for (long d = 2; d * d <= m; d++) if (m % d == 0) { trial = false; break; }
            Assert.Equal(trial, ModP.IsPrime(m));
        }
    }

    // Deliberately NOT the same peel-off loop retyped, which would only catch a typo: the question here
    // is the membership one. Over every prime q <= m, q divides m exactly when PrimeFactors reports it.
    [Fact]
    public void PrimeFactors_ReportsExactlyThePrimesThatDivide()
    {
        for (int m = 2; m <= 500; m++)
        {
            var reported = ModP.PrimeFactors(m).ToList();
            Assert.Equal(reported.OrderBy(q => q).ToList(), reported);      // ascending, as documented
            Assert.Equal(reported.Distinct().Count(), reported.Count);      // distinct, as documented
            for (int q = 2; q <= m; q++)
            {
                bool qIsPrime = true;
                for (int d = 2; d * d <= q; d++) if (q % d == 0) { qIsPrime = false; break; }
                if (!qIsPrime) continue;
                Assert.Equal(m % q == 0, reported.Contains(q));
            }
        }
    }

    // ---- a root of exact multiplicative order: the cyclotomic atom Seed and LevelCollision share ----

    [Theory]
    [InlineData(8)]
    [InlineData(12)]
    [InlineData(20)]
    [InlineData(24)]
    public void RootOfOrder_HasExactlyThatOrder(int order)
    {
        // a prime with order | p-1, found the way the callers find it
        long p = 0;
        for (long cand = 2L * order + 1; cand < 10_000_000; cand += order)
            if ((cand - 1) % order == 0 && ModP.IsPrime(cand)) { p = cand; break; }
        Assert.True(p > 0);

        long z = ModP.RootOfOrder(order, p);
        Assert.NotEqual(0, z);
        Assert.Equal(1, (long)BigInteger.ModPow(z, order, p));               // z^order = 1
        // the order is EXACTLY this, read without ModP.PrimeFactors: a shared miss there would make
        // implementation and oracle skip the same divisor. Every proper divisor is tried instead.
        for (int d = 1; d < order; d++)
            if (order % d == 0)
                Assert.NotEqual(1, (long)BigInteger.ModPow(z, d, p));
    }

    // ---- the rank, against ranks known by construction ----

    [Fact]
    public void Rank_OfTheIdentity_IsItsSize()
    {
        for (int d = 1; d <= 6; d++)
        {
            var rows = Enumerable.Range(0, d)
                .Select(i => Enumerable.Range(0, d).Select(j => i == j ? 1L : 0L).ToArray()).ToList();
            Assert.Equal(d, ModP.Rank(rows, Big));
        }
    }

    // the anti-vacuity partner of the line above: a deficient matrix must come back deficient.
    // (Most ways of breaking an elimination return the FULL size or zero, so both ends are pinned.)
    [Fact]
    public void Rank_OfADependentSet_IsTheIndependentCount()
    {
        var rows = new List<long[]>
        {
            new[] { 1L, 2L, 3L },
            new[] { 2L, 4L, 6L },      // twice the first
            new[] { 0L, 1L, 1L },
            new[] { 1L, 3L, 4L },      // first + third
        };
        Assert.Equal(2, ModP.Rank(rows, Big));
    }

    // the rank must be the rank OVER GF(p), not over Q: these rows are independent over the integers
    // and dependent mod the larger prime, which is exactly why the callers take the MAXIMUM over two.
    [Fact]
    public void Rank_IsTakenOverTheField_NotOverTheIntegers()
    {
        long p = Big;
        var rows = new List<long[]> { new[] { 1L, 1L }, new[] { 1L, 1L + p } };
        Assert.Equal(1, ModP.Rank(rows, p));                 // mod p the two rows coincide
        Assert.Equal(2, ModP.Rank(rows, ModP.Primes.Min())); // at the other prime they do not
        Assert.Equal(2, ModP.Rank(rows));                    // and the max over the list recovers it
    }

    // the top-residue face of the rank: entries just below p, where an int truncation of any
    // intermediate product silently changes the answer. The expected value is the BigInteger
    // determinant, computed here, not a number asserted from the shape of the matrix.
    [Fact]
    public void Rank_IsExact_AtTopResidueEntries()
    {
        long p = Big;
        var rows = new List<long[]>
        {
            new[] { p - 1, p - 2, p - 3 },
            new[] { p - 2, p - 3, p - 1 },
            new[] { p - 3, p - 1, p - 2 },
        };
        var m = new BigInteger[3, 3];
        for (int i = 0; i < 3; i++) for (int j = 0; j < 3; j++) m[i, j] = rows[i][j];
        BigInteger det = m[0, 0] * (m[1, 1] * m[2, 2] - m[1, 2] * m[2, 1])
                       - m[0, 1] * (m[1, 0] * m[2, 2] - m[1, 2] * m[2, 0])
                       + m[0, 2] * (m[1, 0] * m[2, 1] - m[1, 1] * m[2, 0]);
        Assert.NotEqual(BigInteger.Zero, det % p);   // else the line below would assert nothing about rank 3
        Assert.Equal(3, ModP.Rank(rows, p));
    }

    // a rectangular case: more rows than columns, the shape BlindSeat's Krylov matrix has.
    [Fact]
    public void Rank_HandlesMoreRowsThanColumns()
    {
        var rows = new List<long[]>
        {
            new[] { 1L, 0L },
            new[] { 0L, 1L },
            new[] { 3L, 5L },
            new[] { 7L, 11L },
        };
        Assert.Equal(2, ModP.Rank(rows, Big));
    }

    // ---- the arguments themselves: this primitive exists to stop a silent wrong residue, so an
    // input it cannot serve must say so rather than return a plausible number ----

    // The docstring licenses ragged input to nobody: cols comes from the first row, so a short row
    // would fault mid-elimination and a long one would be silently truncated. F161/F162 are meant to
    // be written in this primitive, so the refusal is the contract, not a convenience.
    [Fact]
    public void Rank_Refuses_RaggedRows()
    {
        var shortRow = new List<long[]> { new[] { 1L, 2L, 3L }, new[] { 4L, 5L } };
        var longRow = new List<long[]> { new[] { 1L, 2L }, new[] { 3L, 4L, 5L } };
        Assert.Throws<ArgumentException>(() => ModP.Rank(shortRow, Big));
        Assert.Throws<ArgumentException>(() => ModP.Rank(longRow, Big));
    }

    // A non-positive modulus would reach MulMod as a huge unsigned one and return a plausible residue.
    [Fact]
    public void EveryEntryPoint_Refuses_ANonPositiveModulus()
    {
        foreach (long bad in new[] { 0L, -1L, -998244353L })
        {
            Assert.Throws<ArgumentOutOfRangeException>(() => ModP.Mod(5, bad));
            Assert.Throws<ArgumentOutOfRangeException>(() => ModP.MulMod(5, 7, bad));
            Assert.Throws<ArgumentOutOfRangeException>(() => ModP.ModPow(5, 3, bad));
            Assert.Throws<ArgumentOutOfRangeException>(() => ModP.ModInverse(5, bad));
            Assert.Throws<ArgumentOutOfRangeException>(() => ModP.Rank(new List<long[]> { new[] { 1L } }, bad));
        }
    }

    // Zero has no inverse; Fermat returns 0 for it, which multiplies into 0 and not into 1.
    [Fact]
    public void ModInverse_Refuses_Zero()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => ModP.ModInverse(0, Big));
        Assert.Throws<ArgumentOutOfRangeException>(() => ModP.ModInverse(Big, Big));
        Assert.Throws<ArgumentOutOfRangeException>(() => ModP.ModInverse(-Big, Big));
    }

    // The documented contract is that entries may be unreduced and of either sign. The pivot search is
    // where that bites: an entry equal to p is ZERO in the field, and an unreduced search would take it
    // as a pivot and report rank where there is none.
    [Fact]
    public void Rank_ReadsUnreducedAndNegativeEntries_InTheField()
    {
        long p = Big;
        Assert.Equal(0, ModP.Rank(new List<long[]> { new[] { p, p }, new[] { 0L, p } }, p));
        Assert.Equal(1, ModP.Rank(new List<long[]> { new[] { p, 1L }, new[] { 2 * p, 2L } }, p));
        // and the signed side, which is the shape Seed hands it (a pencil of -1 and +1 entries)
        Assert.Equal(1, ModP.Rank(new List<long[]> { new[] { -1L, -2L }, new[] { -2L, -4L } }, p));
        Assert.Equal(2, ModP.Rank(new List<long[]> { new[] { -1L, -2L }, new[] { -2L, -3L } }, p));
    }

    // M-4's partner: the deficiency end read at top-residue entries, not only at entries below ten.
    // Without this, every top-residue rank assertion is a FULL rank, which a broken elimination returns
    // for free.
    [Fact]
    public void Rank_SeesDeficiency_AtTopResidueEntries()
    {
        long p = Big;
        var rows = new List<long[]>
        {
            new[] { p - 1, p - 2, p - 3 },
            new[] { p - 2, p - 4, p - 6 },   // twice the first, mod p
            new[] { p - 3, p - 6, p - 9 },   // three times the first, mod p
        };
        Assert.Equal(1, ModP.Rank(rows, p));
    }

    [Fact]
    public void Rank_OfNothing_IsZero()
    {
        Assert.Equal(0, ModP.Rank(new List<long[]>(), Big));
        Assert.Equal(0, ModP.Rank(new List<long[]> { new[] { 0L, 0L }, new[] { 0L, 0L } }, Big));
    }
}
