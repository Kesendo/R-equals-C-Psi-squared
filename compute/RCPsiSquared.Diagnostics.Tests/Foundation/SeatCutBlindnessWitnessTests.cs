using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Diagnostics.Foundation;

namespace RCPsiSquared.Diagnostics.Tests.Foundation;

/// <summary>From-below pins for <see cref="SeatCutBlindnessWitness"/>, the F157 live lab (proof
/// <c>docs/proofs/PROOF_BLIND_SEAT_SPAN_AND_NODE_LEMMA.md</c>, gates
/// <c>simulations/blind_seat_span_proof.py</c> and <c>simulations/seat_cut_blindness.py</c>).
///
/// <para>The witness is an independent construction, not a port: it rebuilds the integer
/// single-excitation matrix in C# and ranks it over GF(p), so agreeing with the committed Python
/// run is two computations meeting rather than one repeated. The expected values below are read
/// from that run and from the proof file's tables, never recomputed from the closed form the
/// witness is being checked against.</para>
///
/// <para>Two-sided throughout: every positive count sits beside a zero the same routine produces,
/// and the span identity is pinned where it HOLDS and where it BREAKS, since a test that only ever
/// saw it hold would pass equally for a witness that returned 1 + blind by construction.</para></summary>
public class SeatCutBlindnessWitnessTests
{
    [Fact]
    public void UniformHeisenbergN7_TheLiveCountMeetsTheClosedForm()
    {
        var w = new SeatCutBlindnessWitness(7);
        Assert.True(w.IsUniform);
        Assert.Equal(7, w.ClosedFormAgreements());
        Assert.Equal(new[] { 3 }, w.BlindSeats());
        Assert.Equal(3, w.Blind(3));
    }

    [Fact]
    public void UniformXyN7_BlindsThreeSeatsWhereHeisenbergBlindsOne()
    {
        var xy = new SeatCutBlindnessWitness(7, SeatCutBook.Xy);
        Assert.Equal(new[] { 1, 3, 5 }, xy.BlindSeats());
        Assert.Equal(1, xy.Blind(1));
        Assert.Equal(3, xy.Blind(3));
        Assert.Equal(1, xy.Blind(5));
        Assert.Equal(7, xy.ClosedFormAgreements());
    }

    [Fact]
    public void TheEndSeatsAreZero_AndAGenericChainHasNoBlindSeatAtAll()
    {
        // The falsifiers. Lemma J1 forces both ends to 0 on any zero-free profile, and on the
        // Heisenberg book an irregular zero-free chain generically blinds nobody: a routine that
        // returned N, or that mistook the Krylov rank for its complement, would fail here.
        var irregular = new SeatCutBlindnessWitness(6, SeatCutBook.Heisenberg, new long[] { 1, -2, 3, -1, 5 });
        Assert.Empty(irregular.BlindSeats());
        for (int n = 3; n <= 12; n++)
        {
            var w = new SeatCutBlindnessWitness(n);
            Assert.Equal(0, w.Blind(0));
            Assert.Equal(0, w.Blind(n - 1));
        }
    }

    [Fact]
    public void ParityForcesEveryOddSeatOfAnOddXyChain_AtSignedAndIrregularProfiles()
    {
        // The third kind: not a coincidence the uniform law counts, but a parity fact that survives
        // any zero-free profile. Heisenberg is the control and stays blind-free on the same bonds.
        var irregular = new[]
        {
            new long[] { 1, -2, 3, -1, 2, 5, -3, 1 },
            new long[] { 7, 7, 2, -9, 4, 4, -1, 6 },
        };
        foreach (var bonds in irregular)
        {
            var xy = new SeatCutBlindnessWitness(9, SeatCutBook.Xy, bonds);
            // Every odd seat is forced, at a signed and irregular profile where no gcd law applies.
            Assert.Equal(new[] { 1, 3, 5, 7 }, xy.BlindSeats());
            foreach (int seat in new[] { 1, 3, 5, 7 })
                Assert.Equal(1, xy.Blind(seat));
            // The paired zero, and it is the whole content: the SAME profile on the Heisenberg book
            // blinds nobody, so what fires here is the XY book's zero diagonal and not the disorder.
            Assert.Empty(new SeatCutBlindnessWitness(9, SeatCutBook.Heisenberg, bonds).BlindSeats());
        }

        // On the UNIFORM chain an even seat can be blind too, by coincidence rather than by parity:
        // seat 4 carries gcd(5, 10) - 1 = 4. Asserting "even seats are zero" would have been false
        // here, and the difference between the forced kind and the met kind is exactly this.
        var uniform = new SeatCutBlindnessWitness(9, SeatCutBook.Xy);
        Assert.Equal(new[] { 1, 3, 4, 5, 7 }, uniform.BlindSeats());
        Assert.Equal(4, uniform.Blind(4));
    }

    [Fact]
    public void OnAZeroFreeChain_TheSpanIdentityHoldsAtEverySeat_BothBooks()
    {
        // Corollary B, live: Lemma J makes the spectrum simple, so 1 + blind IS the kernel.
        foreach (var book in new[] { SeatCutBook.Heisenberg, SeatCutBook.Xy })
            foreach (int n in new[] { 3, 4, 5, 6, 7, 8 })
            {
                var w = new SeatCutBlindnessWitness(n, book, Enumerable.Range(1, n - 1).Select(i => (long)i).ToArray());
                for (int seat = 0; seat < n; seat++)
                    Assert.Equal(1 + w.Blind(seat), w.Span(seat));
            }
    }

    [Fact]
    public void TheZeroBondPathBreaksTheSpanAtTwoSeats_AndTheOtherZeroBondPathDoesNot()
    {
        // The discriminating pair from the proof file's Corollary C, and the whole reason the fence
        // was moved off the zero bond and onto the spectrum of H on the Krylov complement.
        var broken = new SeatCutBlindnessWitness(6, SeatCutBook.Heisenberg, new long[] { 1, 1, 0, 1, 1 });
        Assert.Equal(new[] { 3, 4, 3, 3, 4, 3 },
            Enumerable.Range(0, 6).Select(broken.Blind).ToArray());
        Assert.Equal(new int?[] { 4, 7, 4, 4, 7, 4 },
            Enumerable.Range(0, 6).Select(broken.Span).ToArray());
        Assert.Equal(2, Enumerable.Range(0, 6).Count(s => broken.Span(s) != 1 + broken.Blind(s)));

        // [1,0,1] is cut and degenerate too, and holds at every seat.
        var holds = new SeatCutBlindnessWitness(4, SeatCutBook.Heisenberg, new long[] { 1, 0, 1 });
        for (int seat = 0; seat < 4; seat++)
            Assert.Equal(1 + holds.Blind(seat), holds.Span(seat));
    }

    [Fact]
    public void TheCountWalksPastTheWall()
    {
        // The count is an N x N elimination, so the law is checked where no eigendecomposition goes.
        // N = 60 and 120 here rather than the guard's 200: the constructor counts every seat, so the
        // cost is O(N^4) in the suite and N = 200 alone took most of a minute. MirrorWorld's
        // BlindSeatTests pin N = 200, which is where that reading belongs.
        foreach (int n in new[] { 60, 120 })
        {
            var w = new SeatCutBlindnessWitness(n);
            Assert.Equal(n, w.ClosedFormAgreements());
            var xy = new SeatCutBlindnessWitness(n, SeatCutBook.Xy);
            Assert.Equal(n, xy.ClosedFormAgreements());
        }
    }

    [Fact]
    public void AboveTheSpanGuard_TheSpanIsNotComputedRatherThanGuessed()
    {
        var w = new SeatCutBlindnessWitness(SeatCutBlindnessWitness.MaxSpanN + 1);
        Assert.Null(w.Span(0));
        Assert.NotNull(new SeatCutBlindnessWitness(SeatCutBlindnessWitness.MaxSpanN).Span(0));
    }

    [Fact]
    public void OffAUniformProfile_TheClosedFormIsReportedInapplicableRatherThanWrong()
    {
        var w = new SeatCutBlindnessWitness(6, SeatCutBook.Heisenberg, new long[] { 1, 2, 3, 4, 5 });
        Assert.False(w.IsUniform);
        Assert.Null(w.ClosedFormAgreements());
        Assert.Null(w.ClosedForm(2));
        // A zero bond is not a uniform profile either, even when every nonzero entry agrees.
        Assert.False(new SeatCutBlindnessWitness(6, SeatCutBook.Heisenberg, new long[] { 1, 1, 0, 1, 1 }).IsUniform);
    }

    [Fact]
    public void ACouplingDivisibleByARankingPrime_DoesNotFakeBlindness()
    {
        // A rank over GF(p) is blind to a factor of p in every entry, so a coupling divisible by a
        // ranking prime would report the whole space blind IF THAT PRIME WERE THE ONLY ONE. The
        // defence is that there are two and the larger rank is kept: p2 does not divide a multiple
        // of p1, so the reading must come out identical to the plain uniform chain.
        const long p1 = 2147483647L;
        var w = new SeatCutBlindnessWitness(6, SeatCutBook.Heisenberg,
            new[] { p1, p1, p1, p1, p1 });
        var plain = new SeatCutBlindnessWitness(6);
        Assert.Equal(Enumerable.Range(0, 6).Select(plain.Blind).ToArray(),
                     Enumerable.Range(0, 6).Select(w.Blind).ToArray());
        Assert.Equal(new[] { 1, 4 }, w.BlindSeats());

        // A non-uniform profile has no closed form to catch it, so the reduction is the only guard
        // there. Same profile scaled by the prime, same answer.
        var scaled = new SeatCutBlindnessWitness(6, SeatCutBook.Heisenberg,
            new[] { p1, 2 * p1, 3 * p1, 4 * p1, 5 * p1 });
        var unscaled = new SeatCutBlindnessWitness(6, SeatCutBook.Heisenberg, new long[] { 1, 2, 3, 4, 5 });
        Assert.Equal(Enumerable.Range(0, 6).Select(unscaled.Blind).ToArray(),
                     Enumerable.Range(0, 6).Select(scaled.Blind).ToArray());

        // And the case two primes could NOT survive, divisibility by both at once, needs a coupling
        // of 2147483647 * 999999937; the magnitude guard refuses that outright, so the two defences
        // are exhaustive rather than probabilistic.
        Assert.True(p1 * 999999937L > SeatCutBlindnessWitness.MaxCoupling);
    }

    [Fact]
    public void AnOversizedCoupling_ThrowsRatherThanWrappingSilently()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => new SeatCutBlindnessWitness(
            4, SeatCutBook.Heisenberg, new[] { SeatCutBlindnessWitness.MaxCoupling + 1, 1L, 1L }));
    }

    [Theory]
    [InlineData(1)]
    [InlineData(0)]
    [InlineData(SeatCutBlindnessWitness.MaxN + 1)]
    public void OutOfRangeN_Throws(int n)
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => new SeatCutBlindnessWitness(n));
    }

    [Fact]
    public void AWrongLengthProfile_Throws()
    {
        Assert.Throws<ArgumentException>(() => new SeatCutBlindnessWitness(6, SeatCutBook.Heisenberg, new long[] { 1, 1 }));
    }

    [Fact]
    public void OutOfRangeSeat_Throws()
    {
        var w = new SeatCutBlindnessWitness(5);
        Assert.Throws<ArgumentOutOfRangeException>(() => w.Blind(-1));
        Assert.Throws<ArgumentOutOfRangeException>(() => w.Blind(5));
        Assert.Throws<ArgumentOutOfRangeException>(() => w.Span(5));
    }

    [Fact]
    public void TheInspectTree_IsLiveAndCarriesTheFalsifiers()
    {
        var w = new SeatCutBlindnessWitness(7);
        var children = w.Children.ToList();
        Assert.Contains(children, c => c.DisplayName.Contains("falsifier"));
        Assert.Contains(children, c => c.DisplayName.Contains("span identity"));
        Assert.Contains(children, c => c.Provenance == NodeProvenance.Live);
        Assert.Contains("MATCH", children.Single(c => c.DisplayName.Contains("vs the closed form")).Summary);
    }
}
