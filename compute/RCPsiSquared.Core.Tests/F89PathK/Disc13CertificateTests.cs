using System;
using RCPsiSquared.Core.F89PathK;
using Xunit;
using Xunit.Abstractions;

namespace RCPsiSquared.Core.Tests.F89PathK;

/// <summary>The certificate-grade D-device on the weight-general path
/// (<see cref="FoldResultantCertificate.WeightCertifyDiscMultiplicity"/>): the (1,3)@N=6 scout's
/// deg/v/layers upgraded from two-prime agreement to the lc-divisor-bound certification, the same
/// proof shape as the committed (1,2) device. Controls: the weight path must reproduce the (1,2)
/// certificate values at N=4 (52/24/[20,4], ~90 ms, Category EVEN_DISC13_SCOUT) and N=6
/// (926/536/[124,133], 499 primes, Category SLOW_DISC13_CERT). Value pins measured 2026-08-11 and
/// unchanged since; the (1,2)@N=6 runtime of that day, ~11 s per parity, is the reference the
/// cost paragraph below calibrates against. Run: dotnet test compute/RCPsiSquared.Core.Tests
/// -c Release --filter "Category=SLOW_DISC13_CERT|Category=EVEN_DISC13_SCOUT" (~30 s since the
/// 2026-08-13 speedups, against ~2 min before them; the scout category also carries
/// Disc13LayerScoutTests' fast members, a SEPARATE xunit collection, so that filter runs it
/// concurrently and its cores come out of the certificate's own parallelism: time this class with
/// --filter FullyQualifiedName~Disc13CertificateTests, which is how the figures below were taken).
/// <para>Since 2026-08-13 the two FOLD-FIXED blocks here, (1,2)@N=4 and (1,3)@N=6, are certified
/// on the q-even halved sampling path (the checkerboard corollary: D = disc_Λ(F_res) is even in
/// q, so D(q) = E(q²) and E needs half the nodes). Every certified value is bit-identical to the
/// dense path, which is what N4_QEvenHalving_ReproducesTheDenseCertificateExactly gates at N=4 and
/// what the (1,3)@N=6 members below demonstrate at full size, since they now reach the values
/// committed from the dense path on 2026-08-11 by the halved route. The same day the prime loop
/// went PARALLEL, which is the larger of the two levers and an independent one.
/// <para>Measured 2026-08-13, one machine (24 logical cores, 12 physical), one session, three code
/// states of this file run in turn, so the two levers compose in a known order. They do not
/// SEPARATE: the parallel state was measured only with the halving already in, so no interaction
/// term exists. (1,3)@N=6, R-odd then R-even: serial and dense 109.2 s / 50.4 s; serial and halved
/// 62.0 s / 28.6 s (ratios 0.568 and 0.567, so 1.76×, short of the naive 2× because the
/// arbitrary-node Newton costs more per node than the unit-spaced one and the pre-loop is
/// untouched apart from the Taylor shift the licence read adds); parallel and halved 13.9 s /
/// 7.1 s, i.e. 4.5× and 4.0× on top, 7.9× and 7.1× overall. Well short of 24× because at N=6 the
/// SERIAL pre-loop is already a real share of the total. How large a share is DERIVED, not
/// measured, and the assumed core count moves it: solving 62.0 = s + (62.0 − s)/f for the R-odd
/// pair gives s ≈ 11.8 s at f = 24 and s ≈ 9.5 s at f = 12, so read it as 10-12 s. No core count
/// touches that part; skipLeadingFormSharpening is the lever that does, and it was NOT in any of
/// the three states. Expect the parallel factor to improve at a large block, where sampling
/// dominates, and never quote 24×. The dense (1,2)@N=6 sibling ran ~12 s per parity in the first
/// two states against its committed ~11 s, which is what makes the states comparable; it is 3.6 s
/// in the third, from the same parallelism.</para></summary>
public class Disc13CertificateTests
{
    private readonly ITestOutputHelper _out;
    public Disc13CertificateTests(ITestOutputHelper output) => _out = output;

    [Theory]
    [InlineData(4, 1, 2, false)]
    [InlineData(4, 1, 2, true)]
    [Trait("Category", "EVEN_DISC13_SCOUT")]
    public void N4_Control_WeightCertificateMatchesCommitted(int n, int wKet, int wBra, bool rOdd)
    {
        var sw = System.Diagnostics.Stopwatch.StartNew();
        var r = FoldResultantCertificate.WeightCertifyDiscMultiplicity(n, wKet, wBra, rOdd);
        _out.WriteLine($"N{n} ({wKet},{wBra}) rOdd={rOdd}: deg={r.TrueDiscriminantDegree} v={r.TrueQValuationD} " +
                       $"layers=[{string.Join(", ", r.DiscLayerDegrees)}] mD={r.InfinityRepeatedD} dBound={r.DiscriminantDegreeBound} " +
                       $"lcBound={r.LcDivisorBoundD} sampled={r.PrimesSampled} layerPrime={r.LayerPrime} " +
                       $"certified={r.DiscLayersCertified}  [{sw.Elapsed}]");
        Assert.True(r.DiscLayersCertified);
        Assert.Equal(52, r.TrueDiscriminantDegree);
        Assert.Equal(24, r.TrueQValuationD);
        Assert.Equal(new[] { 20, 4 }, r.DiscLayerDegrees);
        // THE PRIME-ORDER PIN for the parallel loop, and the only assertion in this file that
        // catches a fold in completion order: 1073741833 is the FIRST split prime above 2^30, and
        // the device reads its layer profile at the FIRST entry attaining both trueDegD and the
        // minimal v. Many primes attain both, so this value is not about the mathematics at all;
        // it says the results were folded in prime order. Scramble the order and some other
        // prime surfaces here while every other assertion in this file still passes.
        Assert.Equal(1073741833L, r.LayerPrime);
    }

    [Theory]
    [InlineData(6, 1, 2, false)]
    [InlineData(6, 1, 2, true)]
    [Trait("Category", "SLOW_DISC13_CERT")]
    public void N6_Control_WeightCertificateMatchesCommitted12(int n, int wKet, int wBra, bool rOdd)
    {
        var sw = System.Diagnostics.Stopwatch.StartNew();
        var r = FoldResultantCertificate.WeightCertifyDiscMultiplicity(n, wKet, wBra, rOdd);
        _out.WriteLine($"N{n} ({wKet},{wBra}) rOdd={rOdd}: deg={r.TrueDiscriminantDegree} v={r.TrueQValuationD} " +
                       $"layers=[{string.Join(", ", r.DiscLayerDegrees)}] mD={r.InfinityRepeatedD} dBound={r.DiscriminantDegreeBound} " +
                       $"lcBound={r.LcDivisorBoundD} sampled={r.PrimesSampled} layerPrime={r.LayerPrime} " +
                       $"certified={r.DiscLayersCertified}  [{sw.Elapsed}]");
        Assert.True(r.DiscLayersCertified);
        Assert.Equal(926, r.TrueDiscriminantDegree);
        Assert.Equal(536, r.TrueQValuationD);
        Assert.Equal(new[] { 124, 133 }, r.DiscLayerDegrees);
        // The halving's negative control by FACT and not by predicate: this block is not
        // fold-fixed, its exact support-parity read declines, and the odd layer degree 133 is
        // why it must (an even-in-q D has its non-zero roots in ± pairs).
        Assert.False(r.QEvenHalvingUsed);
    }

    /// <summary>The (1,3)@N=6 certificate, both parities (R-even 924 sampled primes past the
    /// lc-divisor bound 923; R-odd 1210 past 1209; runtimes in the class doc, halved since
    /// 2026-08-13, the prime counts unchanged by that). The scout's
    /// two-prime deg/v/layers upgraded to certified invariants: TrueDiscriminantDegree and
    /// TrueQValuationD are what the arc's Sturm closing step will lean on, exactly as the (1,2)
    /// gcd certificate leaned on the (1,2) D-device's. The leading-form m_D (51 / 42) sharpens the
    /// interpolation bound to 1620 / 2172 (the scouts use the crude 1722 / 2256); the certified
    /// layer profile reproduces the scout exactly, and MaxDiscMultiplicity = 4 bounds the true
    /// maximum from above (the mod-p lift is one-way upward; "Klein-forced ceiling" is the census
    /// sections' interpretation of the value, not something this method certifies). NOTE the
    /// scope the (1,2) device's consumers must not import: its MaxDiscMultiplicity ≤ 2 IS the
    /// β-exotic exclusion, and at (1,3) that exclusion does NOT fire, a Puiseux-3/2 point being
    /// hideable inside the multiplicity-4 layer.</summary>
    [Theory]
    [InlineData(6, 1, 3, false)]
    [InlineData(6, 1, 3, true)]
    [Trait("Category", "SLOW_DISC13_CERT")]
    public void N6_Certificate_Disc13(int n, int wKet, int wBra, bool rOdd)
    {
        var sw = System.Diagnostics.Stopwatch.StartNew();
        var r = FoldResultantCertificate.WeightCertifyDiscMultiplicity(n, wKet, wBra, rOdd,
            log: s => _out.WriteLine("  " + s));
        _out.WriteLine($"N{n} ({wKet},{wBra}) rOdd={rOdd}: deg={r.TrueDiscriminantDegree} v={r.TrueQValuationD} " +
                       $"layers=[{string.Join(", ", r.DiscLayerDegrees)}] mD={r.InfinityRepeatedD} dBound={r.DiscriminantDegreeBound} " +
                       $"lcBound={r.LcDivisorBoundD} sampled={r.PrimesSampled} layerPrime={r.LayerPrime} " +
                       $"certified={r.DiscLayersCertified}  [{sw.Elapsed}]");
        Assert.True(r.DiscLayersCertified);
        Assert.Equal(rOdd ? 2124 : 1572, r.TrueDiscriminantDegree);
        Assert.Equal(rOdd ? 1188 : 876, r.TrueQValuationD);
        Assert.Equal(rOdd ? new[] { 0, 132, 0, 168 } : new[] { 0, 108, 0, 120 }, r.DiscLayerDegrees);
        Assert.Equal(rOdd ? 42 : 51, r.InfinityRepeatedD);
        Assert.Equal(rOdd ? 2172 : 1620, r.DiscriminantDegreeBound);
        Assert.Equal(4, r.MaxDiscMultiplicity);
        Assert.True(r.PrimesSampled > r.LcDivisorBoundD);
        // THE PRIME-ORDER PIN for the parallel loop, and the only assertion in this file that
        // catches a fold in completion order: 1073741833 is the FIRST split prime above 2^30, and
        // the device reads its layer profile at the FIRST entry attaining both trueDegD and the
        // minimal v. Many primes attain both, so this value is not about the mathematics at all;
        // it says the results were folded in prime order. Scramble the order and some other
        // prime surfaces here while every other assertion in this file still passes.
        Assert.Equal(1073741833L, r.LayerPrime);
        // Fold-fixed (bra weight 3 = N/2), so the q-even halving carries this certificate: half
        // the sample points per prime, every certified value above unchanged.
        Assert.True(r.QEvenHalvingUsed);
        Assert.Equal((r.DiscriminantDegreeBound / 2) + 1 + 24, r.SamplesPerPrime);
    }

    /// <summary>Skipping the leading form reaches the SAME certificate through a slightly cruder
    /// bound. The sharpened bound subtracts 2·m_D from resDeg·(resDeg−1), and m_D costs a
    /// Qc SquarefreeLayers call whose unnormalized BigRational growth was measured at ~8 min for
    /// degree 80 and extrapolates to hours or days at degree 224. What it buys there is under 1%
    /// of the sample points, so at a large block the sharpening is the expensive half of the
    /// pre-loop and the cheap half of the answer. The scout methods of this file already skip it
    /// for exactly this reason; this gate makes the CERTIFIED path able to, and pins that the
    /// certified invariants do not move when it does: only DiscriminantDegreeBound rises and
    /// InfinityRepeatedD goes to its −1 sentinel. At N=4 the crude bound is 8·7 = 56 against the
    /// sharpened 52, i.e. 2 extra sample points, which is why this control is nearly free.</summary>
    [Theory]
    [InlineData(4, 1, 2, false)]
    [InlineData(4, 1, 2, true)]
    [Trait("Category", "EVEN_DISC13_SCOUT")]
    public void N4_CrudeBound_ReachesTheSameCertificate(int n, int wKet, int wBra, bool rOdd)
    {
        var sharp = FoldResultantCertificate.WeightCertifyDiscMultiplicity(n, wKet, wBra, rOdd);
        var crude = FoldResultantCertificate.WeightCertifyDiscMultiplicity(
            n, wKet, wBra, rOdd, skipLeadingFormSharpening: true);
        _out.WriteLine($"N{n} ({wKet},{wBra}) rOdd={rOdd}: sharpened dBound={sharp.DiscriminantDegreeBound} " +
                       $"(m_D={sharp.InfinityRepeatedD}, {sharp.SamplesPerPrime} samples/prime) vs crude " +
                       $"dBound={crude.DiscriminantDegreeBound} (m_D={crude.InfinityRepeatedD}, " +
                       $"{crude.SamplesPerPrime} samples/prime)");

        Assert.False(sharp.SharpeningSkipped);
        Assert.True(crude.SharpeningSkipped);
        Assert.Equal(-1, crude.InfinityRepeatedD);
        Assert.Equal(56, crude.DiscriminantDegreeBound);              // resDeg·(resDeg−1) = 8·7, resDeg = 8
        Assert.True(crude.DiscriminantDegreeBound > sharp.DiscriminantDegreeBound);

        Assert.True(crude.DiscLayersCertified);
        Assert.Equal(sharp.TrueDiscriminantDegree, crude.TrueDiscriminantDegree);
        Assert.Equal(sharp.TrueQValuationD, crude.TrueQValuationD);
        Assert.Equal(sharp.DiscLayerDegrees, crude.DiscLayerDegrees);
        Assert.Equal(sharp.MaxDiscMultiplicity, crude.MaxDiscMultiplicity);
        Assert.Equal(sharp.LcDivisorBoundD, crude.LcDivisorBoundD);   // the cap is a norm read, not a degree read
        Assert.Equal(sharp.PrimesSampled, crude.PrimesSampled);
        Assert.Equal(sharp.LayerPrime, crude.LayerPrime);
    }

    /// <summary>The fold-fixed predicate: a block is fixed by a fold leg exactly when N is even and
    /// one weight side sits at N/2 (f_P fixes bra weight N/2, f_Q fixes ket weight N/2). It is NOT
    /// the halving's precondition, and the distinction is the point of pinning it: the device
    /// decides on its own exact ℤ[i] support-parity read, so a block that is not fold-fixed would
    /// halve if it passed that read, while a fold-fixed block that fails it throws. The predicate
    /// is the reason to EXPECT the read to pass, and the reason a failed read is a finding.
    /// Negative controls: (1,2)@N=6 and (1,2)@N=8, where the halving must not fire, and does not,
    /// as a FACT about D and not only as policy. At (1,2)@N=6 the certified layer profile
    /// [124, 133] carries an odd degree, which an even-in-q D cannot produce (its non-zero roots
    /// come in ± pairs, so every layer degree would be even); at (1,2)@N=8 the committed gcd
    /// certificate has deg Im D = 6085, odd for the same reason.</summary>
    [Theory]
    [InlineData(4, 1, 2, true)]         // bra weight 2 = N/2
    [InlineData(4, 2, 1, true)]         // ket weight 2 = N/2
    [InlineData(6, 1, 3, true)]         // the census block
    [InlineData(8, 1, 4, true)]         // the payoff block
    [InlineData(6, 1, 2, false)]        // negative control, layer degree 133 is odd
    [InlineData(8, 1, 2, false)]        // negative control, the headline device run
    [InlineData(7, 1, 3, false)]        // odd N has no N/2 seat at all
    [Trait("Category", "EVEN_DISC13_SCOUT")]
    public void FoldFixedPredicate_IsWhyTheHalvingIsExpected(int n, int wKet, int wBra, bool expected)
        => Assert.Equal(expected, FoldResultantCertificate.IsFoldFixedBlock(n, wKet, wBra));

    /// <summary>The q-even halving reproduces the dense device EXACTLY, at the fast fold-fixed
    /// control (1,2)@N=4, both parities. This is the whole safety argument for the speedup: it is
    /// not a new certificate, it is the same certificate reached with half the sample points, so
    /// every certified field must agree bit for bit with the dense path, and the halved path must
    /// actually have halved (SamplesPerPrime, not a claim about wall clock). Where the soundness
    /// sits: on the exact ℤ[i] support-parity read plus the pre-existing proven dBound, which
    /// together determine E from dBound/2 + 1 values. NOT on the 24 verification nodes, which
    /// over-determine the dense path but under-determine the halved one and are a smoke test for
    /// the interpolation code rather than a falsifier of the premise. Read the margin here before
    /// trusting this gate too far, and read it in FREE ROOTS, the only unambiguous unit: at N=4,
    /// dBound = 52, so after the 27 interpolation nodes the difference D − dp still has
    /// dBound/2 − 1 = 25 free roots against 24 verification points. It is therefore NOT pinned by
    /// node count, which is what makes this A/B a real test of the halving. Had dBound been ≤ 48,
    /// the 24 points would have determined D outright, the halved path would have been correct
    /// for a reason having nothing to do with evenness, and the gate would have proved nothing.
    /// The margin is one root wide. The (1,3)@N=6 members carry the full-size version of the same
    /// comparison, with far more room, against values committed from the dense path.</summary>
    [Theory]
    [InlineData(4, 1, 2, false)]
    [InlineData(4, 1, 2, true)]
    [Trait("Category", "EVEN_DISC13_SCOUT")]
    public void N4_QEvenHalving_ReproducesTheDenseCertificateExactly(int n, int wKet, int wBra, bool rOdd)
    {
        var dense = FoldResultantCertificate.WeightCertifyDiscMultiplicity(
            n, wKet, wBra, rOdd, forceDenseSampling: true);
        var halved = FoldResultantCertificate.WeightCertifyDiscMultiplicity(n, wKet, wBra, rOdd);
        _out.WriteLine($"N{n} ({wKet},{wBra}) rOdd={rOdd}: dense samples/prime = {dense.SamplesPerPrime} " +
                       $"(halving={dense.QEvenHalvingUsed}), halved = {halved.SamplesPerPrime} " +
                       $"(halving={halved.QEvenHalvingUsed}), dBound={halved.DiscriminantDegreeBound}");

        Assert.False(dense.QEvenHalvingUsed);
        Assert.True(halved.QEvenHalvingUsed);
        Assert.Equal((halved.DiscriminantDegreeBound / 2) + 1 + 24, halved.SamplesPerPrime);
        Assert.Equal(halved.DiscriminantDegreeBound + 1 + 24, dense.SamplesPerPrime);

        Assert.True(dense.DiscLayersCertified);
        Assert.True(halved.DiscLayersCertified);
        Assert.Equal(dense.InfinityRepeatedD, halved.InfinityRepeatedD);
        Assert.Equal(dense.DiscriminantDegreeBound, halved.DiscriminantDegreeBound);
        Assert.Equal(dense.TrueDiscriminantDegree, halved.TrueDiscriminantDegree);
        Assert.Equal(dense.TrueQValuationD, halved.TrueQValuationD);
        Assert.Equal(dense.DiscLayerDegrees, halved.DiscLayerDegrees);
        Assert.Equal(dense.MaxDiscMultiplicity, halved.MaxDiscMultiplicity);
        Assert.Equal(dense.LayerPrime, halved.LayerPrime);
        Assert.Equal(dense.PrimesSampled, halved.PrimesSampled);
        Assert.Equal(dense.LcDivisorBoundD, halved.LcDivisorBoundD);
        // The batching pin for the parallel prime loop: the device stops at exactly cap + 1
        // sampled primes, the count the sequential loop stopped at, so a batch that overshot or
        // undershot moves this number. Note what this A/B canNOT gate, since both sides run the
        // same parallel loop: a fold in completion order rather than prime order would scramble
        // both identically and pass. The guard against THAT is the literal LayerPrime pin in the
        // value-pinning controls, which is prime-order-dependent output.
        Assert.Equal(halved.LcDivisorBoundD + 1, halved.PrimesSampled);
    }
}
