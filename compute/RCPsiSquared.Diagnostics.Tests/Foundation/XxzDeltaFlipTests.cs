using System;
using System.Linq;
using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.F89PathK;
using RCPsiSquared.Core.Numerics;
using RCPsiSquared.Diagnostics.Foundation;
using Xunit;

namespace RCPsiSquared.Diagnostics.Tests.Foundation;

/// <summary>The F89 Delta proposal/certification tests: validate the (q,Δ) XxzCoherenceBlock against
/// the existing F89 block at Δ=0 and certify the N=4 Delta=0 control. Every nonzero-Delta proposal must pass
/// independent full-block coincidence, seed correspondence and isolated-pair character checks.
/// See docs/superpowers/plans/2026-06-27-f89-path4-delta-test.md.</summary>
public class XxzDeltaFlipTests
{
    [Theory]
    [Trait("Category", "TASK8_DELTA_CERTIFICATION")]
    [InlineData(0.0, false, .001)]
    [InlineData(0.1, false, .001)]
    [InlineData(0.0, true, .001)]
    [InlineData(0.1, true, .001)]
    [InlineData(0.0, false, 1.0)]
    [InlineData(0.1, false, 1.0)]
    [InlineData(0.0, true, 1.0)]
    [InlineData(0.1, true, 1.0)]
    public void Path6_UntrustedGlobalOrTrackedLocator_HasNoDefinitiveCharacter(double delta, bool residualOnly, double searchTolerance)
    {
        var result = XxzCoherenceBlock.TrackDiabolicUnderDelta(7, new Complex(.6788, 0),
            new Complex(-4.557, 0), delta, residualOnly: residualOnly, coalesceTol: searchTolerance);
        Console.WriteLine($"N7 unsafe locator delta={delta:R}, residual={residualOnly}, verdict={result.Verdict}, alg={result.Algebraic}, gap={result.Gap:R}");
        Assert.Equal(XxzCoherenceBlock.DeltaFlipVerdict.Uncertified, result.Verdict);
        Assert.Equal(0, result.Algebraic);
        Assert.Equal(0, result.Geometric);
        Assert.Null(result.Survived);
    }

    [Theory]
    [Trait("Category", "TASK8_DELTA_CERTIFICATION")]
    [InlineData(.001)]
    [InlineData(1.0)]
    public void Path5_TrackedNearSplit_IsUncertifiedRegardlessOfSearchTolerance(double searchTolerance)
    {
        var result = XxzCoherenceBlock.TrackDiabolicUnderDelta(6, new Complex(.7581, .260),
            new Complex(-5.392, 1.653), .02, residualOnly: true, coalesceTol: searchTolerance);
        Console.WriteLine($"N6 tracked near split tol={searchTolerance:R}, verdict={result.Verdict}, gap={result.Gap:R}");
        Assert.True(result.Gap > 1e-6 && result.Gap < .001);
        Assert.Equal(XxzCoherenceBlock.DeltaFlipVerdict.Uncertified, result.Verdict);
        Assert.Null(result.Survived);
    }

    [Fact]
    [Trait("Category", "TASK8_DELTA_CERTIFICATION")]
    public void CharacterAtDiabolicNear_RejectsUnrelatedLambdaSeed()
    {
        var result = XxzCoherenceBlock.CharacterAtDiabolicNear(4, 0,
            new Complex(GaloisMonodromyWitness.QEp, 0), new Complex(100, 0));
        Assert.Equal(0, result.Algebraic); // Unavailable character, not a borrowed unrelated pair.
        Assert.Equal(0, result.Geometric);
        Assert.True(double.IsNaN(result.Departure));
        Assert.Equal(XxzCoherenceBlock.DeltaFlipVerdict.Uncertified, result.Verdict);
        Assert.Null(result.Survived);
    }

    [Fact]
    [Trait("Category", "TASK8_DELTA_CERTIFICATION")]
    public void DefaultTracker_UnrelatedLambdaSeed_CannotBorrowGlobalMinimumCharacter()
    {
        var result = XxzCoherenceBlock.TrackDiabolicUnderDelta(5, new Complex(.6407, .180),
            new Complex(100, 100), 0);
        Assert.Equal(XxzCoherenceBlock.DeltaFlipVerdict.Uncertified, result.Verdict);
        Assert.Equal(0, result.Algebraic);
    }

    private static Complex[] SpectrumOf(Matrix<Complex> m)
        => m.Evd().EigenValues.ToArray();

    // multiset agreement (every a has a b within tol, and vice versa), robust to ordering, unlike sort-by-Re
    // which is unstable across the Re=-4 cluster (many eigenvalues share Re to float noise).
    private static void AssertSameSpectrum(Complex[] a, Complex[] b, double tol)
    {
        Assert.Equal(a.Length, b.Length);
        foreach (var x in a) Assert.True(b.Min(y => (x - y).Magnitude) < tol, $"{x} has no match in b");
        foreach (var y in b) Assert.True(a.Min(x => (x - y).Magnitude) < tol, $"{y} has no match in a");
    }

    // ΔTask 1: the (q,Δ) builder at Δ=0 IS the XY (SE,DE) block; at N=4 its R=+1 symmetric-sector spectrum
    // must equal F89Path3OcticBlock.BuildSeDeSymBlock(q, γ=1) to machine precision (the trusted anchor).
    [Fact]
    public void XxzBlock_Delta0_MatchesF89SymBlockSpectrum_AtQ2()
    {
        var xxz = XxzCoherenceBlock.SeDeSymSpectrum(4, new Complex(2.0, 0), 0.0);
        var f89 = SpectrumOf(F89Path3OcticBlock.BuildSeDeSymBlock(2.0, 1.0));   // (j=2, γ=1)
        AssertSameSpectrum(xxz, f89, 1e-9);
    }

    // The Delta=0 N4 control passes the full certificate. The two nonzero-Delta proposals
    // remain split above the independent coincidence bound, so neither Jordan character nor
    // departure is available. This is a locator boundary, not a non-existence or lifting result.
    [Fact]
    public void N4_DeltaProposals_RequireStrictFullPairCoincidence()
    {
        double qEp = GaloisMonodromyWitness.QEp;                 // sqrt((-1+sqrt13)/6) ~ 0.658983
        var lamEp = new Complex(-4, 2 * qEp);

        var r0 = XxzCoherenceBlock.CharacterAtDiabolicNear(4, 0.0, new Complex(qEp, 0), lamEp);
        Assert.Equal(2, r0.Algebraic);
        Assert.Equal(2, r0.Geometric);                          // geo==alg => DIABOLIC (semisimple)
        Assert.Equal(XxzCoherenceBlock.DeltaFlipVerdict.Diabolic, r0.Verdict);
        Assert.True(r0.Gap <= XxzCoherenceBlock.FullBlockCoincidenceTolerance);
        Assert.True(r0.Departure < 1e-6, $"d=0 departure {r0.Departure} should be ~0");

        var r2 = XxzCoherenceBlock.CharacterAtDiabolicNear(4, 0.02, new Complex(qEp, 0), lamEp);
        var r10 = XxzCoherenceBlock.CharacterAtDiabolicNear(4, 0.10, new Complex(qEp, 0), lamEp);
        Console.WriteLine($"N4 strict Delta=0: {r0}");
        Console.WriteLine($"N4 strict Delta=.02: {r2}");
        Console.WriteLine($"N4 strict Delta=.10: {r10}");
        foreach (var (delta, result) in new[] { (.02, r2), (.10, r10) })
        {
            AssertUncertifiedSplit(4, delta, result);
            Assert.Null(result.Survived);
            Assert.False(result.IsCertifiedDiabolic);
        }
    }

    // Residual-only Δ-test (the N>=6 fix): the full sym spectrum floods on AT-locked degeneracies, so the box
    // scan in TrackDiabolicUnderDelta captures AT crossings at N=6 (q_candidate jumps, false LIFTs). ResidualRootsTrackedXxz
    // uses base-label nearest-neighbour continuation from (q0=2, Δ=0), where the labels equal the F89 residual.
    // This is proposal tracking, not an invariant-set proof; full-block character remains required.
    // The base comparison pins the XxzCoherenceBlock-vs-F89 spectrum convention at path-5.
    [Fact]
    public void ResidualRootsTrackedXxz_Path5_MatchesLocatorResidual_AtBase()
    {
        var xxz = XxzCoherenceBlock.ResidualRootsTrackedXxz(5, new Complex(2, 0), 0.0)
            .OrderBy(z => z.Real).ThenBy(z => z.Imaginary).ToArray();
        var loc = PathKMonodromyScout.ResidualRootsAt(5, new Complex(2, 0))
            .OrderBy(z => z.Real).ThenBy(z => z.Imaginary).ToArray();
        Assert.Equal(32, xxz.Length);                            // F_d degree for path-5 (not 45 = full sym block)
        Assert.Equal(loc.Length, xxz.Length);
        for (int i = 0; i < xxz.Length; i++)
            Assert.True((xxz[i] - loc[i]).Magnitude < 1e-7, $"strand {i}: xxz {xxz[i]} vs locator {loc[i]}");
    }

    // N6 tracked proposals retain the Delta=0 diabolic controls. The positive-Delta search
    // does not reach strict full-block coincidence; those readings are unknown, not defect/lift evidence.
    [Fact]
    public void Path5_ResidualProposals_RequireIndependentFullPairCertification()
    {
        var diabolics = new[]
        {
            (q: new Complex(0.7090, -0.219), lam: new Complex(-4.151, 1.615)),   // rung-near (the full-block Δ-test captured AT here)
            (q: new Complex(0.7581, 0.260), lam: new Complex(-5.392, 1.653)),    // rung-far
            (q: new Complex(1.0561, 0.238), lam: new Complex(-5.187, 1.441)),    // rung-far
        };
        foreach (var (q, lam) in diabolics)
        {
            var d0 = XxzCoherenceBlock.TrackDiabolicUnderDelta(6, q, lam, 0.0, residualOnly: true);
            Assert.Equal(XxzCoherenceBlock.DeltaFlipVerdict.Diabolic, d0.Verdict);              // Gate 0
            Assert.True((d0.QCandidate - q).Magnitude < 0.03, $"Δ=0 q_candidate={d0.QCandidate} must stay at the residual diabolic {q}, not jump to an AT crossing");
            Assert.True(d0.Gap > 1e-13, $"the residual diabolic has a finite gap (got {d0.Gap:E2}); an AT capture is exactly 0");
            var d = XxzCoherenceBlock.TrackDiabolicUnderDelta(6, q, lam, 0.1, residualOnly: true);
            AssertUncertifiedSplit(6, .1, d);
            Console.WriteLine($"N6 tracked qSeed={q}, Delta=.1, gap={d.Gap:R}, verdict={d.Verdict}");
        }

        // The rung-far proposal is also split: no Jordan character is supplied by a near pair.
        var clean = XxzCoherenceBlock.TrackDiabolicUnderDelta(6, new Complex(0.7581, 0.260), new Complex(-5.392, 1.653), 0.1, residualOnly: true);
        AssertUncertifiedSplit(6, .1, clean);
    }

    // N5 Delta=0 diabolics are positive controls. Neither their positive-Delta proposals nor
    // the approximate defective-control seeds reach strict full-pair coincidence with this locator.
    [Fact]
    public void Path4_ProposalsAndDefectiveControl_RequireIndependentFullPairCertification()
    {
        var diabolics = new[]
        {
            (q: new Complex(0.6407, 0.180), lam: new Complex(-4.077, -1.115)),   // clean
            (q: new Complex(0.7654, 0.024), lam: new Complex(-4.371, -2.056)),   // near-real
            (q: new Complex(1.9447, 1.217), lam: new Complex(-2.455, -3.473)),   // far
        };
        foreach (var (q, lam) in diabolics)
        {
            var d0 = XxzCoherenceBlock.TrackDiabolicUnderDelta(5, q, lam, 0.0);
            Assert.Equal(XxzCoherenceBlock.DeltaFlipVerdict.Diabolic, d0.Verdict);   // Gate 0
            var d = XxzCoherenceBlock.TrackDiabolicUnderDelta(5, q, lam, 0.05);
            AssertUncertifiedSplit(5, .05, d);
            Console.WriteLine($"N5 default qSeed={q}, Delta=.05, gap={d.Gap:R}, verdict={d.Verdict}");
        }

        // The smaller Delta proposal is still split at the strict full-block tolerance.
        var clean02 = XxzCoherenceBlock.TrackDiabolicUnderDelta(5, new Complex(0.6407, 0.180), new Complex(-4.077, -1.115), 0.02);
        AssertUncertifiedSplit(5, .02, clean02);
        Console.WriteLine($"N5 default clean Delta=.02, gap={clean02.Gap:R}, verdict={clean02.Verdict}");

        // Approximate defective seeds are no substitute for independently reached coincidence.
        var ctrlQ = new Complex(0.9938, 0.183);
        var ctrlLam = new Complex(-4.712, 0.824);
        AssertUncertifiedSplit(5, 0, XxzCoherenceBlock.TrackDiabolicUnderDelta(5, ctrlQ, ctrlLam, 0));
        AssertUncertifiedSplit(5, .05, XxzCoherenceBlock.TrackDiabolicUnderDelta(5, ctrlQ, ctrlLam, .05));
    }

    private static void AssertUncertifiedSplit(int n, double delta, XxzCoherenceBlock.DeltaTrackResult result)
    {
        var pair = XxzCoherenceBlock.SeDeSymSpectrum(n, result.QCandidate, delta)
            .OrderBy(z => (z - result.LambdaCandidate).Magnitude).Take(2).ToArray();
        Assert.True((pair[0] - pair[1]).Magnitude > 1e-6);
        Assert.True(Math.Abs((pair[0] - pair[1]).Magnitude - result.Gap) < 1e-10);
        Assert.Equal(XxzCoherenceBlock.DeltaFlipVerdict.Uncertified, result.Verdict);
        Assert.Equal(0, result.Algebraic);
        Assert.Equal(0, result.Geometric);
        Assert.True(double.IsNaN(result.Departure));
        Assert.Null(result.Survived);
    }

    // ---- Fixed-complement compressed XXZ proposals ----
    // The tracked residual path (ResidualRootsTrackedXxz) breaks at k=6 (the F_53 AT-degeneracy flood: leaked
    // AT strands wearing residual labels, min-gap zero across the box), exactly as the XY tracker did before
    // the Delta=0 exact restriction. The XXZ port supplies compressed proposals, not residual eigenvalues:
    // M_xxz(q,Δ) = (A + qC + qΔ·G)/2 in F89's mirror basis (G the ZZ-
    // frequency generator, diagonal −2i·zzDiag), compressed onto the q-independent AT complement U_res. These
    // gates validate the matrix port; full-block coincidence/correspondence and character remain required.

    // Foundation: the XXZ block built in F89's mirror basis (A + qC + qΔ·G) has the SAME spectrum as the
    // independently-constructed XxzCoherenceBlock.BuildSym at Δ≠0. Pins the ZZ generator (the −2i·zzDiag
    // diagonal and its ×2-cleared scaling: a reflection-invariant diagonal contributes 2·value to 2M for both
    // orbit lengths) against the trusted block, with NO compression involved.
    [Theory]
    [InlineData(5)]   // N=6
    [InlineData(6)]   // N=7
    public void AllRootsXxz_F89Basis_MatchesXxzBlockSpectrum_AtDelta(int k)
    {
        var q = new Complex(2.0, 0.3);
        const double delta = 0.1;
        var f89 = PathKMonodromyScout.AllRootsXxz(k, q, delta);
        var xxz = XxzCoherenceBlock.SeDeSymSpectrum(k + 1, q, delta);
        AssertSameSpectrum(f89, xxz, 1e-7);
    }

    // At Δ=0 the exact XXZ residual roots ARE the XY exact residual roots (G drops out): a wiring guard that
    // the port does not perturb the established Δ=0 science.
    [Theory]
    [InlineData(5)]
    [InlineData(6)]
    public void ResidualRootsCompressedXxz_Delta0_EqualsXyExactResidual(int k)
    {
        var q = new Complex(1.7, 0.0);
        var xxz = PathKMonodromyScout.ResidualRootsCompressedXxz(k, q, 0.0);
        var xy = PathKMonodromyScout.ResidualRootsExact(k, q);
        AssertSameSpectrum(xxz, xy, 1e-9);
    }

    // Compression gives 32/53 distinct proposal roots in this generic sample. It is not an invariant
    // residual subspace at Delta!=0 and does not certify a subset of the full spectrum or a perturbation bound.
    [Theory]
    [InlineData(5, 32)]
    [InlineData(6, 53)]
    public void ResidualRootsCompressedXxz_HasDistinctProposalRoots_AtDelta(int k, int fdDegree)
    {
        var q = new Complex(2.0, 0.3);
        const double delta = 0.1;
        var res = PathKMonodromyScout.ResidualRootsCompressedXxz(k, q, delta);
        Assert.Equal(fdDegree, res.Length);
        Assert.True(PathKMonodromyScout.MinGap(res) > 1e-6,
            $"the sampled compressed roots must be distinct; got min gap {PathKMonodromyScout.MinGap(res):E2}");
    }

    // Delta=0 remains calibrated. At Delta!=0 a compressed proposal is not a replacement for
    // the tracked full-spectrum locator, even at N=6; its own full-pair gap decides certification.
    [Fact]
    public void Path5_CompressedProposals_RequireFullPairBeforeCharacter()
    {
        var diabolics = new[]
        {
            (q: new Complex(0.7090, -0.219), lam: new Complex(-4.151, 1.615)),
            (q: new Complex(0.7581, 0.260), lam: new Complex(-5.392, 1.653)),
            (q: new Complex(1.0561, 0.238), lam: new Complex(-5.187, 1.441)),
        };
        foreach (var (q, lam) in diabolics)
        {
            var d0 = XxzCoherenceBlock.TrackDiabolicUnderDelta(6, q, lam, 0.0, exact: true);
            Assert.Equal(XxzCoherenceBlock.DeltaFlipVerdict.Diabolic, d0.Verdict);
            Assert.True((d0.QCandidate - q).Magnitude < 0.03, $"Δ=0 q_candidate={d0.QCandidate} must stay at the diabolic {q}, not jump");
            Assert.True(d0.Gap > 1e-13, $"the residual diabolic has a finite gap (got {d0.Gap:E2}); an AT capture is exactly 0");
            var d = XxzCoherenceBlock.TrackDiabolicUnderDelta(6, q, lam, 0.1, exact: true);
            var pair = XxzCoherenceBlock.SeDeSymSpectrum(6, d.QCandidate, 0.1)
                .OrderBy(z => (z - d.LambdaCandidate).Magnitude).Take(2).ToArray();
            double fullGap = (pair[0] - pair[1]).Magnitude;
            Assert.True(Math.Abs(d.Gap - fullGap) < 1e-10);
            Assert.Equal(fullGap > 1e-3, d.Verdict == XxzCoherenceBlock.DeltaFlipVerdict.Uncertified);
            Console.WriteLine($"N6 compressed proposal: qSeed={q}, fullGap={fullGap:R}, verdict={d.Verdict}");
        }

        // Wrong-input control: the tracked locator proposes this point, but the strict full-block check leaves it Uncertified.
        var clean = XxzCoherenceBlock.TrackDiabolicUnderDelta(6, new Complex(0.7581, 0.260), new Complex(-5.392, 1.653), 0.1, exact: true);
        Assert.Equal(XxzCoherenceBlock.DeltaFlipVerdict.Uncertified, clean.Verdict);
        Assert.True(clean.Gap > 1e-3);
        Assert.Equal(0, clean.Geometric);
        Assert.True(double.IsNaN(clean.Departure));
    }

    [Theory]
    [Trait("Category", "TASK8_DELTA_CERTIFICATION")]
    [InlineData(0.0, 0.0, true)]
    [InlineData(0.0001, 0.0, false)] // Close but split full eigenvalues.
    [InlineData(0.0, 0.01, false)] // Coincident full pair unrelated to proposed midpoint.
    public void CompressedCertification_RequiresCoincidenceAndCorrespondence(double split, double midpoint, bool certified)
    {
        var matrix = MathNet.Numerics.LinearAlgebra.Matrix<Complex>.Build.DenseDiagonal(3, 3,
            i => i == 2 ? 1 : i * split);
        var result = XxzCoherenceBlock.CertifyFullBlockProposal(matrix, Complex.Zero, new Complex(midpoint, 0));
        Assert.Equal(certified, result.IsCertifiedDiabolic);
        Assert.Equal(certified ? XxzCoherenceBlock.DeltaFlipVerdict.Diabolic : XxzCoherenceBlock.DeltaFlipVerdict.Uncertified, result.Verdict);
        Assert.True(Math.Abs(result.Gap - split) < 1e-12);
    }

    [Theory]
    [Trait("Category", "TASK8_DELTA_CERTIFICATION")]
    [InlineData(RCPsiSquared.Core.Numerics.EpCharacter.EpKind.Diabolic, XxzCoherenceBlock.DeltaFlipVerdict.Diabolic)]
    [InlineData(RCPsiSquared.Core.Numerics.EpCharacter.EpKind.Defective, XxzCoherenceBlock.DeltaFlipVerdict.Defective)]
    [InlineData(RCPsiSquared.Core.Numerics.EpCharacter.EpKind.Normal, XxzCoherenceBlock.DeltaFlipVerdict.Uncertified)]
    [InlineData(RCPsiSquared.Core.Numerics.EpCharacter.EpKind.NearEp, XxzCoherenceBlock.DeltaFlipVerdict.Uncertified)]
    [InlineData((RCPsiSquared.Core.Numerics.EpCharacter.EpKind)99, XxzCoherenceBlock.DeltaFlipVerdict.Uncertified)]
    public void CompressedCharacterMapping_RejectsNoncoalescentKinds(RCPsiSquared.Core.Numerics.EpCharacter.EpKind kind,
        XxzCoherenceBlock.DeltaFlipVerdict expected)
        => Assert.Equal(expected, XxzCoherenceBlock.CertifiedCharacterVerdict(kind));

    [Fact]
    [Trait("Category", "TASK8_DELTA_CERTIFICATION")]
    public void UncertifiedSurvival_IsUnknown_NotFalse()
    {
        var unknown = new XxzCoherenceBlock.DeltaTrackResult(
            XxzCoherenceBlock.DeltaFlipVerdict.Uncertified, 0, 0, double.NaN, Complex.Zero, Complex.Zero, 0.02);
        Assert.Null(unknown.Survived);
        Assert.False(unknown.IsCertifiedDiabolic);
        Assert.True((unknown with { Verdict = XxzCoherenceBlock.DeltaFlipVerdict.Diabolic }).Survived);
        Assert.False((unknown with { Verdict = XxzCoherenceBlock.DeltaFlipVerdict.Defective }).Survived);
        Assert.False((unknown with { Verdict = XxzCoherenceBlock.DeltaFlipVerdict.Lifted }).Survived);
    }

    [Theory]
    [Trait("Category", "TASK8_DELTA_CERTIFICATION")]
    [InlineData(0.0, 0.001)]
    [InlineData(0.0, 1000.0)]
    [InlineData(0.1, 1000.0)]
    public void CompressedProposal_RejectsWrongLambdaSeed(double delta, double proposalTolerance)
    {
        var q = new Complex(0.6788, 0);
        var result = XxzCoherenceBlock.TrackDiabolicUnderDelta(7, q, new Complex(100, 0), delta,
            coalesceTol: proposalTolerance, exact: true);
        Assert.Equal(XxzCoherenceBlock.DeltaFlipVerdict.Uncertified, result.Verdict);
        Assert.Equal(0, result.Algebraic);
        Assert.Equal(0, result.Geometric);
        Assert.Null(result.Survived);
        Assert.True(double.IsNaN(result.Departure));

        // Same public path and q with the actual seed must retain its calibrated Delta=0 pair.
        var valid = XxzCoherenceBlock.TrackDiabolicUnderDelta(7, q, new Complex(-4.557, 0), 0,
            coalesceTol: proposalTolerance, exact: true);
        Assert.Equal(XxzCoherenceBlock.DeltaFlipVerdict.Diabolic, valid.Verdict);
        Assert.Equal(2, valid.Algebraic);
        Assert.Equal(2, valid.Geometric);
        Assert.True(valid.Gap <= XxzCoherenceBlock.FullBlockCoincidenceTolerance);
    }

    [Theory]
    [Trait("Category", "TASK8_DELTA_CERTIFICATION")]
    [InlineData(0.02, 1.0)] // Inflating proposal tolerance must not authorize a split full pair.
    [InlineData(0.0001, 0.001)] // Even a near pair below the old search tolerance is not coincident.
    public void CompressedProposal_SearchToleranceCannotCertifySplitPair(double delta, double proposalTolerance)
    {
        var d = XxzCoherenceBlock.TrackDiabolicUnderDelta(7, new Complex(0.6788, 0),
            new Complex(-4.557, 0), delta, coalesceTol: proposalTolerance, exact: true);
        Console.WriteLine($"strict control Delta={delta:R}, gap={d.Gap:R}, verdict={d.Verdict}");
        Assert.True(d.Gap > 1e-6);
        Assert.True(d.Gap < proposalTolerance);
        Assert.Equal(XxzCoherenceBlock.DeltaFlipVerdict.Uncertified, d.Verdict);
        Assert.Equal(0, d.Algebraic);
        Assert.Equal(0, d.Geometric);
        Assert.True(double.IsNaN(d.Departure));
    }

    // The Delta=0 control is a genuine full-block diabolic. Its compressed Delta proposal
    // must not acquire a Jordan/lift verdict when the nearest full-block pair remains split.
    [Theory]
    [Trait("Category", "TASK8_DELTA_CERTIFICATION")]
    [InlineData(1.1264, -4.942, 0.10)]
    [InlineData(1.3038, -5.171, 0.10)]
    [InlineData(2.6280, -4.343, 0.10)]
    [InlineData(0.6788, -4.557, 0.10)]
    [InlineData(0.6788, -4.557, 0.02)]
    public void Path6_CompressedDeltaProposals_AreUncertifiedWhenFullPairIsSplit(double realQ, double realLambda, double delta)
    {
        var q = new Complex(realQ, 0); var lam = new Complex(realLambda, 0);
        const double coalesceTol = 1e-3;
        var d0 = XxzCoherenceBlock.TrackDiabolicUnderDelta(7, q, lam, 0.0, exact: true);
        Assert.Equal(XxzCoherenceBlock.DeltaFlipVerdict.Diabolic, d0.Verdict);
        Assert.Equal(2, d0.Algebraic);
        Assert.Equal(2, d0.Geometric);
        Assert.True(d0.Gap < coalesceTol);
        var d = XxzCoherenceBlock.TrackDiabolicUnderDelta(7, q, lam, delta, exact: true);
        Assert.Equal("Uncertified", d.Verdict.ToString());
        var fullPair = XxzCoherenceBlock.SeDeSymSpectrum(7, d.QCandidate, delta)
            .OrderBy(z => (z - d.LambdaCandidate).Magnitude).Take(2).ToArray();
        double actualGap = (fullPair[0] - fullPair[1]).Magnitude;
        Assert.True(actualGap > coalesceTol, $"q={realQ}, Delta={delta}, full gap={actualGap:R}");
        Assert.True(Math.Abs(d.Gap - actualGap) < 1e-10);
        Assert.Equal(0, d.Algebraic); // unavailable, not a character reading
        Assert.Equal(0, d.Geometric);
        Assert.True(double.IsNaN(d.Departure));
        Console.WriteLine($"N7 compressed proposal: qSeed={realQ:R}, Delta={delta:R}, qCandidate={d.QCandidate}, fullGap={d.Gap:R}, verdict={d.Verdict}");
    }
}
