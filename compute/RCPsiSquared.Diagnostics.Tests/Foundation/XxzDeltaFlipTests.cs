using System;
using System.Linq;
using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.F89PathK;
using RCPsiSquared.Core.Numerics;
using RCPsiSquared.Diagnostics.Foundation;
using Xunit;

namespace RCPsiSquared.Diagnostics.Tests.Foundation;

/// <summary>The F89 path-4 Δ-test: does XXZ anisotropy Δ≠0 flip the N=5 diabolics (defective / lifted), as
/// it flips the N=4 diabolic defective? Validates the (q,Δ) XxzCoherenceBlock against the existing F89 block
/// at Δ=0, reproduces the committed N=4 Δ-flip table (Gate A), then tracks path-4 diabolics under Δ (Gate B).
/// See docs/superpowers/plans/2026-06-27-f89-path4-delta-test.md.</summary>
public class XxzDeltaFlipTests
{
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

    // ΔTask 2 / Gate A: the C# (q,Δ) tooling reproduces the committed N=4 Δ-flip table
    // (DIABOLIC_BY_INTEGRABILITY.md:32-39 / f89_zz_break_gate.py): Δ=0 diabolic (geo=alg=2, dep≈0) -> Δ>0
    // defective (geo=1<alg=2, dep growing ~linearly). The discriminant is geo vs alg (R-3), NOT EpCharacter.Kind
    // (which labels a small-departure Jordan block as "Normal"); the quantitative dep-bands guard the ZZ
    // convention (a 4x error flips qualitatively but misses the calibration).
    [Fact]
    public void N4_DeltaFlip_ReproducesCommittedTable()
    {
        double qEp = GaloisMonodromyWitness.QEp;                 // sqrt((-1+sqrt13)/6) ~ 0.658983
        var lamEp = new Complex(-4, 2 * qEp);

        var r0 = XxzCoherenceBlock.CharacterAtDiabolicNear(4, 0.0, new Complex(qEp, 0), lamEp);
        Assert.Equal(2, r0.Algebraic);
        Assert.Equal(2, r0.Geometric);                          // geo==alg => DIABOLIC (semisimple)
        Assert.True(r0.Departure < 1e-6, $"d=0 departure {r0.Departure} should be ~0");

        var r2 = XxzCoherenceBlock.CharacterAtDiabolicNear(4, 0.02, new Complex(qEp, 0), lamEp);
        Assert.Equal(2, r2.Algebraic);
        Assert.Equal(1, r2.Geometric);                          // geo<alg => DEFECTIVE (Jordan)
        Assert.InRange(r2.Departure, 0.012, 0.030);             // table 0.022 (mine 0.018, locator-shifted q*)

        var r10 = XxzCoherenceBlock.CharacterAtDiabolicNear(4, 0.10, new Complex(qEp, 0), lamEp);
        Assert.Equal(1, r10.Geometric);
        Assert.InRange(r10.Departure, 0.08, 0.15);              // table 0.112
    }

    // Residual-only Δ-test (the N>=6 fix): the full sym spectrum floods on AT-locked degeneracies, so the box
    // scan in TrackDiabolicUnderDelta captures AT crossings at N=6 (q* jumps, false LIFTs). ResidualRootsTrackedXxz
    // tracks the residual SET from the base (q0=2, Δ=0) - where it equals the locator's F89 residual - and the SET
    // is Δ-stable (ZZ Hermitian => AT rate Re λ = −2γ⟨n_XY⟩ is Δ-independent). At the base it must reproduce the
    // locator's residual exactly (this also pins the XxzCoherenceBlock-vs-F89 spectrum convention at path-5).
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

    // Gate B at N=6 (the residual-only Δ-test): the full-block box scan captures AT crossings at N=6 (q* jumps,
    // gap exactly 0, false LIFTs); residualOnly tracks the residual SET, so each path-5 diabolic reads cleanly
    // DIABOLIC at Δ=0 (q* at the seed, finite gap) and flips DEFECTIVE / LIFTS at Δ>0 - integrability-protection
    // confirmed off N=4 at N=6.
    [Fact]
    public void Path5_Diabolics_DieUnderDelta_ResidualOnly()
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
            Assert.True((d0.QStar - q).Magnitude < 0.03, $"Δ=0 q*={d0.QStar} must stay at the residual diabolic {q}, not jump to an AT crossing");
            Assert.True(d0.Gap > 1e-13, $"the residual diabolic has a finite gap (got {d0.Gap:E2}); an AT capture is exactly 0");
            var d = XxzCoherenceBlock.TrackDiabolicUnderDelta(6, q, lam, 0.1, residualOnly: true);
            Assert.False(d.Survived, $"diabolic {q} must defect or lift at Δ=0.1 (got {d.Verdict}, dep={d.Departure})");
        }

        // the clean rung-far diabolic's Jordan flip is explicit: geo 2->1, departure on (the integrability signature).
        var clean = XxzCoherenceBlock.TrackDiabolicUnderDelta(6, new Complex(0.7581, 0.260), new Complex(-5.392, 1.653), 0.1, residualOnly: true);
        Assert.Equal(XxzCoherenceBlock.DeltaFlipVerdict.Defective, clean.Verdict);
        Assert.Equal(1, clean.Geometric);
        Assert.True(clean.Departure > 0.001, $"the Jordan departure should turn on under Δ (got {clean.Departure})");
    }

    // ΔTask 3 / Gate B: the experiment. Each path-4 diabolic (complex-q) IS a diabolic at Δ=0 (Gate 0) and
    // DOES NOT SURVIVE at Δ>0 (defects to a Jordan EP, or lifts entirely) - the falsifiable synthesis claim.
    // A path-4 DEFECTIVE EP (control) is non-diabolic at Δ=0 and Δ does not make it diabolic: Δ perturbs
    // everything but kills the diabolic CHARACTER specifically. Confirms integrability-protection generalizes
    // to N=5 (DIABOLIC_BY_INTEGRABILITY's gate, off N=4).
    [Fact]
    public void Path4_Diabolics_DieUnderDelta_ControlStaysPut()
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
            Assert.False(d.Survived, $"diabolic {q} must defect or lift at Δ=0.05 (got {d.Verdict}, dep={d.Departure})");
        }

        // the clean diabolic's character flip is explicit: geo 2->1, dep turns on (the Jordan off-diagonal).
        var clean02 = XxzCoherenceBlock.TrackDiabolicUnderDelta(5, new Complex(0.6407, 0.180), new Complex(-4.077, -1.115), 0.02);
        Assert.Equal(XxzCoherenceBlock.DeltaFlipVerdict.Defective, clean02.Verdict);
        Assert.Equal(1, clean02.Geometric);
        Assert.True(clean02.Departure > 0.005, $"the Jordan departure should turn on at Δ=0.02 (got {clean02.Departure})");

        // control: a path-4 DEFECTIVE EP stays non-diabolic at Δ=0 and Δ=0.05 (Δ does not create a diabolic).
        var ctrlQ = new Complex(0.9938, 0.183);
        var ctrlLam = new Complex(-4.712, 0.824);
        Assert.NotEqual(XxzCoherenceBlock.DeltaFlipVerdict.Diabolic,
            XxzCoherenceBlock.TrackDiabolicUnderDelta(5, ctrlQ, ctrlLam, 0.0).Verdict);
        Assert.NotEqual(XxzCoherenceBlock.DeltaFlipVerdict.Diabolic,
            XxzCoherenceBlock.TrackDiabolicUnderDelta(5, ctrlQ, ctrlLam, 0.05).Verdict);
    }
}
