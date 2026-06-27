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

    // multiset agreement (every a has a b within tol, and vice versa) — robust to ordering, unlike sort-by-Re
    // which is unstable across the Re=−4 cluster (many eigenvalues share Re to float noise).
    private static void AssertSameSpectrum(Complex[] a, Complex[] b, double tol)
    {
        Assert.Equal(a.Length, b.Length);
        foreach (var x in a) Assert.True(b.Min(y => (x - y).Magnitude) < tol, $"{x} has no match in b");
        foreach (var y in b) Assert.True(a.Min(x => (x - y).Magnitude) < tol, $"{y} has no match in a");
    }

    // ΔTask 1 — the (q,Δ) builder at Δ=0 IS the XY (SE,DE) block: at N=4 its R=+1 symmetric-sector spectrum
    // must equal F89Path3OcticBlock.BuildSeDeSymBlock(q, γ=1) to machine precision (the trusted anchor).
    [Fact]
    public void XxzBlock_Delta0_MatchesF89SymBlockSpectrum_AtQ2()
    {
        var xxz = XxzCoherenceBlock.SeDeSymSpectrum(4, new Complex(2.0, 0), 0.0);
        var f89 = SpectrumOf(F89Path3OcticBlock.BuildSeDeSymBlock(2.0, 1.0));   // (j=2, γ=1)
        AssertSameSpectrum(xxz, f89, 1e-9);
    }

    // ΔTask 2 / Gate A — the C# (q,Δ) tooling MUST reproduce the committed N=4 Δ-flip table
    // (DIABOLIC_BY_INTEGRABILITY.md:32-39 / f89_zz_break_gate.py): Δ=0 diabolic (g1=2, dep≈0) → Δ>0 defective
    // (g1=1, dep≈0.022 at Δ=0.02, ≈0.112 at Δ=0.10). The quantitative dep-bands are the only guard on the ZZ
    // convention (R-3) — keep them (a 4× ZZ error still flips qualitatively but misses the dep calibration).
    [Fact]
    public void N4_DeltaFlip_ReproducesCommittedTable()
    {
        double qEp = GaloisMonodromyWitness.QEp;                 // = sqrt((-1+sqrt13)/6) ≈ 0.658983
        var lamEp = new Complex(-4, 2 * qEp);                    // −4γ + 2iJ
        // The load-bearing discriminant is geo vs alg (R-3), NOT EpCharacter.Kind: Kind labels a small-dep
        // Jordan block (geo<alg but dep/‖A‖<1e-2) as "Normal", masking the defect. geo<alg IS defective.
        var r0 = XxzCoherenceBlock.CharacterAtDiabolicNear(4, 0.0, new Complex(qEp, 0), lamEp);
        Assert.Equal(2, r0.Algebraic);
        Assert.Equal(2, r0.Geometric);                          // geo==alg ⟹ DIABOLIC (semisimple)
        Assert.True(r0.Departure < 1e-6, $"Δ=0 departure {r0.Departure} should be ≈0");

        var r2 = XxzCoherenceBlock.CharacterAtDiabolicNear(4, 0.02, new Complex(qEp, 0), lamEp);
        Assert.Equal(2, r2.Algebraic);
        Assert.Equal(1, r2.Geometric);                          // geo<alg ⟹ DEFECTIVE (Jordan)
        Assert.InRange(r2.Departure, 0.012, 0.030);             // table: 0.022 (mine 0.018, locator-shifted q*)

        var r10 = XxzCoherenceBlock.CharacterAtDiabolicNear(4, 0.10, new Complex(qEp, 0), lamEp);
        Assert.Equal(1, r10.Geometric);
        Assert.InRange(r10.Departure, 0.08, 0.15);              // table: 0.112
    }
}
