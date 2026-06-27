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
}
