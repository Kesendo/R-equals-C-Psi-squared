using System;
using System.Collections.Generic;
using System.Linq;
using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.F89PathK;
using RCPsiSquared.Core.Numerics;
using RCPsiSquared.Diagnostics.Foundation;
using Xunit;
using Xunit.Abstractions;

namespace RCPsiSquared.Diagnostics.Tests.Foundation;

/// <summary>
/// The MEMBER instrument (Task 4 of the sectorbraid large-N program): the spectator W transports
/// the (1,2) near-defective 2-plane up the |p−q̃|=1 band, the §7 fold gauge maps it to the μ blocks,
/// and a Householder 2-plane residual gives a certified FROM-ABOVE upper bound on σ_min at every
/// diamond member, with two sparse matvecs and no solves. These tests pin the four public methods of
/// <see cref="SectorWitnessTransport"/> from below: the decisive one (Test 0) composes the new
/// matrix-free ApplyW with the Task-1 census CSR operators, closing the seam between the two code
/// lineages (BuildW's big-endian JW basis and WeightCoherenceBlock's little-endian Configs basis).
/// </summary>
public sealed class SectorWitnessTransportTests
{
    private readonly ITestOutputHelper _out;
    public SectorWitnessTransportTests(ITestOutputHelper o) => _out = o;

    // ---------------------------------------------------------------------------------------------
    // Test (0): THE decisive from-below gate [F3] — ApplyW intertwines the census's OWN CSR pencil.
    // BuildW and WeightCoherenceBlock have never been composed in this codebase; different lineages
    // and endianness. If this passes, the translation is correct against the Task-1 operators.
    // ---------------------------------------------------------------------------------------------
    [Theory]
    [InlineData(4, 1, 2)]
    [InlineData(5, 1, 2)]
    [InlineData(5, 2, 3)]
    public void ApplyW_Intertwines_TheCensusCsrPencil(int n, int p, int w)
    {
        var q = new Complex(0.7, 0.0);
        var src = WeightCoherenceSectorCsr.BuildFull(n, p, w, q);
        var dst = WeightCoherenceSectorCsr.BuildFull(n, p + 1, w + 1, q);
        var rng = new Random(9);
        var v = new Complex[src.Dim];
        for (int i = 0; i < v.Length; i++) v[i] = new Complex(rng.NextDouble() - 0.5, rng.NextDouble() - 0.5);
        var lv = new Complex[src.Dim];
        CsrOps.Multiply(src, v, lv);
        var wLv = SectorWitnessTransport.ApplyW(n, p, w, lv);
        var wv = SectorWitnessTransport.ApplyW(n, p, w, v);
        var lWv = new Complex[dst.Dim];
        CsrOps.Multiply(dst, wv, lWv);
        double diff = 0, scale = 0;
        for (int i = 0; i < dst.Dim; i++) { diff += (lWv[i] - wLv[i]).Magnitude; scale += lWv[i].Magnitude; }
        _out.WriteLine($"N={n} ({p},{w}): intertwining residual {diff:E3} scale {scale:E3}");
        Assert.True(diff <= 1e-12 * (1 + scale), $"intertwining residual {diff}");
    }

    // ---------------------------------------------------------------------------------------------
    // Test (i): ApplyW vs the dense BuildW THROUGH the explicit endianness bridge. Configs is
    // little-endian (bit s = site s); PopcountStates is big-endian (site 0 = MSB). The two enumerate
    // the SAME ascending masks, so the bridge is the per-index bit-reversal permutation. We assert
    // equality up to ONE global per-block sign and RECORD that sign (measured: +1).
    // ---------------------------------------------------------------------------------------------
    [Theory]
    [InlineData(4, 1, 2)]
    [InlineData(5, 2, 3)]
    public void ApplyW_MatchesBuildW_ThroughEndiannessBridge(int n, int p, int w)
    {
        var rng = new Random(17);
        var kets = WeightCoherenceBlock.Configs(n, p);
        var bras = WeightCoherenceBlock.Configs(n, w);
        int nBra = bras.Count;
        var v = new Complex[kets.Count * nBra];
        for (int i = 0; i < v.Length; i++) v[i] = new Complex(rng.NextDouble() - 0.5, rng.NextDouble() - 0.5);

        // little-endian Configs -> big-endian PopcountStates (same list, reversed bit meaning)
        var vBE = BridgeToBigEndian(n, p, w, v);
        var wMat = SpectatorIntertwiner.BuildW(n, p, w);
        var wBE = (wMat * MathNet.Numerics.LinearAlgebra.Vector<Complex>.Build.DenseOfArray(vBE)).ToArray();
        var bridgedBack = BridgeFromBigEndian(n, p + 1, w + 1, wBE);

        var applyW = SectorWitnessTransport.ApplyW(n, p, w, v);
        Assert.Equal(bridgedBack.Length, applyW.Length);

        // measure the global sign at the largest component, assert global consistency, record it
        int imax = 0; double amax = 0;
        for (int i = 0; i < applyW.Length; i++) if (applyW[i].Magnitude > amax) { amax = applyW[i].Magnitude; imax = i; }
        Complex ratio = bridgedBack[imax] / applyW[imax];
        int sign = (int)Math.Round(ratio.Real);
        _out.WriteLine($"N={n} ({p},{w}): measured global bridge sign = {sign} (ratio {ratio.Real:F6}{ratio.Imaginary:+0.0e+0;-0.0e+0}i)");
        Assert.True(Math.Abs(ratio.Imaginary) < 1e-10 && Math.Abs(Math.Abs(ratio.Real) - 1) < 1e-9,
            $"the bridge relation must be a real global sign; got {ratio}");
        double diff = 0, scale = 0;
        for (int i = 0; i < applyW.Length; i++)
        { diff += (bridgedBack[i] - sign * applyW[i]).Magnitude; scale += applyW[i].Magnitude; }
        Assert.True(diff <= 1e-12 * (1 + scale), $"bridged BuildW must equal sign*ApplyW; residual {diff}");
        Assert.Equal(1, sign);   // the recorded value: the two constructions coincide exactly (sign +1)
    }

    // ---------------------------------------------------------------------------------------------
    // Test (ii): the §7 fold gauge preserves the shifted residual norm with s ↦ −s−2N. This checks
    // the exact P·𝒟 convention against L(p,N−w) = −P𝒟·L(p,w)·𝒟Pᵀ − 2N·I (holds at every q, Δ=0).
    // ---------------------------------------------------------------------------------------------
    public static IEnumerable<object[]> FoldCases()
    {
        foreach (var q in new[] { new Complex(0.620878, 0), new Complex(0.9, 0.3) })
            foreach (var pw in new[] { (5, 1, 2), (5, 2, 3) })
                yield return new object[] { pw.Item1, pw.Item2, pw.Item3, q };
    }

    [Theory]
    [MemberData(nameof(FoldCases))]
    public void FoldGauge_PreservesShiftedResidual(int n, int p, int w, Complex q)
    {
        var s = new Complex(-3.9, 0);
        var foldShift = -s - 2.0 * n;
        var lpw = WeightCoherenceSectorCsr.BuildFull(n, p, w, q);
        var rng = new Random(31);
        var v = new Complex[lpw.Dim];
        for (int i = 0; i < v.Length; i++) v[i] = new Complex(rng.NextDouble() - 0.5, rng.NextDouble() - 0.5);

        double lhs = ShiftedResidualRatio(lpw, s, v);

        var folded = SectorWitnessTransport.ApplyFoldGauge(n, p, w, v);
        var lpNw = WeightCoherenceSectorCsr.BuildFull(n, p, n - w, q);
        double rhs = ShiftedResidualRatio(lpNw, foldShift, folded);

        _out.WriteLine($"N={n} ({p},{w}) q={q}: lhs {lhs:E12} rhs {rhs:E12} |diff| {Math.Abs(lhs - rhs):E3}");
        Assert.True(Math.Abs(lhs - rhs) <= 1e-12 * (1 + lhs), $"fold gauge must preserve the shifted residual; lhs {lhs} rhs {rhs}");
    }

    // ---------------------------------------------------------------------------------------------
    // Test (iii): end-to-end at the N=5 seed q*=0.620878. MemberUpperBounds returns exactly the
    // ExpectedMembers(5) keys; each bound is a deep FROM-ABOVE reading (< memberTol/100) and stays
    // within a factor 2 of the dense reference (>= 0.5·dense). Ratios recorded, never gating.
    // ---------------------------------------------------------------------------------------------
    [Fact]
    public void MemberUpperBounds_N5_EndToEnd()
    {
        var seed = RealDefectiveSeeds.ForN(5).Single(s => Math.Abs(s.QStar - 0.620878) < 1e-6);
        var (qRefined, lambda, pairGap) = SectorShellCensus.RefineSeed(seed);
        var q = new Complex(qRefined, 0);
        var mu = -lambda - 2.0 * 5;
        double memberTol = 10.0 * pairGap;

        var bounds = SectorWitnessTransport.MemberUpperBounds(seed, m => _out.WriteLine(m));
        var expected = SectorShellCensus.ExpectedMembers(5);

        Assert.True(bounds.Keys.ToHashSet().SetEquals(expected),
            $"keys must equal ExpectedMembers(5). got [{string.Join(",", bounds.Keys)}] expected [{string.Join(",", expected)}]");

        foreach (var (key, val) in bounds)
        {
            var shift = key.Shift == "lambdaA" ? lambda : mu;
            double dense = DenseSectorSigma(5, key.P, key.W, val.CarriedOdd, q, shift);
            double ratio = dense > 0 ? val.Bound / dense : double.NaN;
            _out.WriteLine($"cell ({key.P},{key.W}) {key.Shift} carriedOdd={val.CarriedOdd}: bound {val.Bound:E4} dense {dense:E4} ratio {ratio:F3} memberTol/100 {memberTol / 100:E3}");
            Assert.True(val.Bound < memberTol / 100,
                $"({key.P},{key.W}) {key.Shift}: bound {val.Bound:E4} must be < memberTol/100 {memberTol / 100:E4}");
            Assert.True(val.Bound >= 0.5 * dense,
                $"({key.P},{key.W}) {key.Shift}: bound {val.Bound:E4} must be >= 0.5*dense {0.5 * dense:E4}");
        }
    }

    // ---------------------------------------------------------------------------------------------
    // Test (iv): the alternation-law parity landing. R·W = (−1)^{p+q̃}·W·R, so one W-step off the
    // odd-(p+q̃) band block (1,2) FLIPS the carried R-parity; the fold preserves it. Run both seed
    // parities at N=5.
    // ---------------------------------------------------------------------------------------------
    [Theory]
    [InlineData(0.620878, +1)]   // R-even seed: (1,2)->(2,3) lands R-odd
    [InlineData(0.643037, -1)]   // R-odd  seed: (1,2)->(2,3) lands R-even
    public void AlternationLaw_ParityLanding(double qStar, int rParity)
    {
        const int n = 5;
        var seed = RealDefectiveSeeds.ForN(n).Single(s => Math.Abs(s.QStar - qStar) < 1e-6 && s.RParity == rParity);
        var (v1, _, _, _, _) = SectorWitnessTransport.WitnessPlane(seed);

        // predicted carried parity after one W-step to (2,3): seedOdd XOR ((p-1) mod 2), p=2 => flip
        bool seedOdd = seed.RParity < 0;
        bool predictOdd = seedOdd ^ true;

        var at23 = SectorWitnessTransport.ApplyW(n, 1, 2, v1);
        double frac23 = ParityFraction(n, 2, 3, at23, predictOdd);
        _out.WriteLine($"seed q*={qStar} R={rParity}: (2,3) fraction in predicted parity(odd={predictOdd}) = {frac23:F6}");
        Assert.True(frac23 >= 0.999, $"W-step must land >=99.9% in the alternation-predicted parity; got {frac23:F6}");

        // the fold (2,3)->(2,2) preserves the measured parity
        var at22 = SectorWitnessTransport.ApplyFoldGauge(n, 2, 3, at23);
        double frac22 = ParityFraction(n, 2, 2, at22, predictOdd);
        _out.WriteLine($"seed q*={qStar}: (2,2) fraction in same parity(odd={predictOdd}) = {frac22:F6}");
        Assert.True(frac22 >= 0.999, $"the fold must preserve the parity; got {frac22:F6}");
    }

    // ---------------------------------------------------------------------------------------------
    // Test (v): TwoPlaneResidual sanity on a hand-checkable near-defective 2x2 Jordan block
    // [[0,1],[ε²,0]] (σ_min = ε²). The Householder route must resolve the quadratic depth; the Gram
    // route would read ~sqrt of a lost-to-rounding difference (~3e-8) and this test kills any
    // silent regression to it.
    // ---------------------------------------------------------------------------------------------
    [Fact]
    public void TwoPlaneResidual_ResolvesQuadraticDepth()
    {
        double eps = 1e-6;
        double eps2 = eps * eps;   // 1e-12
        var jordan = new WeightCoherenceSectorCsr.Csr(
            2, new[] { 0, 1, 2 }, new[] { 1, 0 }, new[] { Complex.One, new Complex(eps2, 0) });
        var e1 = new Complex[] { Complex.One, Complex.Zero };
        var e2 = new Complex[] { Complex.Zero, Complex.One };
        double val = SectorWitnessTransport.TwoPlaneResidual(jordan, Complex.Zero, e1, e2);
        _out.WriteLine($"Jordan [[0,1],[e2,0]] eps2={eps2:E3}: TwoPlaneResidual = {val:E6}");
        Assert.True(Math.Abs(val - eps2) <= 1e-3 * eps2, $"expected sigma_min = eps^2 = {eps2:E3}, got {val:E3}");
    }

    // ---- helpers ----

    private static double ShiftedResidualRatio(WeightCoherenceSectorCsr.Csr m, Complex shift, Complex[] v)
    {
        var y = new Complex[m.Dim];
        CsrOps.MultiplyShifted(m, shift, v, y);
        double ny = 0, nv = 0;
        for (int i = 0; i < m.Dim; i++) { ny += y[i].Real * y[i].Real + y[i].Imaginary * y[i].Imaginary; }
        for (int i = 0; i < v.Length; i++) { nv += v[i].Real * v[i].Real + v[i].Imaginary * v[i].Imaginary; }
        return Math.Sqrt(ny) / Math.Sqrt(nv);
    }

    private static double DenseSectorSigma(int n, int p, int w, bool odd, Complex q, Complex shift)
    {
        var (a, d) = WeightCoherenceBlock.BuildReflectionSectorColumnMajor(n, p, w, q, odd);
        for (int i = 0; i < d; i++) a[(long)i * d + i] -= shift;
        return ShiftedSigmaMin.EstimateColumnMajor(a, d).SigmaMin;
    }

    /// <summary>Fraction of ‖v‖² sitting in the requested R-parity (odd = −1 eigenspace of the
    /// reflection). even part = (v + Rv)/2, odd part = (v − Rv)/2.</summary>
    private static double ParityFraction(int n, int p, int w, Complex[] v, bool odd)
    {
        var perm = WeightCoherenceBlock.ReflectionPermutation(n, p, w);
        double num = 0, den = 0;
        for (int i = 0; i < v.Length; i++)
        {
            Complex rv = v[perm[i]];                       // (Rv)[i] = v[perm[i]] (perm is an involution)
            Complex part = odd ? (v[i] - rv) / 2.0 : (v[i] + rv) / 2.0;
            num += part.Real * part.Real + part.Imaginary * part.Imaginary;
            den += v[i].Real * v[i].Real + v[i].Imaginary * v[i].Imaginary;
        }
        return num / den;
    }

    private static int RevBits(int c, int n)
    {
        int r = 0;
        for (int s = 0; s < n; s++) if (((c >> s) & 1) != 0) r |= 1 << (n - 1 - s);
        return r;
    }

    private static Complex[] BridgeToBigEndian(int n, int p, int w, Complex[] v)
    {
        var kets = WeightCoherenceBlock.Configs(n, p);
        var bras = WeightCoherenceBlock.Configs(n, w);
        int nBra = bras.Count;
        var ketPos = new Dictionary<int, int>(); for (int i = 0; i < kets.Count; i++) ketPos[kets[i]] = i;
        var braPos = new Dictionary<int, int>(); for (int i = 0; i < bras.Count; i++) braPos[bras[i]] = i;
        var vBE = new Complex[v.Length];
        for (int t = 0; t < v.Length; t++)
        {
            int kc = kets[t / nBra], bc = bras[t % nBra];
            int beK = ketPos[RevBits(kc, n)], beB = braPos[RevBits(bc, n)];
            vBE[beK * nBra + beB] = v[t];
        }
        return vBE;
    }

    private static Complex[] BridgeFromBigEndian(int n, int p, int w, Complex[] vBE)
    {
        var kets = WeightCoherenceBlock.Configs(n, p);
        var bras = WeightCoherenceBlock.Configs(n, w);
        int nBra = bras.Count;
        var ketPos = new Dictionary<int, int>(); for (int i = 0; i < kets.Count; i++) ketPos[kets[i]] = i;
        var braPos = new Dictionary<int, int>(); for (int i = 0; i < bras.Count; i++) braPos[bras[i]] = i;
        var v = new Complex[vBE.Length];
        for (int t = 0; t < v.Length; t++)
        {
            int kc = kets[t / nBra], bc = bras[t % nBra];
            int beK = ketPos[RevBits(kc, n)], beB = braPos[RevBits(bc, n)];
            v[t] = vBE[beK * nBra + beB];
        }
        return v;
    }
}
