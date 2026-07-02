using System;
using System.Collections.Generic;
using System.Linq;
using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.F89PathK;
using RCPsiSquared.Diagnostics.Foundation;
using Xunit;
using Xunit.Abstractions;

namespace RCPsiSquared.Diagnostics.Tests.Foundation;

/// <summary>The R-ODD deep-loci probe (the sectorbraid arc's immediate next action, premise 2's open
/// half). The fold-resultant certificate Res_λ(F_18, F_corner(−λ−2N)) covers the R-EVEN residual only
/// (F_18 = the ×2-cleared S₂-symmetric sector's H_B-mixed factor); the known deep loci
/// q* = 1.8141 ± 0.3666i are R-even and covered. This gate builds the R-ODD sector of the full (SE,DE)
/// block (<see cref="F89PathKSeDeBlock.ROddBasis"/> + the exact R-odd AT/residual split), validates the
/// split (spectrum bookkeeping, AT q-independence, the known real R-odd defective EP at q ≈ 2.805
/// reproduced), then scans the R-odd residual over the SAME window the R-even census used
/// ((0.2, 3.0) × (−1.5, 1.5), cell 0.05) for complex defective EPs and reads each against the moved
/// rate window of the (4,4) corner (the ComplexQRateWindowTests instrument).
///
/// <para>THE DECISION RULE (pre-committed): an R-odd complex defective EP with Re λ_A &lt; −6 (deep, the
/// moved window silent) would seed a corner-leaking diamond of its own and require a SECOND resultant on
/// the R-odd residual factor; if none exists in the window, the deep-loci set is complete on the R-even
/// side and ONE resultant on F_18 suffices. Run:
/// <c>dotnet test "compute/RCPsiSquared.Diagnostics.Tests" --filter "Category=RODDPROBE"
/// --logger "console;verbosity=detailed"</c>.</para></summary>
public sealed class ROddDeepLociProbeTests
{
    private readonly ITestOutputHelper _out;
    public ROddDeepLociProbeTests(ITestOutputHelper o) => _out = o;

    private const int N = 5;                 // the deciding case: k = N−1 = 4, nBlock = 5
    private const int K = N - 1;

    /// <summary>Symmetric sorted multiset distance (the witness's MatchDist idiom).</summary>
    private static double MatchDist(Complex[] a, Complex[] b)
    {
        if (a.Length != b.Length) return double.NaN;
        var sa = a.OrderBy(z => Math.Round(z.Real, 9)).ThenBy(z => Math.Round(z.Imaginary, 9)).ToArray();
        var sb = b.OrderBy(z => Math.Round(z.Real, 9)).ThenBy(z => Math.Round(z.Imaginary, 9)).ToArray();
        double m = 0;
        for (int i = 0; i < sa.Length; i++) m = Math.Max(m, (sa[i] - sb[i]).Magnitude);
        return m;
    }

    [Fact]
    [Trait("Category", "RODDPROBE")]
    public void N5_ROddSector_SplitBookkeeping_And_AtQIndependence()
    {
        var u = F89PathKSeDeBlock.ROddBasis(N);
        Assert.Equal(50, u.GetLength(0));                                 // nBlock·C(nBlock,2)
        Assert.Equal(24, u.GetLength(1));                                 // (50 − 2 fixed singletons)/2

        var q = new Complex(0.7, 0.3);                                    // generic complex q
        var atOdd = PathKMonodromyScout.AtRootsExactROdd(K, q);
        var resOdd = PathKMonodromyScout.ResidualRootsExactROdd(K, q);
        _out.WriteLine($"R-odd sector at N={N}: dim 24 = AT {atOdd.Length} ⊎ residual {resOdd.Length}");
        Assert.Equal(24, atOdd.Length + resOdd.Length);

        // Full (SE,DE) block spectrum = R-even (the ×2-cleared sym block's λ) ⊎ R-odd (AT ⊎ residual).
        var (a, c) = PathKMonodromyScout.BuildLinear(N);
        var even = PathKMonodromyScout.AllRootsAt(a, c, q);               // 26 R-even eigenvalues
        var fullM = Matrix<Complex>.Build.DenseOfArray(F89PathKSeDeBlock.BuildFullBlock(N, q));
        var full = fullM.Evd().EigenValues.ToArray();                     // all 50
        var union = even.Concat(atOdd).Concat(resOdd).ToArray();
        double split = MatchDist(full, union);
        _out.WriteLine($"full spectrum (50) vs R-even (26) ⊎ R-odd (24): match dist {split:E2}");
        Assert.True(split < 1e-8, $"the R split must reproduce the full spectrum; got {split:E2}");

        // The R-odd AT-locked strands are rate-locked and q-LINEAR: λ(q) = r0 + q·(i·f), r0 ∈ {−2,−6} the
        // pinned rung, f real (read at real q=1). At real q the rate sits on the rung; the frequency scales
        // with q (the R-7 scope note's drift is this linearity, not a bug).
        if (atOdd.Length > 0)
        {
            var atReal = PathKMonodromyScout.AtRootsExactROdd(K, new Complex(1, 0));
            double offRung = atReal.Select(z => Math.Min(Math.Abs(z.Real + 2), Math.Abs(z.Real + 6))).Max();
            Assert.True(offRung < 1e-9, $"at real q the R-odd AT rates must sit on Re ∈ {{−2,−6}}; off by {offRung:E2}");

            var q2 = new Complex(2.1, -0.4);
            var predicted = atReal.Select(z =>
                new Complex(Math.Abs(z.Real + 2) < Math.Abs(z.Real + 6) ? -2 : -6, 0)
                + Complex.ImaginaryOne * z.Imaginary * q2).ToArray();
            double lin = MatchDist(predicted, PathKMonodromyScout.AtRootsExactROdd(K, q2));
            _out.WriteLine($"R-odd AT rate-lock (real q): off-rung {offRung:E2}; q-linearity λ(q)=r0+q·(i·f) at q={q2}: {lin:E2}");
            Assert.True(lin < 1e-8, $"R-odd AT roots must be q-linear off the rungs; deviation {lin:E2}");
        }
    }

    [Fact]
    [Trait("Category", "RODDPROBE")]
    public void N5_ROdd_KnownRealDefectiveEP_Reproduced()
    {
        // The third real N=5 defective EP (q ≈ 2.805, λ ≈ −4.488) is the known R-ODD one
        // (PROOF_CODIM1_BY_ADDITIVITY §6; outside the octic sub-story). The R-odd scan must find it.
        var pts = PathKMonodromyScout.FindDiabolicsExactROdd(
            k: K, reLo: 2.5, reHi: 3.1, imLo: -0.15, imHi: 0.15, cell: 0.02);
        foreach (var p in pts)
            _out.WriteLine($"q*={p.QValue.Real:F6}{(p.QValue.Imaginary >= 0 ? "+" : "-")}{Math.Abs(p.QValue.Imaginary):F6}i  " +
                           $"λ={p.MergeLambda.Real:F4}{(p.MergeLambda.Imaginary >= 0 ? "+" : "-")}{Math.Abs(p.MergeLambda.Imaginary):F4}i  " +
                           $"{(p.IsSemisimple ? "diabolic" : "DEFECTIVE")}  gap={p.Gap:E2}  expo={p.GapScalingExponent:F3}");

        var anchor = pts.FirstOrDefault(p =>
            !p.IsSemisimple && Math.Abs(p.QValue.Imaginary) < 1e-3 && Math.Abs(p.QValue.Real - 2.805) < 0.05);
        Assert.True(anchor is not null,
            "the known real R-odd defective EP at q ≈ 2.805 must be reproduced by the R-odd residual scan");
        Assert.True(Math.Abs(anchor!.MergeLambda.Real + 4.488) < 0.05 && Math.Abs(anchor.MergeLambda.Imaginary) < 1e-3,
            $"its eigenvalue must be λ ≈ −4.488 (real); got {anchor.MergeLambda}");
    }

    [Fact]
    [Trait("Category", "RODDPROBE")]
    public void N5_ROdd_ComplexDefectiveLoci_TheDeepLociProbe()
    {
        const int m = 4;   // corner (p_c+1,p_c+1) = (4,4), as in COMPLEXQWINDOW

        var a = Matrix<Complex>.Build.DenseOfArray(WeightCoherenceBlock.Build(N, m, m, Complex.Zero));
        var c = Matrix<Complex>.Build.DenseOfArray(WeightCoherenceBlock.Build(N, m, m, Complex.One)) - a;
        double kNorm = c.L2Norm();
        _out.WriteLine($"corner ({m},{m}) N={N}: dim={a.RowCount}, ‖K‖₂={kNorm:F4}\n");

        // The SAME window and cell the R-even census used (ComplexQRateWindowTests / FindDiabolicsExact).
        var pts = PathKMonodromyScout.FindDiabolicsExactROdd(
            k: K, reLo: 0.2, reHi: 3.0, imLo: -1.5, imHi: 1.5, cell: 0.05);
        var defective = pts.Where(p => !p.IsSemisimple && Math.Abs(p.QValue.Imaginary) > 1e-3)
                           .OrderBy(p => Math.Abs(p.QValue.Imaginary)).ToList();
        _out.WriteLine($"FindDiabolicsExactROdd(k={K}): {pts.Count} coalescences, {defective.Count} complex defective EPs");
        foreach (var p in pts.OrderBy(x => x.QValue.Real).ThenBy(x => x.QValue.Imaginary))   // the full R-odd census, for the record
            _out.WriteLine($"  census: q*={p.QValue.Real,7:F4}{(p.QValue.Imaginary >= 0 ? "+" : "-")}{Math.Abs(p.QValue.Imaginary):F4}i  " +
                           $"λ={p.MergeLambda.Real,8:F4}{(p.MergeLambda.Imaginary >= 0 ? "+" : "-")}{Math.Abs(p.MergeLambda.Imaginary):F4}i  " +
                           $"{(p.IsSemisimple ? "diabolic " : "DEFECTIVE")}  gap={p.Gap:E1}  expo={p.GapScalingExponent:F3}  loopId={p.LoopIsIdentity}");
        _out.WriteLine("");

        // The R-parity refinement: W = Σ c_l†ρc_l is site-reflection-symmetric, so it COMMUTES with R and an
        // R-odd (1,2) defective chain arrives in the corner still R-odd. Its kernel death therefore needs μ
        // absent only from the R-ODD SECTOR of the corner, whose moved window is its own (smaller) Bendixson
        // box. Build the corner's reflection permutation ((a,b) ↦ (rev a, rev b)) and its R-odd columns.
        var permC = CornerReflectionPermutation(N, m);
        var oddCols = new List<(int T, int T2)>();
        for (int t = 0; t < permC.Length; t++)
            if (permC[t] > t) oddCols.Add((t, permC[t]));
        double inv2 = 1.0 / Math.Sqrt(2.0);
        var uOdd = Matrix<Complex>.Build.Dense(permC.Length, oddCols.Count, (r, s) =>
            r == oddCols[s].T ? new Complex(inv2, 0) : r == oddCols[s].T2 ? new Complex(-inv2, 0) : Complex.Zero);
        _out.WriteLine($"corner R split: dim {permC.Length} = R-even {permC.Length - oddCols.Count} ⊎ R-odd {oddCols.Count}\n");

        int deep = 0, movedClosed = 0, parityClosed = 0, open = 0;
        foreach (var d in defective)
        {
            Complex qs = d.QValue, lamA = d.MergeLambda;
            Complex mu = -lamA - 2.0 * N;                                // the holomorphic fold value

            var corner = Matrix<Complex>.Build.DenseOfArray(WeightCoherenceBlock.Build(N, m, m, qs));
            double commuteR = (Permute(corner, permC) - corner).L2Norm();  // R commutes with the corner block
            Assert.True(commuteR < 1e-12, $"R must commute with the corner block; ‖RMR−M‖ = {commuteR:E2}");
            double cornerExcludesMu = corner.Evd().EigenValues.Select(z => (z - mu).Magnitude).Min();

            var h = (corner + corner.ConjugateTranspose()) / 2.0;        // Bendixson Hermitian part at THIS q
            double hBottom = h.Evd().EigenValues.Select(z => z.Real).Min();
            double widen = Math.Abs(qs.Imaginary) * kNorm;
            bool movedExcluded = mu.Real < hBottom - 1e-9;

            var cornerOdd = uOdd.ConjugateTranspose() * corner * uOdd;   // the R-ODD corner sector
            var hOdd = (cornerOdd + cornerOdd.ConjugateTranspose()) / 2.0;
            double hOddBottom = hOdd.Evd().EigenValues.Select(z => z.Real).Min();
            double oddExcludesMu = cornerOdd.Evd().EigenValues.Select(z => (z - mu).Magnitude).Min();
            bool parityExcluded = !movedExcluded && mu.Real < hOddBottom - 1e-9;

            bool isDeep = lamA.Real < -6.0;
            if (isDeep) deep++;
            else if (movedExcluded) movedClosed++;
            else if (parityExcluded) parityClosed++;
            else open++;

            _out.WriteLine(
                $"q*={qs.Real,7:F4}{(qs.Imaginary >= 0 ? "+" : "-")}{Math.Abs(qs.Imaginary):F4}i  " +
                $"λ_A={lamA.Real,8:F4}{(lamA.Imaginary >= 0 ? "+" : "-")}{Math.Abs(lamA.Imaginary):F4}i  Reμ={mu.Real,8:F4}");
            _out.WriteLine(
                $"      movedWin bottom={hBottom,7:F3}  widen={widen,5:F3}(margin {lamA.Real + 6.0,6:F3})  " +
                $"corner min|spec−μ|={cornerExcludesMu:E2}");
            _out.WriteLine(
                $"      R-ODD corner sector: movedWin bottom={hOddBottom,7:F3}  min|spec−μ|={oddExcludesMu:E2}  →  " +
                $"{(isDeep ? "DEEP (Re λ_A < −6)"
                   : movedExcluded ? "MOVED-WINDOW CLOSED (full corner)"
                   : parityExcluded ? "R-PARITY MOVED-WINDOW CLOSED (R-odd sector)"
                   : "OPEN (window silent both ways; absence census-only)")}");
        }

        _out.WriteLine("");
        _out.WriteLine(
            $"VERDICT N={N}: {defective.Count} R-odd complex defective loci; {deep} DEEP (Re λ_A < −6); " +
            $"{movedClosed} closed by the full-corner moved window; {parityClosed} closed by the R-parity-refined " +
            $"(R-odd corner sector) moved window; {open} OPEN.");
        _out.WriteLine(deep == 0 && open == 0
            ? "⟹ NO R-odd locus reaches the resultant: the deep-loci set is complete on the R-even side alone, " +
              "ONE resultant on F_18 suffices."
            : deep == 0
              ? $"⟹ no R-odd DEEP loci, but {open} loci escape both windows: their corner absence is census-only " +
                "and needs either a sharpened window argument or the R-odd residual factor in the resultant."
              : $"⟹ {deep} R-odd deep loci seed corner-leaking diamonds; a SECOND resultant on the R-odd residual " +
                "factor is needed.");

        // Structural sanity, pre-committed: the scan must at least rediscover the known real R-odd EP
        // (q ≈ 2.805 lies inside the window), so an empty result means the instrument failed, not "no loci".
        Assert.Contains(pts, p => !p.IsSemisimple && Math.Abs(p.QValue.Imaginary) < 1e-3
                                  && Math.Abs(p.QValue.Real - 2.805) < 0.06);
    }

    /// <summary>The site-reflection permutation on the (w,w) weight-coherence basis: |a⟩⟨b| ↦ |rev a⟩⟨rev b|
    /// (bit reversal over n sites), in <see cref="WeightCoherenceBlock.Build"/>'s ordering (kets outer, bras
    /// inner, both ascending-mask). An involution commuting with the block (the chain's bond set is
    /// reversal-symmetric and n_diff is reversal-invariant).</summary>
    private static int[] CornerReflectionPermutation(int n, int w)
    {
        var configs = WeightCoherenceBlock.Configs(n, w);
        var index = new Dictionary<int, int>();
        for (int i = 0; i < configs.Count; i++) index[configs[i]] = i;
        int Rev(int mask)
        {
            int r = 0;
            for (int b = 0; b < n; b++) if (((mask >> b) & 1) != 0) r |= 1 << (n - 1 - b);
            return r;
        }
        int nc = configs.Count;
        var perm = new int[nc * nc];
        for (int ki = 0; ki < nc; ki++)
            for (int bi = 0; bi < nc; bi++)
                perm[ki * nc + bi] = index[Rev(configs[ki])] * nc + index[Rev(configs[bi])];
        return perm;
    }

    /// <summary>R·M·R for a basis permutation R (perm[t] = the image index; an involution).</summary>
    private static Matrix<Complex> Permute(Matrix<Complex> mIn, int[] perm)
    {
        int d = mIn.RowCount;
        var outM = Matrix<Complex>.Build.Dense(d, d);
        for (int r = 0; r < d; r++)
            for (int s = 0; s < d; s++)
                outM[perm[r], perm[s]] = mIn[r, s];
        return outM;
    }
}
