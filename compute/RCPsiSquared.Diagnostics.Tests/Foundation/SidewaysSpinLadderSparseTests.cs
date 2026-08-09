using System;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;
using System.Numerics;
using RCPsiSquared.Diagnostics.Foundation;
using Xunit;
using Xunit.Abstractions;

namespace RCPsiSquared.Diagnostics.Tests.Foundation;

/// <summary>Gates for <see cref="SidewaysSpinLadderSparse"/>, in the order the trust is built:
/// (1) the sparse builders are entry-identical to the dense ones at N=5 (exact, integer test vectors, so
/// the comparison is a case-1 exact zero); (2) the whole sparse walk reproduces the dense walk's measured
/// N=7 numbers (CONTROL exactly 0.0, CG norms 2, √6, √6, 2); (3) only then is N=9 read, the falsifiable
/// next case of <see cref="RCPsiSquared.Core.Symmetry.SidewaysSpinLadderClaim"/>: ℓ = 3, predicted S⁺
/// transport norms √6, √10, √12, √12, √10, √6 on both chains, seven interior sectors each, with the
/// middle blocks (10584², 15876²) reached by inverse iteration instead of dense zgeev.</summary>
public class SidewaysSpinLadderSparseTests
{
    private static readonly CultureInfo Inv = CultureInfo.InvariantCulture;
    private readonly ITestOutputHelper _out;
    public SidewaysSpinLadderSparseTests(ITestOutputHelper o) => _out = o;

    private static string F(double v) => v.ToString("E1", Inv);

    [Fact(DisplayName = "sparse builders ≡ dense builders at N=5 (CSR entries and ladder applies, exact)")]
    [Trait("Category", "SIDEWAYS")]
    public void SparseBuilders_MatchDense_N5()
    {
        const int n = 5;
        var (qStar, _) = SidewaysSpinLadderSparse.Seed(n);
        double j = 2.0 * qStar;
        foreach (var (p, q) in new[] { (1, 3), (2, 2), (2, 4), (3, 3), (3, 1) })
        {
            // CSR vs dense block: every dense entry appears in the CSR row (and vice versa by count).
            var dense = SidewaysSpinLadderChain.BuildBlock(n, p, q, qStar).ToArray();
            var csr = SidewaysSpinLadderSparse.BuildCsr(n, p, q, j);
            int dim = dense.GetLength(0);
            Assert.Equal(dim, csr.Dim);
            int denseNnz = 0;
            for (int r = 0; r < dim; r++)
                for (int c = 0; c < dim; c++)
                    if (dense[r, c] != Complex.Zero) denseNnz++;
            int csrNnzNonZero = 0;
            for (int r = 0; r < dim; r++)
                for (int k = csr.RowPtr[r]; k < csr.RowPtr[r + 1]; k++)
                {
                    Assert.True(dense[r, csr.ColIdx[k]] == csr.Values[k],
                        $"({p},{q}) CSR entry [{r},{csr.ColIdx[k]}] = {csr.Values[k]} != dense {dense[r, csr.ColIdx[k]]}");
                    if (csr.Values[k] != Complex.Zero) csrNnzNonZero++;
                }
            Assert.Equal(denseNnz, csrNnzNonZero);

            // Ladder applies vs dense ladder matrices, on an integer test vector (entries in {−2..2}),
            // so both routes sum the same integers and the comparison is exactly 0.0.
            var rng = new Random(97 + p * 10 + q);
            var x = new Complex[dim];
            for (int i = 0; i < dim; i++) x[i] = new Complex(rng.Next(-2, 3), rng.Next(-2, 3));
            foreach (bool phi in new[] { true, false })
            {
                int tq = q + (phi ? 1 : -1);
                if (tq < 0 || tq > n || p + 1 > n) continue;
                var mDense = (phi ? SpectatorIntertwiner.BuildW(n, p, q)
                                  : SpectatorIntertwiner.BuildSPlus(n, p, q)).ToArray();
                int tdim = mDense.GetLength(0);
                var yDense = new Complex[tdim];
                for (int r = 0; r < tdim; r++)
                {
                    Complex acc = Complex.Zero;
                    for (int c = 0; c < dim; c++) acc += mDense[r, c] * x[c];
                    yDense[r] = acc;
                }
                var ySparse = new Complex[tdim];
                SidewaysSpinLadderSparse.ApplyLadder(n, p, q, phi, x, ySparse);
                for (int r = 0; r < tdim; r++)
                    Assert.True(yDense[r] == ySparse[r],
                        $"({p},{q}) {(phi ? "Φ" : "S⁺")} apply row {r}: dense {yDense[r]} != sparse {ySparse[r]}");
            }
        }
    }

    [Fact(DisplayName = "inverse iteration finds the fold eigenpair on the N=5 (1,3) block (λ and ‖S⁺v‖ = √2)")]
    [Trait("Category", "SIDEWAYS")]
    public void InverseIteration_FindsFoldPair_N5()
    {
        const int n = 5;
        var (qStar, lamA) = SidewaysSpinLadderSparse.Seed(n);
        double fold = -lamA - 2.0 * n;
        var block = SidewaysSpinLadderSparse.BuildCsr(n, 1, 3, 2.0 * qStar);
        var (lambda, v, eigRes, outer) = SidewaysSpinLadderSparse.EigNearShift(block, fold, seed: 1,
            log: s => _out.WriteLine("  " + s));
        _out.WriteLine($"λ = {lambda.Real.ToString("F6", Inv)} vs fold {fold.ToString("F6", Inv)}; " +
                       $"eig residual {F(eigRes)} in {outer} outer iterations");
        Assert.True((lambda - fold).Magnitude < 1e-2, $"inverse iteration missed the fold: λ = {lambda}");
        // The floor is the LU solve accuracy plus the pair split; the norm check below is the criterion
        // that matters, the eigen residual only has to be deep enough that the out-of-pair-plane error
        // cannot move the norm at 1e-6.
        Assert.True(eigRes < 1e-9, $"eigen residual {eigRes:E1} did not converge");
        var img = new Complex[SpectatorIntertwiner.BuildSPlus(n, 1, 3).RowCount];
        SidewaysSpinLadderSparse.ApplyLadder(n, 1, 3, phi: false, v, img);
        double nrm = Math.Sqrt(img.Sum(z => z.Real * z.Real + z.Imaginary * z.Imaginary));
        _out.WriteLine($"‖S⁺v‖ = {nrm.ToString("0.000000", Inv)} (CG √2 = {Math.Sqrt(2.0).ToString("0.000000", Inv)})");
        Assert.True(Math.Abs(nrm - Math.Sqrt(2.0)) < 1e-6, $"‖S⁺v‖ = {nrm} != √2");
    }

    [Fact(DisplayName = "the N=9 seed is the recorded census entry (RealDefectiveSeeds, first R-parity +1 at N=9)")]
    [Trait("Category", "SIDEWAYS")]
    public void SeedN9_IsTheRecordedCensusEntry()
    {
        var (qStar, lamA) = SidewaysSpinLadderSparse.Seed(9);
        Assert.Contains(RCPsiSquared.Core.F89PathK.RealDefectiveSeeds.All,
            s => s.N == 9 && s.QStar == qStar && s.LambdaA == lamA && s.RParity == +1);
        Assert.Equal(qStar,
            RCPsiSquared.Core.F89PathK.RealDefectiveSeeds.ForN(9).First(s => s.RParity == +1).QStar);
    }

    // denseCutoff 500 forces the (2,4)/(3,3)-class blocks (735, 1225) through the LU shift-invert path,
    // so the LU eigen path is validated where the dense walk's answer (2, √6, √6, 2) is already gated,
    // and the walk is DIFFED rung by rung against a dense SidewaysSpinLadderChain.Run(7) in the same test.
    [Fact(DisplayName = "the sparse walk reproduces the dense walk at N=7 (rung-by-rung diff, LU path forced)")]
    [Trait("Category", "SLOW_SIDEWAYS")]
    public void SparseWalk_MatchesDense_N7() => RunAndGate(7, denseCutoff: 500, diffAgainstDense: true);

    // The falsifiable next case of SidewaysSpinLadderClaim: ℓ = 3, predicted norms
    // √6, √10, √12, √12, √10, √6, seven interior sectors per chain, 28 = 4·9−8 with the band's half.
    // Seed: the first R-parity +1 census entry at N=9 (RealDefectiveSeeds), q* = 0.591760, λ_A = −5.1206.
    // First run 2026-08-09: CONFIRMED on both chains to six decimals.
    [Fact(DisplayName = "N=9, ℓ=3: the falsifiable CG prediction √6,√10,√12,√12,√10,√6 (LU shift-invert past the dense wall)")]
    [Trait("Category", "SLOW_SIDEWAYS9")]
    public void Gate_N9() => RunAndGate(9);

    private void RunAndGate(int n, int denseCutoff = SidewaysSpinLadderSparse.DenseCutoff,
        bool diffAgainstDense = false)
    {
        var run = SidewaysSpinLadderSparse.Run(n, s => _out.WriteLine("  " + s), denseCutoff);
        _out.WriteLine($"N = {n}, q* = {run.QStar.ToString(Inv)}, λ_A = {run.LambdaA.ToString(Inv)}, " +
                       $"fold = {run.Fold.ToString("0.####", Inv)}");
        _out.WriteLine($"CONTROL worst residual: {F(run.WorstControlResidual)}");
        Assert.True(run.WorstControlResidual == 0.0,
            $"CONTROL: sparse intertwining residual not exactly 0.0: {run.WorstControlResidual:E3}");

        foreach (var chain in run.Chains)
        {
            var norms = chain.Rungs.Select(r => r.Norm).ToList();
            _out.WriteLine($"chain p+q̃ = {chain.Total}: norms [{string.Join(", ", norms.Select(x => x.ToString("0.000000", Inv)))}]"
                           + $" vs CG [{string.Join(", ", run.CgPredicted.Select(x => x.ToString("0.000000", Inv)))}]; "
                           + $"terminal {F(chain.Terminal)}, ratio {F(chain.TerminalRatio)}, "
                           + $"survivors {chain.Survivors}/{chain.TerminalSourceDim} (target dim {chain.TerminalTargetDim})");

            // Every walked interior source must actually carry the fold (the census-verified SHAPE,
            // re-checked here per source). The bound is 1e-2 and its error model is measured, not
            // wished: the fold PAIR members sit ≤ ~1.5e-3 from the recorded 4-decimal fold (N=9 members
            // −12.878060 / −12.880829, split 2.77e-3; the recorded λ_A carries a ~1e-3 member-noise
            // floor per the RealDefectiveSeeds docstring), while the nearest UNRELATED eigenvalue
            // measured is 4.6e-2 away (the 3024-dim N=9 block; 1.7e-1 at 324; ~2.9e-1 at N≤7). The
            // margin shrinks with block dimension (~dim^−0.4), extrapolating to ≈2e-2 at the 15876
            // block: still above the gate, recorded here as the thinnest measured margin.
            // run.FoldTolerance is informational (circular as a gate, see the Run comment).
            foreach (var rung in chain.Rungs)
                Assert.True((rung.Lambda - new Complex(run.Fold, 0)).Magnitude < 1e-2,
                    $"({rung.P},{rung.Q}) λ = {rung.Lambda} misses the fold {run.Fold} beyond 1e-2");

            // The CG prediction, no free parameter (the N=9 falsification when n = 9).
            Assert.True(norms.Count == run.CgPredicted.Count
                        && norms.Zip(run.CgPredicted).All(t => Math.Abs(t.First - t.Second) < 1e-6),
                $"chain p+q̃={chain.Total}: norms are NOT the CG coefficients for ℓ = {(n - 3) / 2.0}: "
                + $"[{string.Join(", ", norms)}] vs [{string.Join(", ", run.CgPredicted)}]");

            // Image eigenvector residual: one law per eigen path, each measured (case 2). Dense zgeev
            // rungs: residual/(eps·√targetDim) stays O(10..100), the dense gate's law. LU shift-invert
            // rungs: the floor is the inverse-iteration pair-mixing decay, not eps·√dim — measured
            // 4.3e-13..1.2e-11 at N=7 (forced) and 4.3e-12..1.4e-10 at N=9; gate < 1e-8, two decades
            // above the worst measurement. (The first N=9 run showed 1.2e-6 here and dressed it as a
            // "√eps near-defective floor"; an independent review reproduction showed that number was an
            // unconverged in-plane MIX of the pair kept by a buggy pair-split acceptance guard — the
            // members reach 8.6e-14 — so that two-regime story is deliberately NOT kept.)
            foreach (var rung in chain.Rungs)
            {
                double eps = Math.Pow(2, -52);
                double scaled = rung.EigResidual / (eps * Math.Sqrt(rung.TargetDim));
                if (rung.ViaInverseIteration)
                    Assert.True(rung.EigResidual < 1e-8,
                        $"({rung.P},{rung.Q}) image residual {rung.EigResidual:E1} above the LU "
                        + "inverse-iteration floor law (measured 4e-13..1.4e-10 at N=7,9; gate 1e-8)");
                else
                    Assert.True(scaled < 100.0,
                        $"({rung.P},{rung.Q}) image residual {rung.EigResidual:E1} = {scaled:0.0}·eps·√dim "
                        + "(law: O(10), same as the dense gate)");
            }

            Assert.True(chain.Survivors >= chain.TerminalTargetDim,
                $"chain p+q̃={chain.Total}: survivors {chain.Survivors} < target dim {chain.TerminalTargetDim} (negative control failed)");
            Assert.True(chain.TerminalRatio < 1e-12,
                $"chain p+q̃={chain.Total}: terminal ratio {chain.TerminalRatio:E1} not at the floor (highest weight S⁺|ℓ,ℓ⟩ = 0)");
        }

        if (diffAgainstDense)
        {
            // The actual dense diff (the display name's promise): the dense walk at the same N, rung by
            // rung against the sparse walk's norms; both eigensolver-limited, so 1e-8 is generous.
            var dense = SidewaysSpinLadderChain.Run(n);
            Assert.Equal(dense.Chains.Count, run.Chains.Count);
            foreach (var (dc, sc) in dense.Chains.Zip(run.Chains))
            {
                Assert.Equal(dc.Total, sc.Total);
                var dn = dc.Rungs.Select(r => r.Norm).ToList();
                var sn = sc.Rungs.Select(r => r.Norm).ToList();
                _out.WriteLine($"dense-diff p+q̃={dc.Total}: max |Δnorm| = "
                               + $"{(dn.Count == sn.Count ? dn.Zip(sn).Max(t => Math.Abs(t.First - t.Second)) : double.NaN):E1}");
                Assert.True(dn.Count == sn.Count && dn.Zip(sn).All(t => Math.Abs(t.First - t.Second) < 1e-8),
                    $"p+q̃={dc.Total}: sparse norms [{string.Join(", ", sn)}] deviate from dense [{string.Join(", ", dn)}]");
            }
        }
    }
}
