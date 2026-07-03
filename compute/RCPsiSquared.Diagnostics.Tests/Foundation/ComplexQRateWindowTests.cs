using System;
using System.Linq;
using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.F89PathK;
using RCPsiSquared.Diagnostics.Foundation;
using Xunit;
using Xunit.Abstractions;

namespace RCPsiSquared.Diagnostics.Tests.Foundation;

/// <summary>The COMPLEX-q rate window (the moved Bendixson), extending §6 of PROOF_CODIM1_BY_ADDITIVITY to
/// complex loci (remainder R1, the arc's ONE open item). The rate-window lemma excludes the diamond's fold
/// value from the corner block (p_c+1,p_c+1) at REAL loci; the proof stops there ("Bendixson needs real q").
///
/// <para>The extension: Bendixson holds at ANY q, Re λ ∈ [λ_min(H), λ_max(H)] with H = (L+L†)/2. For the
/// corner pencil L(q) = A + qC (A = −2·diag(n_diff) real, C anti-Hermitian) the Hermitian part is
/// H(q) = A + Im(q)·K, K = iC Hermitian, so the rate window [−2(N−3), 0] does not vanish at complex q, it
/// MOVES (each edge by ≤ |Im q|·‖K‖₂, Weyl). The value the core (p_c,p_c) carries and W maps to the corner is
/// the HOLOMORPHIC fold μ = −λ_A − 2N (the §7 cross-fold ket-flip (p_c,p_c−1)→(p_c,p_c) at q̄ composed with
/// the F1 form λ_A(q̄) = conj(λ_A(q)); NOT the antiunitary −conj(λ_A)−2N, which coincides only at real q).
/// So kernel death survives at complex q wherever
/// <code>Re λ_A(q*) + 6  &gt;  |Im q*|·‖K‖₂    (N-uniform, 2N − 2(N−3) = 6),</code>
/// equivalently Re μ &lt; λ_min(H(q*)). This closes the complex defective loci where λ_A stays inside the (1,2)
/// block's real-q Bendixson floor (Re λ_A &gt; −6); the DEEP loci where λ_A dives below −6 at complex q push μ
/// up into the corner window and are left to the holomorphic resultant Res_λ(F_18, F_corner(−λ−2N)).</para>
///
/// <para>Asserts (N=5): the holomorphic-fold READING (core carries μ, numerically to locus precision; the
/// exact identity is now DERIVED, a full-spectrum §7 corollary, gate HolomorphicFoldIdentityTests), the absence (corner excludes μ,
/// numerical census), and that the moved window closes ≥1 complex locus. Run:
/// <c>dotnet test "compute/RCPsiSquared.Diagnostics.Tests" --filter "Category=COMPLEXQWINDOW"
/// --logger "console;verbosity=detailed"</c>.</para></summary>
public sealed class ComplexQRateWindowTests
{
    private readonly ITestOutputHelper _out;
    public ComplexQRateWindowTests(ITestOutputHelper o) => _out = o;

    [Fact]
    [Trait("Category", "COMPLEXQWINDOW")]
    public void N5_ComplexDefectiveLoci_MovedRateWindow()
    {
        const int N = 5, pc = 3, m = 4;   // core (p_c,p_c) = (3,3); corner (p_c+1,p_c+1) = (4,4)

        // The corner block's linear parts L(q) = A + qC: A = L(0) (dephasing diagonal), C = L(1) − L(0).
        var a = Matrix<Complex>.Build.DenseOfArray(WeightCoherenceBlock.Build(N, m, m, Complex.Zero));
        var c = Matrix<Complex>.Build.DenseOfArray(WeightCoherenceBlock.Build(N, m, m, Complex.One)) - a;
        double kNorm = c.L2Norm();                                   // ‖K‖₂ = ‖iC‖₂ = ‖C‖₂
        double winBottom = a.Diagonal().Select(z => z.Real).Min();   // −2(N−3) = −4
        double winTop    = a.Diagonal().Select(z => z.Real).Max();   // 0
        _out.WriteLine($"corner ({m},{m}) N={N}: dim={a.RowCount}, real rate window [{winBottom:F3},{winTop:F3}], ‖K‖₂={kNorm:F4}\n");

        // The (1,2) residual's DEFECTIVE EPs (P₁₀ zeros), right half-plane, complex only (|Im q| > 1e-3).
        var pts = PathKMonodromyScout.FindDiabolicsExact(
            k: N - 1, reLo: 0.2, reHi: 3.0, imLo: -1.5, imHi: 1.5, cell: 0.05);
        var defective = pts.Where(p => !p.IsSemisimple && Math.Abs(p.QValue.Imaginary) > 1e-3)
                           .OrderBy(p => Math.Abs(p.QValue.Imaginary)).ToList();
        _out.WriteLine($"FindDiabolicsExact(k={N - 1}): {pts.Count} coalescences, {defective.Count} complex defective EPs\n");

        int movedClosed = 0;
        foreach (var d in defective)
        {
            Complex qs = d.QValue, lamA = d.MergeLambda;
            Complex mu     = -lamA - 2.0 * N;                        // the HOLOMORPHIC fold (what the core carries)
            Complex muConj = -Complex.Conjugate(lamA) - 2.0 * N;     // the antiunitary value (coincides only at real q)

            var core   = Matrix<Complex>.Build.DenseOfArray(WeightCoherenceBlock.Build(N, pc, pc, qs));
            var corner = Matrix<Complex>.Build.DenseOfArray(WeightCoherenceBlock.Build(N, m, m, qs));
            var coreEv   = core.Evd().EigenValues;
            var cornerEv = corner.Evd().EigenValues;
            double coreCarriesMu   = coreEv.Select(z => (z - mu).Magnitude).Min();       // holomorphic fold identity
            double coreCarriesConj = coreEv.Select(z => (z - muConj).Magnitude).Min();   // the conj: NOT carried
            double cornerExcludesMu = cornerEv.Select(z => (z - mu).Magnitude).Min();     // absence (census)

            var h = (corner + corner.ConjugateTranspose()) / 2.0;   // Bendixson Hermitian part at THIS complex q
            double hBottom = h.Evd().EigenValues.Select(z => z.Real).Min();
            double widen = Math.Abs(qs.Imaginary) * kNorm;
            bool movedExcluded = mu.Real < hBottom - 1e-9;          // sharp moved-window exclusion
            bool uniform       = (lamA.Real + 6.0) > widen;          // N-uniform ‖K‖ sufficient condition
            bool insideBlockFloor = lamA.Real > -6.0;                // λ_A inside the (1,2) block's real-q window
            if (movedExcluded) movedClosed++;

            // Structural facts, asserted:
            Assert.True(coreCarriesMu < 1e-2,
                $"core (3,3) must carry the holomorphic fold μ=−λ_A−2N; got {coreCarriesMu:E2} at q={qs}");
            Assert.True(cornerExcludesMu > 1e-2,
                $"corner (4,4) must exclude μ (the absence, census); got {cornerExcludesMu:E2} at q={qs}");

            _out.WriteLine(
                $"q*={qs.Real,7:F4}{(qs.Imaginary >= 0 ? "+" : "-")}{Math.Abs(qs.Imaginary):F4}i  " +
                $"λ_A={lamA.Real,8:F4}{(lamA.Imaginary >= 0 ? "+" : "-")}{Math.Abs(lamA.Imaginary):F4}i  Reμ={mu.Real,8:F4}");
            _out.WriteLine(
                $"      core carries μ=−λ_A−2N: {coreCarriesMu:E2}  (conj −conj(λ_A)−2N: {coreCarriesConj:E2})  |  " +
                $"corner excludes μ: {cornerExcludesMu:E2}");
            _out.WriteLine(
                $"      movedWin bottom={hBottom,7:F3}  widen={widen,5:F3}(margin {lamA.Real + 6.0,6:F3})  " +
                $"{(movedExcluded ? "MOVED-WINDOW CLOSED" : "needs resultant")}" +
                $"  [λ_A {(insideBlockFloor ? ">" : "<")} −6{(uniform ? ", uniform✓" : "")}]");
        }

        _out.WriteLine("");
        _out.WriteLine($"SUMMARY N={N}: {movedClosed}/{defective.Count} complex defective loci CLOSED by the moved rate window " +
                       $"(Re λ_A > −6); {defective.Count - movedClosed} deep loci (Re λ_A < −6) reduced to the holomorphic " +
                       $"resultant Res_λ(F_18, F_corner(−λ−2N)). Absence (core carries μ, corner excludes it) holds at all {defective.Count}.");

        Assert.NotEmpty(defective);
        Assert.True(movedClosed >= 1,
            $"the moved rate window must structurally close at least one complex defective locus; got {movedClosed}");
    }
}
