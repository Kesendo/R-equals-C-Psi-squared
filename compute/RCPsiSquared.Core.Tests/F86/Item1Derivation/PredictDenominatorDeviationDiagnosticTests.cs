using System.Linq;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.CoherenceBlocks;
using RCPsiSquared.Core.F86.Item1Derivation;
using RCPsiSquared.Core.Symmetry;
using Xunit;
using Xunit.Abstractions;

namespace RCPsiSquared.Core.Tests.F86.Item1Derivation;

/// <summary>Diagnostic for the k=31,32 failure of
/// <c>C2FullBlockSigmaAnatomyTests.PredictDenominator_AtKHigherStretch_MatchesExtractedFromAnatomy</c>.
///
/// <para>Finding: the small deviation (~2-4e-4) at k=31,32 is NOT a bug, it is Vandermonde
/// extraction conditioning. The verification test extracts polynomial coefficients via a
/// double-precision Vandermonde solve, then checks <c>|coef·D − integer| &lt; 1e-4</c> with a
/// FIXED ABSOLUTE tolerance. The integrality residual scales as <c>cond(V)·ε_machine</c>: as
/// the Vandermonde size n = ⌊(k+1)/2⌋ grows (13 → 15 → 16), cond(V) grows (~7e7 → ~2e9 → ~8e9)
/// and the unscaled 1e-4 tolerance is eventually exceeded — at n=16, i.e. k=31.</para>
///
/// <para>Model: <c>observed ≈ C · cond(V) · ε_machine</c> with C ~ O(10²) — empirically
/// <c>observed/(cond(V)·ε)</c> stays in a bounded band across all k. The residual does NOT
/// scale with D (k=31 has D=7.87e6 &gt; k=32's D=2.10e6, yet it tracks cond(V), not D); the
/// worst-case <c>D·‖coef‖·cond·ε</c> bound over-shoots the actual residual by ~5 orders.</para>
///
/// <para>Structural context: <c>docs/proofs/PROOF_F86B_OBSTRUCTION.md</c> — the F86b
/// obstruction proof (g_eff / D_k admit no closed form by the six explored routes).
/// The k=31,32 verification rows are kept as a deliberate red signal; this diagnostic
/// characterises that deviation, it does not suppress it.</para>
/// </summary>
public class PredictDenominatorDeviationDiagnosticTests
{
    private readonly ITestOutputHelper _out;
    public PredictDenominatorDeviationDiagnosticTests(ITestOutputHelper @out) => _out = @out;

    private static CoherenceBlock C2Block(int N) =>
        new CoherenceBlock(N: N, n: 1, gammaZero: 0.05);

    [Fact]
    public void Deviation_IsExplainedByVandermondeConditioning()
    {
        const double epsMachine = 2.220446049250313e-16;
        int[] ks = { 25, 30, 31, 32 };
        var rows = new List<(int k, double observed, double condScale)>();

        foreach (int k in ks)
        {
            int N = k + 1;
            var block = C2Block(N);
            var anatomy = C2FullBlockSigmaAnatomy.BuildFaOnly(block);

            // Replicate ExtractRawPolynomialCoefficients EXACTLY:
            // F_a witnesses sorted by Bloch index n.
            var fa = anatomy.SigmaSpectrum
                .Where(w => w.BlochIndexN.HasValue)
                .OrderBy(w => w.BlochIndexN!.Value)
                .ToList();

            int n = fa.Count;                 // = floor(N/2) = floor((k+1)/2)
            int chainN = N;
            double[] yVals = new double[n];
            double[] rhs = new double[n];
            for (int i = 0; i < n; i++)
            {
                int blochN = fa[i].BlochIndexN!.Value;
                yVals[i] = F89PathKAtLockMechanismClaim.BlochEigenvalueY(chainN, blochN);
                rhs[i] = fa[i].Sigma * chainN * chainN * (chainN - 1);
            }

            // Vandermonde V[i,j] = y_i^j (low-to-high powers), same as the source.
            var V = Matrix<double>.Build.Dense(n, n, (i, j) => System.Math.Pow(yVals[i], j));
            var rhsVec = Vector<double>.Build.DenseOfArray(rhs);

            double cond = V.ConditionNumber();
            var coef = V.Solve(rhsVec);
            int D = F89UnifiedFaClosedFormClaim.PredictDenominator(k);

            double coefL2 = System.Math.Sqrt(coef.Sum(c => c * c));
            double observed = Enumerable.Range(0, n)
                .Max(i => System.Math.Abs(coef[i] * D - System.Math.Round(coef[i] * D)));

            // The honest conditioning scale: a degree-(n-1) Vandermonde solve amplifies
            // machine epsilon by cond(V). The integrality residual empirically scales as
            // cond(V)·ε — not with D, not with ‖coef‖.
            double condScale = cond * epsMachine;

            _out.WriteLine(
                $"k={k} n={n} cond(V)={cond:E3} D={D} ‖coef‖={coefL2:E3} " +
                $"observed={observed:E3} cond·ε={condScale:E3} " +
                $"observed/(cond·ε)={observed / condScale:F1} " +
                $"observed/(D·ε)={observed / (D * epsMachine):E2}");

            rows.Add((k, observed, condScale));
        }

        // --- The proof: the deviation is Vandermonde conditioning, not a bug ---

        // (1) observed = Θ(cond(V)·ε): bounded on BOTH sides. The upper bound (< 1000×)
        //     rules out a pathological blow-up — the residual never exceeds the conditioning
        //     scale by more than an O(10²-10³) constant. The lower bound (> 10×) shows the
        //     residual really IS at the conditioning scale, not accidentally orders below.
        //     Together: the deviation is exactly the Vandermonde solve amplifying ε_machine.
        foreach (var (k, observed, condScale) in rows)
        {
            Assert.InRange(observed / condScale, 10.0, 1000.0);
        }

        // (2) cond(V)·ε is what crosses the fixed 1e-4 tolerance. cond(V) grows with the
        //     Vandermonde size n = ⌊(k+1)/2⌋; the original test's 1e-4 tolerance does not
        //     scale with it. observed crosses 1e-4 between k=30 (n=15) and k=31 (n=16) —
        //     exactly where the Vandermonde grows — and does so while D at k=31 (7.87e6)
        //     EXCEEDS D at k=32 (2.10e6), proving D is not the driver.
        double CondScale(int k) => rows.First(r => r.k == k).condScale;
        double Obs(int k) => rows.First(r => r.k == k).observed;

        Assert.True(CondScale(25) < CondScale(30) && CondScale(30) < CondScale(31),
            "cond(V)·ε must grow with the Vandermonde size from k=25 through k=31.");
        Assert.True(Obs(25) < 1e-4 && Obs(30) < 1e-4,
            "k=25,30 (n≤15): observed integrality residual below the fixed 1e-4 tolerance — " +
            "the verification test passes here.");
        Assert.True(Obs(31) > 1e-4 && Obs(32) > 1e-4,
            "k=31,32 (n=16): observed integrality residual above the fixed 1e-4 tolerance — " +
            "the verification test fails here, purely from the larger Vandermonde's " +
            "conditioning amplifying machine epsilon, not a bug.");

        // PROVEN: the k=31,32 deviation is degree-15 → degree-16 Vandermonde extraction
        // conditioning measured against an unscaled fixed 1e-4 tolerance. Plain algebra.
    }
}
