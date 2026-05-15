using System.Linq;
using RCPsiSquared.Core.CoherenceBlocks;
using RCPsiSquared.Core.F86.Item1Derivation;
using RCPsiSquared.Core.Symmetry;
using Xunit;
using Xunit.Abstractions;

namespace RCPsiSquared.Core.Tests.Symmetry;

/// <summary>Tests for <see cref="F89AmplitudeLayerClaim"/>'s Angle A decomposition.
/// Verifies the path-3 algebraic anchor (33+14√5)/9 numerically against the
/// C2FullBlockSigmaAnatomy F_a witnesses, plus the structural identities
/// |S_c(2)|² = (10+4√5)/3 and ‖Mv(2)‖² = (25+4√5)/15 from
/// PROOF_F89_PATH_D_CLOSED_FORM.md Angle A.</summary>
public class F89AmplitudeLayerClaimTests
{
    private readonly ITestOutputHelper _out;
    public F89AmplitudeLayerClaimTests(ITestOutputHelper @out) => _out = @out;

    private static CoherenceBlock C2Block(int N) =>
        new CoherenceBlock(N: N, n: 1, gammaZero: 0.05);

    /// <summary>Path-3 closed-form values from
    /// PROOF_F89_PATH_D_CLOSED_FORM.md, evaluated numerically.</summary>
    private static class Path3 {
        public const double Sqrt5 = 2.2360679774997896964091736687747626778;
        public static double ScSquaredN2 => (10.0 + 4.0 * Sqrt5) / 3.0;
        public static double ScSquaredN4 => (10.0 - 4.0 * Sqrt5) / 3.0;
        public static double MvSquaredN2 => (25.0 + 4.0 * Sqrt5) / 15.0;
        public static double MvSquaredN4 => (25.0 - 4.0 * Sqrt5) / 15.0;
        public static double PnN2 => (33.0 + 14.0 * Sqrt5) / 9.0;
        public static double PnN4 => (33.0 - 14.0 * Sqrt5) / 9.0;
    }

    [Fact]
    public void ComputePn_AtPath3N2_FromAnatomy_MatchesAnchor()
    {
        // Anatomy directly produces sigma_2 at path-3 (chainN=4) via inverse iteration.
        var anatomy = C2FullBlockSigmaAnatomy.BuildFaOnly(C2Block(N: 4));
        var witness = anatomy.SigmaSpectrum.FirstOrDefault(w => w.BlochIndexN == 2);
        Assert.NotNull(witness);
        double pn = F89AmplitudeLayerClaim.ComputePn(witness!.Sigma, chainN: 4);
        Assert.Equal(Path3.PnN2, pn, precision: 8);
    }

    [Fact]
    public void Path3AnchorPn_ReturnsTuple_MatchingProofClosedForm()
    {
        var (rational2, sqrtCoef2, denom2) = F89AmplitudeLayerClaim.Path3AnchorPn(2);
        Assert.Equal(33, rational2);
        Assert.Equal(14, sqrtCoef2);
        Assert.Equal(9, denom2);

        var (rational4, sqrtCoef4, denom4) = F89AmplitudeLayerClaim.Path3AnchorPn(4);
        Assert.Equal(33, rational4);
        Assert.Equal(-14, sqrtCoef4);
        Assert.Equal(9, denom4);
    }

    [Fact]
    public void DecomposeAngleA_AtPath3FromAnatomy_RecoversScAndMvAnchors()
    {
        // Build C2 anatomy at path-3 (chainN=4)
        var anatomy = C2FullBlockSigmaAnatomy.BuildFaOnly(C2Block(N: 4));
        var faModes = anatomy.SigmaSpectrum.Where(w => w.BlochIndexN.HasValue).ToList();

        Assert.Equal(2, faModes.Count); // {n=2, n=4}

        foreach (var witness in faModes)
        {
            int n = witness.BlochIndexN!.Value;
            var (scSquared, mvSquared) = F89AmplitudeLayerClaim.DecomposeAngleA(
                witness.ProbeOverlapSquared,
                witness.SKernelDiagonal,
                chainN: 4);

            _out.WriteLine($"n={n}: |S_c|²={scSquared:G10}, ‖Mv‖²={mvSquared:G10}");

            double scExpected = n == 2 ? Path3.ScSquaredN2 : Path3.ScSquaredN4;
            double mvExpected = n == 2 ? Path3.MvSquaredN2 : Path3.MvSquaredN4;

            // Numerical tolerance: inverse-iteration extracts at ~1e-12 typical
            Assert.Equal(scExpected, scSquared, precision: 8);
            Assert.Equal(mvExpected, mvSquared, precision: 8);
        }
    }

    [Fact]
    public void VerifyAngleA_AtPath3FromAnatomy_ResidualBelowTolerance()
    {
        var anatomy = C2FullBlockSigmaAnatomy.BuildFaOnly(C2Block(N: 4));
        var faModes = anatomy.SigmaSpectrum.Where(w => w.BlochIndexN.HasValue).ToList();

        foreach (var witness in faModes)
        {
            var (scSquared, mvSquared) = F89AmplitudeLayerClaim.DecomposeAngleA(
                witness.ProbeOverlapSquared,
                witness.SKernelDiagonal,
                chainN: 4);
            double residual = F89AmplitudeLayerClaim.VerifyAngleA(
                witness.Sigma, chainN: 4, scSquared, mvSquared);
            _out.WriteLine($"n={witness.BlochIndexN}: residual={residual:E2}");
            // By construction the identity holds; this checks no FP drift exceeds tolerance.
            Assert.True(residual < F89AmplitudeLayerClaim.AngleATolerance,
                $"Angle A residual {residual:E2} exceeds tolerance {F89AmplitudeLayerClaim.AngleATolerance:E2}");
        }
    }

    [Theory]
    [InlineData(5, 2)]   // path-4: chainN=5, 2 F_a modes
    [InlineData(6, 3)]   // path-5: chainN=6, 3 F_a modes
    [InlineData(7, 3)]   // path-6: chainN=7, 3 F_a modes
    public void DecomposeAngleA_AtHigherPath_PreservesIdentityToHighPrecision(
        int chainN, int faCount)
    {
        var anatomy = C2FullBlockSigmaAnatomy.BuildFaOnly(C2Block(chainN));
        var faModes = anatomy.SigmaSpectrum.Where(w => w.BlochIndexN.HasValue).ToList();
        Assert.Equal(faCount, faModes.Count);

        foreach (var witness in faModes)
        {
            var (scSquared, mvSquared) = F89AmplitudeLayerClaim.DecomposeAngleA(
                witness.ProbeOverlapSquared,
                witness.SKernelDiagonal,
                chainN);
            double residual = F89AmplitudeLayerClaim.VerifyAngleA(
                witness.Sigma, chainN, scSquared, mvSquared);
            Assert.True(residual < F89AmplitudeLayerClaim.AngleATolerance,
                $"Angle A residual {residual:E2} at chainN={chainN}, n={witness.BlochIndexN} exceeds tolerance");
        }
    }
}
