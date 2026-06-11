using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.F1;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;
using Xunit;
using ComplexMatrix = MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex>;

namespace RCPsiSquared.Core.Tests.Symmetry;

public class MomentTowerPumpChannelClaimTests
{
    private static (LindbladBitBPiBreakMagnitude F113, F84ThermalAmplitudeDampingPi2Inheritance F84) MakeParents()
    {
        var part2 = new F108Part2Pi2XEvenAlwaysPalindromic();
        var f112 = new LindbladBitBPiBalance(
            new F108Part1Pi2EvenAlwaysPalindromic(part2),
            new LindbladBitAPiBalance(part2));
        var f113 = new LindbladBitBPiBreakMagnitude(f112);

        var ladder = new Pi2DyadicLadderClaim();
        var f81 = new F81Pi2Inheritance(
            new F1PalindromeIdentity(),
            ladder,
            new Pi2OperatorSpaceMirrorClaim(),
            new Pi2I4MemoryLoopClaim());
        var f84 = new F84ThermalAmplitudeDampingPi2Inheritance(
            ladder,
            new F82T1AmplitudeDampingPi2Inheritance(ladder, f81));

        return (f113, f84);
    }

    private static MomentTowerPumpChannelClaim MakeClaim()
    {
        var (f113, f84) = MakeParents();
        return new MomentTowerPumpChannelClaim(f113, f84);
    }

    // ------------------------------------------------------------------
    // Claim metadata
    // ------------------------------------------------------------------

    [Fact]
    public void Tier_IsTier1Derived()
    {
        var claim = MakeClaim();
        Assert.Equal(Tier.Tier1Derived, claim.Tier);
    }

    [Fact]
    public void Anchor_ReferencesProofAndVerifier()
    {
        var claim = MakeClaim();
        Assert.Contains("PROOF_MOMENT_TOWER_PUMP_CHANNEL.md", claim.Anchor);
        Assert.Contains("moment_tower_pump_channel.py", claim.Anchor);
    }

    [Fact]
    public void TypedParents_AreWired()
    {
        var claim = MakeClaim();
        Assert.NotNull(claim.F113);
        Assert.NotNull(claim.F84);
        // The F113 parent carries its own F112 parent (the typed bit_b chain stays intact).
        Assert.NotNull(claim.F113.Parent);
    }

    [Fact]
    public void Constructor_RejectsNullParents()
    {
        var (f113, f84) = MakeParents();
        Assert.Throws<ArgumentNullException>(() => new MomentTowerPumpChannelClaim(null!, f84));
        Assert.Throws<ArgumentNullException>(() => new MomentTowerPumpChannelClaim(f113, null!));
    }

    [Fact]
    public void NotAnIZ2AxisClaim_CubeMapCountsUnchanged()
    {
        // Pinned invariant: the channel is cross-axis structural like AntilinearTriangleClaim
        // and MirrorGroupD4Claim; PolarityCubeMap counts stay untouched.
        Assert.False(typeof(IZ2AxisClaim).IsAssignableFrom(typeof(MomentTowerPumpChannelClaim)));
    }

    // ------------------------------------------------------------------
    // Self-check battery (N = 2 and N = 3 dense operators, built in the ctor)
    // ------------------------------------------------------------------

    [Fact]
    public void Battery_AllCasesPass()
    {
        var claim = MakeClaim();
        Assert.Equal(10, claim.Cases.Count);
        foreach (var c in claim.Cases)
            Assert.True(c.Passes, $"battery case '{c.Name}' failed: expected {c.Expected}, got {c.Actual}");
        Assert.Equal(claim.Cases.Count, claim.PassCount);
    }

    [Fact]
    public void Battery_PumpDirections_AreExact()
    {
        var claim = MakeClaim();
        var pump = claim.Cases.Single(c => c.Name.StartsWith("pump directions", StringComparison.Ordinal));
        Assert.True(pump.Passes, $"pump directions out of tolerance: {pump.Actual}");
    }

    [Fact]
    public void Battery_F113Bridge_UsesTypedParentClosedForm()
    {
        var claim = MakeClaim();
        var bridge = claim.Cases.Single(c => c.Name.StartsWith("F113 bridge", StringComparison.Ordinal));
        Assert.True(bridge.Passes, $"F113 bridge out of tolerance: {bridge.Actual}");
    }

    [Fact]
    public void Battery_CurvatureFingerprint_VisibleAndBlind()
    {
        var claim = MakeClaim();
        var visible = claim.Cases.Single(c => c.Name.StartsWith("curvature fingerprint", StringComparison.Ordinal));
        var blind = claim.Cases.Single(c => c.Name.StartsWith("curvature blind spot", StringComparison.Ordinal));
        Assert.True(visible.Passes, $"Y₀ fingerprint failed: {visible.Actual}");
        Assert.True(blind.Passes, $"Z₀ invisibility failed: {blind.Actual}");
    }

    [Fact]
    public void Battery_GirthTwoWitness_SilentThenFiring()
    {
        var claim = MakeClaim();
        var witness = claim.Cases.Single(c => c.Name.StartsWith("girth-2 witness", StringComparison.Ordinal));
        Assert.True(witness.Passes, $"girth-2 witness failed: {witness.Actual}");
    }

    [Fact]
    public void Summary_CarriesTheLawAndTheBridge()
    {
        var claim = MakeClaim();
        Assert.Contains("Δγ_l", claim.Summary);
        Assert.Contains("−4^N·slope⟨H⟩", claim.Summary);
        Assert.Contains($"{claim.Cases.Count}/{claim.Cases.Count} battery PASS", claim.Summary);
    }

    // ------------------------------------------------------------------
    // Static helper contracts
    // ------------------------------------------------------------------

    [Fact]
    public void PredictPumpSlope_ValidatesArguments()
    {
        var h = ComplexMatrix.Build.DenseIdentity(4);
        Assert.Throws<ArgumentNullException>(() =>
            MomentTowerPumpChannelClaim.PredictPumpSlope(null!, 1, new[] { 0.1, 0.2 }));
        Assert.Throws<ArgumentNullException>(() =>
            MomentTowerPumpChannelClaim.PredictPumpSlope(h, 1, null!));
        Assert.Throws<ArgumentOutOfRangeException>(() =>
            MomentTowerPumpChannelClaim.PredictPumpSlope(h, 0, new[] { 0.1, 0.2 }));
        // 4×4 matrix against N = 3 rates: dimension mismatch.
        Assert.Throws<ArgumentException>(() =>
            MomentTowerPumpChannelClaim.PredictPumpSlope(h, 1, new[] { 0.1, 0.2, 0.3 }));
    }

    [Fact]
    public void MomentTower_ReadsSingleZDrive()
    {
        // H = Z₀ at N = 2: t₁(0) = Tr(Z₀·Z₀) = 4, t₁(1) = Tr(Z₁·Z₀) = 0.
        var h = ComplexMatrix.Build.Dense(4, 4);
        h[0, 0] = Complex.One; h[1, 1] = Complex.One;
        h[2, 2] = -Complex.One; h[3, 3] = -Complex.One;
        var t1 = MomentTowerPumpChannelClaim.MomentTower(h, 1, 2);
        Assert.Equal(4.0, t1[0], 12);
        Assert.Equal(0.0, t1[1], 12);
    }

    // ------------------------------------------------------------------
    // Direct mathematical spot-check, independent of the claim's battery:
    // the slope law recomputed by applying L directly to dense matrices
    // (no superoperator, different H and rates than the battery uses).
    // ------------------------------------------------------------------

    [Fact]
    public void SlopeLaw_DirectMatrixCheck()
    {
        // Single-site operators (site 0 = leftmost tensor factor).
        var i2 = ComplexMatrix.Build.DenseIdentity(2);
        var x = ComplexMatrix.Build.Dense(2, 2);
        x[0, 1] = Complex.One; x[1, 0] = Complex.One;
        var y = ComplexMatrix.Build.Dense(2, 2);
        y[0, 1] = new Complex(0, -1); y[1, 0] = new Complex(0, 1);
        var z = ComplexMatrix.Build.Dense(2, 2);
        z[0, 0] = Complex.One; z[1, 1] = -Complex.One;
        var sm = ComplexMatrix.Build.Dense(2, 2);
        sm[0, 1] = Complex.One;                       // σ⁻ = [[0, 1], [0, 0]]
        var sp = sm.ConjugateTranspose();             // σ⁺

        ComplexMatrix At0(ComplexMatrix op) => op.KroneckerProduct(i2);
        ComplexMatrix At1(ComplexMatrix op) => i2.KroneckerProduct(op);

        var h = At0(z).Multiply(new Complex(0.4, 0))
              + At0(x) * At1(x) * new Complex(0.9, 0)
              + At0(y) * At1(y) * new Complex(0.6, 0)
              + At1(z).Multiply(new Complex(0.25, 0));

        double[] gDeph = { 0.11, 0.06 };
        double[] gDown = { 0.03, 0.015 };
        double[] gUp = { 0.007, 0.002 };
        double[] delta = { gDown[0] - gUp[0], gDown[1] - gUp[1] };

        ComplexMatrix Diss(ComplexMatrix c, ComplexMatrix rho)
        {
            var cdc = c.ConjugateTranspose() * c;
            return c * rho * c.ConjugateTranspose()
                 - (cdc * rho + rho * cdc).Multiply(new Complex(0.5, 0));
        }

        ComplexMatrix L(ComplexMatrix rho)
        {
            var minusI = new Complex(0, -1);
            var result = (h * rho - rho * h).Multiply(minusI);
            var ops0 = new[] { (gDeph[0], At0(z)), (gDown[0], At0(sm)), (gUp[0], At0(sp)) };
            var ops1 = new[] { (gDeph[1], At1(z)), (gDown[1], At1(sm)), (gUp[1], At1(sp)) };
            foreach (var (rate, c) in ops0.Concat(ops1))
                result += Diss(c, rho).Multiply(new Complex(rate, 0));
            return result;
        }

        var mixed = ComplexMatrix.Build.DenseIdentity(4).Multiply(new Complex(0.25, 0));
        var rhoDot = L(mixed);

        var power = ComplexMatrix.Build.DenseIdentity(4);
        for (int j = 1; j <= 3; j++)
        {
            power *= h;
            double slopeDirect = (power * rhoDot).Trace().Real;
            double slopePredicted = MomentTowerPumpChannelClaim.PredictPumpSlope(h, j, delta);
            Assert.True(Math.Abs(slopeDirect - slopePredicted) <= 1e-12,
                $"slope law mismatch at j = {j}: direct {slopeDirect:E17} vs predicted {slopePredicted:E17}");
        }

        // The law is load-bearing: zeroing the net pump kills the slope on the same data.
        double balanced = MomentTowerPumpChannelClaim.PredictPumpSlope(h, 1, new[] { 0.0, 0.0 });
        Assert.Equal(0.0, balanced, 15);
    }
}
