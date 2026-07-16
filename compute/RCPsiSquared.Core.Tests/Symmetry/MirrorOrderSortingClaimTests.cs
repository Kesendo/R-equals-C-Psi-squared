using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.BlockSpectrum;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Core.SymmetryFamily;
using Xunit;
using ComplexMatrix = MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex>;

namespace RCPsiSquared.Core.Tests.Symmetry;

/// <summary>Tests for <see cref="MirrorOrderSortingClaim"/> (F131): mirror conjugation
/// reflects a parameter scan, M·G(x₀ + s·δ)·M⁻¹ = G(x₀ + σ_eff·s·δ) with
/// σ_eff = σ_op·χ_M, and for a readout of definite mirror parity q the response orders
/// sort by q·σ_eff into four cells (generic / EVEN / ODD / IDENTICALLY ZERO). The claim
/// carries Theorem A (the unitary column, F71 site reversal, unconditional) in its
/// self-check battery at N = 3, moment-level: the exact conjugation identity gives the
/// per-order parity Tr[O·L(−t)^k·ρ₀] = q·Tr[O·L(+t)^k·ρ₀] for ALL k (all times by
/// analyticity); the battery spot-checks k = 0..6 (no eigensolver and no expm).</summary>
public class MirrorOrderSortingClaimTests
{
    private static MirrorOrderSortingClaim MakeClaim()
    {
        var kleinV4 = new Pi2KleinV4DephaseSwapGroup();
        var f114 = new CommutatorDConjugationSign(kleinV4);
        var mirror = new MirrorGroupD4Claim(
            new KleinEightCellClaim(new KleinFourCellClaim()),
            f114,
            kleinV4);
        var part2 = new F108Part2Pi2XEvenAlwaysPalindromic();
        var f112 = new LindbladBitBPiBalance(
            new F108Part1Pi2EvenAlwaysPalindromic(part2),
            new LindbladBitAPiBalance(part2));
        var triangle = new AntilinearTriangleClaim(mirror, f114, f112);

        var sectors = new JointPopcountSectors();
        var f71 = new F71MirrorBlockRefinement(sectors);
        var inventory = new SymmetryFamilyInventory();
        return new MirrorOrderSortingClaim(
            new ChiralKClaim(),
            triangle,
            new F71AntiPalindromicGammaSpectralInvariance(sectors, f71),
            new F92BondAntiPalindromicJSpectralInvariance(sectors, f71, inventory),
            new F93DetuningAntiPalindromicSpectralInvariance(sectors, f71, inventory));
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
    public void Anchor_ReferencesProofGateAndAdoption()
    {
        var claim = MakeClaim();
        Assert.Contains("PROOF_MIRROR_ORDER_SORTING.md", claim.Anchor);
        Assert.Contains("PROOF_ZETA2_ANTI_PROTECTION.md", claim.Anchor);
        Assert.Contains("mirror_order_sorting.py", claim.Anchor);
        Assert.Contains("OrderSorting.cs", claim.Anchor);
    }

    [Fact]
    public void TypedParents_AreWired()
    {
        var claim = MakeClaim();
        Assert.NotNull(claim.ChiralK);
        Assert.NotNull(claim.Triangle);
        Assert.NotNull(claim.F91);
        Assert.NotNull(claim.F92);
        Assert.NotNull(claim.F93);
    }

    [Fact]
    public void Constructor_RejectsNullParents()
    {
        var kleinV4 = new Pi2KleinV4DephaseSwapGroup();
        var f114 = new CommutatorDConjugationSign(kleinV4);
        var mirror = new MirrorGroupD4Claim(
            new KleinEightCellClaim(new KleinFourCellClaim()), f114, kleinV4);
        var part2 = new F108Part2Pi2XEvenAlwaysPalindromic();
        var f112 = new LindbladBitBPiBalance(
            new F108Part1Pi2EvenAlwaysPalindromic(part2),
            new LindbladBitAPiBalance(part2));
        var triangle = new AntilinearTriangleClaim(mirror, f114, f112);
        var sectors = new JointPopcountSectors();
        var f71 = new F71MirrorBlockRefinement(sectors);
        var inventory = new SymmetryFamilyInventory();
        var chiral = new ChiralKClaim();
        var f91 = new F71AntiPalindromicGammaSpectralInvariance(sectors, f71);
        var f92 = new F92BondAntiPalindromicJSpectralInvariance(sectors, f71, inventory);
        var f93 = new F93DetuningAntiPalindromicSpectralInvariance(sectors, f71, inventory);

        Assert.Throws<ArgumentNullException>(() => new MirrorOrderSortingClaim(null!, triangle, f91, f92, f93));
        Assert.Throws<ArgumentNullException>(() => new MirrorOrderSortingClaim(chiral, null!, f91, f92, f93));
        Assert.Throws<ArgumentNullException>(() => new MirrorOrderSortingClaim(chiral, triangle, null!, f92, f93));
        Assert.Throws<ArgumentNullException>(() => new MirrorOrderSortingClaim(chiral, triangle, f91, null!, f93));
        Assert.Throws<ArgumentNullException>(() => new MirrorOrderSortingClaim(chiral, triangle, f91, f92, null!));
    }

    [Fact]
    public void NotAnIZ2AxisClaim_CrossAxisStructural()
    {
        Assert.False(typeof(IZ2AxisClaim).IsAssignableFrom(typeof(MirrorOrderSortingClaim)));
    }

    // ------------------------------------------------------------------
    // Self-check battery (N = 3, moment-level, built in the ctor)
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
    public void Battery_ConjugationIdentity_AllThreeAxes()
    {
        var claim = MakeClaim();
        foreach (var axis in new[] { "gamma axis", "J axis", "h axis" })
        {
            var c = claim.Cases.Single(x => x.Name.Contains(axis));
            Assert.True(c.Passes, $"conjugation identity on the {axis} failed: {c.Actual}");
        }
    }

    [Fact]
    public void Battery_TheFence_RejectsAMixedScan()
    {
        var claim = MakeClaim();
        var fence = claim.Cases.Single(c => c.Name.Contains("fence"));
        Assert.True(fence.Passes, $"the fence must reject at O(1): {fence.Actual}");
    }

    [Fact]
    public void Battery_ResponseCells_SortByReadoutParity()
    {
        var claim = MakeClaim();
        Assert.True(claim.Cases.Single(c => c.Name.Contains("even cell")).Passes);
        Assert.True(claim.Cases.Single(c => c.Name.Contains("odd cell")).Passes);
        Assert.True(claim.Cases.Single(c => c.Name.Contains("zero cell")).Passes);
    }

    [Fact]
    public void Battery_Leak_IsExactlyAffine()
    {
        var claim = MakeClaim();
        var leak = claim.Cases.Single(c => c.Name.Contains("leak"));
        Assert.True(leak.Passes, $"the eps-leak must halve exactly: {leak.Actual}");
    }

    [Fact]
    public void Summary_CarriesTheFourCellsAndTheBattery()
    {
        var claim = MakeClaim();
        Assert.Contains("q", claim.Summary);
        Assert.Contains("σ_eff", claim.Summary);
        Assert.Contains($"{claim.Cases.Count}/{claim.Cases.Count} battery PASS", claim.Summary);
    }

    // ------------------------------------------------------------------
    // Direct mathematical spot-check, independent of the claim's battery
    // code: the third-moment odd-cell parity at N = 3 with hand-built
    // operators. k = 3 is the lowest non-vacuous order for a diagonal
    // readout (the dissipator acts only off-diagonal, so pure-D chains
    // read nothing on the diagonal; the first scan-carrying chain is
    // L_H·D·L_H). N = 2 would be degenerate here: (|e₀⟩ + |e₁⟩)/√2 is an
    // eigenvector of the two-site XX+YY bond, so the whole moment tower
    // is stationary there.
    // ------------------------------------------------------------------

    [Fact]
    public void ThirdMoment_OddReadout_IsOddInTheScan_DirectCheck()
    {
        var sx = Dense(new Complex[,] { { 0, 1 }, { 1, 0 } });
        var sy = Dense(new Complex[,] { { 0, new Complex(0, -1) }, { new Complex(0, 1), 0 } });
        var sz = Dense(new Complex[,] { { 1, 0 }, { 0, -1 } });
        var id = ComplexMatrix.Build.DenseIdentity(2);
        int d = 8;

        ComplexMatrix Site(ComplexMatrix single, int site)
        {
            var f0 = site == 0 ? single : id;
            var f1 = site == 1 ? single : id;
            var f2 = site == 2 ? single : id;
            return f0.KroneckerProduct(f1).KroneckerProduct(f2);
        }

        var h = Site(sx, 0) * Site(sx, 1) + Site(sy, 0) * Site(sy, 1)
              + Site(sx, 1) * Site(sx, 2) + Site(sy, 1) * Site(sy, 2);
        var z0 = Site(sz, 0);
        var z2 = Site(sz, 2);
        var idD = ComplexMatrix.Build.DenseIdentity(d);
        var idD2 = ComplexMatrix.Build.DenseIdentity(d * d);
        var minusI = new Complex(0, -1);

        ComplexMatrix Superop(double g0, double g2) =>
            (h.KroneckerProduct(idD) - idD.KroneckerProduct(h.Transpose())).Multiply(minusI)
            + (z0.KroneckerProduct(z0.Transpose()) - idD2).Multiply(new Complex(g0, 0))
            + (Site(sz, 1).KroneckerProduct(Site(sz, 1).Transpose()) - idD2).Multiply(new Complex(0.08, 0))
            + (z2.KroneckerProduct(z2.Transpose()) - idD2).Multiply(new Complex(g2, 0));

        // psi = (|e0> + |e2>)/sqrt2; with site 0 the leftmost kron factor, e0 = 4 and e2 = 1.
        var vec = MathNet.Numerics.LinearAlgebra.Vector<Complex>.Build.Dense(d * d);
        int u = 4, v = 1;
        vec[u * d + u] = new Complex(0.5, 0);
        vec[v * d + v] = new Complex(0.5, 0);
        vec[u * d + v] = new Complex(0.5, 0);
        vec[v * d + u] = new Complex(0.5, 0);

        var oOdd = z0 - z2;
        double Moment(double g0, double g2)
        {
            var sup = Superop(g0, g2);
            var w = sup * (sup * (sup * vec));
            double m = 0;
            for (int i = 0; i < d; i++)
                for (int j = 0; j < d; j++)
                    m += (oOdd[i, j] * w[j * d + i]).Real;
            return m;
        }

        const double baseG = 0.06, dt = 0.03;
        double plus = Moment(baseG - dt, baseG + dt);
        double minus = Moment(baseG + dt, baseG - dt);
        Assert.True(Math.Abs(plus) > 1e-6, $"odd third moment must be non-vacuous: {plus}");
        Assert.True(Math.Abs(plus + minus) < 1e-12, $"odd third moment must be odd in t: {plus} vs {minus}");
    }

    private static ComplexMatrix Dense(Complex[,] cells) => ComplexMatrix.Build.DenseOfArray(cells);
}
