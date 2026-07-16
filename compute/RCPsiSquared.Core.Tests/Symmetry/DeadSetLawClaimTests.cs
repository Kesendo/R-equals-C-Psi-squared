using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.BlockSpectrum;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Core.SymmetryFamily;
using Xunit;
using ComplexMatrix = MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex>;

namespace RCPsiSquared.Core.Tests.Symmetry;

/// <summary>Tests for <see cref="DeadSetLawClaim"/> (F132): the dead-set law. The
/// mirror-composition antiunitaries V_g = (U_g·X^N)∘conj are exact symmetries of the
/// one-sided-watched XY+field Lindblad flow at every fixed h, their kill sign is the
/// N-free mod-4 function of the Majorana degree (ε_odd = (−1)^(d(d−1)/2),
/// ε_even = (−1)^(d(d+1)/2)), and together with the popcount blocks and the conserved
/// degree the identically-zero-at-every-h readout set closes into one line:
/// alive ⟺ (K_pop ∧ d ≡ 0 mod 4) ∨ (coherence ∧ K_coh ∧ d = N). The claim carries an
/// N = 3 battery (matrix + moment level, no eigensolver): the flipped Hamiltonian, the
/// flow symmetry on a generic matrix, the mod-4 identity on all 64 strings three ways,
/// the automatic coherence signs, degree conservation of the adjoint flow, the full
/// 63-string moment census for BOTH preps, and the zz fence (free-world scope).</summary>
public class DeadSetLawClaimTests
{
    private static DeadSetLawClaim MakeClaim()
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
        var chiral = new ChiralKClaim();
        var f131 = new MirrorOrderSortingClaim(
            chiral,
            triangle,
            new F71AntiPalindromicGammaSpectralInvariance(sectors, f71),
            new F92BondAntiPalindromicJSpectralInvariance(sectors, f71, inventory),
            new F93DetuningAntiPalindromicSpectralInvariance(sectors, f71, inventory));
        return new DeadSetLawClaim(chiral, triangle, f131);
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
    public void Anchor_ReferencesOwnerGateAndRegistry()
    {
        var claim = MakeClaim();
        Assert.Contains("LATTICE_DEAD_SET_RULE.md", claim.Anchor);
        Assert.Contains("lattice_dead_set_rule.py", claim.Anchor);
        Assert.Contains("ANALYTICAL_FORMULAS.md", claim.Anchor);
    }

    [Fact]
    public void TypedParents_AreWired()
    {
        var claim = MakeClaim();
        Assert.NotNull(claim.ChiralK);
        Assert.NotNull(claim.Triangle);
        Assert.NotNull(claim.OrderSorting);
    }

    [Fact]
    public void Constructor_RejectsNullParents()
    {
        var full = MakeClaim();
        Assert.Throws<ArgumentNullException>(() => new DeadSetLawClaim(null!, full.Triangle, full.OrderSorting));
        Assert.Throws<ArgumentNullException>(() => new DeadSetLawClaim(full.ChiralK, null!, full.OrderSorting));
        Assert.Throws<ArgumentNullException>(() => new DeadSetLawClaim(full.ChiralK, full.Triangle, null!));
    }

    [Fact]
    public void NotAnIZ2AxisClaim_CrossAxisStructural()
    {
        Assert.False(typeof(IZ2AxisClaim).IsAssignableFrom(typeof(DeadSetLawClaim)));
    }

    // ------------------------------------------------------------------
    // Self-check battery (N = 3, matrix + moment level, built in the ctor)
    // ------------------------------------------------------------------

    [Fact]
    public void Battery_AllCasesPass()
    {
        var claim = MakeClaim();
        Assert.Equal(8, claim.Cases.Count);
        foreach (var c in claim.Cases)
            Assert.True(c.Passes, $"battery case '{c.Name}' failed: expected {c.Expected}, got {c.Actual}");
        Assert.Equal(claim.Cases.Count, claim.PassCount);
    }

    [Fact]
    public void Battery_ComposedMirror_FlipsTheHamiltonian()
    {
        var claim = MakeClaim();
        Assert.True(claim.Cases.Single(c => c.Name.Contains("flips H")).Passes);
    }

    [Fact]
    public void Battery_FlowSymmetry_HoldsForBothGauges()
    {
        var claim = MakeClaim();
        Assert.True(claim.Cases.Single(c => c.Name.Contains("flow symmetry")).Passes);
    }

    [Fact]
    public void Battery_Mod4Identity_AllStringsThreeWays()
    {
        var claim = MakeClaim();
        Assert.True(claim.Cases.Single(c => c.Name.Contains("mod-4 identity")).Passes);
    }

    [Fact]
    public void Battery_AutomaticCoherenceSigns()
    {
        var claim = MakeClaim();
        Assert.True(claim.Cases.Single(c => c.Name.Contains("automatic")).Passes);
    }

    [Fact]
    public void Battery_DegreeConservation()
    {
        var claim = MakeClaim();
        Assert.True(claim.Cases.Single(c => c.Name.Contains("degree")
            && c.Name.Contains("preserves")).Passes);
    }

    [Fact]
    public void Battery_MomentCensus_BothPreps()
    {
        var claim = MakeClaim();
        Assert.True(claim.Cases.Single(c => c.Name.Contains("population prep")).Passes);
        Assert.True(claim.Cases.Single(c => c.Name.Contains("coherence prep")).Passes);
    }

    [Fact]
    public void Battery_ZzFence_FreeWorldScope()
    {
        var claim = MakeClaim();
        Assert.True(claim.Cases.Single(c => c.Name.Contains("zz fence")).Passes);
    }

    [Fact]
    public void Summary_CarriesTheOneLinerAndTheBattery()
    {
        var claim = MakeClaim();
        Assert.Contains("d ≡ 0 mod 4", claim.Summary);
        Assert.Contains($"{claim.Cases.Count}/{claim.Cases.Count} battery PASS", claim.Summary);
    }

    // ------------------------------------------------------------------
    // Direct mathematical spot-checks, independent of the claim's battery
    // code.
    // ------------------------------------------------------------------

    /// <summary>The mod-4 formulas at the four residues, hand-evaluated: the odd-sites
    /// gauge kills d ≡ 2, 3 (mod 4); the even-sites gauge kills d ≡ 1, 2 (mod 4).</summary>
    [Fact]
    public void EpsDegree_KillClasses_HandChecked()
    {
        // eps_odd = (−1)^(d(d−1)/2): d = 0,1 → +1; d = 2,3 → −1 (period 4).
        // eps_even = (−1)^(d(d+1)/2): d = 0,3 → +1; d = 1,2 → −1 (period 4).
        int[] expectedOdd = { +1, +1, -1, -1 };
        int[] expectedEven = { +1, -1, -1, +1 };
        for (int d = 0; d <= 11; d++)
        {
            Assert.Equal(expectedOdd[d % 4], DeadSetLawClaim.EpsDegree(d, evenGauge: false));
            Assert.Equal(expectedEven[d % 4], DeadSetLawClaim.EpsDegree(d, evenGauge: true));
        }
    }

    /// <summary>The Majorana-degree recursion on hand-computable strings: a bare X at
    /// site l has degree 2l + 1 (l Z-pairs below the cap); the Jordan-Wigner string
    /// Z…ZX (tail left) has degree 1; pure-Z strings have even degree 2·|T| when no
    /// caps sit above; the all-Z string is the full parity, degree 2N.</summary>
    [Fact]
    public void MajoranaDegree_HandChecked()
    {
        Assert.Equal(1, DeadSetLawClaim.MajoranaDegree("XII"));
        Assert.Equal(3, DeadSetLawClaim.MajoranaDegree("IXI"));
        Assert.Equal(5, DeadSetLawClaim.MajoranaDegree("IIX"));
        Assert.Equal(1, DeadSetLawClaim.MajoranaDegree("ZZX"));   // the JW string a_k
        Assert.Equal(2, DeadSetLawClaim.MajoranaDegree("ZII"));
        Assert.Equal(6, DeadSetLawClaim.MajoranaDegree("ZZZ"));   // Γ = the full parity
        Assert.Equal(3, DeadSetLawClaim.MajoranaDegree("XXX"));   // one cap per site
        Assert.Equal(4, DeadSetLawClaim.MajoranaDegree("XIX"));
    }

    /// <summary>The single-site Z on an unwatched site is a doubly-mirrored zero: at
    /// N = 3 with the population prep, every moment Tr[Z₀·L(h)^k·ρ₀] vanishes although
    /// its popcount blocks are populated and its degree sector is fed (d = 2). Direct
    /// hand-built check, independent of the battery's census code.</summary>
    [Fact]
    public void SingleSiteZ_AllMomentsVanish_DirectCheck()
    {
        var sx = Dense(new Complex[,] { { 0, 1 }, { 1, 0 } });
        var sy = Dense(new Complex[,] { { 0, new Complex(0, -1) }, { new Complex(0, 1), 0 } });
        var sz = Dense(new Complex[,] { { 1, 0 }, { 0, -1 } });
        var id = ComplexMatrix.Build.DenseIdentity(2);
        const int d = 8;

        ComplexMatrix Site(ComplexMatrix single, int site)
        {
            var f0 = site == 0 ? single : id;
            var f1 = site == 1 ? single : id;
            var f2 = site == 2 ? single : id;
            return f0.KroneckerProduct(f1).KroneckerProduct(f2);
        }

        var h = Site(sx, 0) * Site(sx, 1) + Site(sy, 0) * Site(sy, 1)
              + Site(sx, 1) * Site(sx, 2) + Site(sy, 1) * Site(sy, 2)
              + Site(sz, 0).Multiply(new Complex(0.4, 0))
              + Site(sz, 1).Multiply(new Complex(-0.24, 0))
              + Site(sz, 2).Multiply(new Complex(0.12, 0));
        var z2 = Site(sz, 2);
        double gamma = 0.5;

        ComplexMatrix Rhs(ComplexMatrix rho)
        {
            var comm = (h * rho - rho * h).Multiply(new Complex(0, -1));
            // one-sided watching on the last site: D[Z₂]ρ = γ(Z₂ρZ₂ − ρ)
            return comm + (z2 * rho * z2 - rho).Multiply(new Complex(gamma, 0));
        }

        var rho0 = ComplexMatrix.Build.Dense(d, d);
        rho0[1, 1] = new Complex(0.5, 0);
        rho0[6, 6] = new Complex(0.5, 0);

        var z0 = Site(sz, 0);
        var w = rho0;
        for (int k = 0; k <= 6; k++)
        {
            Complex tr = Complex.Zero;
            for (int i = 0; i < d; i++)
                for (int j = 0; j < d; j++)
                    tr += z0[i, j] * w[j, i];
            Assert.True(tr.Magnitude < 1e-10 * Math.Max(1.0, w.FrobeniusNorm()),
                $"moment k={k} of Z₀ must vanish, got {tr.Magnitude}");
            w = Rhs(w);
        }
    }

    private static ComplexMatrix Dense(Complex[,] cells) => ComplexMatrix.Build.DenseOfArray(cells);
}
