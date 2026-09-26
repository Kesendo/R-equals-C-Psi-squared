using RCPsiSquared.Core.BlockSpectrum;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Numerics;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Diagnostics.Foundation;
using RCPsiSquared.Diagnostics.Knowledge;

namespace RCPsiSquared.Diagnostics.Tests.Foundation;

public class CompressedDensityLocusClaimTests
{
    [Fact]
    public void ClaimStandsOnTheOwnersOfItsPremisesAndEndpoints()
    {
        var registry = KnowledgeRegistryFactory.BuildDefault();
        var claim = registry.Get<CompressedDensityLocusClaim>();

        Assert.Equal(Tier.Tier1Derived, claim.Tier);
        Assert.Same(registry.Get<AbsorptionTheoremClaim>(), claim.Absorption);
        Assert.Same(registry.Get<F71AntiPalindromicGammaSpectralInvariance>(), claim.Locus);
        Assert.Same(registry.Get<JointPopcountSectors>(), claim.Sectors);
        Assert.Same(registry.Get<SeedRungGramClaim>(), claim.SeedRungGram);
        Assert.Same(registry.Get<FrozenDivisorClaim>(), claim.FrozenDivisor);
        Assert.Contains("CompressedDensityN11Witness", claim.Anchor);
    }

    [Fact]
    public void PhysicalN11CellsExposeTheMixedParityContrastAndIdentityFailure()
    {
        var r = new CompressedDensityN11Witness().Reading;

        Assert.Equal(4, r.FrequencyMultiplicity);
        Assert.Equal(0, r.FrequencyMembershipMismatchCount);
        Assert.Equal(new BigRational(-1, 144) * Sqrt23.Sqrt2, r.LeftLocalCross);
        Assert.Equal(new BigRational(1, 144) * Sqrt23.Sqrt2, r.RightLocalCross);
        Assert.Equal(new BigRational(-1, 72) * Sqrt23.Sqrt2, r.ContrastCross);
        Assert.Equal(Sqrt23.Rational(1, 324), r.ContrastTraceSquare);
        Assert.Equal(new BigRational(1, 36) * Sqrt23.Sqrt2, r.BalancedPhysicalCross);
        Assert.Equal(Sqrt23.Zero, r.BalancedConditionalPredictionCross);
        Assert.Equal(
            new[]
            {
                Sqrt23.Rational(14392, 81), Sqrt23.Rational(15836, 81), Sqrt23.Rational(13031, 162),
                Sqrt23.Rational(44, 3), Sqrt23.One,
            },
            r.BalancedRoomCharacteristicPolynomial);
        Assert.All(r.Checks, check => Assert.True(check.Passes, check.Detail));
    }

    [Fact]
    public void ZeroFrequencyRoomIsF143sGramAndReachesBothEnds()
    {
        var r = new CompressedDensityN11Witness().Reading;

        Assert.Equal(11, r.ZeroFrequencyMultiplicity);
        Assert.Equal(0, r.ZeroFrequencyMembershipMismatchCount);
        Assert.Equal(Sqrt23.Zero, r.IdentityRayleigh);
        Assert.Equal(Sqrt23.Rational(-4), r.ChiralRayleigh);
        Assert.Contains(r.Checks, c => c.Name.StartsWith("zero-frequency room is F143") && c.Passes);
    }

    [Fact]
    public void EveryN11RoomCarriesTheLowerEndpointOncePerChiralTransposeOrbit()
    {
        var r = new CompressedDensityN11Witness().Reading;
        Assert.Equal(55, r.RoomCount);
        Assert.Equal(45, r.RoomsReachingLowerEndpoint);
        Assert.True(r.Checks.Single(c => c.Name.StartsWith("every room, balanced profile")).Passes);
        Assert.True(r.Checks.Single(c => c.Name.StartsWith("0 lies in the zero-frequency room only")).Passes);
    }

    [Fact]
    public void OffLocusRayleighIsExactAndBelowTheNominalInterval()
    {
        var r = new CompressedDensityN11Witness().Reading;
        Assert.Equal(new BigRational(-1, 72) * (Sqrt23.Rational(22) + new BigRational(5) * Sqrt23.Sqrt3),
            r.OffLocusRayleigh);
        Assert.Equal(-1, (r.OffLocusRayleigh + Sqrt23.Rational(4, 11)).Sign);
    }

    [Fact]
    public void ExactSignDecidesNearCancellations()
    {
        // 1 + √2 − √3 − √6·(1/10): positive; √2 + √3 − √6·(13/10): negative (√2+√3 ≈ 3.146, 1.3√6 ≈ 3.184).
        Assert.Equal(1, (Sqrt23.One + Sqrt23.Sqrt2 - Sqrt23.Sqrt3 - new BigRational(1, 10) * Sqrt23.Sqrt6).Sign);
        Assert.Equal(-1, (Sqrt23.Sqrt2 + Sqrt23.Sqrt3 - new BigRational(13, 10) * Sqrt23.Sqrt6).Sign);
        Assert.Equal(0, (Sqrt23.Sqrt2 * Sqrt23.Sqrt3 - Sqrt23.Sqrt6).Sign);
        // (√2+√3)² = 5 + 2√6 exactly; a 10⁻³⁰ nudge either way is read through the √3-level branch.
        var square = (Sqrt23.Sqrt2 + Sqrt23.Sqrt3) * (Sqrt23.Sqrt2 + Sqrt23.Sqrt3);
        var exact = square - Sqrt23.Rational(5) - new BigRational(2) * Sqrt23.Sqrt6;
        Assert.True(exact.IsZero);
        var tiny = new BigRational(1, System.Numerics.BigInteger.Pow(10, 30));
        Assert.Equal(1, (Sqrt23.Sqrt6 - Sqrt23.Sqrt2 * Sqrt23.Sqrt3 + tiny * Sqrt23.Sqrt3 - tiny * Sqrt23.Sqrt2).Sign);
        Assert.Equal(-1, (Sqrt23.Sqrt2 - Sqrt23.Sqrt3 + tiny * Sqrt23.Sqrt6).Sign);
        var x = Sqrt23.One + Sqrt23.Sqrt2 - new BigRational(3, 7) * Sqrt23.Sqrt3 + new BigRational(2, 5) * Sqrt23.Sqrt6;
        Assert.Equal(Sqrt23.One, x * x.Inverse());
    }

    [Fact]
    public void KetOnlyMutantLosesTheContrastAndBothEndpoints()
    {
        var physical = new CompressedDensityN11Witness().Reading;
        var mutant = new CompressedDensityN11Witness((ket, bra, site) => ket == site ? 1 : 0).Reading;

        Assert.True(physical.AllPass);
        Assert.Equal(Sqrt23.Zero, mutant.ContrastCross);
        Assert.NotEqual(Sqrt23.Zero, mutant.IdentityRayleigh);
        Assert.NotEqual(Sqrt23.Rational(-4), mutant.ChiralRayleigh);
        Assert.False(mutant.AllPass);
        foreach (var name in new[]
                 {
                     "identity reaches the upper endpoint", "chiral difference reaches the lower endpoint",
                     "zero-frequency room is F143's −4γbar(I − G), balanced profile",
                 })
            Assert.False(mutant.Checks.Single(c => c.Name == name).Passes, name);
    }

    [Fact]
    public void OrRuleMutantBreaksTheEndpointHalfOnItsOwn()
    {
        var mutant = new CompressedDensityN11Witness((ket, bra, site) => ket == site || bra == site ? 1 : 0).Reading;
        Assert.False(mutant.Checks.Single(c => c.Name == "identity reaches the upper endpoint").Passes);
        Assert.False(mutant.Checks.Single(c => c.Name == "mirror contrast cross").Passes);
    }
}
