using System.Reflection;
using RCPsiSquared.Core.BlockSpectrum;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Diagnostics.Foundation;

namespace RCPsiSquared.Diagnostics.Tests.Foundation;

public class MissingPhaseRelaxationScaleClaimTests
{
    private static (JDefectLightMigrationClaim Light, SeatCutBlindnessClaim Seat, JointPopcountSectors Sectors)
        Parents() =>
        (
            new JDefectLightMigrationClaim(),
            new SeatCutBlindnessClaim(new F4KernelDimensionByComponentsClaim()),
            new JointPopcountSectors()
        );

    [Fact]
    public void Constructor_HoldsTheThreeTypedParents()
    {
        var lightMigration = new JDefectLightMigrationClaim();
        var seatBlindness = new SeatCutBlindnessClaim(new F4KernelDimensionByComponentsClaim());
        var sectors = new JointPopcountSectors();

        var claim = new MissingPhaseRelaxationScaleClaim(lightMigration, seatBlindness, sectors);

        Assert.Same(lightMigration, claim.LightMigration);
        Assert.Same(seatBlindness, claim.SeatBlindness);
        Assert.Same(sectors, claim.JointSectors);
    }

    [Fact]
    public void Constructor_RejectsEveryMissingParent()
    {
        var (light, seat, sectors) = Parents();

        Assert.Equal("lightMigration", Assert.Throws<ArgumentNullException>(
            () => new MissingPhaseRelaxationScaleClaim(null!, seat, sectors)).ParamName);
        Assert.Equal("seatBlindness", Assert.Throws<ArgumentNullException>(
            () => new MissingPhaseRelaxationScaleClaim(light, null!, sectors)).ParamName);
        Assert.Equal("jointSectors", Assert.Throws<ArgumentNullException>(
            () => new MissingPhaseRelaxationScaleClaim(light, seat, null!)).ParamName);
    }

    [Fact]
    public void Claim_OwnsOnlyTheExactN7LocalGermAndCompletePeripheralSplit()
    {
        var (light, seat, sectors) = Parents();
        var claim = new MissingPhaseRelaxationScaleClaim(light, seat, sectors);

        Assert.Equal(Tier.Tier1Derived, claim.Tier);
        Assert.Equal(7, MissingPhaseRelaxationScaleClaim.SiteCount);
        Assert.Equal(10, MissingPhaseRelaxationScaleClaim.PeripheralDimension);
        Assert.Equal(2, MissingPhaseRelaxationScaleClaim.PerFrequencySlowClusterRank);
        Assert.Equal(4, MissingPhaseRelaxationScaleClaim.SlowRateMultiplicity);
        Assert.Equal(2, MissingPhaseRelaxationScaleClaim.PuncturedKernelDimension);
        Assert.Equal(0.5, MissingPhaseRelaxationScaleClaim.QuadraticCoefficientInGamma);
        Assert.Equal(0.5, MissingPhaseRelaxationScaleClaim.CubicCoefficientInGamma);
        Assert.Contains("Delta_A(epsilon) = (gamma/2) epsilon^2 + (gamma/2) epsilon^3 + O_gamma(epsilon^4)",
            claim.GapGerm);
        Assert.Contains("w_c = epsilon^2/4 + epsilon^3/4 + O_gamma(epsilon^4)", claim.CentreLightGerm);
        Assert.Contains("two conjugate rank-2 clusters", claim.PeripheralSplit);
        Assert.Contains("four directions", claim.PeripheralSplit);
        Assert.Contains("10-dimensional", claim.PeripheralSplit);
        Assert.Contains("2-dimensional", claim.PeripheralSplit);
        Assert.Contains("locally slowest", claim.PeripheralSplit);
    }

    [Fact]
    public void Claim_ExposesAllThreeTypedParentsAsInspectableChildren()
    {
        var (light, seat, sectors) = Parents();
        var claim = new MissingPhaseRelaxationScaleClaim(light, seat, sectors);
        var children = claim.Children.ToList();

        Assert.Contains(children, child => ReferenceEquals(light, child));
        Assert.Contains(children, child => ReferenceEquals(seat, child));
        Assert.Contains(children, child => ReferenceEquals(sectors, child));
    }

    [Fact]
    public void Narrative_PinsTheLocalScopeAndTheExistingGeneralMechanism()
    {
        var (light, seat, sectors) = Parents();
        var claim = new MissingPhaseRelaxationScaleClaim(light, seat, sectors);
        string narrative = string.Join("\n", new[]
        {
            claim.Name,
            claim.Summary,
            claim.Scope,
            claim.GeneralMechanism,
            claim.Anchor,
            string.Join("\n", claim.Children.Select(child => child.Summary))
        });

        Assert.Contains("N = 7", narrative);
        Assert.Contains("A block", narrative);
        Assert.Contains("joint-popcount (1,1)", narrative);
        Assert.Contains("XY path", narrative);
        Assert.Contains("c = 3", narrative);
        Assert.Contains("h_01 = 2(1 + epsilon)", narrative);
        Assert.Contains("h_j,j+1 = 2 for j = 1,...,5", narrative);
        Assert.Contains("centre-only gamma > 0", narrative);
        Assert.Contains("fixed-gamma punctured epsilon neighbourhood", narrative);
        Assert.Contains("not an all-odd-N theorem", narrative);
        Assert.Contains("not the full 4^7 Liouvillian", narrative);
        Assert.Contains("not an observable lifetime", narrative);
        Assert.Contains(nameof(JDefectLightMigrationClaim), narrative);
        Assert.Contains("general Re-drift/light mechanism", narrative);
        Assert.Contains("docs/proofs/PROOF_MISSING_PHASE_RELAXATION_SCALE.md", claim.Anchor);
        Assert.Contains("simulations/missing_phase_relaxation_scale.py", claim.Anchor);
        Assert.Contains("simulations/results/missing_phase_relaxation_scale.json", claim.Anchor);
    }

    [Fact]
    public void Claim_HasNoPublicMethodThatCanMisstateTheCubicGermAsAnExactFiniteEpsilonGap()
    {
        var publicOperations = typeof(MissingPhaseRelaxationScaleClaim)
            .GetMethods(BindingFlags.Public | BindingFlags.Instance | BindingFlags.Static | BindingFlags.DeclaredOnly)
            .Where(method => !method.IsSpecialName)
            .ToList();

        Assert.Empty(publicOperations);
    }
}
