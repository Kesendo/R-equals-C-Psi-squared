using System.Globalization;
using System.Text.Json;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Diagnostics.Foundation;

namespace RCPsiSquared.Diagnostics.Tests.Foundation;

/// <summary>Every number the claim states is read back out of the committed ball certificate,
/// so a producer rerun that moves a radius, an enclosure or a remainder turns this gate red
/// rather than leaving the claim's prose ahead of its evidence. The certificate's own premise
/// bindings are checked too: a premise that moves under a frozen certificate is the failure the
/// producer records those hashes to catch, and it went unnoticed once because nothing read them.</summary>
[Trait("Category", "ROUTE_B_A2_N6_REMAINDER")]
public class RouteBN6RemainderBoundClaimTests
{
    private const string Artifact = "simulations/results/route_b_n6_remainder_ball.json";

    private static string RepoRoot()
    {
        foreach (string start in new[] { AppContext.BaseDirectory, Directory.GetCurrentDirectory() })
        for (var directory = new DirectoryInfo(start); directory != null; directory = directory.Parent)
            if (File.Exists(Path.Combine(directory.FullName, Artifact)))
                return directory.FullName;
        throw new FileNotFoundException($"Cannot locate {Artifact}.");
    }

    private static JsonDocument Certificate() =>
        JsonDocument.Parse(File.ReadAllText(Path.Combine(RepoRoot(), Artifact)));

    /// <summary>The one canonical digest the F163 certificate consumer already owns. A gate
    /// that exists to catch two implementations disagreeing must not add a third.</summary>
    private static string DigestOf(string path) =>
        RouteBN6UnfoldingCertificate.CanonicalUtf8LfSha256(File.ReadAllText(path));

    private static double Number(string ball)
    {
        Assert.DoesNotContain("+/-", ball);
        return double.Parse(ball, CultureInfo.InvariantCulture);
    }

    private sealed record Attempt(double RadiusEta, int EtaOrder, double VariableRadius,
        double EpsilonDisplay, double HalfRadiusRemainder, IReadOnlyList<double> NewtonImages);

    /// <summary>A branch is certified on the LARGEST radius at which the contraction succeeded;
    /// a profile's common radius is then the smallest of its branches' certified radii, because
    /// holomorphy on a disk restricts to every smaller disk.</summary>
    private static IReadOnlyList<Attempt> CertifiedPerBranch(JsonElement profiles, string profile)
    {
        var certified = new List<Attempt>();
        foreach (var branch in profiles.GetProperty(profile).EnumerateArray())
        {
            Attempt? best = null;
            foreach (var attempt in branch.GetProperty("attempts").EnumerateArray())
            {
                if (!attempt.GetProperty("success").GetBoolean()) continue;
                var read = new Attempt(
                    Number(attempt.GetProperty("radius_eta").GetString()!),
                    attempt.GetProperty("eta_order").GetInt32(),
                    Number(attempt.GetProperty("radius_variables").GetString()!),
                    Number(attempt.GetProperty("epsilon_radius_display").GetString()!),
                    Number(attempt.GetProperty("half_radius_relative_remainder_upper").GetString()!
                        .TrimStart('[').Split(" +/- ")[0]),
                    attempt.GetProperty("newton_image_bounds").EnumerateArray()
                        .Select(b => Number(b.GetString()!)).ToList());
                if (best is null || read.RadiusEta > best.RadiusEta) best = read;
            }
            Assert.True(best is not null, $"{profile}: a branch carries no successful attempt");
            certified.Add(best!);
        }
        Assert.Equal(2, certified.Count);
        return certified;
    }

    private static double CommonRadius(IReadOnlyList<Attempt> branches) => branches.Min(b => b.RadiusEta);

    [Fact]
    public void CertificateCoversTheCrossingTheClaimNames_AtTheClaimedPrecision()
    {
        using var certificate = Certificate();
        var root = certificate.RootElement;
        Assert.True(root.GetProperty("all_branches_certified").GetBoolean());
        Assert.Equal(RouteBN6RemainderBoundClaim.Crossing, root.GetProperty("seed_id").GetString());
        Assert.Equal(RouteBN6RemainderBoundClaim.PrecisionBits, root.GetProperty("precision").GetInt32());
    }

    /// <summary>The certificate binds its producers and premises by digest precisely so that a
    /// premise moving underneath it is visible. Reading them back is the only thing that makes
    /// that binding do any work.</summary>
    [Fact]
    public void EveryProducerAndPremiseDigest_StillMatchesTheFileOnDisk()
    {
        string root = RepoRoot();
        using var certificate = Certificate();
        var document = certificate.RootElement;
        Assert.Equal(document.GetProperty("script_sha256").GetString(),
            DigestOf(Path.Combine(root, "simulations/route_b_n6_remainder_ball.py")));
        Assert.Equal(document.GetProperty("base_script_sha256").GetString(),
            DigestOf(Path.Combine(root, "simulations/route_b_n6_ball_base.py")));
        int premises = 0;
        foreach (var premise in document.GetProperty("premises").EnumerateObject())
        {
            Assert.Equal(premise.Value.GetString(), DigestOf(Path.Combine(root, premise.Name)));
            premises++;
        }
        Assert.Equal(3, premises);
    }

    [Theory]
    [InlineData("one", RouteBN6RemainderBoundClaim.OneCommonRadiusExponent, 1)]
    [InlineData("even", RouteBN6RemainderBoundClaim.EvenCommonRadiusExponent, 1)]
    [InlineData("odd", RouteBN6RemainderBoundClaim.OddCommonRadiusExponent, 2)]
    public void ClaimedCommonRadiusAndEtaOrder_AreTheCertifiedOnes(string profile, int exponent, int etaOrder)
    {
        using var certificate = Certificate();
        var branches = CertifiedPerBranch(certificate.RootElement.GetProperty("profiles"), profile);
        Assert.All(branches, b => Assert.Equal(etaOrder, b.EtaOrder));
        Assert.Equal(Math.Pow(2.0, -exponent), CommonRadius(branches));
    }

    /// <summary>rho is the yardstick every Newton image is measured against, so the claim's rho
    /// must be the one the producer actually ran with, not a number chosen beside it.</summary>
    [Fact]
    public void ClaimedRho_IsTheVariableRadiusTheProducerRanWith()
    {
        using var certificate = Certificate();
        double rho = Math.Pow(2.0, -RouteBN6RemainderBoundClaim.ImageRadiusExponent);
        int seen = 0;
        foreach (var profile in new[] { "one", "even", "odd" })
        foreach (var branch in CertifiedPerBranch(certificate.RootElement.GetProperty("profiles"), profile))
        {
            Assert.Equal(rho, branch.VariableRadius);
            seen++;
        }
        Assert.Equal(6, seen);
    }

    /// <summary>eta = epsilon^2 on the odd profile, so the reported epsilon radius must be the
    /// exact square root of the certified eta radius. Both are dyadic, so this is exact.</summary>
    [Fact]
    public void OddProfileEpsilonRadius_IsTheSquareRootOfItsEtaRadius()
    {
        using var certificate = Certificate();
        var branches = CertifiedPerBranch(certificate.RootElement.GetProperty("profiles"), "odd");
        foreach (var branch in branches)
        {
            Assert.Equal(2, branch.EtaOrder);
            Assert.Equal(branch.RadiusEta, branch.EpsilonDisplay * branch.EpsilonDisplay);
            Assert.Equal(Math.Pow(2.0, -RouteBN6RemainderBoundClaim.OddEpsilonRadiusExponent), branch.EpsilonDisplay);
        }
        Assert.Equal(2 * RouteBN6RemainderBoundClaim.OddEpsilonRadiusExponent,
            RouteBN6RemainderBoundClaim.OddCommonRadiusExponent);
    }

    /// <summary>Every certified Newton image lands strictly inside rho, and the largest of them
    /// is pinned exactly rather than inside a window wide enough to hide its last digit.</summary>
    [Fact]
    public void EveryCertifiedNewtonImage_LiesStrictlyInsideRho()
    {
        using var certificate = Certificate();
        double rho = Math.Pow(2.0, -RouteBN6RemainderBoundClaim.ImageRadiusExponent);
        double worst = 0.0;
        int seen = 0;
        foreach (var profile in new[] { "one", "even", "odd" })
        foreach (var branch in CertifiedPerBranch(certificate.RootElement.GetProperty("profiles"), profile))
        foreach (double image in branch.NewtonImages)
        {
            Assert.True(image < rho, $"{profile}: Newton image {image} is not inside rho={rho}");
            worst = Math.Max(worst, image);
            seen++;
        }
        Assert.Equal(12, seen);
        Assert.Equal(RouteBN6RemainderBoundClaim.LargestCertifiedNewtonImage, worst);
    }

    /// <summary>The claim's three percentages are the worst branch at half the profile's COMMON
    /// radius. Which branch that is cannot be argued from s alone: the wider branch sits at a
    /// smaller s but carries a larger prefactor, so the gate reads the winner off the certificate
    /// and pins that it is the branch certified exactly on the common radius.</summary>
    [Theory]
    [InlineData("one", RouteBN6RemainderBoundClaim.OneWorstBranchHalfRadiusPercent)]
    [InlineData("even", RouteBN6RemainderBoundClaim.EvenWorstBranchHalfRadiusPercent)]
    [InlineData("odd", RouteBN6RemainderBoundClaim.OddWorstBranchHalfRadiusPercent)]
    public void ClaimedWorstBranchRemainder_IsTheTightestFiveDecimalUpperBound(string profile, double claimed)
    {
        using var certificate = Certificate();
        var branches = CertifiedPerBranch(certificate.RootElement.GetProperty("profiles"), profile);
        double common = CommonRadius(branches);
        double worst = 0.0;
        double onTheCommonRadius = 0.0;
        foreach (var branch in branches)
        {
            // half_radius_relative_remainder_upper is taken at half the branch's OWN radius,
            // where s = 1/2. At half the COMMON radius s = common/(2*own), and the stored value
            // rescales by f(s)/f(1/2) with f(s) = s^2/(1-s).
            double s = common / (2.0 * branch.RadiusEta);
            double atHalfCommon = branch.HalfRadiusRemainder * (s * s / (1.0 - s)) / (0.25 / 0.5);
            if (branch.RadiusEta == common)
            {
                // s is exactly 1/2 here, so the rescaling factor is exactly one.
                Assert.Equal(branch.HalfRadiusRemainder, atHalfCommon);
                onTheCommonRadius = Math.Max(onTheCommonRadius, atHalfCommon);
            }
            worst = Math.Max(worst, atHalfCommon);
        }
        // The claim reads the worst branch off the certificate rather than arguing it from s:
        // a branch certified on a wider disk sits at a smaller s but carries a larger prefactor.
        Assert.Equal(onTheCommonRadius, worst);
        double percent = worst * 100.0;
        Assert.True(percent <= claimed, $"{profile}: claimed {claimed}% is not an upper bound for {percent}%");
        Assert.True(percent > claimed - 1e-5, $"{profile}: claimed {claimed}% is looser than the tightest 5-decimal bound");
    }

    [Fact]
    public void ClaimFencesTheResultToOneCrossingAndCarriesF163AsItsPremise()
    {
        var registry = RCPsiSquared.Diagnostics.Knowledge.KnowledgeRegistryFactory.BuildDefault();
        var claim = Assert.IsType<RouteBN6RemainderBoundClaim>(
            registry.All().Single(c => c is RouteBN6RemainderBoundClaim));
        // The tier and the single parent edge are gated in Runtime.Tests, which owns registration;
        // this test owns only the fencing words.
        Assert.Contains(RouteBN6RemainderBoundClaim.Crossing, claim.Name);
        Assert.Contains("BRANCH-SPECIFIC", claim.Name);
        Assert.Contains("WORST-BRANCH", claim.Name);
        Assert.Contains("no uniform family radius", claim.Name);
        Assert.Contains("no physical-real-q operating point", claim.Name);
        Assert.Contains("sufficient", claim.Name);
        Assert.Contains("266-locus family", claim.Summary);
        Assert.NotNull(claim.Unfolding);
    }
}
