using System.Numerics;
using System.Text.Json.Nodes;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Diagnostics.Foundation;

namespace RCPsiSquared.Diagnostics.Tests.Foundation;

public class RouteBN6A2UnfoldingWitnessTests
{
    private static RouteBN6A2UnfoldingWitness NewWitness(string? path = null) => new(path);
    private static readonly Lazy<RouteBN6A2UnfoldingWitness> Shared = new(() => NewWitness());

    [Fact]
    public void ProvenanceFingerprintUsesCanonicalUtf8LfForLfCrLfAndLoneCr()
    {
        const string expected = "e49c81e2d2f84e259d40e2fb8192f3bcd198b355184845d76d8f58807d0d78ee";
        Assert.Equal(expected, RouteBN6UnfoldingCertificate.CanonicalUtf8LfSha256("alpha\nbeta\n"));
        Assert.Equal(expected, RouteBN6UnfoldingCertificate.CanonicalUtf8LfSha256("alpha\r\nbeta\r\n"));
        Assert.Equal(expected, RouteBN6UnfoldingCertificate.CanonicalUtf8LfSha256("alpha\rbeta\r"));
    }

    [Fact]
    public void ConstructionIsLazyAndDisplayMetadataDoesNotReconstruct()
    {
        var witness = NewWitness();
        Assert.Equal(0, (int)witness.ReconstructionCount);
        var root = (IInspectable)witness;
        _ = root.DisplayName;
        _ = root.Payload;
        Assert.Equal(0, (int)witness.ReconstructionCount);
    }

    [Fact]
    public void RecomputedChildrenAreMarkedLiveAndConsumedProofIsMarkedStored()
    {
        var children = Shared.Value.Children.ToDictionary(c => c.DisplayName);
        Assert.Equal(NodeProvenance.Stored, children["source and proof boundary"].Provenance);
        foreach (string name in new[] { "exact pencil and inventory", "live irreducibility modulo 367",
            "live modular plane", "live local response invariants", "controls through the effective-matrix reading" })
            Assert.Equal(NodeProvenance.Live, children[name].Provenance);
    }

    [Fact]
    public void TheLiveRootDefinesTheThreeEndBondEpsilonConventions()
    {
        string prose = string.Join("\n", Shared.Value.Children.Select(c => c.Summary));
        Assert.Contains("(1+epsilon,1)", prose);
        Assert.Contains("(1+epsilon/2,1+epsilon/2)", prose);
        Assert.Contains("(1+epsilon/2,1-epsilon/2)", prose);
    }

    [Fact]
    public void ExactPencilsJoinTheCompleteInventoryAndStateTheConsumedCertificate()
    {
        var read = Shared.Value.Reading;
        Assert.Equal(90, (int)read.FullDimension);
        Assert.Equal(45, (int)read.EvenDimension);
        Assert.Equal(45, (int)read.OddDimension);
        Assert.Equal(32, (int)read.ResidualDegree);
        Assert.Equal(13, (int)read.AtDegree);
        Assert.Equal(133, (int)read.EvenLoci);
        Assert.Equal(133, (int)read.OddLoci);
        Assert.True((bool)read.ParityTransportExact);
        Assert.True((bool)read.HermitianRealTExact);
        Assert.True((bool)read.SourcePencilMatches);
        Assert.Contains("consumed", (string)read.CertificateSource);
        Assert.Contains("bounded-CRT", (string)read.CertificateSource);
    }

    [Fact]
    public void RecomputesTheIrreducibleDegree133ReductionAtTheNonSplitPrime367()
    {
        var read = Shared.Value.Reading;
        Assert.Equal(new[] { 133 }, read.CycleType.ToArray());
        var reducible = new BigInteger[134];
        reducible[0] = reducible[133] = BigInteger.One;
        int[] cycle = Shared.Value.CycleTypeModulo367(reducible)!;
        Assert.NotEqual(new[] { 133 }, cycle);
        Assert.Equal(133, cycle.Sum());
        Assert.Contains(1, cycle);
    }

    [Fact]
    public void ExposedCycleTypeCannotBeMutatedBeforeALaterInspect()
    {
        var witness = NewWitness();
        IReadOnlyList<int> exposed = witness.Reading.CycleType;
        var list = Assert.IsAssignableFrom<IList<int>>(exposed);
        Assert.Throws<NotSupportedException>(() => list[0] = 7);
        Assert.Equal(133, witness.Reading.CycleType[0]);
        Assert.Contains("[133]", witness.Children.Single(c =>
            c.DisplayName == "live irreducibility modulo 367").Summary);
    }

    [Fact]
    public void RecomputesTheRepeatedRootAndSemisimplePlaneModulo101()
    {
        var m = Shared.Value.Reading.Mod101;
        Assert.Equal(101, (int)m.Prime);
        Assert.Equal(31, (int)m.T);
        Assert.Equal(0, (int)m.A2Value);
        Assert.Equal(65, (int)m.A2Derivative);
        Assert.Equal(1, (int)m.ResidualGcdDegree);
        Assert.Equal(48, (int)m.LambdaCleared);
        Assert.Equal(24, (int)m.LambdaPhysical);
        Assert.Equal(2, (int)m.Nullity);
        Assert.Equal(2, (int)m.SquaredNullity);
        Assert.Equal(2, (int)m.ProjectorRank);
        Assert.Equal(39, (int)m.ComplementCharacteristicAtLambda);
        Assert.True((bool)m.ProjectorIdempotent);
        Assert.True((bool)m.ProjectorInEvenSector);
        Assert.True((bool)m.ReducedResolventIdentity);
    }

    [Fact]
    public void TheThreeProfilesComputeTheirOwnNonzeroUnfoldingDiscriminants()
    {
        var m = Shared.Value.Reading.Mod101;
        Assert.Equal(new[] { 19, 4, 0, 16 }, (int[])m.OneEnd.ToArray());
        Assert.Equal(new[] { 19, 4, 0, 16 }, (int[])m.EqualEnds.ToArray());
        Assert.Equal(new[] { 19, 52, 18, 23 }, (int[])m.OppositeEnds.ToArray());
        Assert.True((bool)m.OddFirstOrderZero);
        Assert.True((bool)m.OneEqualsEqualFirstOrder);
    }

    [Fact]
    public void ScalarAndCommutingResponsesFailToUnfoldAndEvenAdmixtureReachesThePlane()
    {
        var m = Shared.Value.Reading.Mod101;
        Assert.Equal(0, (int)m.ScalarControlOmega);
        Assert.Equal(0, (int)m.CommutingControlOmega);
        Assert.True((bool)m.EvenAdmixtureNonzero);
    }

    [Fact]
    public void AChangedConsumedA2CoefficientIsRejectedThroughTheSameLoader()
    {
        string source = Find("simulations/results/route_b_n6_exact_unfolding.json");
        var payload = JsonNode.Parse(File.ReadAllText(source))!;
        var coefficients = payload["layer"]!["a2_coefficients_lowest_first"]!.AsArray();
        coefficients[0] = (BigInteger.Parse(coefficients[0]!.GetValue<string>()) + 1).ToString();
        string temporary = Path.Combine(Path.GetTempPath(), $"f163-coefficient-{Guid.NewGuid():N}.json");
        try
        {
            File.WriteAllText(temporary, payload.ToJsonString());
            var witness = NewWitness(temporary);
            Assert.Throws<InvalidDataException>(() => { object ignored = witness.Reading; });
        }
        finally { File.Delete(temporary); }
    }

    [Fact]
    public void AChangedExactProducerFingerprintIsRejectedThroughTheSameLoader()
    {
        string source = Find("simulations/results/route_b_n6_exact_unfolding.json");
        var payload = JsonNode.Parse(File.ReadAllText(source))!;
        payload["provenance"]!["source_producer_sha256"] = new string('0', 64);
        string temporary = Path.Combine(Path.GetTempPath(), $"f163-source-producer-{Guid.NewGuid():N}.json");
        try
        {
            File.WriteAllText(temporary, payload.ToJsonString());
            var witness = NewWitness(temporary);
            Assert.Throws<InvalidDataException>(() => { object ignored = witness.Reading; });
        }
        finally { File.Delete(temporary); }
    }

    [Fact]
    public void AChangedTransitiveProducerDependencyFingerprintIsRejectedThroughTheSameLoader()
    {
        string source = Find("simulations/results/route_b_n6_exact_unfolding.json");
        var payload = JsonNode.Parse(File.ReadAllText(source))!;
        var dependencies = payload["provenance"]!["source_producer_dependencies"]!.AsObject();
        Assert.Contains("simulations/o2b_gcd_certificate.py", dependencies.Select(pair => pair.Key));
        dependencies["simulations/o2b_gcd_certificate.py"] = new string('0', 64);
        string temporary = Path.Combine(Path.GetTempPath(), $"f163-transitive-dependency-{Guid.NewGuid():N}.json");
        try
        {
            File.WriteAllText(temporary, payload.ToJsonString());
            var witness = NewWitness(temporary);
            Assert.Throws<InvalidDataException>(() => { object ignored = witness.Reading; });
        }
        finally { File.Delete(temporary); }
    }

    [Theory]
    [InlineData("script_sha256", "wrapper")]
    [InlineData("fixture_sha256", "fixture")]
    public void AChangedWrapperOrFixtureFingerprintIsRejectedThroughTheSameLoader(string field, string label)
    {
        string source = Find("simulations/results/route_b_n6_exact_unfolding.json");
        var payload = JsonNode.Parse(File.ReadAllText(source))!;
        payload["provenance"]![field] = new string('0', 64);
        string temporary = Path.Combine(Path.GetTempPath(), $"f163-{label}-{Guid.NewGuid():N}.json");
        try
        {
            File.WriteAllText(temporary, payload.ToJsonString());
            var witness = NewWitness(temporary);
            Assert.Throws<InvalidDataException>(() => { object ignored = witness.Reading; });
        }
        finally { File.Delete(temporary); }
    }

    [Fact]
    public void DuplicateRootLayerIsRejectedEvenWhenTheValidValueComesLast()
    {
        string json = File.ReadAllText(Find("simulations/results/route_b_n6_exact_unfolding.json"));
        int rootStart = json.IndexOf('{') + 1;
        string duplicated = json.Insert(rootStart, "\n  \"layer\": null,");
        AssertInvalidThroughPublicConstructor(duplicated);
    }

    [Fact]
    public void DuplicateNestedA2DegreeIsRejectedEvenWhenTheValidValueComesLast()
    {
        string json = File.ReadAllText(Find("simulations/results/route_b_n6_exact_unfolding.json"));
        const string valid = "\"a2_degree\": 133";
        Assert.Equal(1, json.Split(valid, StringSplitOptions.None).Length - 1);
        string duplicated = json.Replace(valid, "\"a2_degree\": 999,\n    \"a2_degree\": 133", StringComparison.Ordinal);
        AssertInvalidThroughPublicConstructor(duplicated);
    }

    [Fact]
    public void OverrideCertificatePathReportsTheActualLoadedPath()
    {
        string source = Find("simulations/results/route_b_n6_exact_unfolding.json");
        string temporary = Path.Combine(Path.GetTempPath(), $"f163-override-{Guid.NewGuid():N}.json");
        try
        {
            File.Copy(source, temporary);
            var witness = NewWitness(temporary);
            Assert.Contains(Path.GetFullPath(temporary), witness.Reading.CertificateSource, StringComparison.Ordinal);
        }
        finally { File.Delete(temporary); }
    }

    [Fact]
    public void EveryChildRendersAndOneInstanceReconstructsOnlyOnce()
    {
        var witness = Shared.Value;
        object first = witness.Reading;
        var root = (IInspectable)witness;
        Assert.Contains("F163", root.DisplayName);
        Assert.Contains("complex-q", root.Summary);
        var children = root.Children.ToArray();
        Assert.NotEmpty(children);
        Assert.All(children, child =>
        {
            Assert.False(string.IsNullOrWhiteSpace(child.DisplayName));
            Assert.False(string.IsNullOrWhiteSpace(child.Summary));
            _ = child.Children.ToArray();
            _ = child.Payload;
        });
        _ = root.Summary;
        _ = root.Children.ToArray();
        Assert.Same(first, (object)witness.Reading);
        Assert.Equal(1, (int)witness.ReconstructionCount);
    }

    private static string Find(string relative)
    {
        for (var directory = new DirectoryInfo(AppContext.BaseDirectory); directory != null; directory = directory.Parent)
        {
            string path = Path.Combine(directory.FullName, relative);
            if (File.Exists(path)) return path;
        }
        throw new FileNotFoundException(relative);
    }

    private static void AssertInvalidThroughPublicConstructor(string json)
    {
        string temporary = Path.Combine(Path.GetTempPath(), $"f163-duplicate-{Guid.NewGuid():N}.json");
        try
        {
            File.WriteAllText(temporary, json);
            var witness = new RouteBN6A2UnfoldingWitness(temporary);
            var error = Assert.Throws<InvalidDataException>(() => { object ignored = witness.Reading; });
            Assert.Contains("Duplicate JSON property", error.ToString(), StringComparison.Ordinal);
        }
        finally { File.Delete(temporary); }
    }
}
