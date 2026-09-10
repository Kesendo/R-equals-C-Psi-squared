using System.Text.Json.Nodes;
using RCPsiSquared.Diagnostics.Foundation;

namespace RCPsiSquared.Diagnostics.Tests.Foundation;

/// <summary>From-below pins for <see cref="RouteBN5RealQSemisimpleWitness"/> (F164, proof
/// docs/proofs/PROOF_N5_REAL_Q_DIABOLIC.md). Every number here is read off the witness, and the two
/// structural residuals are compared to 0.0 exactly because there is an exact route to them.</summary>
public class RouteBN5RealQSemisimpleWitnessTests
{
    private static readonly RouteBN5RealQSemisimpleWitness W = new();

    [Fact]
    public void SectorDimensions_AreTheShapeTheProofConsumes()
        => Assert.Equal((50, 26, 24), W.SectorDimensions);

    [Fact]
    public void TheCouplingBookIsPinnedBySpectrum_NotAssumed()
        => Assert.DoesNotContain("NONE", W.PinnedBook, StringComparison.Ordinal);

    [Fact]
    public void TheCouplingBookIsUnambiguous()
        => Assert.DoesNotContain("AMBIGUOUS", W.PinnedBook, StringComparison.Ordinal);

    [Fact]
    public void N6ExactRationalTBoxesExcludeTheRealQAxis()
    {
        var reading = W.N6RealQExclusion;
        Assert.Equal(266, reading.LocusCount);
        Assert.Equal(0, reading.RealAxisTouchCount);
        Assert.True(reading.AllExcluded);
    }

    [Fact]
    public void SummarySeparatesExactFromNumericalEvidenceAndNamesParameterVersusSpectralResult()
    {
        Assert.Contains("q is the settable parameter", W.Summary);
        Assert.Contains("lambda is the resulting spectral eigenvalue", W.Summary);
        Assert.DoesNotContain("observed eigenvalue", W.Summary);
        Assert.Contains("F164 certifies all 26 R-odd loci", W.Summary);
        Assert.Contains("38 distinct loci are exact-certified", W.Summary);
        Assert.Contains("20 nonreal-q R-even loci are numerical-only", W.Summary);
        Assert.Contains("classifier-local provenance", W.Summary);
    }

    [Fact]
    public void N6ExactRationalAxisGateRejectsABoxThatTouchesZero()
    {
        var json = JsonNode.Parse(File.ReadAllText(N6AtlasPath()))!;
        var real = json["loci"]![0]!["tBox"]!["real"]!;
        real["lower"] = RationalNode(0, 1);
        real["upper"] = RationalNode(0, 1);

        Assert.Throws<InvalidDataException>(() =>
            RouteBN5RealQSemisimpleWitness.ReadN6RealQExclusion(json.ToJsonString()));
        var reading = new N6RealQExclusionReading(266, 1);
        Assert.False(reading.AllExcluded);
        Assert.Contains("265 of 266", reading.Summary);
        Assert.DoesNotContain("266 of 266", reading.Summary);
    }

    [Fact]
    public void N6ReaderRejectsAMissingLocus()
    {
        var json = JsonNode.Parse(File.ReadAllText(N6AtlasPath()))!;
        json["loci"]!.AsArray().RemoveAt(0);
        Assert.Throws<InvalidDataException>(() =>
            RouteBN5RealQSemisimpleWitness.ReadN6RealQExclusion(json.ToJsonString()));
    }

    [Fact]
    public void N6ReaderRejectsAFabricatedOneLocusDocument()
    {
        var source = JsonNode.Parse(File.ReadAllText(N6AtlasPath()))!;
        var one = new JsonObject
        {
            ["schemaVersion"] = 3,
            ["n"] = 6,
            ["a2Degrees"] = new JsonObject { ["E"] = 133, ["O"] = 133 },
            ["loci"] = new JsonArray(source["loci"]![0]!.DeepClone()),
        };
        Assert.False(new N6RealQExclusionReading(1, 0).AllExcluded);
        Assert.Throws<InvalidDataException>(() =>
            RouteBN5RealQSemisimpleWitness.ReadN6RealQExclusion(one.ToJsonString()));
    }

    [Fact]
    public void N6ReaderRejectsDuplicateIdsAndWrongSchema()
    {
        var duplicate = JsonNode.Parse(File.ReadAllText(N6AtlasPath()))!;
        duplicate["loci"]![1]!["id"] = duplicate["loci"]![0]!["id"]!.GetValue<string>();
        Assert.Throws<InvalidDataException>(() =>
            RouteBN5RealQSemisimpleWitness.ReadN6RealQExclusion(duplicate.ToJsonString()));

        var wrongSchema = JsonNode.Parse(File.ReadAllText(N6AtlasPath()))!;
        wrongSchema["schemaVersion"] = 2;
        Assert.Throws<InvalidDataException>(() =>
            RouteBN5RealQSemisimpleWitness.ReadN6RealQExclusion(wrongSchema.ToJsonString()));

        var wrongN = JsonNode.Parse(File.ReadAllText(N6AtlasPath()))!;
        wrongN["n"] = 5;
        Assert.Throws<InvalidDataException>(() =>
            RouteBN5RealQSemisimpleWitness.ReadN6RealQExclusion(wrongN.ToJsonString()));
    }

    [Fact]
    public void N6ReaderRejectsAParameterConventionThatChangesTheRealQAxis()
    {
        var changed = JsonNode.Parse(File.ReadAllText(N6AtlasPath()))!;
        changed["conventions"]!["parameter"] = "t=qCSharp";
        Assert.Throws<InvalidDataException>(() =>
            RouteBN5RealQSemisimpleWitness.ReadN6RealQExclusion(changed.ToJsonString()));
    }

    [Fact]
    public void N6ReaderRejectsBoxesDetachedFromTheSemanticallyCertifiedCarrier()
    {
        var changed = JsonNode.Parse(File.ReadAllText(N6AtlasPath()))!;
        foreach (JsonNode? locus in changed["loci"]!.AsArray())
        {
            locus!["tBox"]!["real"] = new JsonObject
            {
                ["lower"] = RationalNode(1, 1),
                ["upper"] = RationalNode(2, 1),
            };
        }
        Assert.Throws<InvalidDataException>(() =>
            RouteBN5RealQSemisimpleWitness.ReadN6RealQExclusion(changed.ToJsonString()));
    }

    [Fact]
    public void N6ReaderRejectsWrongParityCountsAndUncertifiedMultiplicity()
    {
        var parity = JsonNode.Parse(File.ReadAllText(N6AtlasPath()))!;
        parity["loci"]![0]!["parity"] = "O";
        Assert.Throws<InvalidDataException>(() =>
            RouteBN5RealQSemisimpleWitness.ReadN6RealQExclusion(parity.ToJsonString()));

        var multiplicity = JsonNode.Parse(File.ReadAllText(N6AtlasPath()))!;
        multiplicity["loci"]![0]!["algebraicMultiplicity"] = 1;
        Assert.Throws<InvalidDataException>(() =>
            RouteBN5RealQSemisimpleWitness.ReadN6RealQExclusion(multiplicity.ToJsonString()));
    }

    /// <summary>Step 2: at real t the R-odd block is real symmetric. Exact route, so exact comparison.</summary>
    [Fact]
    public void OnTheHermitianAxis_TheOddBlockIsExactlyRealSymmetric()
        => Assert.True(W.HermitianAxisResidual == 0.0,
            $"real-symmetry residual {W.HermitianAxisResidual:E3} is not exactly zero");

    /// <summary>The control the axis check needs: OFF the axis the same residual must be large, or the
    /// check above would pass for a reason that has nothing to do with t being real.</summary>
    [Fact]
    public void OffTheHermitianAxis_TheSameResidualIsLarge()
    {
        var odd = RouteBN5RealQSemisimpleWitness.OddSectorColumns();
        var full = MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex>.Build.DenseOfArray(
            RouteBN5RealQSemisimpleWitness.BuildFull(
                new System.Numerics.Complex(1.0, 0.0), new[] { 1.0, 1.0, 1.0, 1.0 }));
        var block = odd.ConjugateTranspose() * full * odd;
        double worst = 0;
        for (int r = 0; r < block.RowCount; r++)
            for (int c = 0; c < block.ColumnCount; c++)
                worst = Math.Max(worst, Math.Abs(block[r, c].Imaginary));
        Assert.True(worst > 1.0, $"a real coupling should leave a large imaginary part, got {worst:E3}");
    }

    [Fact]
    public void AtThePhysicalPoint_TheReadingIsNullityTwo()
    {
        var (s1, s2, s3) = W.SingularTriple;
        Assert.True(s1 < 1e-10, $"s1 = {s1:E3}");
        Assert.True(s2 < 1e-10, $"s2 = {s2:E3}");
        Assert.True(s3 > 1e-3, $"s3 = {s3:E3} should be clear of the pair");
    }

    [Theory]
    [InlineData(1.0, 1.0, 1.0, 1.0)]
    [InlineData(1.0, 0.0, 0.0, 1.0)]
    [InlineData(1.0, -1.0, -1.0, 1.0)]
    public void PalindromicDetuneDirections_ReachTheOddSector(double a, double b, double c, double d)
    {
        var pattern = new[] { a, b, c, d };
        Assert.True(RouteBN5RealQSemisimpleWitness.IsPalindromic(pattern));
        Assert.False(RouteBN5RealQSemisimpleWitness.IsAntiPalindromic(pattern));
        Assert.True(RouteBN5RealQSemisimpleWitness.DirectionReach(pattern) > 0.5,
            $"reach {RouteBN5RealQSemisimpleWitness.DirectionReach(pattern):E3}");
    }

    /// <summary>The half that refuses the unscoped flatness claim: an anti-palindromic direction cannot
    /// move the R-odd sector at all, so no linear response law applies along it.</summary>
    [Theory]
    [InlineData(1.0, -1.0, 1.0, -1.0)]
    [InlineData(1.0, 0.0, 0.0, -1.0)]
    [InlineData(0.0, 1.0, -1.0, 0.0)]
    public void AntiPalindromicDetuneDirections_CannotReachTheOddSectorAtAll(double a, double b, double c, double d)
    {
        var pattern = new[] { a, b, c, d };
        Assert.True(RouteBN5RealQSemisimpleWitness.IsAntiPalindromic(pattern));
        // Exact route (the +-1 orbit basis), so exact comparison, not a tolerance.
        Assert.True(RouteBN5RealQSemisimpleWitness.DirectionReach(pattern) == 0.0,
            $"reach {RouteBN5RealQSemisimpleWitness.DirectionReach(pattern):E3} is not exactly zero");
    }

    /// <summary>The theorem is a property of the DELTA, not of the base: on a non-palindromic base an
    /// anti-palindromic delta still projects to exactly zero, even though the reflection no longer
    /// commutes there. What the base's symmetry buys is the sector's invariance, gated separately.</summary>
    [Fact]
    public void OnANonPalindromicBase_TheAntiPalindromicDeltaStillProjectsToZero()
    {
        var q = new System.Numerics.Complex(1.1292509708747671, 0.0);
        var bad = new[] { 1.0, 1.3, 0.7, 1.9 };
        var build = MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex>.Build;
        var lBad = build.DenseOfArray(RouteBN5RealQSemisimpleWitness.BuildFull(q, bad));
        var moved = build.DenseOfArray(RouteBN5RealQSemisimpleWitness.BuildFull(
            q, bad.Zip(new[] { 1.0, -1.0, 1.0, -1.0 }, (b, w) => b + w).ToArray()));
        var u = RouteBN5RealQSemisimpleWitness.OddSectorColumnsExact();
        double reach = ((u.ConjugateTranspose() * (moved - lBad) * u)).Enumerate().Max(z => z.Magnitude);
        Assert.True(reach == 0.0, $"reach {reach:E3} on an asymmetric base is not exactly zero");
    }

    /// <summary>A vanishing perturbation must not come back as 0.0, because 0.0 is exactly the value that
    /// means "this direction misses the sector". A sentinel that cannot be told from the theorem's answer
    /// is not a guard.</summary>
    [Fact]
    public void AVanishingDetuneIsRefused_NotReportedAsAMissedSector()
        => Assert.Throws<ArgumentException>(
            () => RouteBN5RealQSemisimpleWitness.DirectionReach(new[] { 0.0, 0.0, 0.0, 0.0 }));

    /// <summary>The two named classes are not a dichotomy: a direction that is neither palindromic nor
    /// anti-palindromic is outside the theorem and does reach the sector.</summary>
    [Fact]
    public void ADirectionThatIsNeitherClassStillReachesTheSector()
    {
        var pattern = new[] { 1.0, 2.0, 0.0, 0.0 };
        Assert.False(RouteBN5RealQSemisimpleWitness.IsPalindromic(pattern));
        Assert.False(RouteBN5RealQSemisimpleWitness.IsAntiPalindromic(pattern));
        Assert.True(RouteBN5RealQSemisimpleWitness.DirectionReach(pattern) > 0.1);
    }

    /// <summary>The positive control the restriction needs: on the uniform profile the reflection commutes
    /// with the block exactly, which is what makes the R-odd restriction a sub-spectrum rather than a
    /// compression. Exact route, so exact comparison.</summary>
    [Fact]
    public void OnTheUniformProfile_TheReflectionCommutesExactly()
        => Assert.True(RouteBN5RealQSemisimpleWitness.ReflectionCommutator(null) == 0.0,
            $"commutator {RouteBN5RealQSemisimpleWitness.ReflectionCommutator(null):E3} is not exactly zero");

    [Fact]
    public void OnANonPalindromicProfile_TheReflectionDoesNotCommute()
        => Assert.True(RouteBN5RealQSemisimpleWitness.ReflectionCommutator(new[] { 1.0, 1.3, 0.7, 1.9 }) > 1.0);

    [Fact]
    public void TheWitnessRendersItsWholeSubtreeWithoutThrowing()
    {
        Assert.NotEmpty(W.Summary);
        Assert.All(W.Children.ToList(), child => Assert.NotEmpty(child.Summary));
    }

    private static JsonObject RationalNode(long numerator, long denominator) => new()
    {
        ["numerator"] = numerator.ToString(System.Globalization.CultureInfo.InvariantCulture),
        ["denominator"] = denominator.ToString(System.Globalization.CultureInfo.InvariantCulture),
    };

    private static string N6AtlasPath()
    {
        DirectoryInfo? cursor = new(AppContext.BaseDirectory);
        while (cursor is not null)
        {
            string candidate = Path.Combine(cursor.FullName, "simulations", "results", "route_b_a2_n6.json");
            if (File.Exists(candidate)) return candidate;
            cursor = cursor.Parent;
        }
        throw new FileNotFoundException("Could not locate simulations/results/route_b_a2_n6.json");
    }
}
