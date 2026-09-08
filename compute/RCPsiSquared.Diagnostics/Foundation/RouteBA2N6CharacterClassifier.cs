using System.Numerics;
using System.Globalization;
using System.Text.Json;
using System.Text.Json.Serialization;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.Numerics;

namespace RCPsiSquared.Diagnostics.Foundation;

public sealed record RouteBA2N6ReconciliationReport(
    int TotalLoci, int ConsumedLoci, int UnresolvedLoci,
    IReadOnlyDictionary<string, int> ByAlgebraicSource,
    IReadOnlyDictionary<string, int> ByCharacterSource,
    IReadOnlyDictionary<EpCharacter.EpKind, int> ByVerdict,
    [property: JsonConverter(typeof(RouteBA2N6MarginJsonConverter))] double MinimumIsolationMargin);

/// <summary>Output encoding of the absence of numerical character evidence.</summary>
internal sealed class RouteBA2N6MarginJsonConverter : JsonConverter<double>
{
    public override double Read(ref Utf8JsonReader reader, Type typeToConvert, JsonSerializerOptions options) =>
        reader.TokenType == JsonTokenType.String && reader.GetString() == "not-applicable"
            ? double.PositiveInfinity : reader.GetDouble();

    public override void Write(Utf8JsonWriter writer, double value, JsonSerializerOptions options)
    {
        if (double.IsPositiveInfinity(value)) writer.WriteStringValue("not-applicable");
        else writer.WriteNumberValue(value);
    }
}

/// <summary>Character at each exact-artifact N=6 seed, read on its full 45-dimensional parity sector.
/// No exact-rank fallback is consulted: unresolved character remains an ambiguity.</summary>
public sealed class RouteBA2N6CharacterClassifier
{
    private readonly RouteBA2N6Inventory inventory;
    private readonly Func<Matrix<Complex>, Complex, double, EpCharacter.Reading> readCharacter;

    public RouteBA2N6CharacterClassifier(RouteBA2N6Inventory inventory)
        : this(inventory, (block, lambda, radius) => EpCharacter.Characterize(block, lambda, radius)) { }

    internal RouteBA2N6CharacterClassifier(RouteBA2N6Inventory inventory,
        Func<Matrix<Complex>, Complex, double, EpCharacter.Reading> readCharacter)
    {
        this.inventory = inventory ?? throw new ArgumentNullException(nameof(inventory));
        this.readCharacter = readCharacter ?? throw new ArgumentNullException(nameof(readCharacter));
    }

    internal static bool IsHermitianEligible(RouteBA2N6Locus locus) =>
        locus.TBox.ImLo.Numerator.IsZero && locus.TBox.ImHi.Numerator.IsZero;

    public A2CharacterReading Classify(RouteBA2N6Locus locus) => ClassifyRadii(locus)[1];

    /// <summary>Consume every certified direct-t locus once, in ordinal ID order.
    /// Numerical character is stable evidence, never an exact-rank certificate.
    /// The margin is the minimum returned middle-contour margin among numerical
    /// character readings; +Infinity means no numerical character was needed.</summary>
    public RouteBA2N6ReconciliationReport ReconcileAll() => BuildAtlasManifest().Reconciliation;

    /// <summary>Classify every certified locus once and bind the accepted readings to the
    /// exact inventory geometry used by the N=6 discovery atlas.</summary>
    public RouteBA2N6AtlasManifest BuildAtlasManifest()
    {
        var classified = inventory.Loci.OrderBy(locus => locus.Id, StringComparer.Ordinal)
            .Select(locus => (Locus: locus, Contours: ClassifyRadii(locus)))
            .ToArray();
        RouteBA2N6ReconciliationReport reconciliation = ReconcileReadings(
            classified.Select(item => item.Contours[1]));
        return RouteBA2N6AtlasManifestBuilder.Build(inventory, classified, reconciliation);
    }

    internal RouteBA2N6ReconciliationReport ReconcileReadings(IEnumerable<A2CharacterReading> readings)
    {
        static ExactComplexBox Conjugate(ExactComplexBox box) =>
            new(box.ReLo, box.ReHi, box.ImHi.Negate(), box.ImLo.Negate());
        static void Require(bool condition, string reason)
        {
            if (!condition) throw new InvalidOperationException($"N6 reconciliation: {reason}");
        }

        Require(inventory.Loci.Select(l => l.Id).Distinct(StringComparer.Ordinal).Count() == inventory.Loci.Count,
            "duplicate inventory ID");
        Require(inventory.Loci.Count(l => l.Parity == A2Parity.Even) == inventory.A2Degrees.E
            && inventory.Loci.Count(l => l.Parity == A2Parity.Odd) == inventory.A2Degrees.O
            && inventory.Loci.Count == inventory.A2Degrees.E + inventory.A2Degrees.O,
            "inventory totals do not equal certified parity A2 degrees");
        var byId = inventory.Loci.ToDictionary(l => l.Id, StringComparer.Ordinal);
        foreach (var locus in inventory.Loci)
        {
            Require(byId.TryGetValue(locus.ConjugationPartnerId, out var conjugate)
                && conjugate.ConjugationPartnerId == locus.Id && conjugate.Parity == locus.Parity
                && conjugate.TBox == Conjugate(locus.TBox)
                && conjugate.LambdaClearedBox == Conjugate(locus.LambdaClearedBox),
                $"conjugation partner map fails at {locus.Id}");
            Require(byId.TryGetValue(locus.ParityPartnerId, out var partner)
                && partner.ParityPartnerId == locus.Id && partner.Parity != locus.Parity
                && partner.TBox == locus.TBox.Negate()
                && partner.LambdaClearedBox == locus.LambdaClearedBox,
                $"parity partner map fails at {locus.Id}");
        }

        var consumed = new Dictionary<string, A2CharacterReading>(StringComparer.Ordinal);
        foreach (var reading in readings)
        {
            Require(byId.ContainsKey(reading.LocusId), $"foreign reading {reading.LocusId}");
            Require(consumed.TryAdd(reading.LocusId, reading), $"duplicate consumption {reading.LocusId}");
            var locus = byId[reading.LocusId];
            Require(reading.Algebraic == locus.AlgebraicMultiplicity && reading.Algebraic == 2
                && (reading.Kind == EpCharacter.EpKind.Diabolic && reading.Geometric == 2
                    || reading.Kind == EpCharacter.EpKind.Defective && reading.Geometric == 1),
                $"unresolved character at {locus.Id}");
            switch (reading.Source)
            {
                case A2CharacterSource.HermitianAxis:
                    Require(IsHermitianEligible(locus) && reading.Kind == EpCharacter.EpKind.Diabolic
                        && reading.FullBlockHermiticityResidual is double residual && double.IsFinite(residual)
                        && residual < 1e-12, $"HermitianAxis evidence missing at {locus.Id}");
                    break;
                case A2CharacterSource.EpCharacter:
                    Require(!IsHermitianEligible(locus)
                        && reading.IsolationMargin is double margin && double.IsFinite(margin) && margin > 0
                        && reading.RelativeDeparture is double departure && double.IsFinite(departure)
                        && (reading.Kind == EpCharacter.EpKind.Diabolic && departure < 1e-6
                            || reading.Kind == EpCharacter.EpKind.Defective && departure > 5e-2),
                        $"EpCharacterStable evidence missing at {locus.Id}");
                    break;
                default:
                    // N=6 has no executed exact-rank route. A certificate label cannot
                    // turn a numerical reading into one, even if an artifact is supplied.
                    Require(false, $"ExactRankExecuted is unavailable at {locus.Id}");
                    break;
            }
        }
        Require(consumed.Count == byId.Count, "missing locus consumption");
        foreach (var locus in inventory.Loci)
        {
            var reading = consumed[locus.Id];
            Require(reading.Kind == consumed[locus.ConjugationPartnerId].Kind
                && reading.Kind == consumed[locus.ParityPartnerId].Kind,
                $"partner character disagreement at {locus.Id}");
        }
        var accepted = consumed.Values.ToArray();
        var sources = new Dictionary<string, int>(StringComparer.Ordinal)
        {
            ["HermitianAxis"] = accepted.Count(r => r.Source == A2CharacterSource.HermitianAxis),
            ["EpCharacterStable"] = accepted.Count(r => r.Source == A2CharacterSource.EpCharacter),
            ["ExactRankExecuted"] = accepted.Count(r => r.Source == A2CharacterSource.ExactRank)
        };
        return new(byId.Count, accepted.Length, byId.Count - accepted.Length,
            new Dictionary<string, int>(StringComparer.Ordinal)
            { ["ExactAlgebraic"] = accepted.Count(r => r.Algebraic == byId[r.LocusId].AlgebraicMultiplicity) },
            sources, accepted.GroupBy(r => r.Kind).ToDictionary(g => g.Key, g => g.Count()),
            accepted.Where(r => r.Source == A2CharacterSource.EpCharacter)
                .Select(r => r.IsolationMargin!.Value).DefaultIfEmpty(double.PositiveInfinity).Min());
    }

    public IReadOnlyList<A2CharacterReading> ClassifyRadii(RouteBA2N6Locus locus)
    {
        if (!inventory.Loci.Contains(locus))
            throw new A2CharacterUncertifiedException(locus.Id, "Locus does not belong to this inventory", []);
        try
        {
            // The loader already checked q=-it and physical lambda=Lambda/2 exactly.
            // Convert the authoritative decimal seeds once; box midpoints are not seeds.
            Complex q = ParseSeed(locus.QPhysicalCSharpSeed);
            Complex lambda = ParseSeed(locus.LambdaPhysicalSeed);
            var block = RouteBA2CharacterClassifier.FullParityBlock(inventory.N, locus.Parity, q);
            double? hermiticity = null;
            if (IsHermitianEligible(locus))
            {
                hermiticity = RouteBA2CharacterClassifier.HermiticityResidual(block);
                if (q.Real != 0 || lambda.Imaginary != 0
                    || !double.IsFinite(hermiticity.Value) || hermiticity.Value >= 1e-12)
                    throw new A2CharacterUncertifiedException(locus.Id,
                        $"Real-t full-sector Hermiticity failed: residual={hermiticity:G17}", []);
            }
            double[] distances = RouteBA2CharacterKernel.PairDistances(locus.Id,
                block.Evd().EigenValues, lambda);
            if (!hermiticity.HasValue)
                return RouteBA2CharacterKernel.Read(locus.Id, distances,
                    radius => readCharacter(block, lambda, radius));

            // Exact artifact multiplicity plus executed Hermiticity determines semisimplicity.
            // All contours still exclude a third eigenvalue, including every AT direction.
            return new[] { 0.25, 0.5, 0.75 }.Select(fraction =>
            {
                double radius = RouteBA2CharacterKernel.Radius(locus.Id, distances, fraction);
                return new A2CharacterReading(locus.Id, EpCharacter.EpKind.Diabolic,
                    locus.AlgebraicMultiplicity, 2, null, radius, distances[2] - radius,
                    A2CharacterSource.HermitianAxis) { FullBlockHermiticityResidual = hermiticity };
            }).ToArray();
        }
        catch (A2CharacterUncertifiedException) { throw; }
        catch (Exception error) when (error is ArithmeticException or ArgumentException or InvalidOperationException
            or FormatException)
        {
            throw new A2CharacterUncertifiedException(locus.Id,
                $"Numerical characterization failed: {error.Message}", [], error);
        }
    }

    private static Complex ParseSeed(N6SeedText seed)
    {
        double real = double.Parse(seed.Real, NumberStyles.Float, CultureInfo.InvariantCulture);
        double imag = double.Parse(seed.Imag, NumberStyles.Float, CultureInfo.InvariantCulture);
        if (!double.IsFinite(real) || !double.IsFinite(imag))
            throw new ArithmeticException("Seed cannot be represented as a finite complex double.");
        return new(real, imag);
    }
}
