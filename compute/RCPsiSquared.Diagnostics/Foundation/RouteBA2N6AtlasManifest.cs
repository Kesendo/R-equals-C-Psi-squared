using System.Globalization;

namespace RCPsiSquared.Diagnostics.Foundation;

public sealed record RouteBA2N6AtlasContour(
    double Radius, double IsolationMargin, double? RelativeDeparture);

public sealed record RouteBA2N6AtlasLocus(
    string Id,
    string Parity,
    double TReal,
    double TImag,
    string TRealExactMidpoint,
    string TImagExactMidpoint,
    string ConjugationPartnerId,
    string ParityPartnerId,
    int AlgebraicMultiplicity,
    int GeometricMultiplicity,
    string Verdict,
    string CharacterSource,
    double? FullBlockHermiticityResidual,
    IReadOnlyList<RouteBA2N6AtlasContour> Contours);

public sealed record RouteBA2N6AtlasOrbit(string Id, IReadOnlyList<string> MemberIds);

public sealed record RouteBA2N6AtlasManifest(
    int SchemaVersion,
    int N,
    string SourceArtifact,
    string SourceSha256,
    string CoordinateMeaning,
    RouteBA2N6ReconciliationReport Reconciliation,
    IReadOnlyList<RouteBA2N6AtlasLocus> Loci,
    IReadOnlyList<RouteBA2N6AtlasOrbit> Orbits);

internal static class RouteBA2N6AtlasManifestBuilder
{
    internal static RouteBA2N6AtlasManifest Build(
        RouteBA2N6Inventory inventory,
        IReadOnlyList<(RouteBA2N6Locus Locus, IReadOnlyList<A2CharacterReading> Contours)> classified,
        RouteBA2N6ReconciliationReport reconciliation)
    {
        if (string.IsNullOrWhiteSpace(inventory.SourcePayloadSha256))
            throw new InvalidOperationException("N6 atlas: source payload digest is unavailable.");

        var loci = classified.Select(item =>
        {
            A2CharacterReading middle = item.Contours[1];
            var t = item.Locus.TBox.Midpoint;
            return new RouteBA2N6AtlasLocus(
                item.Locus.Id,
                item.Locus.Parity == A2Parity.Even ? "E" : "O",
                t.Real,
                t.Imaginary,
                MidpointText(item.Locus.TBox.ReLo, item.Locus.TBox.ReHi),
                MidpointText(item.Locus.TBox.ImLo, item.Locus.TBox.ImHi),
                item.Locus.ConjugationPartnerId,
                item.Locus.ParityPartnerId,
                middle.Algebraic,
                middle.Geometric,
                middle.Kind.ToString(),
                middle.Source switch
                {
                    A2CharacterSource.HermitianAxis => "HermitianAxis",
                    A2CharacterSource.EpCharacter => "EpCharacterStable",
                    A2CharacterSource.ExactRank => "ExactRankExecuted",
                    _ => throw new InvalidOperationException(
                        $"N6 atlas: unknown character source at {item.Locus.Id}.")
                },
                middle.FullBlockHermiticityResidual,
                item.Contours.Select(reading => new RouteBA2N6AtlasContour(
                    reading.Radius ?? throw new InvalidOperationException(
                        $"N6 atlas: missing contour radius at {reading.LocusId}."),
                    reading.IsolationMargin ?? throw new InvalidOperationException(
                        $"N6 atlas: missing isolation margin at {reading.LocusId}."),
                    reading.RelativeDeparture)).ToArray());
        }).ToArray();

        var byId = inventory.Loci.ToDictionary(locus => locus.Id, StringComparer.Ordinal);
        var unseen = inventory.Loci.Select(locus => locus.Id).ToHashSet(StringComparer.Ordinal);
        var orbits = new List<RouteBA2N6AtlasOrbit>();
        while (unseen.Count > 0)
        {
            string seed = unseen.Min()!;
            var members = new HashSet<string>(StringComparer.Ordinal) { seed };
            var frontier = new Queue<string>();
            frontier.Enqueue(seed);
            while (frontier.Count > 0)
            {
                RouteBA2N6Locus locus = byId[frontier.Dequeue()];
                foreach (string partner in new[] { locus.ConjugationPartnerId, locus.ParityPartnerId })
                    if (members.Add(partner)) frontier.Enqueue(partner);
            }
            string[] ordered = members.Order(StringComparer.Ordinal).ToArray();
            foreach (string id in ordered) unseen.Remove(id);
            orbits.Add(new($"N6-A2-ORBIT-{orbits.Count:D3}", ordered));
        }

        return new RouteBA2N6AtlasManifest(
            1,
            inventory.N,
            "simulations/results/route_b_a2_n6.json",
            inventory.SourcePayloadSha256,
            "exact t-box midpoint converted to display double; not an algebraic root",
            reconciliation,
            loci,
            orbits);
    }

    private static string MidpointText(ExactRational lower, ExactRational upper)
    {
        ExactRational midpoint = (lower + upper).Half();
        return midpoint.Denominator.IsOne
            ? midpoint.Numerator.ToString(CultureInfo.InvariantCulture)
            : string.Create(CultureInfo.InvariantCulture,
                $"{midpoint.Numerator}/{midpoint.Denominator}");
    }
}
