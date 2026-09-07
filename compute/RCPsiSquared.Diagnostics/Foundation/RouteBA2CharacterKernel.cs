using System.Globalization;
using System.Numerics;
using RCPsiSquared.Core.Numerics;

namespace RCPsiSquared.Diagnostics.Foundation;

/// <summary>The fixed Route B pair-isolation and three-contour character contract.</summary>
internal static class RouteBA2CharacterKernel
{
    internal static double[] PairDistances(string locusId, IEnumerable<Complex> roots, Complex lambda)
    {
        double[] distances = roots.Select(value => (value - lambda).Magnitude).Order().ToArray();
        if (distances.Length < 3 || distances.Any(d => !double.IsFinite(d))
            || !(distances[2] > 100 * distances[1]))
            throw new A2CharacterUncertifiedException(locusId,
                $"Pair is not isolated; nearest distances=[{string.Join(", ", distances.Take(3).Select(d => d.ToString("G17", CultureInfo.InvariantCulture)))}]", []);
        return distances;
    }

    internal static double Radius(string locusId, double[] distances, double fraction)
    {
        double radius = distances[1] + fraction * (distances[2] - distances[1]);
        if (!double.IsFinite(radius) || !(radius > distances[1]) || !(distances[2] > radius)
            || distances.Count(d => d < radius) != 2)
            throw new A2CharacterUncertifiedException(locusId,
                $"Radius {radius:G17} does not contain exactly two roots with positive margins", []);
        return radius;
    }

    internal static IReadOnlyList<A2CharacterReading> Read(string locusId, IEnumerable<Complex> roots,
        Complex lambda, Func<double, EpCharacter.Reading> readCharacter) =>
        Read(locusId, PairDistances(locusId, roots, lambda), readCharacter);

    internal static IReadOnlyList<A2CharacterReading> Read(string locusId, double[] distances,
        Func<double, EpCharacter.Reading> readCharacter)
    {
        var readings = new List<A2CharacterReading>();
        try
        {
            foreach (double fraction in new[] { 0.25, 0.5, 0.75 })
            {
                double radius = Radius(locusId, distances, fraction);
                var character = readCharacter(radius);
                // A non-finite norm must not turn a finite departure into an accepted zero.
                double relative = double.IsFinite(character.CompressionNorm) && double.IsFinite(character.Departure)
                    ? character.Departure / Math.Max(1, character.CompressionNorm) : double.NaN;
                readings.Add(new(locusId, character.Kind, character.Algebraic, character.Geometric,
                    relative, radius, distances[2] - radius, A2CharacterSource.EpCharacter));
            }
            bool consistent = readings.All(r => r.Algebraic == 2 && r.Kind == readings[0].Kind
                && r.Geometric == readings[0].Geometric && r.RelativeDeparture.HasValue
                && double.IsFinite(r.RelativeDeparture.Value)
                && (r.Kind == EpCharacter.EpKind.Diabolic && r.Geometric == 2 && r.RelativeDeparture < 1e-6
                    || r.Kind == EpCharacter.EpKind.Defective && r.Geometric == 1 && r.RelativeDeparture > 5e-2));
            if (!consistent)
                throw new A2CharacterUncertifiedException(locusId,
                    "Full-sector character fails the fixed three-radius contract", readings.AsReadOnly());
            return readings.AsReadOnly();
        }
        catch (A2CharacterUncertifiedException) { throw; }
        catch (Exception error) when (error is ArithmeticException or ArgumentException or InvalidOperationException)
        {
            throw new A2CharacterUncertifiedException(locusId,
                $"Numerical characterization failed: {error.Message}", readings.AsReadOnly(), error);
        }
    }
}
