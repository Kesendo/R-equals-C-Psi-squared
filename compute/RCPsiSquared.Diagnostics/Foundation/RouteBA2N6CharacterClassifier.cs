using System.Numerics;
using System.Globalization;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.Numerics;

namespace RCPsiSquared.Diagnostics.Foundation;

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
