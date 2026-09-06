using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.F89PathK;
using RCPsiSquared.Core.Numerics;

namespace RCPsiSquared.Diagnostics.Foundation;

public enum A2CharacterSource { EpCharacter, HermitianAxis, ExactRank }

public sealed record A2CharacterReading(string LocusId, EpCharacter.EpKind Kind, int Algebraic,
    int Geometric, double? RelativeDeparture, double? Radius, double? IsolationMargin, A2CharacterSource Source)
{
    public double? FullBlockHermiticityResidual { get; init; }
}

public sealed class A2CharacterUncertifiedException : InvalidOperationException
{
    public string LocusId { get; }
    public IReadOnlyList<A2CharacterReading> Readings { get; }

    public A2CharacterUncertifiedException(string locusId, string reason, IReadOnlyList<A2CharacterReading> readings,
        Exception? innerException = null)
        : base($"{locusId}: {reason}; readings=[{string.Join("; ", readings)}]", innerException)
    {
        LocusId = locusId;
        Readings = readings;
    }
}

/// <summary>Reads full-sector character at the exported seeds, without searching for a nearby coalescence.</summary>
public sealed class RouteBA2CharacterClassifier
{
    private readonly RouteBA2Inventory inventory;
    private readonly IA2ExactRankFallback fallback;
    private readonly Func<A2Root, A2QLocus, double, EpCharacter.Reading> readCharacter;

    public RouteBA2CharacterClassifier(RouteBA2Inventory inventory) : this(inventory, new NoExactRankFallback()) { }

    public RouteBA2CharacterClassifier(RouteBA2Inventory inventory, IA2ExactRankFallback fallback)
        : this(inventory, fallback, (root, locus, radius) => root.Parity == A2Parity.Even
            ? PathKMonodromyScout.CharacterizeAt(inventory.N - 1, locus.QPhysicalCSharpSeed, root.LambdaSeed, radius)
            : PathKMonodromyScout.CharacterizeAtROdd(inventory.N - 1, locus.QPhysicalCSharpSeed, root.LambdaSeed, radius)) { }

    private RouteBA2CharacterClassifier(RouteBA2Inventory inventory, IA2ExactRankFallback fallback,
        Func<A2Root, A2QLocus, double, EpCharacter.Reading> readCharacter)
    {
        this.inventory = inventory ?? throw new ArgumentNullException(nameof(inventory));
        this.fallback = fallback ?? throw new ArgumentNullException(nameof(fallback));
        this.readCharacter = readCharacter;
    }

    internal static RouteBA2CharacterClassifier ForcedAmbiguousForTest(RouteBA2Inventory inventory,
        IA2ExactRankFallback fallback) => new(inventory, fallback, (_, _, _) =>
            new(EpCharacter.EpKind.NearEp, 20, 1, 1, 1, 0, 1, [], 0));

    public IReadOnlyList<A2CharacterReading> ClassifyNonreal() => inventory.Roots
        .Where(root => root.RootKind == A2RootKind.Nonreal)
        .SelectMany(root => root.QLoci.Select(locus => Classify(root, locus))).ToArray();

    public IReadOnlyList<A2CharacterReading> ClassifyAll() => inventory.Roots
        .SelectMany(root => root.QLoci.Select(locus => Classify(root, locus))).ToArray();

    // Full parity sector in HS-orthonormal coordinates, including the AT strands.
    // The R-even orbit-sum matrix used by BuildLinear has a nontrivial metric at N=5;
    // build from the raw coherence block so Hermiticity has its ordinary meaning.
    internal static Matrix<Complex> FullParityBlock(int n, A2Parity parity, Complex q)
    {
        var reflection = F89PathKSeDeBlock.ReflectionPermutation(n);
        var representatives = Enumerable.Range(0, reflection.Length)
            .Where(i => parity == A2Parity.Even ? reflection[i] >= i : reflection[i] > i).ToArray();
        var basis = Matrix<Complex>.Build.Dense(reflection.Length, representatives.Length);
        for (int column = 0; column < representatives.Length; column++)
        {
            int i = representatives[column], partner = reflection[i];
            double weight = i == partner ? 1 : 1 / Math.Sqrt(2);
            basis[i, column] = weight;
            if (i != partner) basis[partner, column] = parity == A2Parity.Even ? weight : -weight;
        }
        var full = Matrix<Complex>.Build.DenseOfArray(F89PathKSeDeBlock.BuildFullBlock(n, q));
        return basis.ConjugateTranspose() * full * basis;
    }

    internal static double HermiticityResidual(Matrix<Complex> block) =>
        (block - block.ConjugateTranspose()).FrobeniusNorm() / Math.Max(1, block.FrobeniusNorm());

    public A2CharacterReading Classify(A2Root root, A2QLocus locus)
    {
        if (!inventory.Roots.Contains(root) || !root.QLoci.Contains(locus))
            throw new A2CharacterUncertifiedException(locus.Id, "Root/locus does not belong to this inventory", []);
        try { return ClassifyNumerical(root, locus); }
        catch (A2CharacterUncertifiedException error) when (root.RootKind != A2RootKind.NegativeReal)
        {
            if (!fallback.TryGet(locus.Id, out var certificate))
                throw new A2CharacterUncertifiedException(locus.Id,
                    $"{error.Message}; exact rank certificate absent; tolerances not loosened", error.Readings, error);
            RouteBA2Inventory.ValidateCertificate(certificate, inventory.Roots);
            if (certificate.LocusId != locus.Id)
                throw new InvalidDataException($"Exact rank certificate ID mismatch: {locus.Id}.");
            int geometric = certificate.MatrixDimension - certificate.ExactRank;
            return new(locus.Id, geometric == 1 ? EpCharacter.EpKind.Defective : EpCharacter.EpKind.Diabolic,
                2, geometric, null, null, null, A2CharacterSource.ExactRank);
        }
    }

    /// <summary>All positive-real and nonreal q lifts requiring future exact rank production.
    /// Deliberately bypasses the fallback so existing certificates cannot hide numerical ambiguity.</summary>
    public IReadOnlyList<string> FindAmbiguousNonHermitianLocusIds()
    {
        var ids = new List<string>();
        foreach (var root in inventory.Roots.Where(r => r.RootKind != A2RootKind.NegativeReal))
        foreach (var locus in root.QLoci)
        {
            try { _ = ClassifyNumerical(root, locus); }
            catch (A2CharacterUncertifiedException) { ids.Add(locus.Id); }
        }
        return ids.Order(StringComparer.Ordinal).ToArray();
    }

    private A2CharacterReading ClassifyNumerical(A2Root root, A2QLocus locus)
    {
        var readings = new List<A2CharacterReading>();
        try
        {
            if (!inventory.Roots.Contains(root) || !root.QLoci.Contains(locus))
                throw new A2CharacterUncertifiedException(locus.Id, "Root/locus does not belong to this inventory", readings);
            int k = inventory.N - 1;
            var q = locus.QPhysicalCSharpSeed;
            var lambda = root.LambdaSeed;
            double? hermiticity = null;
            System.Numerics.Complex[] roots;
            if (root.RootKind == A2RootKind.NegativeReal)
            {
                var block = FullParityBlock(inventory.N, root.Parity, q);
                hermiticity = HermiticityResidual(block);
                if (q.Real != 0 || lambda.Imaginary != 0 || !double.IsFinite(hermiticity.Value)
                    || hermiticity.Value >= 1e-12)
                    throw new A2CharacterUncertifiedException(locus.Id,
                        $"Imaginary-q full-sector Hermiticity failed: residual={hermiticity:G17}", readings);
                roots = block.Evd().EigenValues.ToArray();
            }
            else if (root.Parity == A2Parity.Even)
            {
                var (a, c) = PathKMonodromyScout.BuildLinear(inventory.N);
                roots = PathKMonodromyScout.AllRootsAt(a, c, q);
            }
            else roots = PathKMonodromyScout.AllRootsROdd(k, q);
            double[] distances = roots.Select(value => (value - lambda).Magnitude).Order().ToArray();
            if (distances.Length < 3 || distances.Any(d => !double.IsFinite(d)) || !(distances[2] > 100 * distances[1]))
                throw new A2CharacterUncertifiedException(locus.Id,
                    $"Pair is not isolated; nearest distances=[{string.Join(", ", distances.Take(3).Select(d => d.ToString("G17", System.Globalization.CultureInfo.InvariantCulture)))}]", readings);
            double d1 = distances[1], d2 = distances[2];
            if (hermiticity.HasValue)
            {
                // Multiplicity two comes from PSC1/S1 at the exact A2 root (the Python
                // producer boundary), not from doubling the w count. This full-sector
                // spectrum pins that pair to the exported lambda and excludes AT overlap.
                double radius = d1 + 0.5 * (d2 - d1);
                if (distances.Count(d => d < radius) != 2)
                    throw new A2CharacterUncertifiedException(locus.Id, "Hermitian pair is not unique", readings);
                return new(locus.Id, EpCharacter.EpKind.Diabolic, 2, 2, null, radius,
                    d2 - radius, A2CharacterSource.HermitianAxis)
                    { FullBlockHermiticityResidual = hermiticity };
            }
            foreach (double fraction in new[] { 0.25, 0.5, 0.75 })
            {
                double radius = d1 + fraction * (d2 - d1);
                if (distances.Count(d => d < radius) != 2)
                    throw new A2CharacterUncertifiedException(locus.Id, $"Radius {radius:G17} does not contain exactly two roots", readings);
                var character = readCharacter(root, locus, radius);
                double relative = character.Departure / Math.Max(1, character.CompressionNorm);
                readings.Add(new(locus.Id, character.Kind, character.Algebraic, character.Geometric,
                    relative, radius, d2 - radius, A2CharacterSource.EpCharacter));
            }
            bool consistent = readings.All(r => r.Algebraic == 2 && r.Kind == readings[0].Kind
                && r.Geometric == readings[0].Geometric && r.RelativeDeparture.HasValue
                && double.IsFinite(r.RelativeDeparture.Value)
                && (r.Kind == EpCharacter.EpKind.Diabolic && r.Geometric == 2 && r.RelativeDeparture < 1e-6
                    || r.Kind == EpCharacter.EpKind.Defective && r.Geometric == 1 && r.RelativeDeparture > 5e-2));
            if (!consistent)
                throw new A2CharacterUncertifiedException(locus.Id, "Full-sector character fails the fixed three-radius contract", readings.AsReadOnly());
            return readings[1];
        }
        catch (A2CharacterUncertifiedException) { throw; }
        catch (Exception error) when (error is ArithmeticException or ArgumentException or InvalidOperationException)
        {
            throw new A2CharacterUncertifiedException(locus.Id, $"Numerical characterization failed: {error.Message}", readings.AsReadOnly(), error);
        }
    }
}
