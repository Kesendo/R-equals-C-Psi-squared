using System.Globalization;
using System.Numerics;
using System.Text.Json;
using System.Text.Json.Serialization;

namespace RCPsiSquared.Diagnostics.Foundation;

public enum A2Parity { Even, Odd }
public enum A2RootKind { NegativeReal, PositiveReal, Nonreal }

/// <summary>A reduced rational with a strictly positive denominator.</summary>
public readonly record struct ExactRational : IComparable<ExactRational>
{
    public BigInteger Numerator { get; }
    public BigInteger Denominator { get; }

    public ExactRational(BigInteger numerator, BigInteger denominator)
    {
        if (denominator <= 0) throw new InvalidDataException("A rational denominator must be positive.");
        BigInteger gcd = BigInteger.GreatestCommonDivisor(numerator, denominator);
        Numerator = numerator / gcd;
        Denominator = denominator / gcd;
    }

    public int CompareTo(ExactRational other) =>
        (Numerator * other.Denominator).CompareTo(other.Numerator * Denominator);

    public ExactRational Half() => new(Numerator, 2 * Denominator);
    public ExactRational Negate() => new(-Numerator, Denominator);
    public static ExactRational operator +(ExactRational a, ExactRational b) =>
        new(a.Numerator * b.Denominator + b.Numerator * a.Denominator, a.Denominator * b.Denominator);

    /// <summary>Round the exact ratio to binary64, avoiding two independently rounded integer casts.</summary>
    public double ToDouble()
    {
        if (Numerator.IsZero) return 0;
        BigInteger numerator = BigInteger.Abs(Numerator);
        int exponent = checked((int)(numerator.GetBitLength() - Denominator.GetBitLength()));
        bool belowPower = exponent >= 0
            ? numerator < (Denominator << exponent)
            : (numerator << -exponent) < Denominator;
        if (belowPower) exponent--;
        int shift = Math.Min(1074, 52 - exponent);
        BigInteger dividend = shift >= 0 ? numerator << shift : numerator;
        BigInteger divisor = shift >= 0 ? Denominator : Denominator << -shift;
        BigInteger mantissa = BigInteger.DivRem(dividend, divisor, out BigInteger remainder);
        int halfway = (2 * remainder).CompareTo(divisor);
        if (halfway > 0 || halfway == 0 && !mantissa.IsEven) mantissa++;
        return Numerator.Sign * Math.ScaleB((double)mantissa, -shift);
    }
}

public sealed record ExactComplexBox
{
    public ExactRational ReLo { get; }
    public ExactRational ReHi { get; }
    public ExactRational ImLo { get; }
    public ExactRational ImHi { get; }

    public ExactComplexBox(ExactRational reLo, ExactRational reHi, ExactRational imLo, ExactRational imHi)
    {
        if (reLo.CompareTo(reHi) > 0 || imLo.CompareTo(imHi) > 0)
            throw new InvalidDataException("Inverted complex box.");
        ReLo = reLo; ReHi = reHi; ImLo = imLo; ImHi = imHi;
    }

    public Complex Midpoint => new((ReLo + ReHi).Half().ToDouble(), (ImLo + ImHi).Half().ToDouble());
    public ExactComplexBox Half() => new(ReLo.Half(), ReHi.Half(), ImLo.Half(), ImHi.Half());
    public ExactComplexBox Negate() => new(ReHi.Negate(), ReLo.Negate(), ImHi.Negate(), ImLo.Negate());

    // Decimal seed text has fewer digits than the q bounds. The numerical contract is containment
    // after correctly rounding exact endpoints to binary64, with no fitted tolerance.
    public bool ContainsRounded(Complex seed) => double.IsFinite(seed.Real) && double.IsFinite(seed.Imaginary)
        && seed.Real >= ReLo.ToDouble() && seed.Real <= ReHi.ToDouble()
        && seed.Imaginary >= ImLo.ToDouble() && seed.Imaginary <= ImHi.ToDouble();
}

public sealed record A2QLocus(string Id, ExactComplexBox QUnitHopBox,
    ExactComplexBox QPhysicalCSharpBox, Complex QPhysicalCSharpSeed);

public sealed record A2Root(string Id, A2Parity Parity, A2RootKind RootKind,
    ExactComplexBox WBox, ExactComplexBox LambdaBox, Complex LambdaSeed, IReadOnlyList<A2QLocus> QLoci);

/// <summary>The exact N=5 A2 inventory. This validates its carrier schema, not the producer's algebraic certificate.</summary>
public sealed class RouteBA2Inventory
{
    public int SchemaVersion { get; }
    public int N { get; }
    public IReadOnlyList<A2Root> Roots { get; }
    public IReadOnlyList<A2ExactRankCertificate> ExactRankCertificates { get; }

    private RouteBA2Inventory(int schemaVersion, int n, IReadOnlyList<A2Root> roots,
        IReadOnlyList<A2ExactRankCertificate> exactRankCertificates)
    {
        SchemaVersion = schemaVersion;
        N = n;
        Roots = roots;
        ExactRankCertificates = exactRankCertificates;
    }

    public static RouteBA2Inventory LoadDefault()
    {
        foreach (string start in new[] { AppContext.BaseDirectory, Directory.GetCurrentDirectory() })
        for (var directory = new DirectoryInfo(start); directory != null; directory = directory.Parent)
        {
            string path = Path.Combine(directory.FullName, "simulations", "results", "route_b_a2_n5.json");
            if (File.Exists(path)) return Load(path);
        }
        throw new FileNotFoundException("Cannot locate simulations/results/route_b_a2_n5.json.");
    }

    public static RouteBA2Inventory Load(string path)
    {
        try
        {
            var dto = JsonSerializer.Deserialize<InventoryDto>(File.ReadAllText(path), new JsonSerializerOptions
            {
                PropertyNamingPolicy = JsonNamingPolicy.CamelCase,
                UnmappedMemberHandling = JsonUnmappedMemberHandling.Disallow
            }) ?? throw new InvalidDataException("Empty A2 inventory.");
            Require(dto.SchemaVersion == 1 && dto.N == 5, "Expected schemaVersion=1 and n=5.");
            Require(dto.Conventions != null && dto.Conventions.W == "qUnitHop^2"
                && dto.Conventions.QPhysicalCSharp == "qUnitHop/2", "Unexpected parameter conventions.");
            Require(dto.Sectors != null && dto.Sectors.Count == 2, "Expected two parity sectors.");
            var roots = new List<A2Root>();
            var ids = new HashSet<string>(StringComparer.Ordinal);
            var parities = new HashSet<A2Parity>();
            foreach (var sector in dto.Sectors!)
            {
                A2Parity parity = ParseParity(sector.Parity);
                Require(parities.Add(parity), "Duplicate parity sector.");
                Require(sector.A2Roots != null && sector.A2Roots.Count == (parity == A2Parity.Even ? 16 : 13),
                    "Expected E=16 and O=13 roots.");
                for (int index = 0; index < sector.A2Roots!.Count; index++)
                {
                    var root = sector.A2Roots[index];
                    string expectedId = $"N5-{sector.Parity}-A2-W-{index:D3}";
                    Require(root.Id == expectedId && ids.Add(root.Id), $"Invalid or duplicate root ID: {root.Id}.");
                    Require(ParseParity(root.Parity) == parity, $"Parity mismatch: {root.Id}.");
                    A2RootKind kind = root.RootKind switch
                    {
                        "negativeReal" => A2RootKind.NegativeReal,
                        "positiveReal" => A2RootKind.PositiveReal,
                        "nonreal" => A2RootKind.Nonreal,
                        _ => throw new InvalidDataException($"Unknown root kind: {root.RootKind}.")
                    };
                    var w = ParseBox(root.WBox);
                    var lambda = ParseBox(root.LambdaBox);
                    ValidateSeed(root.WSeed, w, root.Id + " w");
                    var lambdaSeed = ValidateSeed(root.LambdaSeed, lambda, root.Id + " lambda");
                    bool real = w.ImLo.Numerator.IsZero && w.ImHi.Numerator.IsZero;
                    Require(kind switch
                    {
                        A2RootKind.NegativeReal => real && w.ReHi.Numerator.Sign < 0,
                        A2RootKind.PositiveReal => real && w.ReLo.Numerator.Sign > 0,
                        _ => w.ImLo.Numerator.Sign > 0 || w.ImHi.Numerator.Sign < 0
                    }, $"Root kind disagrees with its exact w box: {root.Id}.");
                    Require(root.QLoci != null && root.QLoci.Count == 2, $"Expected two q loci: {root.Id}.");
                    var loci = new List<A2QLocus>();
                    for (int lift = 0; lift < 2; lift++)
                    {
                        var q = root.QLoci![lift];
                        string suffix = lift == 0 ? "plus" : "minus";
                        Require(q.Id == root.Id + "-Q-" + suffix && ids.Add(q.Id), $"Invalid or duplicate q ID: {q.Id}.");
                        var unit = ParseBox(q.QUnitHopBox);
                        var physical = ParseBox(q.QPhysicalCSharpBox);
                        Require(physical == unit.Half(), $"Exact q-half relation failed: {q.Id}.");
                        ValidateSeed(q.QUnitHopSeed, unit, q.Id + " unit q");
                        var physicalSeed = ValidateSeed(q.QPhysicalCSharpSeed, physical, q.Id + " physical q");
                        loci.Add(new(q.Id, unit, physical, physicalSeed));
                    }
                    Require(loci[1].QUnitHopBox == loci[0].QUnitHopBox.Negate(), $"q lifts are not opposite: {root.Id}.");
                    roots.Add(new(root.Id, parity, kind, w, lambda, lambdaSeed, loci.AsReadOnly()));
                }
            }
            Require(roots.Count(r => r.RootKind == A2RootKind.NegativeReal) == 12
                && roots.Count(r => r.RootKind == A2RootKind.PositiveReal) == 1
                && roots.Count(r => r.RootKind == A2RootKind.Nonreal) == 16, "Expected root-kind counts 12/1/16.");
            Require(dto.ExactRankCertificates != null, "Expected exactRankCertificates array.");
            var certificates = new List<A2ExactRankCertificate>();
            var certificateIds = new HashSet<string>(StringComparer.Ordinal);
            foreach (var c in dto.ExactRankCertificates!)
            {
                var certificate = new A2ExactRankCertificate(c.LocusId, c.MatrixDimension, c.ExactRank,
                    c.AlgebraicMultiplicity, c.Criterion, c.DeletedRow, c.DeletedColumn,
                    c.RemainderPolynomial, c.ZeroRemaindersSha256, c.Transport);
                ValidateCertificate(certificate, roots);
                Require(certificateIds.Add(c.LocusId), $"Duplicate exact rank certificate: {c.LocusId}.");
                certificates.Add(certificate);
            }
            return new(dto.SchemaVersion, dto.N, roots.AsReadOnly(), certificates.AsReadOnly());
        }
        catch (Exception error) when (error is JsonException or FormatException or OverflowException or NullReferenceException)
        {
            throw new InvalidDataException($"Invalid Route B A2 artifact: {error.Message}", error);
        }
    }

    internal static void ValidateCertificate(A2ExactRankCertificate certificate, IReadOnlyList<A2Root> roots)
    {
        Require(certificate != null, "Null exact rank certificate.");
        var c = certificate!;
        var root = roots.SingleOrDefault(r => r.QLoci.Any(q => q.Id == c.LocusId));
        Require(root != null && root.RootKind != A2RootKind.NegativeReal,
            $"Exact rank certificate must name an eligible non-Hermitian locus: {c.LocusId}.");
        // Full N=5 parity sectors, including AT subspaces, have dimensions 26 (E) and 24 (O).
        int dimension = root!.Parity == A2Parity.Even ? 26 : 24;
        Require(c.MatrixDimension == dimension && c.AlgebraicMultiplicity == 2,
            $"Exact rank certificate dimension/multiplicity mismatch: {c.LocusId}.");
        Require(c.ExactRank == dimension - 1 || c.ExactRank == dimension - 2,
            $"Exact rank must be n-1 or n-2: {c.LocusId}.");
        Require(c.Transport == null || c.Transport == "q-negation quotient automorphism",
            $"Unknown exact rank transport: {c.LocusId}.");
        if (c.ExactRank == dimension - 1)
        {
            Require(c.Criterion == "nonzero-codimension-one-minor"
                && c.DeletedRow is >= 0 && c.DeletedRow < dimension
                && c.DeletedColumn is >= 0 && c.DeletedColumn < dimension
                && !string.IsNullOrWhiteSpace(c.RemainderPolynomial)
                && !(BigInteger.TryParse(c.RemainderPolynomial, NumberStyles.Integer,
                    CultureInfo.InvariantCulture, out var constant) && constant.IsZero),
                $"Invalid nonzero minor certificate: {c.LocusId}.");
        }
        else
        {
            Require(c.Criterion == "all-codimension-one-minors-zero"
                && c.DeletedRow == null && c.DeletedColumn == null && c.RemainderPolynomial == null
                && c.ZeroRemaindersSha256 is { Length: 64 } digest && digest.All(Uri.IsHexDigit),
                $"Invalid all-minors-zero certificate: {c.LocusId}.");
        }
    }

    private static void Require(bool condition, string message)
    {
        if (!condition) throw new InvalidDataException(message);
    }

    private static A2Parity ParseParity(string value) => value switch
    {
        "E" => A2Parity.Even,
        "O" => A2Parity.Odd,
        _ => throw new InvalidDataException($"Unknown parity: {value}.")
    };

    private static ExactComplexBox ParseBox(BoxDto dto) => new(
        ParseRational(dto.ReLo), ParseRational(dto.ReHi), ParseRational(dto.ImLo), ParseRational(dto.ImHi));

    private static ExactRational ParseRational(RationalDto dto) => new(
        BigInteger.Parse(dto.Numerator, NumberStyles.AllowLeadingSign, CultureInfo.InvariantCulture),
        BigInteger.Parse(dto.Denominator, NumberStyles.AllowLeadingSign, CultureInfo.InvariantCulture));

    private static Complex ValidateSeed(SeedDto dto, ExactComplexBox box, string id)
    {
        var seed = new Complex(double.Parse(dto.Real, NumberStyles.Float, CultureInfo.InvariantCulture),
            double.Parse(dto.Imag, NumberStyles.Float, CultureInfo.InvariantCulture));
        Require(box.ContainsRounded(seed), $"Seed lies outside its exact box rounded to binary64: {id}.");
        return seed;
    }

    private sealed class InventoryDto
    {
        public required int SchemaVersion { get; init; }
        public required int N { get; init; }
        public required ConventionsDto Conventions { get; init; }
        public required List<SectorDto> Sectors { get; init; }
        public required List<ExactRankCertificateDto> ExactRankCertificates { get; init; }
    }
    private sealed class ExactRankCertificateDto
    {
        public required string LocusId { get; init; }
        public required int MatrixDimension { get; init; }
        public required int ExactRank { get; init; }
        public required int AlgebraicMultiplicity { get; init; }
        public required string Criterion { get; init; }
        public required int? DeletedRow { get; init; }
        public required int? DeletedColumn { get; init; }
        public required string? RemainderPolynomial { get; init; }
        public required string? ZeroRemaindersSha256 { get; init; }
        public required string? Transport { get; init; }
    }
    private sealed class ConventionsDto
    {
        public required string W { get; init; }
        public required string QPhysicalCSharp { get; init; }
    }
    private sealed class SectorDto
    {
        public required string Parity { get; init; }
        public required List<RootDto> A2Roots { get; init; }
    }
    private sealed class RootDto
    {
        public required string Id { get; init; }
        public required string Parity { get; init; }
        public required string RootKind { get; init; }
        public required BoxDto WBox { get; init; }
        public required BoxDto LambdaBox { get; init; }
        public required SeedDto WSeed { get; init; }
        public required SeedDto LambdaSeed { get; init; }
        public required List<QLocusDto> QLoci { get; init; }
    }
    private sealed class QLocusDto
    {
        public required string Id { get; init; }
        public required BoxDto QUnitHopBox { get; init; }
        public required BoxDto QPhysicalCSharpBox { get; init; }
        public required SeedDto QUnitHopSeed { get; init; }
        public required SeedDto QPhysicalCSharpSeed { get; init; }
    }
    private sealed class BoxDto
    {
        public required RationalDto ReLo { get; init; }
        public required RationalDto ReHi { get; init; }
        public required RationalDto ImLo { get; init; }
        public required RationalDto ImHi { get; init; }
    }
    private sealed class RationalDto
    {
        public required string Numerator { get; init; }
        public required string Denominator { get; init; }
    }
    private sealed class SeedDto
    {
        public required string Real { get; init; }
        public required string Imag { get; init; }
    }
}
