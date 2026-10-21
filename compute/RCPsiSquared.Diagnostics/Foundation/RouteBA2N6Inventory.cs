using System.Globalization;
using System.Numerics;
using System.Security.Cryptography;
using System.Text;
using System.Text.Json;
using System.Text.Json.Serialization;

namespace RCPsiSquared.Diagnostics.Foundation;

public sealed record N6SeedText(string Real, string Imag);
public sealed record N6A2Degrees(int E, int O);
public sealed record N6LayerIdentity(
    int DiscriminantDegree, int Valuation, int A1Degree, int A2Degree,
    string Constant, string ProofModulus, string ProofBound);
public sealed record RouteBA2N6Locus(
    string Id, A2Parity Parity, int AlgebraicMultiplicity,
    ExactComplexBox TBox, ExactComplexBox QPhysicalCSharpBox,
    ExactComplexBox LambdaClearedBox, ExactComplexBox LambdaPhysicalBox,
    N6SeedText TSeed, N6SeedText QPhysicalCSharpSeed,
    N6SeedText LambdaClearedSeed, N6SeedText LambdaPhysicalSeed,
    string ConjugationPartnerId, string ParityPartnerId);

public static class ExactComplexBoxN6Extensions
{
    /// <summary>Apply qCSharp=-i*t exactly: (Re q, Im q)=(Im t,-Re t).</summary>
    public static ExactComplexBox RotateTToPhysicalQ(this ExactComplexBox box) =>
        new(box.ImLo, box.ImHi, box.ReHi.Negate(), box.ReLo.Negate());
}

/// <summary>
/// Strict consumer of the N=6 direct-t Route B carrier. It validates the carrier schema and exact
/// maps; it does not reproduce the producer's polynomial or root-count certificates.
/// </summary>
public sealed class RouteBA2N6Inventory
{
    private const string ExpectedModel = "open-uniform-XY-delta0-field0-gamma1-SEket-DEbra";
    private const string ExpectedSourcePencilDigest =
        "cf1549c54a116132373e481d0ce7a7ea03409c07f75e3f6ab5b9fd0f738dc916";
    private const string ExpectedLayerIdentityDigest =
        "3c70aa262dbee5e52f61475b8d600d3328696e33a131f388d1f93528db4c8a10";

    public int SchemaVersion { get; }
    public int N { get; }
    public N6A2Degrees A2Degrees { get; }
    public N6LayerIdentity LayerIdentity { get; }
    public IReadOnlyList<RouteBA2N6Locus> Loci { get; }
    public IReadOnlyList<A2ExactRankCertificate> ExactRankCertificates { get; }

    private RouteBA2N6Inventory(int schemaVersion, int n, N6A2Degrees a2Degrees,
        N6LayerIdentity layerIdentity, IReadOnlyList<RouteBA2N6Locus> loci,
        IReadOnlyList<A2ExactRankCertificate> exactRankCertificates)
    {
        SchemaVersion = schemaVersion;
        N = n;
        A2Degrees = a2Degrees;
        LayerIdentity = layerIdentity;
        Loci = loci;
        ExactRankCertificates = exactRankCertificates;
    }

    public static RouteBA2N6Inventory LoadDefault()
    {
        foreach (string start in new[] { AppContext.BaseDirectory, Directory.GetCurrentDirectory() })
        for (var directory = new DirectoryInfo(start); directory != null; directory = directory.Parent)
        {
            string path = Path.Combine(directory.FullName, "simulations", "results", "route_b_a2_n6.json");
            if (File.Exists(path)) return Load(path);
        }
        throw new FileNotFoundException("Cannot locate simulations/results/route_b_a2_n6.json.");
    }

    public static RouteBA2N6Inventory Load(string path)
    {
        try
        {
            byte[] payload = File.ReadAllBytes(path);
            RejectDuplicateJsonProperties(payload);
            var dto = JsonSerializer.Deserialize<InventoryDto>(payload, new JsonSerializerOptions
            {
                PropertyNamingPolicy = JsonNamingPolicy.CamelCase,
                UnmappedMemberHandling = JsonUnmappedMemberHandling.Disallow
            }) ?? throw new InvalidDataException("Empty N=6 A2 inventory.");

            Require(dto.SchemaVersion == 3, "Expected schemaVersion=3.");
            Require(dto.N == 6, "Expected n=6.");
            Require(dto.Model == ExpectedModel, "Unexpected N=6 model.");
            Require(dto.SourcePencilDigest == ExpectedSourcePencilDigest, "Unexpected source pencil digest.");
            Require(dto.Conventions.Parameter == "t=i*qCSharp;qCSharp=-i*t"
                && dto.Conventions.Eigenvalue == "Lambda=2*lambda"
                && dto.Conventions.CoefficientOrder == "lambda-lowest-first,t-lowest-first",
                "Unexpected N=6 conventions.");
            Require(dto.A2Degrees.E == 133 && dto.A2Degrees.O == 133, "Expected A2 degrees E=133 and O=133.");
            Require(dto.LayerIdentity.DiscriminantDegree == 926
                && dto.LayerIdentity.Valuation == 536
                && dto.LayerIdentity.A1Degree == 124
                && dto.LayerIdentity.A2Degree == 133, "Unexpected doubled-layer identity.");

            BigInteger constant = ParseCanonicalUnsigned(dto.LayerIdentity.Constant, "layer constant");
            BigInteger proofModulus = ParseCanonicalUnsigned(dto.LayerIdentity.ProofModulus, "proof modulus");
            BigInteger proofBound = ParseCanonicalUnsigned(dto.LayerIdentity.ProofBound, "proof bound");
            Require(constant > 0 && proofBound > 0 && proofModulus > 2 * proofBound,
                "Invalid doubled-layer proof integers.");
            Require(LayerIdentityDigest(dto.LayerIdentity) == ExpectedLayerIdentityDigest,
                "LayerIdentity does not match the canonical N=6 LayerIdentity.");

            Require(dto.Loci.Count == dto.A2Degrees.E + dto.A2Degrees.O,
                "Locus count does not equal the two A2 degrees.");
            Require(dto.Loci.Count == 266, "Expected exactly 266 N=6 loci.");
            Require(dto.ExactRankCertificates.Count == 0,
                "Schema 3 requires an empty exactRankCertificates array.");

            var loci = new List<RouteBA2N6Locus>(dto.Loci.Count);
            var ids = new HashSet<string>(StringComparer.Ordinal);
            int evenIndex = 0, oddIndex = 0;
            var rootBoxCounts = new Dictionary<A2Parity, int[]>
            {
                [A2Parity.Even] = new int[3],
                [A2Parity.Odd] = new int[3]
            };
            foreach (var item in dto.Loci)
            {
                A2Parity parity = ParseParity(item.Parity);
                int index = parity == A2Parity.Even ? evenIndex++ : oddIndex++;
                string expectedId = $"N6-{item.Parity}-A2-T-{index:D3}";
                Require(item.Id == expectedId && ids.Add(item.Id),
                    $"Invalid, unstable or duplicate N=6 locus ID: {item.Id}.");
                Require(item.AlgebraicMultiplicity == 2,
                    $"Expected algebraicMultiplicity=2: {item.Id}.");

                ExactComplexBox tBox = ParseBox(item.TBox);
                rootBoxCounts[parity][(int)ClassifyTBox(tBox, item.Id)]++;
                ExactComplexBox qBox = ParseBox(item.QPhysicalCSharpBox);
                ExactComplexBox lambdaClearedBox = ParseBox(item.LambdaClearedBox);
                ExactComplexBox lambdaPhysicalBox = ParseBox(item.LambdaPhysicalBox);
                Require(qBox == tBox.RotateTToPhysicalQ(), $"Exact q=-i*t box relation failed: {item.Id}.");
                Require(lambdaPhysicalBox == lambdaClearedBox.Half(),
                    $"Exact lambda=Lambda/2 box relation failed: {item.Id}.");

                var tSeed = ParseSeed(item.TSeed, tBox, item.Id + " t");
                var qSeed = ParseSeed(item.QPhysicalCSharpSeed, qBox, item.Id + " physical q");
                var lambdaClearedSeed = ParseSeed(item.LambdaClearedSeed, lambdaClearedBox,
                    item.Id + " cleared Lambda");
                var lambdaPhysicalSeed = ParseSeed(item.LambdaPhysicalSeed, lambdaPhysicalBox,
                    item.Id + " physical lambda");
                Require(qSeed == RotateSeed(tSeed), $"Exact q=-i*t seed relation failed: {item.Id}.");
                Require(lambdaPhysicalSeed == HalfSeed(lambdaClearedSeed),
                    $"Exact lambda=Lambda/2 seed relation failed: {item.Id}.");

                loci.Add(new(item.Id, parity, item.AlgebraicMultiplicity,
                    tBox, qBox, lambdaClearedBox, lambdaPhysicalBox,
                    new(item.TSeed.Real, item.TSeed.Imag),
                    new(item.QPhysicalCSharpSeed.Real, item.QPhysicalCSharpSeed.Imag),
                    new(item.LambdaClearedSeed.Real, item.LambdaClearedSeed.Imag),
                    new(item.LambdaPhysicalSeed.Real, item.LambdaPhysicalSeed.Imag),
                    item.ConjugationPartnerId, item.ParityPartnerId));
            }
            Require(evenIndex == 133 && oddIndex == 133, "Expected exactly 133 loci per parity.");
            foreach (A2Parity parity in Enum.GetValues<A2Parity>())
            {
                Require(rootBoxCounts[parity].SequenceEqual(new[] { 59, 37, 37 }),
                    $"Expected 59 real, 37 upper and 37 lower t boxes for parity {parity}.");
                RouteBA2N6Locus[] ordered = loci.Where(x => x.Parity == parity).ToArray();
                for (int index = 1; index < ordered.Length; index++)
                    Require(CompareTBoxes(ordered[index - 1].TBox, ordered[index].TBox) < 0,
                        $"Locus IDs do not follow lexicographic exact t-box order for parity {parity}.");
            }
            Require(loci.Take(133).All(x => x.Parity == A2Parity.Even)
                && loci.Skip(133).All(x => x.Parity == A2Parity.Odd),
                "Expected the stable E-then-O locus ordering.");

            ValidateSameParityDisjointness(loci);
            ValidatePartners(loci);

            var degrees = new N6A2Degrees(dto.A2Degrees.E, dto.A2Degrees.O);
            var identity = new N6LayerIdentity(dto.LayerIdentity.DiscriminantDegree,
                dto.LayerIdentity.Valuation, dto.LayerIdentity.A1Degree, dto.LayerIdentity.A2Degree,
                dto.LayerIdentity.Constant, dto.LayerIdentity.ProofModulus, dto.LayerIdentity.ProofBound);
            return new(dto.SchemaVersion, dto.N, degrees, identity, loci.AsReadOnly(),
                Array.Empty<A2ExactRankCertificate>());
        }
        catch (Exception error) when (error is JsonException or FormatException or OverflowException
            or ArgumentException or NullReferenceException)
        {
            throw new InvalidDataException($"Invalid N=6 Route B A2 artifact: {error.Message}", error);
        }
    }

    private static void ValidateSameParityDisjointness(IReadOnlyList<RouteBA2N6Locus> loci)
    {
        foreach (A2Parity parity in Enum.GetValues<A2Parity>())
        {
            RouteBA2N6Locus[] same = loci.Where(x => x.Parity == parity).ToArray();
            for (int i = 0; i < same.Length; i++)
            for (int j = i + 1; j < same.Length; j++)
                Require(!ClosuresOverlap(same[i].TBox, same[j].TBox),
                    $"Same-parity t boxes overlap: {same[i].Id}, {same[j].Id}.");
        }
    }

    private static void RejectDuplicateJsonProperties(ReadOnlySpan<byte> utf8Json)
    {
        var reader = new Utf8JsonReader(utf8Json);
        var containers = new Stack<HashSet<string>?>();
        while (reader.Read())
        {
            switch (reader.TokenType)
            {
                case JsonTokenType.StartObject:
                    containers.Push(new HashSet<string>(StringComparer.Ordinal));
                    break;
                case JsonTokenType.EndObject:
                    containers.Pop();
                    break;
                case JsonTokenType.StartArray:
                    containers.Push(null);
                    break;
                case JsonTokenType.EndArray:
                    containers.Pop();
                    break;
                case JsonTokenType.PropertyName:
                    string property = reader.GetString()
                        ?? throw new InvalidDataException("Null JSON property name.");
                    HashSet<string>? names = containers.Peek();
                    Require(names != null && names.Add(property), $"Duplicate JSON property: {property}.");
                    break;
            }
        }
        Require(containers.Count == 0, "Incomplete JSON container.");
    }

    // Immutable carrier binding only: UTF-8 of the seven named fields below, joined by LF with no
    // trailing LF. This does not re-run or claim to verify the producer's modular proof.
    private static string LayerIdentityDigest(LayerIdentityDto identity)
    {
        string canonical = string.Join('\n',
            $"discriminantDegree={identity.DiscriminantDegree.ToString(CultureInfo.InvariantCulture)}",
            $"valuation={identity.Valuation.ToString(CultureInfo.InvariantCulture)}",
            $"a1Degree={identity.A1Degree.ToString(CultureInfo.InvariantCulture)}",
            $"a2Degree={identity.A2Degree.ToString(CultureInfo.InvariantCulture)}",
            $"constant={identity.Constant}",
            $"proofModulus={identity.ProofModulus}",
            $"proofBound={identity.ProofBound}");
        return Convert.ToHexString(SHA256.HashData(Encoding.UTF8.GetBytes(canonical))).ToLowerInvariant();
    }

    private static void ValidatePartners(IReadOnlyList<RouteBA2N6Locus> loci)
    {
        var byId = loci.ToDictionary(x => x.Id, StringComparer.Ordinal);
        foreach (RouteBA2N6Locus locus in loci)
        {
            Require(byId.TryGetValue(locus.ConjugationPartnerId, out var conjugate),
                $"Missing conjugation partner: {locus.Id}.");
            Require(conjugate!.ConjugationPartnerId == locus.Id && conjugate.Parity == locus.Parity,
                $"Conjugation partner is not a same-parity involution: {locus.Id}.");
            Require(conjugate.TBox == Conjugate(locus.TBox)
                && conjugate.QPhysicalCSharpBox == Conjugate(locus.TBox).RotateTToPhysicalQ()
                && conjugate.LambdaClearedBox == Conjugate(locus.LambdaClearedBox)
                && conjugate.LambdaPhysicalBox == Conjugate(locus.LambdaPhysicalBox),
                $"Conjugation partner boxes disagree: {locus.Id}.");
            Require(byId.TryGetValue(locus.ParityPartnerId, out var parityPartner),
                $"Missing parity partner: {locus.Id}.");
            Require(parityPartner!.ParityPartnerId == locus.Id && parityPartner.Parity != locus.Parity,
                $"Parity partner is not an opposite-parity involution: {locus.Id}.");
            Require(parityPartner.TBox == locus.TBox.Negate()
                && parityPartner.QPhysicalCSharpBox == locus.QPhysicalCSharpBox.Negate()
                && parityPartner.LambdaClearedBox == locus.LambdaClearedBox
                && parityPartner.LambdaPhysicalBox == locus.LambdaPhysicalBox,
                $"Parity transport boxes disagree: {locus.Id}.");
        }
    }

    private static bool ClosuresOverlap(ExactComplexBox a, ExactComplexBox b) =>
        a.ReLo.CompareTo(b.ReHi) <= 0 && b.ReLo.CompareTo(a.ReHi) <= 0
        && a.ImLo.CompareTo(b.ImHi) <= 0 && b.ImLo.CompareTo(a.ImHi) <= 0;

    private static int CompareTBoxes(ExactComplexBox left, ExactComplexBox right)
    {
        int comparison = left.ReLo.CompareTo(right.ReLo);
        if (comparison != 0) return comparison;
        comparison = left.ReHi.CompareTo(right.ReHi);
        if (comparison != 0) return comparison;
        comparison = left.ImLo.CompareTo(right.ImLo);
        return comparison != 0 ? comparison : left.ImHi.CompareTo(right.ImHi);
    }

    private static TBoxKind ClassifyTBox(ExactComplexBox box, string id)
    {
        var zero = new ExactRational(BigInteger.Zero, BigInteger.One);
        bool containsZero = box.ReLo.CompareTo(zero) <= 0 && zero.CompareTo(box.ReHi) <= 0
            && box.ImLo.CompareTo(zero) <= 0 && zero.CompareTo(box.ImHi) <= 0;
        Require(!containsZero, $"The t box must exclude t=0: {id}.");
        if (box.ImLo == zero && box.ImHi == zero) return TBoxKind.Real;
        if (box.ImLo.CompareTo(zero) > 0) return TBoxKind.Upper;
        if (box.ImHi.CompareTo(zero) < 0) return TBoxKind.Lower;
        throw new InvalidDataException($"Invalid t root-box form (axis contact or straddling): {id}.");
    }

    private static ExactComplexBox Conjugate(ExactComplexBox box) =>
        new(box.ReLo, box.ReHi, box.ImHi.Negate(), box.ImLo.Negate());

    private static ExactComplexBox ParseBox(BoxDto box) =>
        new(ParseRational(box.Real.Lower), ParseRational(box.Real.Upper),
            ParseRational(box.Imag.Lower), ParseRational(box.Imag.Upper));

    private static ExactRational ParseRational(RationalDto dto)
    {
        BigInteger numerator = ParseCanonicalInteger(dto.Numerator, "rational numerator");
        BigInteger denominator = ParseCanonicalUnsigned(dto.Denominator, "rational denominator");
        Require(denominator > 0, "A rational denominator must be positive.");
        Require(BigInteger.GreatestCommonDivisor(BigInteger.Abs(numerator), denominator).IsOne,
            "A rational endpoint must be reduced.");
        return new(numerator, denominator);
    }

    private static ExactSeed ParseSeed(SeedDto dto, ExactComplexBox box, string id)
    {
        var seed = new ExactSeed(ParseDecimal(dto.Real), ParseDecimal(dto.Imag));
        Require(Contains(box, seed), $"Seed lies outside its exact box: {id}.");
        return seed;
    }

    private static ExactSeed RotateSeed(ExactSeed seed) => new(seed.Imag, seed.Real.Negate());
    private static ExactSeed HalfSeed(ExactSeed seed) => new(seed.Real.Half(), seed.Imag.Half());

    private static bool Contains(ExactComplexBox box, ExactSeed seed) =>
        box.ReLo.CompareTo(seed.Real) <= 0 && seed.Real.CompareTo(box.ReHi) <= 0
        && box.ImLo.CompareTo(seed.Imag) <= 0 && seed.Imag.CompareTo(box.ImHi) <= 0;

    private static ExactRational ParseDecimal(string text)
    {
        Require(!string.IsNullOrEmpty(text), "Empty decimal seed.");
        int cursor = 0;
        bool negative = false;
        if (text[cursor] is '+' or '-')
        {
            negative = text[cursor] == '-';
            cursor++;
            Require(cursor < text.Length, "Decimal seed contains only a sign.");
        }

        int integerStart = cursor;
        while (cursor < text.Length && char.IsAsciiDigit(text[cursor])) cursor++;
        Require(cursor > integerStart, "Decimal seed needs integer digits.");
        string integerDigits = text[integerStart..cursor];
        string fractionalDigits = "";
        if (cursor < text.Length && text[cursor] == '.')
        {
            cursor++;
            int fractionalStart = cursor;
            while (cursor < text.Length && char.IsAsciiDigit(text[cursor])) cursor++;
            Require(cursor > fractionalStart, "Decimal point must be followed by digits.");
            fractionalDigits = text[fractionalStart..cursor];
        }

        int exponent = 0;
        if (cursor < text.Length && text[cursor] is 'e' or 'E')
        {
            cursor++;
            bool exponentNegative = false;
            if (cursor < text.Length && text[cursor] is '+' or '-')
            {
                exponentNegative = text[cursor] == '-';
                cursor++;
            }
            int exponentStart = cursor;
            while (cursor < text.Length && char.IsAsciiDigit(text[cursor])) cursor++;
            Require(cursor > exponentStart, "Decimal exponent needs digits.");
            exponent = int.Parse(text[exponentStart..cursor], NumberStyles.None, CultureInfo.InvariantCulture);
            if (exponentNegative) exponent = -exponent;
        }
        Require(cursor == text.Length, "Unexpected character in decimal seed.");
        Require(Math.Abs((long)exponent) <= 10_000, "Decimal exponent is unreasonably large.");

        BigInteger significand = BigInteger.Parse(integerDigits + fractionalDigits,
            NumberStyles.None, CultureInfo.InvariantCulture);
        if (negative) significand = -significand;
        int scale = checked(fractionalDigits.Length - exponent);
        return scale >= 0
            ? new(significand, BigInteger.Pow(10, scale))
            : new(significand * BigInteger.Pow(10, -scale), BigInteger.One);
    }

    private static BigInteger ParseCanonicalInteger(string text, string label)
    {
        Require(!string.IsNullOrEmpty(text), $"Missing {label}.");
        int firstDigit = text[0] == '-' ? 1 : 0;
        Require(firstDigit < text.Length && text[firstDigit..].All(char.IsAsciiDigit),
            $"Invalid {label}.");
        Require(text[firstDigit] != '0' || text.Length - firstDigit == 1,
            $"Non-canonical {label}.");
        BigInteger value = BigInteger.Parse(text, NumberStyles.AllowLeadingSign, CultureInfo.InvariantCulture);
        Require(!value.IsZero || text == "0", $"Non-canonical {label}.");
        return value;
    }

    private static BigInteger ParseCanonicalUnsigned(string text, string label)
    {
        Require(!string.IsNullOrEmpty(text) && text.All(char.IsAsciiDigit), $"Invalid {label}.");
        Require(text[0] != '0' || text.Length == 1, $"Non-canonical {label}.");
        return BigInteger.Parse(text, NumberStyles.None, CultureInfo.InvariantCulture);
    }

    private static A2Parity ParseParity(string value) => value switch
    {
        "E" => A2Parity.Even,
        "O" => A2Parity.Odd,
        _ => throw new InvalidDataException($"Unknown N=6 parity: {value}.")
    };

    private static void Require(bool condition, string message)
    {
        if (!condition) throw new InvalidDataException(message);
    }

    private readonly record struct ExactSeed(ExactRational Real, ExactRational Imag);
    private enum TBoxKind { Real, Upper, Lower }

    private sealed class InventoryDto
    {
        public required int SchemaVersion { get; init; }
        public required int N { get; init; }
        public required string Model { get; init; }
        public required ConventionsDto Conventions { get; init; }
        public required A2DegreesDto A2Degrees { get; init; }
        public required LayerIdentityDto LayerIdentity { get; init; }
        public required List<LocusDto> Loci { get; init; }
        public required List<JsonElement> ExactRankCertificates { get; init; }
        public required string SourcePencilDigest { get; init; }
    }

    private sealed class ConventionsDto
    {
        public required string Parameter { get; init; }
        public required string Eigenvalue { get; init; }
        public required string CoefficientOrder { get; init; }
    }

    private sealed class A2DegreesDto
    {
        [JsonPropertyName("E")] public required int E { get; init; }
        [JsonPropertyName("O")] public required int O { get; init; }
    }

    private sealed class LayerIdentityDto
    {
        public required int DiscriminantDegree { get; init; }
        public required int Valuation { get; init; }
        public required int A1Degree { get; init; }
        public required int A2Degree { get; init; }
        public required string Constant { get; init; }
        public required string ProofModulus { get; init; }
        public required string ProofBound { get; init; }
    }

    private sealed class LocusDto
    {
        public required string Id { get; init; }
        public required string Parity { get; init; }
        public required int AlgebraicMultiplicity { get; init; }
        public required BoxDto TBox { get; init; }
        public required BoxDto QPhysicalCSharpBox { get; init; }
        public required BoxDto LambdaClearedBox { get; init; }
        public required BoxDto LambdaPhysicalBox { get; init; }
        public required SeedDto TSeed { get; init; }
        public required SeedDto QPhysicalCSharpSeed { get; init; }
        public required SeedDto LambdaClearedSeed { get; init; }
        public required SeedDto LambdaPhysicalSeed { get; init; }
        public required string ConjugationPartnerId { get; init; }
        public required string ParityPartnerId { get; init; }
    }

    private sealed class BoxDto
    {
        public required IntervalDto Real { get; init; }
        public required IntervalDto Imag { get; init; }
    }

    private sealed class IntervalDto
    {
        public required RationalDto Lower { get; init; }
        public required RationalDto Upper { get; init; }
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
