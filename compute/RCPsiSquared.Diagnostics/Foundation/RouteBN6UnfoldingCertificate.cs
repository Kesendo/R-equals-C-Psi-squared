using System.Globalization;
using System.Numerics;
using System.Security.Cryptography;
using System.Text;
using System.Text.Json;
using RCPsiSquared.Core.F89PathK;

namespace RCPsiSquared.Diagnostics.Foundation;

/// <summary>Source binding for the consumed complete bounded-CRT factorization.
/// This loader does not redo that proof. Its coefficient fingerprint identifies the
/// committed certificate; its fixture is compared coefficient by coefficient with
/// both freshly reconstructed Core pencils. Producer hashes use LF-normalized source
/// text so a Windows checkout's CRLF conversion does not change its identity.</summary>
internal static class RouteBN6UnfoldingCertificate
{
    private const string CoefficientSha256 = "e0a68ad94112d80d0df4d4cdb57fb5b9313502a5b862bcf8980df18c80e74832";
    private const string SourcePencilDigest = "cf1549c54a116132373e481d0ce7a7ea03409c07f75e3f6ab5b9fd0f738dc916";
    private const string Script = "simulations/route_b_n6_exact_unfolding.py";
    private const string SourceProducer = "simulations/route_b_a2_n6.py";
    private static readonly string[] SourceProducerDependencies =
    {
        "simulations/o2b_gcd_certificate.py",
        "simulations/o2b_krein_sign_law.py",
        "simulations/seed_existence_nullity_check.py",
    };
    private const string Fixture = "simulations/tests/fixtures/route_b_a2_n6_residual.json";
    internal const string Artifact = "simulations/results/route_b_n6_exact_unfolding.json";

    internal sealed record Input(BigInteger[] A2, bool SourcePencilMatches, string Description);

    internal static Input Load(string? path, RouteBN6ExactPencil even, RouteBN6ExactPencil odd,
        RouteBA2N6Inventory inventory)
    {
        try
        {
            string resolvedPath = Path.GetFullPath(path ?? Find(Artifact));
            byte[] payload = File.ReadAllBytes(resolvedPath);
            RejectDuplicateJsonProperties(payload);
            using var document = JsonDocument.Parse(payload);
            var layer = document.RootElement.GetProperty("layer");
            var coefficients = layer.GetProperty("a2_coefficients_lowest_first").EnumerateArray()
                .Select(c => c.GetString() ?? throw new InvalidDataException("Null F163 coefficient.")).ToArray();
            Require(Hash(string.Join("\n", coefficients)) == CoefficientSha256,
                "F163 consumed A2 coefficient fingerprint mismatch.");
            var a2 = coefficients.Select(c => BigInteger.Parse(c, CultureInfo.InvariantCulture)).ToArray();
            Require(a2.Length == 134 && !a2[^1].IsZero, "F163 A2 degree changed.");
            Require(layer.GetProperty("a2_degree").GetInt32() == 133
                && layer.GetProperty("a2_exponent").GetInt32() == 2
                && layer.GetProperty("discriminant_degree").GetInt32() == inventory.LayerIdentity.DiscriminantDegree,
                "F163 discriminant layer shape mismatch.");
            var modulus = BigInteger.Parse(layer.GetProperty("proof_modulus").GetString()!, CultureInfo.InvariantCulture);
            var bound = BigInteger.Parse(layer.GetProperty("proof_bound").GetString()!, CultureInfo.InvariantCulture);
            Require(modulus == BigInteger.Parse(inventory.LayerIdentity.ProofModulus, CultureInfo.InvariantCulture)
                && bound == BigInteger.Parse(inventory.LayerIdentity.ProofBound, CultureInfo.InvariantCulture)
                && modulus > 2 * bound && bound > 0,
                "F163 bounded-CRT input does not match the certified inventory.");
            var provenance = document.RootElement.GetProperty("provenance");
            Require(provenance.GetProperty("source_pencil_digest").GetString() == SourcePencilDigest,
                "F163 source pencil digest mismatch.");
            Require(provenance.GetProperty("script_sha256").GetString() == CanonicalUtf8LfSha256(File.ReadAllText(Find(Script))),
                "F163 producer source fingerprint mismatch.");
            Require(provenance.GetProperty("source_producer_sha256").GetString()
                == CanonicalUtf8LfSha256(File.ReadAllText(Find(SourceProducer))),
                "F163 exact source-producer fingerprint mismatch.");
            var dependencyProperties = provenance.GetProperty("source_producer_dependencies")
                .EnumerateObject().OrderBy(property => property.Name, StringComparer.Ordinal).ToArray();
            Require(dependencyProperties.Select(property => property.Name)
                    .SequenceEqual(SourceProducerDependencies),
                "F163 exact source-producer dependency manifest mismatch.");
            foreach (var dependency in dependencyProperties)
                Require(dependency.Value.GetString()
                        == CanonicalUtf8LfSha256(File.ReadAllText(Find(dependency.Name))),
                    $"F163 exact source-producer dependency fingerprint mismatch: {dependency.Name}.");
            string fixtureText = File.ReadAllText(Find(Fixture));
            Require(provenance.GetProperty("fixture_sha256").GetString() == CanonicalUtf8LfSha256(fixtureText),
                "F163 fixture fingerprint mismatch.");
            using var fixture = JsonDocument.Parse(fixtureText);
            bool matches = Matches(fixture.RootElement.GetProperty("rEven"), even)
                && Matches(fixture.RootElement.GetProperty("rOdd"), odd);
            Require(matches, "F163 certificate fixture differs from the rebuilt Core pencils.");
            return new(a2, matches,
                $"Full bounded-CRT factorization consumed from {resolvedPath}; coefficient SHA256 {CoefficientSha256}; " +
                "bound and modulus joined to RouteBA2N6Inventory, both fixture pencils rebuilt in Core; " +
                "wrapper, exact source-producer, its complete transitive local dependency manifest, and fixture " +
                "text hashes checked after LF normalization. " +
                "The CRT proof is not rerun at inspect time.");
        }
        catch (Exception error) when (error is JsonException or KeyNotFoundException or FormatException or InvalidOperationException)
        {
            throw new InvalidDataException("Invalid F163 consumed certificate.", error);
        }
    }

    private static bool Matches(JsonElement source, RouteBN6ExactPencil pencil) =>
        source.GetProperty("rOdd").GetBoolean() == pencil.ROdd
        && source.GetProperty("sectorDimension").GetInt32() == pencil.SectorDimension
        && RowsMatch(source.GetProperty("residualInT"), pencil.ResidualInT)
        && RowsMatch(source.GetProperty("atFactorInT"), pencil.AtFactorInT);

    private static bool RowsMatch(JsonElement element, BigInteger[][] rows)
    {
        var source = element.EnumerateArray().ToArray();
        return source.Length == rows.Length && source.Select((row, index) =>
            row.EnumerateArray().Select(c => BigInteger.Parse(c.GetString()!, CultureInfo.InvariantCulture))
                .SequenceEqual(rows[index])).All(equal => equal);
    }

    internal static bool ParityTransport(BigInteger[][] even, BigInteger[][] odd) =>
        even.Length == odd.Length && even.Select((row, i) =>
            odd[i].SequenceEqual(row.Select((value, j) => (j & 1) == 0 ? value : -value))).All(equal => equal);

    internal static string Find(string relative)
    {
        foreach (string start in new[] { AppContext.BaseDirectory, Directory.GetCurrentDirectory() })
        for (var directory = new DirectoryInfo(start); directory != null; directory = directory.Parent)
        {
            string path = Path.Combine(directory.FullName, relative);
            if (File.Exists(path)) return path;
        }
        throw new FileNotFoundException($"Cannot locate {relative}.");
    }

    internal static void Require(bool condition, string message)
    {
        if (!condition) throw new InvalidDataException(message);
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
                        ?? throw new InvalidDataException("Null F163 JSON property name.");
                    HashSet<string>? names = containers.Peek();
                    Require(names != null && names.Add(property), $"Duplicate JSON property: {property}.");
                    break;
            }
        }
        Require(containers.Count == 0, "Incomplete F163 JSON container.");
    }

    internal static string CanonicalUtf8LfSha256(string text)
    {
        ArgumentNullException.ThrowIfNull(text);
        string canonical = text.Replace("\r\n", "\n", StringComparison.Ordinal)
                               .Replace('\r', '\n');
        return Hash(canonical);
    }

    private static string Hash(string text) => Convert.ToHexString(SHA256.HashData(Encoding.UTF8.GetBytes(text))).ToLowerInvariant();
}
