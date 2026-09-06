namespace RCPsiSquared.Diagnostics.Foundation;

/// <summary>An artifact producer's rank certificate. Consumers validate its shape, not its algebra.</summary>
public sealed record A2ExactRankCertificate(string LocusId, int MatrixDimension, int ExactRank,
    int AlgebraicMultiplicity, string Criterion, int? DeletedRow, int? DeletedColumn,
    string? RemainderPolynomial, string? ZeroRemaindersSha256, string? Transport);

public interface IA2ExactRankFallback
{
    bool TryGet(string locusId, out A2ExactRankCertificate certificate);
}

public sealed class NoExactRankFallback : IA2ExactRankFallback
{
    public bool TryGet(string locusId, out A2ExactRankCertificate certificate)
    {
        certificate = null!;
        return false;
    }
}

/// <summary>Exact ordinal ID lookup; this consumer neither computes nor verifies polynomial minors.</summary>
public sealed class ArtifactExactRankFallback(RouteBA2Inventory inventory) : IA2ExactRankFallback
{
    private readonly IReadOnlyDictionary<string, A2ExactRankCertificate> certificates =
        inventory.ExactRankCertificates.ToDictionary(c => c.LocusId, StringComparer.Ordinal);

    public bool TryGet(string locusId, out A2ExactRankCertificate certificate) =>
        certificates.TryGetValue(locusId, out certificate!);
}
