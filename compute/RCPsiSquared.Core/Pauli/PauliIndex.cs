namespace RCPsiSquared.Core.Pauli;

/// <summary>Encoding for the 4^N Pauli-string basis.
///
/// A flat index k ∈ [0, 4^N) decodes to an N-letter sequence by extracting 2 bits per
/// site, with site N−1 in the least-significant pair (i.e. site 0 is most-significant
/// in big-endian — same convention as <see cref="PauliString.Build"/>).
///
/// Per site, the (a, b) bit pair packs to <c>a + 2·b</c>:
/// <list type="bullet">
///   <item>0 = I (0, 0)</item>
///   <item>1 = X (1, 0)</item>
///   <item>2 = Z (0, 1)</item>
///   <item>3 = Y (1, 1)</item>
/// </list>
/// </summary>
public static class PauliIndex
{
    public static PauliLetter[] FromFlat(long k, int N)
    {
        if (N <= 0) throw new ArgumentOutOfRangeException(nameof(N));
        long max = 1L << (2 * N);
        if (k < 0 || k >= max)
            throw new ArgumentOutOfRangeException(nameof(k), $"k must be in [0, 4^N = {max}); got {k}");
        var letters = new PauliLetter[N];
        long kk = k;
        for (int site = N - 1; site >= 0; site--)
        {
            int i = (int)(kk & 0b11);
            kk >>= 2;
            letters[site] = (PauliLetter)i;
        }
        return letters;
    }

    public static long ToFlat(IReadOnlyList<PauliLetter> letters)
    {
        long k = 0;
        for (int i = 0; i < letters.Count; i++) k = (k << 2) | (long)letters[i];
        return k;
    }

    /// <summary>XY-weight (Σ bit_a) of the Pauli string.</summary>
    public static int TotalBitA(IReadOnlyList<PauliLetter> letters)
    {
        int s = 0;
        foreach (var L in letters) s += L.BitA();
        return s;
    }

    /// <summary>Π²-parity (Σ bit_b mod 2) of the Pauli string.</summary>
    public static int TotalBitBParity(IReadOnlyList<PauliLetter> letters)
    {
        int s = 0;
        foreach (var L in letters) s += L.BitB();
        return s & 1;
    }
}
