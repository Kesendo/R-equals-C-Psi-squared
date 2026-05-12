namespace RCPsiSquared.Core.BlockSpectrum;

/// <summary>Shared helper: per-Hilbert-side F71 mirror lookup table. For each integer
/// <c>x ∈ [0, 2^N)</c> with bit <c>i</c> (LSB = site 0) holding <c>b_i</c>, returns
/// the table <c>mirror[x]</c> with bit <c>i</c> = bit <c>(N − 1 − i)</c> of x. This
/// lookup is the building block of the Liouville-space spatial-mirror map
/// <c>flat ↔ mirror[flat / d] · d + mirror[flat % d]</c> used by
/// <see cref="F71MirrorBlockRefinement"/>, <see cref="F71BilateralBlockRefinement"/>,
/// and the <c>block-spectrum</c> CLI smoke runs.</summary>
public static class F71MirrorIndexHelper
{
    /// <summary>Build the per-Hilbert-side F71 mirror lookup table of size 2^N.</summary>
    public static int[] BuildHilbertMirrorLookup(int N)
    {
        int d = 1 << N;
        var bits = new int[d];
        for (int x = 0; x < d; x++)
        {
            int m = 0;
            for (int i = 0; i < N; i++)
                if ((x & (1 << i)) != 0) m |= 1 << (N - 1 - i);
            bits[x] = m;
        }
        return bits;
    }
}
