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

    /// <summary>Walk a sector's Liouville flat-indices and partition them into F71 orbits:
    /// fixed points (where <c>mirror(flat) == flat</c>) and size-2 pairs <c>(S, Ps)</c> with
    /// <c>S &lt; Ps</c>. Iteration order is ascending by flat index, so the result is
    /// reproducible across runs.
    ///
    /// <para>Used by F71-refined block-decomposition primitives to build the F71-orbit
    /// basis: <c>(|flat⟩ ± |mirror⟩)/√2</c> for size-2 orbits; <c>|flat⟩</c> alone for fixed
    /// points. The caller supplies the Liouville-space mirror function (typically
    /// <c>flat → mirrorBits[flat/d]·d + mirrorBits[flat%d]</c> built from
    /// <see cref="BuildHilbertMirrorLookup"/>).</para></summary>
    public static (List<int> FixedPoints, List<(int S, int Ps)> Pairs) FindOrbitsInSector(
        IReadOnlyList<int> sectorFlat,
        Func<int, int> mirror)
    {
        var fixedPoints = new List<int>();
        var pairs = new List<(int S, int Ps)>();
        var seen = new HashSet<int>(sectorFlat.Count);
        foreach (int flat in sectorFlat.OrderBy(x => x))
        {
            if (!seen.Add(flat)) continue;
            int m = mirror(flat);
            if (m == flat)
            {
                fixedPoints.Add(flat);
            }
            else
            {
                int s = Math.Min(flat, m);
                int ps = Math.Max(flat, m);
                pairs.Add((s, ps));
                seen.Add(ps);
            }
        }
        return (fixedPoints, pairs);
    }
}
