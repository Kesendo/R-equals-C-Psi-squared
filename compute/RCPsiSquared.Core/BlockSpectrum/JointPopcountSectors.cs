using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.BlockSpectrum;

/// <summary>JointPopcountSectors (Tier 1 derived; 2026-05-11):
///
/// <para>The XY+Z-dephasing Liouvillian L on N qubits is exactly block-diagonal under
/// the joint label (popcount_col, popcount_row) of Liouville-space basis indices.
/// Block (p_c, p_r) has size C(N, p_c) · C(N, p_r); the partition has (N+1)² sectors
/// summing to 4^N.</para>
///
/// <list type="bullet">
///   <item>XY swap |01⟩↔|10⟩ preserves Hilbert-side popcount.</item>
///   <item>Z-dephasing is diagonal in computational basis, trivially preserves popcount.</item>
///   <item>Verified bit-exact off-block-Frobenius = 0 at N=2, N=3, N=4 (probe in
///         JointPopcountSectorsTests).</item>
/// </list>
///
/// <para>Cubic-cost speedup over full diagonalization scales as
/// (4^N)³ / Σ (C(N,p_c)·C(N,p_r))³. At N=8: max block 4900 (vs full 65536); max-block
/// memory ≈ 0.38 GB (vs full 68.7 GB); cubic-cost speedup ≈ 515×.</para>
///
/// <para>Distinct from F86's <c>CoherenceBlock(N, n)</c>, which is one specific
/// sub-sector with |p_c − p_r| ≤ 1; JointPopcountSectors is the full (N+1)² family.</para>
///
/// <para>Anchors: <c>compute/RCPsiSquared.Core/BlockSpectrum/JointPopcountSectorBuilder.cs</c>
/// (basis permutation builder), <c>compute/RCPsiSquared.Core.Tests/BlockSpectrum/JointPopcountSectorsTests.cs</c>
/// (off-block-Frobenius=0 verification at N=2,3,4).</para></summary>
public sealed class JointPopcountSectors : Claim
{
    /// <summary>Number of sectors as a function of N: (N+1)².</summary>
    public static int SectorCount(int N) => (N + 1) * (N + 1);

    /// <summary>Size of the (p_c, p_r) sector: C(N, p_c) · C(N, p_r).</summary>
    public static long SectorSize(int N, int pCol, int pRow)
    {
        if (N < 0 || pCol < 0 || pCol > N || pRow < 0 || pRow > N)
            throw new ArgumentOutOfRangeException();
        return Binomial(N, pCol) * Binomial(N, pRow);
    }

    /// <summary>Max-block size across all (p_c, p_r); achieved at (⌊N/2⌋, ⌊N/2⌋).</summary>
    public static long MaxSectorSize(int N) => SectorSize(N, N / 2, N / 2);

    private static long Binomial(int n, int k)
    {
        if (k < 0 || k > n) return 0;
        if (k == 0 || k == n) return 1;
        k = Math.Min(k, n - k);
        long c = 1;
        for (int i = 0; i < k; i++) c = c * (n - i) / (i + 1);
        return c;
    }

    public JointPopcountSectors()
        : base("JointPopcountSectors: XY+Z-dephasing Liouvillian is exactly block-diagonal in joint (popcount_col, popcount_row) sectors; (N+1)² blocks summing to 4^N; bit-exact verified at N=2,3,4.",
               Tier.Tier1Derived,
               "U(1)×U(1) per-side popcount conservation; XY swap preserves Hilbert-side popcount; Z-dephasing diagonal trivially preserves; verified off-block-Frobenius = 0 at N=2,3,4 in JointPopcountSectorsTests")
    {
    }

    public override string DisplayName =>
        "JointPopcountSectors: (N+1)² block-diagonal decomposition of XY+Z-deph L";

    public override string Summary =>
        $"XY+Z-deph L block-diagonal in (popcount_col, popcount_row); (N+1)² sectors; N=8 max-block 4900 (vs full 65536) ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("Sector count formula", summary: "(N+1)²");
            yield return new InspectableNode("Sector size formula", summary: "C(N, p_c) · C(N, p_r)");
            yield return new InspectableNode("N=2 sectors", summary: $"{SectorCount(2)} sectors, max-block {MaxSectorSize(2)} (full D²=16)");
            yield return new InspectableNode("N=3 sectors", summary: $"{SectorCount(3)} sectors, max-block {MaxSectorSize(3)} (full D²=64)");
            yield return new InspectableNode("N=8 sectors", summary: $"{SectorCount(8)} sectors, max-block {MaxSectorSize(8)} (full D²=65536)");
        }
    }
}
