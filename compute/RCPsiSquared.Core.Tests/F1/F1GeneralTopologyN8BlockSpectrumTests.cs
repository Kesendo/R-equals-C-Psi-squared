using System.Diagnostics;
using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.BlockSpectrum;
using RCPsiSquared.Core.ChainSystems;
using RCPsiSquared.Core.F1;
using RCPsiSquared.Core.Pauli;
using RCPsiSquared.Core.Tests.TestHelpers;
using Xunit;
using Xunit.Abstractions;
using ComplexMatrix = MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex>;

namespace RCPsiSquared.Core.Tests.F1;

/// <summary>F1 general-topology verification at N=8 via
/// <see cref="LiouvillianBlockSpectrum.ComputeSpectrumPerBlock"/>. Dogfood of the
/// cubic-cost block-decomposition speedup the infrastructure was designed for:
/// full L_vec at N=8 is 65536² × 16 B = 68.7 GB, past the .NET 2 GB array limit;
/// the largest joint-popcount block is 4900² × 16 B ≈ 0.38 GB, comfortably in
/// commodity RAM. The block path is the only route to the N=8 spectrum within a
/// test budget.
///
/// <para>Four [Fact] methods, each marked
/// <c>[Trait("Category", "SLOW_N8")]</c> so they are excluded from default
/// <c>dotnet test</c> runs. Opt-in via:</para>
/// <code>dotnet test compute/RCPsiSquared.Core.Tests/RCPsiSquared.Core.Tests.csproj
///        --filter "Category=SLOW_N8" --no-build --nologo</code>
///
/// <para>Each test builds H = Heisenberg (XX+YY+ZZ) on the topology's bonds (J=1,
/// per-site γ = 0.5, so σ = N·γ = 4), computes the full Liouvillian spectrum
/// (65 536 eigenvalues across (N+1)² = 81 joint-popcount sectors) via
/// <see cref="LiouvillianBlockSpectrum.ComputeSpectrumPerBlock"/>, and asserts the
/// F1 palindromic-pairing identity {λ_k} = {−2σ − λ_k} as a multiset to tolerance
/// 1e-6 (relaxed from the N=7 dogfood's 1e-7 envelope to absorb MKL Evd accumulation
/// across the 81 block diagonalisations at sector dims up to 4900² for N=8 vs 64
/// blocks at 1225² for N=7).</para>
///
/// <para>Each test additionally runs <see cref="F1SpectrumStatistics.Compute"/> on the
/// per-block spectrum, logging the five metric groups (wall-time profile, pairing
/// precision histogram, spectrum-structure invariants, block-decomposition cost
/// picture, Hamiltonian + dissipator setup) through <see cref="ITestOutputHelper"/>
/// and persisting the metrics as JSON under
/// <c>simulations/results/f1_n8_n9_metrics/&lt;topology&gt;_N8.json</c>.</para>
///
/// <para>Heisenberg (popcount-conserving by XX+YY swap term + ZZ diagonal term)
/// sits inside the <see cref="JointPopcountSectorBuilder"/> block-infrastructure
/// domain (the F1-truly side of the F87 trichotomy), the same domain the N=7
/// dogfood uses. The companion DEBUG-only
/// <c>LiouvillianBlockSpectrum.DebugAssertPopcountConservingH</c> guard (added in
/// commit b07dd4d) enforces the contract structurally if a non-conserving H ever
/// slips in.</para>
///
/// <para>Anchors: <c>docs/proofs/PROOF_F1_GENERAL_TOPOLOGY.md</c> (synthesis,
/// verification-record table populated from the SLOW_N8 + SLOW_N9 JSON capture),
/// <see cref="F1GeneralTopologyVerifiedClaim"/> (typed Tier-2 claim,
/// <c>VerifiedNValues = {5, 6, 7, 8, 9}</c>; N=9 ran 2026-05-19 via the
/// <c>MklDirect</c> ILP64 bridge in
/// <c>F1GeneralTopologyN9BlockSpectrumChainTests</c>; new frontier at N=10 is
/// memory-pressure rather than the LP64 ceiling, see
/// <see cref="F1GeneralTopologyVerifiedClaim.ScaleFrontierBlockedAtN"/>),
/// <c>F1GeneralTopologyN7BlockSpectrumTests</c> (the N=7 dogfood whose pattern
/// this file mirrors at the next N), <see cref="F1SpectrumStatistics"/>
/// (the shared metrics utility).</para></summary>
public class F1GeneralTopologyN8BlockSpectrumTests
{
    private readonly ITestOutputHelper _out;

    public F1GeneralTopologyN8BlockSpectrumTests(ITestOutputHelper output) => _out = output;

    // ----------------------------------------------------------------------
    // Graph builders (mirror the N=7 file's static methods; identical shape
    // adapted to N=8 bond layout)
    // ----------------------------------------------------------------------

    private static Bond[] ChainBonds(int N) =>
        Enumerable.Range(0, N - 1).Select(i => new Bond(i, i + 1, 1.0)).ToArray();

    private static Bond[] RingBonds(int N) =>
        Enumerable.Range(0, N).Select(i => new Bond(i, (i + 1) % N, 1.0)).ToArray();

    private static Bond[] StarBonds(int N) =>
        Enumerable.Range(1, N - 1).Select(i => new Bond(0, i, 1.0)).ToArray();

    /// <summary>K_4 on sites {0,1,2,3} (6 bonds) plus disjoint 4-chain on
    /// {4,5,6,7} (3 bonds). B = 9. Two connected components, disconnected graph.
    /// Degree sequence on K_4: deg = 3 each (4·9 = 36); on the chain
    /// {4-5-6-7}: deg = (1, 2, 2, 1) → 1+4+4+1 = 10. D2(G) = 36 + 10 = 46.</summary>
    private static Bond[] K4PlusDisjoint4ChainBondsN8()
    {
        var bonds = new List<Bond>();
        for (int i = 0; i < 4; i++)
            for (int j = i + 1; j < 4; j++)
                bonds.Add(new Bond(i, j, 1.0));
        bonds.Add(new Bond(4, 5, 1.0));
        bonds.Add(new Bond(5, 6, 1.0));
        bonds.Add(new Bond(6, 7, 1.0));
        return bonds.ToArray();
    }

    // ----------------------------------------------------------------------
    // Test 1-4: N=8 F1 palindromic-pairing via the block infrastructure.
    // Heisenberg (XX+YY+ZZ) on the topology's bonds + uniform Z-dephasing
    // (γ=0.5). σ = N·γ = 4. Spectrum size = 4^8 = 65 536. Tolerance 1e-6.
    // Each test also captures the full F1SpectrumStatistics.TopologyMetrics
    // payload to simulations/results/f1_n8_n9_metrics/.
    // ----------------------------------------------------------------------

    [Fact]
    [Trait("Category", "SLOW_N8")]
    public void Chain_HeisenbergN8_F1PalindromicPairingViaBlockSpectrum()
    {
        VerifyF1PalindromicPairingAtN8(ChainBonds(8),
            topology: "chain (7 bonds)",
            jsonFileName: "chain_N8.json");
    }

    [Fact]
    [Trait("Category", "SLOW_N8")]
    public void Ring_HeisenbergN8_F1PalindromicPairingViaBlockSpectrum()
    {
        VerifyF1PalindromicPairingAtN8(RingBonds(8),
            topology: "ring (8 bonds, cycle 0-1-...-7-0)",
            jsonFileName: "ring_N8.json");
    }

    [Fact]
    [Trait("Category", "SLOW_N8")]
    public void Star_HeisenbergN8_F1PalindromicPairingViaBlockSpectrum()
    {
        VerifyF1PalindromicPairingAtN8(StarBonds(8),
            topology: "star (7 bonds, hub=0)",
            jsonFileName: "star_N8.json");
    }

    [Fact]
    [Trait("Category", "SLOW_N8")]
    public void K4PlusDisjoint4Chain_HeisenbergN8_F1PalindromicPairingViaBlockSpectrum()
    {
        VerifyF1PalindromicPairingAtN8(K4PlusDisjoint4ChainBondsN8(),
            topology: "K_4 on {0,1,2,3} + disjoint 4-chain on {4,5,6,7} (9 bonds, disconnected, two components)",
            jsonFileName: "k4_plus_disjoint_4chain_N8.json");
    }

    // ----------------------------------------------------------------------
    // Internal helpers
    // ----------------------------------------------------------------------

    /// <summary>Build H = (J/4) Σ_b (X_i X_j + Y_i Y_j + Z_i Z_j) on the given
    /// bond graph (extends <see cref="PauliHamiltonian.HeisenbergChain"/> from
    /// the chain to any graph topology). Popcount-conserving (XX+YY swap term +
    /// ZZ diagonal term), so it stays inside the
    /// <see cref="JointPopcountSectorBuilder"/> block-infrastructure domain.</summary>
    private static ComplexMatrix BuildHeisenbergGraphHamiltonian(int N, IReadOnlyList<Bond> bonds, double J = 1.0)
    {
        var terms = new (PauliLetter, PauliLetter, Complex)[]
        {
            (PauliLetter.X, PauliLetter.X, J / 4.0),
            (PauliLetter.Y, PauliLetter.Y, J / 4.0),
            (PauliLetter.Z, PauliLetter.Z, J / 4.0),
        };
        return PauliHamiltonian.Bilinear(N, bonds, terms).ToMatrix();
    }

    /// <summary>Compute the full Liouvillian spectrum at N=8 via
    /// <see cref="LiouvillianBlockSpectrum.ComputeSpectrumPerBlock"/> on Heisenberg
    /// H over the given bonds + uniform Z-dephasing (γ=0.5), and assert the F1
    /// palindromic pairing {λ_k} = {−2σ − λ_k} as a multiset (σ = N·γ = 4 at
    /// N=8, γ=0.5).
    ///
    /// <para>Records the full <see cref="F1SpectrumStatistics.TopologyMetrics"/>
    /// payload (wall-time profile, pairing-precision histogram, spectrum-structure
    /// invariants, block-decomposition cost picture, Hamiltonian + dissipator
    /// setup), logs each metric through <see cref="ITestOutputHelper"/> for
    /// grep-friendly inspection, and persists the JSON under
    /// <c>simulations/results/f1_n8_n9_metrics/&lt;jsonFileName&gt;</c> for the
    /// proof-document verification table to pull from.</para>
    ///
    /// <para>Tolerance 1e-6: empirically the N=7 dogfood budget 1e-7 sits at the edge
    /// for N=8 (max distance observed 1.184e-7 on K_4 + disjoint-4-chain at the first
    /// run, so a ~10× safety factor over the N=7 envelope absorbs the accumulation
    /// from 81 block diagonalisations at sector dims up to 4900² for N=8 vs 64 blocks
    /// at 1225² for N=7). The relaxed tolerance is &lt; 1 part in 10⁶ of σ = 4 and
    /// stays orders of magnitude tighter than any physical-spectrum gap.</para></summary>
    private void VerifyF1PalindromicPairingAtN8(Bond[] bonds, string topology, string jsonFileName)
    {
        const int N = 8;
        const double J = 1.0;
        const double Gamma = 0.5;
        const double Tolerance = 1e-6;
        double sigma = N * Gamma;  // = 4.0

        var totalSw = Stopwatch.StartNew();

        var H = BuildHeisenbergGraphHamiltonian(N, bonds, J: J);
        var gammaPerSite = Enumerable.Repeat(Gamma, N).ToArray();

        // Group 1 timing: ComputeSpectrumPerBlock wall time isolated from scaffolding.
        var computeSw = Stopwatch.StartNew();
        var spectrum = LiouvillianBlockSpectrum.ComputeSpectrumPerBlock(H, gammaPerSite, N);
        computeSw.Stop();

        // Sanity: 4^8 = 65 536 eigenvalues across (N+1)² = 81 sectors, max block 4900².
        Assert.Equal(1 << (2 * N), spectrum.Length);

        // Multiset equality with nearest-neighbour matching. This drives the pass/fail
        // and is the F1 palindromic-pairing dogfood the test exists to verify.
        var shifted = new Complex[spectrum.Length];
        for (int i = 0; i < spectrum.Length; i++)
            shifted[i] = -2.0 * sigma - spectrum[i];
        MultisetAssert.NearestNeighbourEqual(spectrum, shifted, tolerance: Tolerance,
            context: $"{topology} N=8 Heisenberg F1 palindromic pairing");

        totalSw.Stop();

        // F1SpectrumStatistics: capture every metric group from the same spectrum.
        var bondPairs = bonds.Select(b => (b.Site1, b.Site2)).ToArray();
        var metrics = F1SpectrumStatistics.Compute(
            N: N,
            topologyName: topology,
            hamiltonianClass: "Heisenberg XXX (XX+YY+ZZ) at J=1.0, uniform Z-dephasing γ=0.5",
            jValue: J,
            gammaValue: Gamma,
            bonds: bondPairs,
            spectrum: spectrum,
            totalWallSeconds: totalSw.Elapsed.TotalSeconds,
            computeSpectrumWallSeconds: computeSw.Elapsed.TotalSeconds);

        // Grep-friendly metric log (one line per quantity, "metric=value"); shared
        // helper in F1SpectrumStatistics is reused by the N=9 chain test class.
        F1SpectrumStatistics.LogMetrics(metrics, Tolerance, _out.WriteLine);

        // Persist JSON to simulations/results/f1_n8_n9_metrics/<jsonFileName>.
        string outDir = F1SpectrumStatistics.ResolveMetricsDirectory();
        F1SpectrumStatistics.WriteJson(metrics, outDir, jsonFileName);
        _out.WriteLine($"metrics_json_path={Path.Combine(outDir, jsonFileName)}");
    }
}
