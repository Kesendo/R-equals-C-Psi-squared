using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.BlockSpectrum;
using RCPsiSquared.Core.ChainSystems;
using RCPsiSquared.Core.F1;
using RCPsiSquared.Core.Lindblad;
using RCPsiSquared.Core.Pauli;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Core.Tests.TestHelpers;
using ComplexMatrix = MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex>;

namespace RCPsiSquared.Core.Tests.F1;

/// <summary>F1 general-topology verification: C# half of the closure for the last F1
/// OpenQuestion ("general topology beyond chain/ring/star/K_N").
///
/// <para>The Python half lives in <c>simulations/f1_general_topology_verify.py</c>
/// and covers N=5, 6 across named graphs, random connected Erdős-Rényi graphs,
/// disconnected components, weighted edges, and the single-body class. This C# class
/// adds:</para>
///
/// <list type="bullet">
///   <item><b>Graph-aware C# verification at N=5</b> using
///         <see cref="PalindromeResidualScalingClaim"/> in its graph-aware mode
///         (<c>bondCount</c> + <c>degreeSquaredSum</c> arguments). Observed ‖M‖² from
///         <see cref="PalindromeResidual.Build"/> on chain / ring / star / triangle +
///         disjoint bond, matched against the typed claim's <c>Factor</c>.</item>
///   <item><b>F1 palindromic-pairing verification at N=7</b> via
///         <see cref="LiouvillianBlockSpectrum.ComputeSpectrumPerBlock"/> on its proper
///         domain (XY + Z-dephasing, F1-truly). For each Liouvillian eigenvalue λ_k the
///         partner λ_k + 2σ + λ_pair = 0 must hold, i.e., {λ} = {−2σ − λ}. This is the
///         spectral content of the F1 identity Π·L·Π⁻¹ = −L − 2σ·I; the test asserts it
///         as a multiset identity across the full 4^7 = 16 384 eigenvalues. Multiple
///         graph topologies (chain, ring, star, disconnected) are exercised; the only
///         path that scales to N=7 in a test budget (full L_vec at N=7 is 16384² = 4 GB,
///         past the .NET 2 GB array limit).</item>
///   <item><b>N=5 self-test</b> that the block infrastructure is loss-free in the
///         Frobenius-norm sense: Σ_blocks ‖L_b‖²_F = ‖L_full‖²_F to machine precision.
///         Confirms the per-block path's bookkeeping is consistent with the dense
///         baseline.</item>
/// </list>
///
/// <para><b>Infrastructure note on Hamiltonian choice.</b> The
/// <see cref="JointPopcountSectorBuilder"/> block decomposition is valid only for
/// Hamiltonians that preserve the joint popcount label, i.e., XY + Z-dephasing. Non-
/// truly Hamiltonians (XX+YZ, XY+YX, etc.) generally break this label and are outside
/// the block infrastructure's scope. So the F1 scaling-formula tests on non-truly H
/// are done with full <see cref="PalindromeResidual.Build"/> at N=5 (graph-aware
/// mode), and the N=7 block path is used for the F1 palindromic-pairing identity
/// (which IS testable on the block infrastructure's proper truly-H domain).</para>
///
/// <para>Anchor: <c>docs/proofs/PROOF_F1_GENERAL_TOPOLOGY.md</c> (synthesis with
/// citation to PROOF_CROSS_TERM_FORMULA Lemma 3 Corollary for the universality of the
/// (B, D2) parameterisation).</para></summary>
public class F1GeneralTopologyN7BlockSpectrumTests
{
    // ----------------------------------------------------------------------
    // c_H anchor for XX+YZ at N=2 (matches PalindromeResidualScalingClaim.Verify
    // and simulations/f1_general_topology_verify.py section 1). Static-cached so
    // the ~14 s build cost is paid once across all N=5 graph-aware tests.
    // ----------------------------------------------------------------------

    private static readonly double _cHAnchorXxYz = ComputeCHAnchorXxYzAtN2();

    private static double ComputeCHAnchorXxYzAtN2()
    {
        const int N = 2;
        const double gammaZ = 0.1;
        var bonds = new[] { new Bond(0, 1, 1.0) };
        var H = BuildXxYzGraphHamiltonian(N, bonds);
        var gammaPerSite = Enumerable.Repeat(gammaZ, N).ToArray();
        var L = PauliDephasingDissipator.BuildZ(H, gammaPerSite);
        var M = PalindromeResidual.Build(L, N, N * gammaZ, PauliLetter.Z);
        double frob = M.FrobeniusNorm();
        return frob * frob;
    }

    // ----------------------------------------------------------------------
    // Graph builders
    // ----------------------------------------------------------------------

    private static Bond[] ChainBonds(int N) =>
        Enumerable.Range(0, N - 1).Select(i => new Bond(i, i + 1, 1.0)).ToArray();

    private static Bond[] RingBonds(int N) =>
        Enumerable.Range(0, N).Select(i => new Bond(i, (i + 1) % N, 1.0)).ToArray();

    private static Bond[] StarBonds(int N) =>
        Enumerable.Range(1, N - 1).Select(i => new Bond(0, i, 1.0)).ToArray();

    /// <summary>Two-component disconnected graph at N=5: triangle (3-cycle) on {0,1,2}
    /// plus a single bond on {3,4}. Bonds = (0,1), (1,2), (0,2), (3,4); B = 4.
    /// Degree sequence: deg(0)=deg(1)=deg(2)=2, deg(3)=deg(4)=1. D2 = 3·4 + 2·1 = 14.</summary>
    private static Bond[] TriangleAndBondBondsN5() => new[]
    {
        new Bond(0, 1, 1.0),
        new Bond(1, 2, 1.0),
        new Bond(0, 2, 1.0),
        new Bond(3, 4, 1.0),
    };

    /// <summary>K_4 on sites {0,1,2,3} plus disjoint path (4,5)−(5,6) on {4,5,6}.
    /// B = 8 (6 from K_4 + 2 from path). D2 = 4·9 + 1 + 4 + 1 = 42.</summary>
    private static Bond[] K4PlusIsolatedPathBondsN7()
    {
        var bonds = new List<Bond>();
        for (int i = 0; i < 4; i++)
            for (int j = i + 1; j < 4; j++)
                bonds.Add(new Bond(i, j, 1.0));
        bonds.Add(new Bond(4, 5, 1.0));
        bonds.Add(new Bond(5, 6, 1.0));
        return bonds.ToArray();
    }

    // ----------------------------------------------------------------------
    // Test 1-4: N=5 graph-aware scaling (non-truly H = XX+YZ), full
    // PalindromeResidual.Build path
    // ----------------------------------------------------------------------

    [Fact]
    public void Chain_XxYzN5_PalindromeResidualMatchesGraphAwareFactor()
    {
        const int N = 5;
        VerifyGraphAwareScalingAtN5(ChainBonds(N), expectedB: 4, expectedD2: 14, topology: "chain");
    }

    [Fact]
    public void Ring_XxYzN5_PalindromeResidualMatchesGraphAwareFactor()
    {
        const int N = 5;
        // Ring: every site degree 2 → D2 = N·4 = 20.
        VerifyGraphAwareScalingAtN5(RingBonds(N), expectedB: 5, expectedD2: 20, topology: "ring");
    }

    [Fact]
    public void Star_XxYzN5_PalindromeResidualMatchesGraphAwareFactor()
    {
        const int N = 5;
        // Star: hub has degree 4, leaves have degree 1 each → D2 = 16 + 4·1 = 20.
        VerifyGraphAwareScalingAtN5(StarBonds(N), expectedB: 4, expectedD2: 20, topology: "star");
    }

    [Fact]
    public void DisconnectedTriangleAndBond_XxYzN5_PalindromeResidualMatchesGraphAwareFactor()
    {
        // Disconnected: 3-cycle + 1-bond → B = 4, D2 = 14.
        VerifyGraphAwareScalingAtN5(TriangleAndBondBondsN5(), expectedB: 4, expectedD2: 14,
            topology: "triangle + disjoint bond (disconnected)");
    }

    // ----------------------------------------------------------------------
    // Test 5-8: N=7 F1 palindromic-pairing verification via the block
    // infrastructure (XY + Z-dephasing, F1-truly). Computes the full spectrum
    // via LiouvillianBlockSpectrum.ComputeSpectrumPerBlock and asserts that the
    // multiset {λ_k} matches the F1-shifted multiset {−2σ − λ_k}.
    // ----------------------------------------------------------------------

    [Fact]
    public void Chain_XyZdephN7_F1PalindromicPairingViaBlockSpectrum()
    {
        VerifyF1PalindromicPairingViaBlockSpectrumAtN7(ChainBonds(7), topology: "chain");
    }

    [Fact]
    public void Ring_XyZdephN7_F1PalindromicPairingViaBlockSpectrum()
    {
        VerifyF1PalindromicPairingViaBlockSpectrumAtN7(RingBonds(7), topology: "ring");
    }

    [Fact]
    public void Star_XyZdephN7_F1PalindromicPairingViaBlockSpectrum()
    {
        VerifyF1PalindromicPairingViaBlockSpectrumAtN7(StarBonds(7), topology: "star");
    }

    [Fact]
    public void K4PlusIsolatedPath_XyZdephN7_F1PalindromicPairingViaBlockSpectrum()
    {
        VerifyF1PalindromicPairingViaBlockSpectrumAtN7(K4PlusIsolatedPathBondsN7(),
            topology: "K_4 + disjoint 3-chain (disconnected, two components)");
    }

    // ----------------------------------------------------------------------
    // Test 9: dogfood sanity at N=5, block-decomposed sum equals dense sum
    // for ‖L‖²_F on the infrastructure's proper domain (XY + Z-deph).
    // ----------------------------------------------------------------------

    [Fact]
    public void InfrastructureSelfTest_BlockSumEqualsFullDenseSum_XyZdephN5()
    {
        const int N = 5;
        const double gammaZ = 0.1;

        var H = PauliHamiltonian.XYChain(N, J: 1.0).ToMatrix();
        var gammaPerSite = Enumerable.Repeat(gammaZ, N).ToArray();
        var L = PauliDephasingDissipator.BuildZ(H, gammaPerSite);

        double frobFull = L.FrobeniusNorm();
        double normSqFull = frobFull * frobFull;

        var (normSqBlock, _) = ComputeLNormSquaredAndTraceFromBlocks(H, gammaPerSite, N);

        double rel = Math.Abs(normSqFull - normSqBlock) / normSqFull;
        Assert.True(rel < 1e-12,
            $"InfrastructureSelfTest at N=5 XY+Z-deph: block ‖L‖²_F sum {normSqBlock} " +
            $"disagrees with full-dense {normSqFull} (rel = {rel:E3})");
    }

    // ----------------------------------------------------------------------
    // Internal helpers
    // ----------------------------------------------------------------------

    /// <summary>Build H = Σ_b (X_i X_j + Y_i Z_j) on the given bond set, dense
    /// 2^N × 2^N. The canonical non-truly anchor used by
    /// <see cref="PalindromeResidualScalingClaim.Verify"/>.</summary>
    private static ComplexMatrix BuildXxYzGraphHamiltonian(int N, IReadOnlyList<Bond> bonds)
    {
        var terms = new (PauliLetter, PauliLetter, Complex)[]
        {
            (PauliLetter.X, PauliLetter.X, Complex.One),
            (PauliLetter.Y, PauliLetter.Z, Complex.One),
        };
        return PauliHamiltonian.Bilinear(N, bonds, terms).ToMatrix();
    }

    /// <summary>Build H = (J/2)·Σ_b (X_i X_j + Y_i Y_j) on an arbitrary bond graph
    /// (extends <see cref="PauliHamiltonian.XYChain"/> from the chain to any graph
    /// topology). Preserves the joint popcount label and so stays inside the
    /// <see cref="JointPopcountSectorBuilder"/> domain.</summary>
    private static ComplexMatrix BuildXyGraphHamiltonian(int N, IReadOnlyList<Bond> bonds, double J = 1.0)
    {
        var terms = new (PauliLetter, PauliLetter, Complex)[]
        {
            (PauliLetter.X, PauliLetter.X, J / 2.0),
            (PauliLetter.Y, PauliLetter.Y, J / 2.0),
        };
        return PauliHamiltonian.Bilinear(N, bonds, terms).ToMatrix();
    }

    /// <summary>Block-by-block accumulation of ‖L‖²_F and tr(L) via
    /// <see cref="PerBlockLiouvillianBuilder.BuildBlockZ"/>. No full-L
    /// materialisation; this is the dogfood path that scales past N=6 for
    /// XY+Z-deph Hamiltonians.</summary>
    private static (double NormSquared, Complex Trace) ComputeLNormSquaredAndTraceFromBlocks(
        ComplexMatrix H, IReadOnlyList<double> gammaPerSite, int N)
    {
        var decomp = JointPopcountSectorBuilder.Build(N);
        var perm = decomp.Permutation;

        double normSqAccum = 0.0;
        Complex traceAccum = Complex.Zero;
        foreach (var sector in decomp.SectorRanges)
        {
            int size = sector.Size;
            if (size == 0) continue;

            var flatIndices = new int[size];
            for (int k = 0; k < size; k++)
                flatIndices[k] = perm[sector.Offset + k];

            var block = PerBlockLiouvillianBuilder.BuildBlockZ(H, gammaPerSite, flatIndices);
            double frob = block.FrobeniusNorm();
            normSqAccum += frob * frob;
            traceAccum += block.Trace();
        }
        return (normSqAccum, traceAccum);
    }

    /// <summary>Build H on the given bond graph (XX+YZ non-truly), compute observed
    /// ‖M‖²_F via dense <see cref="PalindromeResidual.Build"/>, and assert it matches
    /// the graph-aware <see cref="PalindromeResidualScalingClaim.Factor"/> · c_H
    /// within tight relative tolerance.</summary>
    private static void VerifyGraphAwareScalingAtN5(Bond[] bonds, int expectedB,
                                                     int expectedD2, string topology)
    {
        const int N = 5;
        const double gammaZ = 0.1;

        // Validate graph invariants (catch bond-builder typos).
        Assert.Equal(expectedB, bonds.Length);
        var deg = new int[N];
        foreach (var b in bonds) { deg[b.Site1]++; deg[b.Site2]++; }
        int observedD2 = deg.Sum(d => d * d);
        Assert.Equal(expectedD2, observedD2);

        // Build L, M, observed ‖M‖²_F (dense path; N=5 → 1024 × 1024, ~16 MB).
        var H = BuildXxYzGraphHamiltonian(N, bonds);
        var gammaPerSite = Enumerable.Repeat(gammaZ, N).ToArray();
        var L = PauliDephasingDissipator.BuildZ(H, gammaPerSite);
        var M = PalindromeResidual.Build(L, N, N * gammaZ, PauliLetter.Z);
        double frob = M.FrobeniusNorm();
        double observed = frob * frob;

        // Predicted: c_H · B · 4^(N−2) via graph-aware claim Factor.
        var claim = new PalindromeResidualScalingClaim(N, HamiltonianClass.Main,
            bondCount: expectedB, degreeSquaredSum: expectedD2);
        double predicted = _cHAnchorXxYz * claim.Factor;

        double rel = Math.Abs(observed - predicted) / predicted;
        // 1e-9 absorbs the floating-point accumulation across the 4^N basis transform.
        Assert.True(rel < 1e-9,
            $"{topology} N=5: observed ‖M‖² = {observed}, predicted c_H·F = " +
            $"{_cHAnchorXxYz}·{claim.Factor} = {predicted} (rel = {rel:E3}). " +
            $"B = {expectedB}, D2 = {observedD2}.");
    }

    /// <summary>Compute the full L spectrum at N=7 via
    /// <see cref="LiouvillianBlockSpectrum.ComputeSpectrumPerBlock"/> and assert the
    /// F1 palindromic pairing: every eigenvalue λ has a partner λ' such that
    /// λ + λ' = −2σ where σ = N·γ. This is the spectral content of the F1 identity
    /// Π·L·Π⁻¹ = −L − 2σ·I.
    ///
    /// <para>Equivalent multiset statement: {λ_k} = {−2σ − λ_k}. We check this via
    /// <see cref="MultisetAssert.NearestNeighbourEqual"/> on the original spectrum
    /// vs its F1-shifted image.</para></summary>
    private static void VerifyF1PalindromicPairingViaBlockSpectrumAtN7(Bond[] bonds, string topology)
    {
        const int N = 7;
        const double gammaZ = 0.1;
        double sigma = N * gammaZ;

        var H = BuildXyGraphHamiltonian(N, bonds);
        var gammaPerSite = Enumerable.Repeat(gammaZ, N).ToArray();

        // ComputeSpectrumPerBlock: per-block eig, no full-L materialisation.
        var spectrum = LiouvillianBlockSpectrum.ComputeSpectrumPerBlock(H, gammaPerSite, N);

        // Sanity: 4^7 = 16 384 eigenvalues across (N+1)² = 64 sectors.
        Assert.Equal(1 << (2 * N), spectrum.Length);

        // F1 palindrome: {λ_k} = {−2σ − λ_k} as a multiset.
        var shifted = new Complex[spectrum.Length];
        for (int i = 0; i < spectrum.Length; i++)
            shifted[i] = -2.0 * sigma - spectrum[i];

        // Multiset equality with nearest-neighbour matching. Tolerance 1e-7 accommodates
        // accumulated MKL Evd noise across (N+1)² = 64 block diagonalisations at sector
        // dims up to 1225² at N=7. Direct dense eig at N=5 typically agrees to 1e-9.
        MultisetAssert.NearestNeighbourEqual(spectrum, shifted, tolerance: 1e-7,
            context: $"{topology} N=7 F1 palindromic pairing");
    }
}
