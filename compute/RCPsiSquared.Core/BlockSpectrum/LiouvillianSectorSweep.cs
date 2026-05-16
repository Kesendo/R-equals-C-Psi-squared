using System.Collections.Concurrent;
using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Pauli;
using RCPsiSquared.Core.SymmetryFamily;
using ComplexMatrix = MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex>;

namespace RCPsiSquared.Core.BlockSpectrum;

/// <summary>Per-(p_c, p_r) sector dense Evd of the chain XY + Z-dephasing Liouvillian,
/// capped at a user-supplied <c>sectorDimCap</c>. Composes
/// <see cref="JointPopcountSectorBuilder"/> (sector enumeration) +
/// <see cref="PerBlockLiouvillianBuilder.BuildBlockZ"/> (per-sector L construction in the
/// computational basis) + <see cref="XGlobalChargeConjugationPairing"/> (X⊗N pairing to
/// halve eigendecomposition work) + dense MathNet <c>Evd</c>, exactly like
/// <see cref="LiouvillianBlockSpectrum.ComputeSpectrumPerBlock"/> — but skips any sector
/// whose dimension exceeds <see cref="SectorDimCap"/>.
///
/// <para>Use case: at N where the largest joint-popcount sector no longer fits dense
/// LAPACK on commodity hardware (N≥9 max sector at half-filling), this primitive yields
/// the union of all sector spectra that DO fit. Combined with the F1 palindrome witness
/// (every λ should pair with −λ − 2Σγ inside the collected set), this gives an honest
/// partial-but-substantial spectrum plus an analytical-symmetry verification on the
/// portion we can compute.</para>
///
/// <para>At small N where the cap covers every sector (<c>sectorDimCap ≥ C(N, ⌊N/2⌋)²</c>),
/// the collected eigenvalue set equals the full 4^N spectrum bit-for-bit (verified at
/// N=3, 4 vs <see cref="Lindblad.PauliDephasingDissipator.BuildZ"/> + direct Evd).</para>
///
/// <para><b>F1 witness expectation — what the Brecher mechanism implies.</b> The witness
/// is expected to land at machine precision (~1e-13 at N=10, γ=0.05) for this primitive's
/// chain XY + Z-dephasing setup: the F-trichotomy classifies XY as "truly", and Z-dephasing
/// is the F1-preserving dephasing axis, so there is no Brecher in play. Conditions that
/// DO break F1 ("Brecher" mechanisms, enumerated in
/// <see cref="F1.F1OpenQuestions"/> and witnessed at small N by
/// <c>PalindromeResidualTests.F1_Palindrome_BreaksFor_T1Dissipator</c>):</para>
/// <list type="bullet">
///   <item><b>T1 amplitude damping</b> (Lindblad operator σ⁻ = (X − iY)/2 carries a Y
///         component, bit_b = 1) — F1 residual jumps from FP-noise to <c>O(γ_T1)</c>.</item>
///   <item><b>Depolarising noise</b> — F1 breaks with residual scaling (2/3)Σγ, linear
///         in γ and N (closed form is an F1OpenQuestions item).</item>
///   <item><b>Transverse-field Hamiltonians</b> h_x·X or h_y·Y at the Hamiltonian level —
///         the Z⊗N-Brecher of <c>hypotheses/THE_POLARITY_LAYER.md</c>, takes the system
///         out of the "truly" class.</item>
/// </list>
/// <para>A residual visibly above ~<c>10⁻¹⁰ · N · max(γ_l)</c> in this primitive's witness
/// therefore indicates either an arithmetic bug in the per-block eig path or an unintended
/// Brecher-introducing change to the L source (e.g. swapping in T1 instead of Z-dephasing).</para>
///
/// <para><b>Tier outcome: <see cref="Tier.Tier1Derived"/>.</b> Pure composition of Tier1
/// primitives (<see cref="JointPopcountSectorBuilder"/>, <see cref="PerBlockLiouvillianBuilder"/>,
/// <see cref="XGlobalChargeConjugationPairing"/>, dense LAPACK Evd) plus an arithmetic F1
/// palindrome residual witness.</para>
///
/// <para>Anchor: <see cref="LiouvillianBlockSpectrum"/> (parent full-coverage Claim);
/// <see cref="F1.F1PalindromeIdentity"/> (the palindrome whose residual is witnessed);
/// <see cref="F1.F1OpenQuestions"/> (enumerated F1-break mechanisms);
/// <c>docs/proofs/MIRROR_SYMMETRY_PROOF.md</c>.</para>
/// </summary>
public sealed class LiouvillianSectorSweep : Claim
{
    public int N { get; }
    public double J { get; }
    public IReadOnlyList<double> GammaPerSite { get; }
    public int SectorDimCap { get; }

    /// <summary>Theoretical full Liouville-space dimension <c>4^N</c>; the long type
    /// guards against the int32 overflow at N=16.</summary>
    public long FullDim { get; }

    public IReadOnlyList<SweepSector> Included { get; }
    public IReadOnlyList<SweepSector> Skipped { get; }
    public Complex[] CollectedEigenvalues { get; }
    public double CoverageFraction => (double)CollectedEigenvalues.LongLength / FullDim;

    /// <summary>Max over collected λ of min distance to its F1 mirror −λ − 2Σγ in the
    /// collected set. At full coverage under uniform γ this is at FP-noise level
    /// (1e-10..1e-13 typically); a substantially larger value either signals F1 violation
    /// (impossible — F1 is Tier1Derived) or that the collected set is not closed under
    /// the F1 involution (some partner sectors were skipped).</summary>
    public double F1PalindromeResidualMax { get; }

    public static LiouvillianSectorSweep Build(int N, IReadOnlyList<double> gammaPerSite,
        int sectorDimCap, double J = 1.0)
    {
        if (N < 1) throw new ArgumentOutOfRangeException(nameof(N), N, "N must be ≥ 1.");
        if (gammaPerSite is null) throw new ArgumentNullException(nameof(gammaPerSite));
        if (gammaPerSite.Count != N)
            throw new ArgumentException(
                $"gammaPerSite length {gammaPerSite.Count} != N {N}", nameof(gammaPerSite));
        if (sectorDimCap < 1)
            throw new ArgumentOutOfRangeException(nameof(sectorDimCap), sectorDimCap, "sectorDimCap must be ≥ 1.");

        var H = PauliHamiltonian.XYChain(N, J).ToMatrix();
        var decomp = JointPopcountSectorBuilder.Build(N);
        int sectorCount = decomp.SectorRanges.Count;

        var included = new List<SweepSector>();
        var skipped = new List<SweepSector>();
        for (int i = 0; i < sectorCount; i++)
        {
            var sec = decomp.SectorRanges[i];
            var sweepSec = new SweepSector(sec.PCol, sec.PRow, sec.Size);
            if (sec.Size <= sectorDimCap) included.Add(sweepSec);
            else skipped.Add(sweepSec);
        }

        Complex[] collectedArr;
        if (skipped.Count == 0)
        {
            // Battle-tested full-coverage path: delegate to the parent primitive.
            collectedArr = LiouvillianBlockSpectrum.ComputeSpectrumPerBlock(H, gammaPerSite, N);
        }
        else
        {
            // Capped path: enumerate primaries that fit, Evd each, collect via X⊗N follower
            // copies for sectors that also fit. Sectors above the cap are dropped entirely.
            var perm = decomp.Permutation;
            var (primarySectorIndices, followerToPrimary) =
                XGlobalChargeConjugationPairing.PartitionByXNPairing(
                    N, decomp.SectorRanges, s => (s.PCol, s.PRow), s => s.Size);

            var includedPrimaries = new List<int>();
            foreach (var pIdx in primarySectorIndices)
                if (decomp.SectorRanges[pIdx].Size <= sectorDimCap)
                    includedPrimaries.Add(pIdx);

            var primaryEigs = new ConcurrentDictionary<int, Complex[]>();
            int outerDop = Math.Max(1, Environment.ProcessorCount / 4);
            var po = new ParallelOptions { MaxDegreeOfParallelism = outerDop };
            Parallel.ForEach(includedPrimaries, po, sIdx =>
            {
                var sector = decomp.SectorRanges[sIdx];
                int size = sector.Size;
                if (size == 0) { primaryEigs[sIdx] = Array.Empty<Complex>(); return; }
                var flatIndices = new int[size];
                for (int k = 0; k < size; k++) flatIndices[k] = perm[sector.Offset + k];
                var block = PerBlockLiouvillianBuilder.BuildBlockZ(H, gammaPerSite, flatIndices);
                var blockEigs = block.Evd().EigenValues;
                var arr = new Complex[size];
                for (int i = 0; i < size; i++) arr[i] = blockEigs[i];
                primaryEigs[sIdx] = arr;
            });

            var collected = new List<Complex>(capacity: included.Sum(s => s.Dim));
            for (int sIdx = 0; sIdx < sectorCount; sIdx++)
            {
                var sec = decomp.SectorRanges[sIdx];
                if (sec.Size > sectorDimCap) continue;
                int sourceIdx = primaryEigs.ContainsKey(sIdx) ? sIdx : followerToPrimary[sIdx];
                if (!primaryEigs.TryGetValue(sourceIdx, out var eigs)) continue;
                collected.AddRange(eigs);
            }
            collectedArr = collected.ToArray();
        }

        double sumGamma = gammaPerSite.Sum();
        double palindromeRes = ComputeF1PalindromeResidual(collectedArr, sumGamma);

        return new LiouvillianSectorSweep(N, J, gammaPerSite.ToArray(), sectorDimCap,
            FullDimOf(N), included, skipped, collectedArr, palindromeRes);
    }

    private static long FullDimOf(int N) => 1L << (2 * N);

    /// <summary>F1 palindrome residual via greedy nearest-neighbour matching: for each λ in
    /// the collected set, find its F1 mirror <c>−λ − 2Σγ</c>'s closest match in the set (not
    /// yet taken by an earlier match). Returns the max distance across all matches. Highly
    /// degenerate Liouvillian spectra can fool a sort-and-zip scan; greedy matching with a
    /// taken-mask is bit-robust at the cost of O(K²) — acceptable for K ≤ a few thousand,
    /// slow at K = 70 k (N=10 with partial coverage). For larger K, this routine remains
    /// correct but becomes the wall-time bottleneck.</summary>
    private static double ComputeF1PalindromeResidual(Complex[] eigs, double sumGamma)
    {
        int K = eigs.Length;
        if (K == 0) return 0.0;
        var taken = new bool[K];
        double maxResidual = 0.0;
        for (int i = 0; i < K; i++)
        {
            var mirror = new Complex(-eigs[i].Real - 2 * sumGamma, -eigs[i].Imaginary);
            int bestJ = -1;
            double bestDist = double.MaxValue;
            for (int j = 0; j < K; j++)
            {
                if (taken[j]) continue;
                double d = (eigs[j] - mirror).Magnitude;
                if (d < bestDist) { bestDist = d; bestJ = j; }
            }
            if (bestJ >= 0) taken[bestJ] = true;
            if (bestDist > maxResidual) maxResidual = bestDist;
        }
        return maxResidual;
    }

    private LiouvillianSectorSweep(int n, double j, double[] gammaPerSite, int sectorDimCap,
        long fullDim,
        IReadOnlyList<SweepSector> included, IReadOnlyList<SweepSector> skipped,
        Complex[] collectedEigenvalues, double f1PalindromeResidualMax)
        : base($"Sector sweep at N={n}, dimCap={sectorDimCap}: {included.Count} sectors included " +
               $"({collectedEigenvalues.LongLength}/{fullDim} eigenvalues, " +
               $"coverage {(double)collectedEigenvalues.LongLength / fullDim:P2}); " +
               $"F1 palindrome residual = {f1PalindromeResidualMax:G3}.",
               Tier.Tier1Derived,
               "compute/RCPsiSquared.Core/BlockSpectrum/LiouvillianBlockSpectrum.cs (parent full-coverage Claim) + " +
               "compute/RCPsiSquared.Core/BlockSpectrum/JointPopcountSectorBuilder.cs (sector enumeration) + " +
               "compute/RCPsiSquared.Core/BlockSpectrum/PerBlockLiouvillianBuilder.cs (computational-basis L block) + " +
               "compute/RCPsiSquared.Core/SymmetryFamily/XGlobalChargeConjugationPairing.cs (X⊗N halving) + " +
               "compute/RCPsiSquared.Core/F1/F1PalindromeIdentity.cs (witnessed palindrome).")
    {
        N = n; J = j;
        GammaPerSite = gammaPerSite;
        SectorDimCap = sectorDimCap;
        FullDim = fullDim;
        Included = included;
        Skipped = skipped;
        CollectedEigenvalues = collectedEigenvalues;
        F1PalindromeResidualMax = f1PalindromeResidualMax;
    }

    public override string DisplayName =>
        $"Sector sweep (N={N}, dimCap={SectorDimCap}): {Included.Count} sectors, " +
        $"{CollectedEigenvalues.LongLength:N0}/{FullDim:N0} eigenvalues";

    public override string Summary =>
        $"coverage={CoverageFraction:P2}, F1 residual={F1PalindromeResidualMax:G3} ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return InspectableNode.RealScalar("N", N);
            yield return InspectableNode.RealScalar("J", J, "G4");
            yield return InspectableNode.RealScalar("sectorDimCap", SectorDimCap);
            yield return InspectableNode.RealScalar("included sectors", Included.Count);
            yield return InspectableNode.RealScalar("skipped sectors", Skipped.Count);
            yield return InspectableNode.RealScalar("collected eigenvalues", CollectedEigenvalues.LongLength);
            yield return InspectableNode.RealScalar("full dimension 4^N", FullDim);
            yield return InspectableNode.RealScalar("coverage fraction", CoverageFraction, "P3");
            yield return InspectableNode.RealScalar("F1 palindrome residual max", F1PalindromeResidualMax, "G3");
        }
    }
}

/// <summary>One (p_c, p_r) sector descriptor for <see cref="LiouvillianSectorSweep"/>:
/// the joint-popcount labels plus the sector's eigendecomposition dimension
/// <c>C(N, p_c) · C(N, p_r)</c>.</summary>
public sealed record SweepSector(int PCol, int PRow, int Dim);
