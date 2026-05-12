using System.Collections.Concurrent;
using System.Diagnostics;
using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.BlockSpectrum;
using RCPsiSquared.Core.Lindblad;
using RCPsiSquared.Core.Pauli;
using RCPsiSquared.Core.SymmetryFamily;
using ComplexMatrix = MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex>;

namespace RCPsiSquared.Cli.Commands;

/// <summary>CLI demonstrator for the BlockSpectrum infrastructure.
///
/// <para>Builds a chain XY+Z-dephasing Liouvillian L at user-specified N, prints sector
/// summary, computes the spectrum either via the full L or via per-block eigendecomposition
/// over the joint-popcount sectors (optionally further refined by the F71 spatial-mirror
/// Z₂), and reports timing + a representative slice of the spectrum.</para>
///
/// <para>Usage: <c>rcpsi block-spectrum --N 6 [--gamma 0.5] [--J 1.0] [--refine f71|none]
/// [--verify | --no-verify] [--top K]</c>.</para>
///
/// <para>The <c>--verify</c> flag additionally diagonalises the full L directly and asserts
/// that the per-block spectrum matches as a multiset to within 1e-9. Default is no-verify
/// (full-L OOMs at N=8). The <c>--no-verify</c> flag is accepted for explicitness but is
/// equivalent to omitting <c>--verify</c>.</para>
///
/// <para>The <c>--top K</c> flag (default 20) controls how many leading eigenvalues are
/// printed (sorted by Re ascending).</para>
///
/// <para>Per-sector timing summary (mean, max, total) is emitted when <c>N ≥ 6</c>, along
/// with peak managed-memory delta over the block-eig computation. The F1 palindrome check
/// (set-equality of the spectrum under <c>λ → −2Σγ − λ</c>) is always emitted as the final
/// structural witness.</para></summary>
public static class BlockSpectrumCommand
{
    public static int Run(string[] args)
    {
        var p = new ArgParser(args);
        p.RequireNoPositional();
        int N = p.RequireInt("N");
        double? gammaScalar = p.OptionalDouble("gamma");
        string? gammaListStr = p.OptionalString("gamma-list");
        double J = p.OptionalDouble("J") ?? 1.0;
        string refineKind = p.OptionalString("refine") ?? "none";
        bool verify = p.HasFlag("verify");
        bool noVerify = p.HasFlag("no-verify");
        int topK = (int)(p.OptionalDouble("top") ?? 20.0);

        if (verify && noVerify)
            throw new ArgumentException("--verify and --no-verify are mutually exclusive.");
        if (noVerify) verify = false;

        if (gammaScalar.HasValue && gammaListStr is not null)
            throw new ArgumentException("--gamma and --gamma-list are mutually exclusive.");

        if (N < 1 || N > 12)
            throw new ArgumentOutOfRangeException(nameof(N), N, "Supported N range: 1..12.");
        if (refineKind != "none" && refineKind != "f71")
            throw new ArgumentException($"unknown --refine value: {refineKind}; expected 'none' or 'f71'.");
        if (topK < 0)
            throw new ArgumentException($"--top must be >= 0; got {topK}.");

        // Build the per-site γ array: scalar broadcasts to N copies, or parse comma-separated list.
        double[] gammaPerSite;
        string gammaDescriptor;
        if (gammaListStr is not null)
        {
            var parts = gammaListStr.Split(',', StringSplitOptions.RemoveEmptyEntries);
            if (parts.Length != N)
                throw new ArgumentException(
                    $"--gamma-list has {parts.Length} entries, expected N={N} (one γ per site).");
            gammaPerSite = parts
                .Select(s => double.Parse(s.Trim(), System.Globalization.CultureInfo.InvariantCulture))
                .ToArray();
            gammaDescriptor = $"[{string.Join(", ", gammaPerSite.Select(g => g.ToString("G6", System.Globalization.CultureInfo.InvariantCulture)))}]";
        }
        else
        {
            double gamma = gammaScalar ?? 0.5;
            gammaPerSite = Enumerable.Repeat(gamma, N).ToArray();
            gammaDescriptor = $"uniform {gamma:G6}";
        }
        double sumGamma = gammaPerSite.Sum();
        double f71AsymNorm = InhomogeneousGammaF71BreakingWitness.F71AsymmetryNorm(gammaPerSite);

        Console.WriteLine($"# block-spectrum: N={N}, J={J:G6}, gamma={gammaDescriptor}, refine={refineKind}, verify={verify}, top={topK}");
        Console.WriteLine($"#   Σγ = {sumGamma:G6}, F71 γ-asymmetry norm sqrt(Σ(γ_l−γ_{{N-1-l}})²) = {f71AsymNorm:E3}");

        // Sector summary (refinement-independent)
        int sectorCount = JointPopcountSectors.SectorCount(N);
        long maxSectorSize = JointPopcountSectors.MaxSectorSize(N);
        long maxBlockMemBytes = maxSectorSize * maxSectorSize * 16L; // complex doubles
        long fullDim = 1L << (2 * N);
        long fullMemBytes = fullDim * fullDim * 16L;

        // Cubic-cost speedup: (4^N)^3 / Σ block_size^3
        double fullCubic = Math.Pow(fullDim, 3);
        double blockCubic = 0;
        for (int pc = 0; pc <= N; pc++)
            for (int pr = 0; pr <= N; pr++)
            {
                long s = JointPopcountSectors.SectorSize(N, pc, pr);
                blockCubic += Math.Pow(s, 3);
            }
        double speedup = fullCubic / blockCubic;

        Console.WriteLine($"# joint-popcount sectors: {sectorCount} = (N+1)^2");
        Console.WriteLine($"# max-block size: {maxSectorSize} (full dim 4^N = {fullDim})");
        Console.WriteLine($"# max-block memory: {FormatBytes(maxBlockMemBytes)} (full L: {FormatBytes(fullMemBytes)})");
        Console.WriteLine($"# cubic-cost speedup vs full eig: ~{speedup:F1}x");

        // Full-L materialisation: required only for --verify (and physically impossible at N≥7
        // because (4^N)² · 16 B exceeds the .NET 2 GB single-array limit). We always need H
        // (Hilbert-space 2^N × 2^N) for the per-block path.
        bool fullLAvailable = N <= 6;
        if (verify && !fullLAvailable)
            throw new ArgumentException($"--verify is not available at N={N}: full-L (4^N × 4^N) exceeds the .NET 2 GB array-size limit at N >= 7. Drop --verify and rely on the F1 palindrome check.");

        // Build the Hilbert-space Hamiltonian once (cheap: 2^N × 2^N).
        var sw = Stopwatch.StartNew();
        var H = PauliHamiltonian.XYChain(N, J).ToMatrix();
        sw.Stop();
        Console.WriteLine($"# built H (Hilbert {1 << N}×{1 << N}) in {sw.ElapsedMilliseconds} ms");

        // Optionally also build full L at small N (≤ 6) for verify / legacy comparison.
        ComplexMatrix? L = null;
        if (fullLAvailable && verify)
        {
            sw.Restart();
            L = PauliDephasingDissipator.BuildZ(H, gammaPerSite);
            sw.Stop();
            Console.WriteLine($"# built full L (Liouville {fullDim}×{fullDim}) in {sw.ElapsedMilliseconds} ms");
        }

        // Compute spectrum via the chosen path. Track managed-memory delta + Gen0/1/2 GC counts
        // around the block-eig phase so we can report peak working-set growth. The per-block
        // path NEVER materialises the full L; only per-block dense matrices of size
        // C(N,p_c)·C(N,p_r) at a time.
        Complex[] spectrum;
        long memBefore = GC.GetTotalMemory(forceFullCollection: true);
        int gen0Before = GC.CollectionCount(0);
        int gen1Before = GC.CollectionCount(1);
        int gen2Before = GC.CollectionCount(2);
        long memDuring;
        sw.Restart();
        double? f71OffBlockNorm = null;
        if (refineKind == "f71")
        {
            var baseDecomp = JointPopcountSectorBuilder.Build(N);
            var refined = F71MirrorBlockRefinement.RefineWithF71(baseDecomp);
            int maxSubBlock = refined.SectorRanges.Max(s => s.Size);
            int nonEmpty = refined.SectorRanges.Count(s => s.Size > 0);
            Console.WriteLine($"# F71 refinement: {nonEmpty} non-empty sub-blocks (of {refined.SectorRanges.Count} total entries)");
            Console.WriteLine($"# F71-refined max sub-block size: {maxSubBlock}");
            spectrum = ComputeSpectrumWithTimingPerBlock(H, gammaPerSite, N, refineF71: true, out var sectorTimings);
            memDuring = GC.GetTotalMemory(forceFullCollection: false);
            EmitSectorTimingSummary(N, sectorTimings);
            // F71 off-block Frobenius: sum across joint-popcount sectors of the Frobenius² of the
            // (even ↔ odd) cross blocks in the rotated union-block. Zero iff γ_l = γ_{N-1-l}.
            f71OffBlockNorm = ComputeF71OffBlockNorm(H, gammaPerSite, N);
        }
        else
        {
            spectrum = ComputeSpectrumWithTimingPerBlock(H, gammaPerSite, N, refineF71: false, out var sectorTimings);
            memDuring = GC.GetTotalMemory(forceFullCollection: false);
            EmitSectorTimingSummary(N, sectorTimings);
        }
        sw.Stop();
        long memAfter = GC.GetTotalMemory(forceFullCollection: true);
        Console.WriteLine($"# computed spectrum in {sw.ElapsedMilliseconds} ms ({spectrum.Length} eigenvalues)");
        Console.WriteLine($"# managed-memory delta during block-eig: peak ≈ {FormatBytes(Math.Max(0, memDuring - memBefore))} (post-collect ≈ {FormatBytes(Math.Max(0, memAfter - memBefore))})");
        Console.WriteLine($"# GC collections during block-eig: gen0={GC.CollectionCount(0) - gen0Before}, gen1={GC.CollectionCount(1) - gen1Before}, gen2={GC.CollectionCount(2) - gen2Before}");

        // Optional verification: full-L direct eig (only available at N ≤ 6).
        if (verify && L is not null)
        {
            sw.Restart();
            var eigsFull = L.Evd().EigenValues.ToArray();
            sw.Stop();
            Console.WriteLine($"# full-L direct eig in {sw.ElapsedMilliseconds} ms ({eigsFull.Length} eigenvalues)");
            bool ok = MultisetMatches(eigsFull, spectrum, tol: 1e-9, out double maxDev);
            Console.WriteLine($"# verify: multiset match = {(ok ? "PASS" : "FAIL")} (max |Δλ| = {maxDev:E3}, tol = 1e-9)");
            if (!ok) return 3;
        }
        else
        {
            Console.WriteLine($"# verify: skipped (no --verify flag{(fullLAvailable ? "" : "; full-L impossible at N>=7 anyway")})");
        }

        // F71 off-block-Frobenius norm in the refined basis (only emitted when --refine f71):
        // 0 iff γ_l = γ_{N-1-l} (γ-distribution palindromic). Nonzero indicates the chain
        // spatial-mirror Z₂ is broken by the γ-asymmetry; magnitude scales with the F71
        // asymmetry norm of the γ-distribution. Witness of InhomogeneousGammaF71BreakingWitness.
        if (f71OffBlockNorm.HasValue)
            Console.WriteLine($"# F71 off-block Frobenius in refined basis: {f71OffBlockNorm.Value:E3} (0 iff γ palindromic; F71 asymmetry norm = {f71AsymNorm:E3})");

        // F1 palindrome check: spectrum invariant under λ → −2·Σγ − λ. Σγ = Σ_l γ_l, computed
        // above from the (possibly inhomogeneous) per-site γ array. Π is γ-blind in that it
        // depends only on the sum across sites, so this check holds regardless of asymmetry.
        bool palindrome = MultisetPalindromeOk(spectrum, sumGamma, tol: 1e-7, out double palMaxDev);
        Console.WriteLine($"# F1 palindrome check (λ → −2Σγ − λ; Σγ = {sumGamma:G6}): {(palindrome ? "PASS" : "FAIL")} (max |Δλ| = {palMaxDev:E3})");

        // First topK eigenvalues sorted by Re
        if (topK > 0)
        {
            int actualTop = Math.Min(topK, spectrum.Length);
            Console.WriteLine($"# first {actualTop} eigenvalues (sorted by Re ascending):");
            var sorted = spectrum.OrderBy(z => z.Real).ThenBy(z => z.Imaginary).Take(actualTop).ToArray();
            Console.WriteLine("#   idx        Re                Im");
            for (int i = 0; i < sorted.Length; i++)
                Console.WriteLine($"  {i,4}  {sorted[i].Real,17:E9}  {sorted[i].Imaginary,17:E9}");
        }

        return 0;
    }

    /// <summary>Per-block eig with per-sector wall-time accumulation. Uses the
    /// <see cref="PerBlockLiouvillianBuilder"/> path which never materialises full L.
    /// Mirrors <see cref="LiouvillianBlockSpectrum.ComputeSpectrumPerBlock"/> and
    /// <see cref="F71MirrorBlockRefinement.ComputeSpectrumPerBlock"/> in structure but
    /// instruments each sector's wall-time so we can summarise where the cost lives.
    /// Memory is O(blockSize²) at any instant — never O(4^N · 4^N).
    ///
    /// <para>H1: per-sector tasks run in <see cref="Parallel.ForEach"/>; per-task wall-time
    /// captured via a <see cref="ConcurrentBag{T}"/>. Outer DOP capped at ProcessorCount/4 to
    /// leave headroom for MKL inside the largest sectors' Evd (BLAS-oversubscription
    /// strategy (c)).</para>
    ///
    /// <para>H3: F71 union-block rotation B' = R^T · unionBlock · R is computed in-place via
    /// <see cref="F71MirrorBlockRefinement.RotateUnionBlockF71InPlace"/> (entry-by-entry from
    /// the Hadamard structure) instead of via two O(n³) matmuls.</para></summary>
    private static Complex[] ComputeSpectrumWithTimingPerBlock(
        ComplexMatrix H, double[] gammaPerSite, int N, bool refineF71,
        out List<(int Size, double Ms)> sectorTimings)
    {
        int liouvilleDim = 1 << (2 * N);
        int d = 1 << N;

        var spectrum = new Complex[liouvilleDim];
        var baseDecomp = JointPopcountSectorBuilder.Build(N);

        // H1: pre-compute per-sector write offsets so each task knows its destination slice.
        int sectorCount = baseDecomp.SectorRanges.Count;
        var writeOffsets = new int[sectorCount];
        int cum = 0;
        for (int i = 0; i < sectorCount; i++)
        {
            writeOffsets[i] = cum;
            cum += baseDecomp.SectorRanges[i].Size;
        }

        // X⊗N pairing (Tier 1, XGlobalChargeConjugationPairing): paired sectors share
        // spectrum exactly; primary = lex-smaller of each pair plus self-paired sectors.
        // Primaries are sorted descending by size so the largest starts first under
        // Parallel.ForEach, overlapping its wall-time with smaller sectors' work.
        var (primarySectorIndices, followerToPrimary) =
            XGlobalChargeConjugationPairing.PartitionByXNPairing(
                N, baseDecomp.SectorRanges, s => (s.PCol, s.PRow), s => s.Size);

        // Per-task timing collector (thread-safe).
        var timingBag = new ConcurrentBag<(int Size, double Ms)>();

        // BLAS-oversubscription strategy (c): cap outer DOP at ProcessorCount/4. Empirically
        // good on a 24-core box: the few largest sectors keep MKL parallelism while many
        // small sectors run concurrently. Outer DOP × MKL threads ≈ ProcessorCount.
        int outerDop = Math.Max(1, Environment.ProcessorCount / 4);
        var po = new ParallelOptions { MaxDegreeOfParallelism = outerDop };

        if (refineF71)
        {
            // Per-side F71 mirror map.
            var mirrorBits = F71MirrorIndexHelper.BuildHilbertMirrorLookup(N);
            int Mirror(int flat) => mirrorBits[flat / d] * d + mirrorBits[flat % d];

            // Phase 2: parallel eig on primary sectors only.
            var primaryEigs = new ConcurrentDictionary<int, (Complex[] Even, Complex[] Odd)>();

            Parallel.ForEach(
                primarySectorIndices, po,
                sIdx =>
                {
                    var sector = baseDecomp.SectorRanges[sIdx];
                    int size = sector.Size;
                    if (size == 0)
                    {
                        primaryEigs[sIdx] = (Array.Empty<Complex>(), Array.Empty<Complex>());
                        return;
                    }
                    var perBlock = Stopwatch.StartNew();

                    // Identify F71 orbits.
                    var sectorFlat = new int[size];
                    for (int k = 0; k < size; k++) sectorFlat[k] = baseDecomp.Permutation[sector.Offset + k];
                    var (fixedPoints, pairs) = F71MirrorIndexHelper.FindOrbitsInSector(sectorFlat, Mirror);
                    int nFix = fixedPoints.Count;
                    int nPairs = pairs.Count;
                    int unionSize = nFix + 2 * nPairs;
                    var unionFlat = new int[unionSize];
                    for (int i = 0; i < nFix; i++) unionFlat[i] = fixedPoints[i];
                    for (int k = 0; k < nPairs; k++) unionFlat[nFix + k] = pairs[k].S;
                    for (int k = 0; k < nPairs; k++) unionFlat[nFix + nPairs + k] = pairs[k].Ps;

                    var unionBlock = PerBlockLiouvillianBuilder.BuildBlockZ(H, gammaPerSite, unionFlat);

                    // H3: in-place rotation (O(n²)) instead of two MathNet matmuls (O(n³)).
                    var rotated = F71MirrorBlockRefinement.RotateUnionBlockF71InPlace(unionBlock, nFix, nPairs);

                    int evenSize = nFix + nPairs;
                    int oddSize = nPairs;
                    Complex[] evenArr = Array.Empty<Complex>();
                    Complex[] oddArr = Array.Empty<Complex>();
                    if (evenSize > 0)
                    {
                        var evenBlock = rotated.SubMatrix(0, evenSize, 0, evenSize);
                        var evenEigs = evenBlock.Evd().EigenValues;
                        evenArr = new Complex[evenSize];
                        for (int i = 0; i < evenSize; i++) evenArr[i] = evenEigs[i];
                        timingBag.Add((evenSize, perBlock.Elapsed.TotalMilliseconds));
                        perBlock.Restart();
                    }
                    if (oddSize > 0)
                    {
                        var oddBlock = rotated.SubMatrix(evenSize, oddSize, evenSize, oddSize);
                        var oddEigs = oddBlock.Evd().EigenValues;
                        oddArr = new Complex[oddSize];
                        for (int i = 0; i < oddSize; i++) oddArr[i] = oddEigs[i];
                        timingBag.Add((oddSize, perBlock.Elapsed.TotalMilliseconds));
                    }
                    primaryEigs[sIdx] = (evenArr, oddArr);
                });

            // Phase 3: sequential write to output array (primaries + followers).
            for (int sIdx = 0; sIdx < sectorCount; sIdx++)
            {
                int sourceIdx = primaryEigs.ContainsKey(sIdx) ? sIdx : followerToPrimary[sIdx];
                var (evenArr, oddArr) = primaryEigs[sourceIdx];
                int write = writeOffsets[sIdx];
                for (int i = 0; i < evenArr.Length; i++) spectrum[write + i] = evenArr[i];
                write += evenArr.Length;
                for (int i = 0; i < oddArr.Length; i++) spectrum[write + i] = oddArr[i];
            }
        }
        else
        {
            // Phase 2: parallel eig on primary sectors only.
            var primaryEigs = new ConcurrentDictionary<int, Complex[]>();
            Parallel.ForEach(
                primarySectorIndices, po,
                sIdx =>
                {
                    var sector = baseDecomp.SectorRanges[sIdx];
                    int size = sector.Size;
                    if (size == 0)
                    {
                        primaryEigs[sIdx] = Array.Empty<Complex>();
                        return;
                    }
                    var perBlock = Stopwatch.StartNew();
                    var flatIndices = new int[size];
                    for (int k = 0; k < size; k++) flatIndices[k] = baseDecomp.Permutation[sector.Offset + k];
                    var block = PerBlockLiouvillianBuilder.BuildBlockZ(H, gammaPerSite, flatIndices);
                    var blockEigs = block.Evd().EigenValues;
                    perBlock.Stop();
                    var arr = new Complex[size];
                    for (int i = 0; i < size; i++) arr[i] = blockEigs[i];
                    primaryEigs[sIdx] = arr;
                    timingBag.Add((size, perBlock.Elapsed.TotalMilliseconds));
                });

            // Phase 3: sequential write to output array (primaries + followers).
            for (int sIdx = 0; sIdx < sectorCount; sIdx++)
            {
                int sourceIdx = primaryEigs.ContainsKey(sIdx) ? sIdx : followerToPrimary[sIdx];
                var eigs = primaryEigs[sourceIdx];
                int write = writeOffsets[sIdx];
                for (int i = 0; i < eigs.Length; i++) spectrum[write + i] = eigs[i];
            }
        }

        sectorTimings = timingBag.ToList();
        return spectrum;
    }

    /// <summary>Per-sector timing summary: mean / max / total wall-time across blocks. Only
    /// emitted at N ≥ 6 where individual block costs become observable. Also lists the top-5
    /// most expensive blocks.
    ///
    /// <para>With X⊗N pairing active, the count reflects only primary sub-blocks actually
    /// eigendecomposed; follower sectors share their X⊗N-partner's spectrum at zero cost
    /// (a memcpy in Phase 3).</para>
    ///
    /// <para>Note: reported total is the sum of per-task wall-clock spans (each task uses
    /// its own Stopwatch); on the parallel path with outer DOP ≈ ProcessorCount/4, total
    /// is roughly DOP × measured wall-time. Use the "computed spectrum in ... ms" line
    /// above for true wall-time; this summary is for relative cost across blocks.</para></summary>
    private static void EmitSectorTimingSummary(int N, List<(int Size, double Ms)> sectorTimings)
    {
        if (N < 6) return;
        if (sectorTimings.Count == 0) return;
        double total = sectorTimings.Sum(t => t.Ms);
        double mean = total / sectorTimings.Count;
        var maxBlock = sectorTimings.OrderByDescending(t => t.Ms).First();
        Console.WriteLine($"# per-sector timing (X⊗N-primary only): count={sectorTimings.Count}, mean={mean:F1} ms, max={maxBlock.Ms:F1} ms (size={maxBlock.Size}), total={total:F1} ms");
        var top5 = sectorTimings.OrderByDescending(t => t.Ms).Take(5).ToArray();
        Console.WriteLine("# top-5 most expensive blocks (size, wall-time ms):");
        for (int i = 0; i < top5.Length; i++)
            Console.WriteLine($"#   #{i + 1}: size={top5[i].Size}, time={top5[i].Ms:F1} ms");
    }

    private static string FormatBytes(long bytes)
    {
        if (bytes < 1024L) return $"{bytes} B";
        if (bytes < 1024L * 1024) return $"{bytes / 1024.0:F2} KB";
        if (bytes < 1024L * 1024 * 1024) return $"{bytes / (1024.0 * 1024):F2} MB";
        return $"{bytes / (1024.0 * 1024 * 1024):F2} GB";
    }

    /// <summary>Greedy nearest-neighbour multiset matching, identical pattern to the
    /// BlockSpectrum tests' AssertMultisetEqual.</summary>
    private static bool MultisetMatches(IReadOnlyList<Complex> expected, IReadOnlyList<Complex> actual, double tol, out double maxDev)
    {
        maxDev = 0;
        if (expected.Count != actual.Count) return false;
        int n = expected.Count;
        var taken = new bool[n];
        for (int i = 0; i < n; i++)
        {
            double bestDist = double.PositiveInfinity;
            int bestJ = -1;
            for (int j = 0; j < n; j++)
            {
                if (taken[j]) continue;
                double dist = (expected[i] - actual[j]).Magnitude;
                if (dist < bestDist) { bestDist = dist; bestJ = j; }
            }
            if (bestJ < 0) return false;
            if (bestDist > maxDev) maxDev = bestDist;
            if (bestDist >= tol) return false;
            taken[bestJ] = true;
        }
        return true;
    }

    /// <summary>F1 palindromic-spectrum check: for each λ in the spectrum, verify that
    /// its mirror image −2·Σγ − λ also appears in the spectrum (within tol). Greedy
    /// nearest-neighbour matching, mirror image computed on Re only (Im invariant under
    /// reflection of Re).</summary>
    private static bool MultisetPalindromeOk(IReadOnlyList<Complex> spectrum, double sumGamma, double tol, out double maxDev)
    {
        maxDev = 0;
        int n = spectrum.Count;
        var mirrored = new Complex[n];
        for (int i = 0; i < n; i++)
            mirrored[i] = new Complex(-2 * sumGamma - spectrum[i].Real, spectrum[i].Imaginary);
        return MultisetMatches(spectrum, mirrored, tol, out maxDev);
    }

    /// <summary>F71 off-block Frobenius in the refined basis (no full-L materialisation).
    ///
    /// <para>For each joint-popcount sector we build the union-block (size = nFix + 2·nPairs)
    /// directly via <see cref="PerBlockLiouvillianBuilder.BuildBlockZ"/>, apply the same
    /// real-orthogonal rotation R as <see cref="F71MirrorBlockRefinement.ComputeSpectrumPerBlock"/>
    /// to project into (F71-even ⊕ F71-odd), and accumulate the Frobenius² of the (even↔odd)
    /// cross blocks. Across joint-popcount sectors L has no cross terms (γ-blind block
    /// structure), so the only F71-off-block contribution lives within each sector.</para>
    ///
    /// <para>Memory: O(unionSize²) at a time, max ≈ MaxSectorSize(N) entries; well below 1 GB
    /// at N=6 and within the F71 refinement budget at N=7,8.</para></summary>
    private static double ComputeF71OffBlockNorm(ComplexMatrix H, double[] gammaPerSite, int N)
    {
        int d = 1 << N;

        // Per-Hilbert-side F71 mirror map (same as ComputeSpectrumWithTimingPerBlock).
        var mirrorBits = F71MirrorIndexHelper.BuildHilbertMirrorLookup(N);
        int Mirror(int flat) => mirrorBits[flat / d] * d + mirrorBits[flat % d];

        var baseDecomp = JointPopcountSectorBuilder.Build(N);
        double offBlockFroSq = 0.0;

        foreach (var sector in baseDecomp.SectorRanges)
        {
            int size = sector.Size;
            if (size == 0) continue;

            // Identify F71 orbits within this sector.
            var sectorFlat = new int[size];
            for (int k = 0; k < size; k++) sectorFlat[k] = baseDecomp.Permutation[sector.Offset + k];
            var (fixedPoints, pairs) = F71MirrorIndexHelper.FindOrbitsInSector(sectorFlat, Mirror);
            int nFix = fixedPoints.Count;
            int nPairs = pairs.Count;
            int unionSize = nFix + 2 * nPairs;
            if (unionSize == 0) continue;

            var unionFlat = new int[unionSize];
            for (int i = 0; i < nFix; i++) unionFlat[i] = fixedPoints[i];
            for (int k = 0; k < nPairs; k++) unionFlat[nFix + k] = pairs[k].S;
            for (int k = 0; k < nPairs; k++) unionFlat[nFix + nPairs + k] = pairs[k].Ps;

            var unionBlock = PerBlockLiouvillianBuilder.BuildBlockZ(H, gammaPerSite, unionFlat);

            // H3: in-place rotation rather than two matmuls.
            var rotated = F71MirrorBlockRefinement.RotateUnionBlockF71InPlace(unionBlock, nFix, nPairs);

            // Accumulate Frobenius² over the two cross blocks (even × odd, odd × even).
            int evenSize = nFix + nPairs;
            int oddSize = nPairs;
            for (int i = 0; i < evenSize; i++)
                for (int j = 0; j < oddSize; j++)
                {
                    var z = rotated[i, evenSize + j];
                    offBlockFroSq += z.Real * z.Real + z.Imaginary * z.Imaginary;
                }
            for (int i = 0; i < oddSize; i++)
                for (int j = 0; j < evenSize; j++)
                {
                    var z = rotated[evenSize + i, j];
                    offBlockFroSq += z.Real * z.Real + z.Imaginary * z.Imaginary;
                }
        }
        return Math.Sqrt(offBlockFroSq);
    }
}
