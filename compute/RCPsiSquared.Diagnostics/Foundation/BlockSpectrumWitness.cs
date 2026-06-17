using System.Globalization;
using System.Numerics;
using System.Text.Json;
using RCPsiSquared.Core.BlockSpectrum;
using RCPsiSquared.Core.ChainSystems;
using RCPsiSquared.Core.F1;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Pauli;
using RCPsiSquared.Core.SymmetryFamily;
using ComplexMatrix = MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex>;

namespace RCPsiSquared.Diagnostics.Foundation;

/// <summary>The joint-popcount block spectrum, browsable live (<c>inspect --root blockspectrum</c>).
/// Surfaces the per-sector Liouvillian structure that the SLOW_N9 test
/// (<see cref="F1GeneralTopologyVerifiedClaim"/>, <c>F1GeneralTopologyN9BlockSpectrumChainTests</c>)
/// banks but never exposes for browsing — the <c>block_spectrum_n9</c> open arc's NextStep.
///
/// <para>The Heisenberg XXX chain H = (J/4)·Σ_b (X_bX_{b+1}+Y_bY_{b+1}+Z_bZ_{b+1}) under uniform
/// Z-dephasing γ (the N=9 banked system: J=1, γ=0.5) is exactly block-diagonal in the joint label
/// (popcount_col, popcount_row): <see cref="JointPopcountSectors"/>, (N+1)² sectors summing to 4^N.
/// What is genuinely recomputed live, cheaply, at any N:</para>
/// <list type="bullet">
///   <item>the decomposition counts — the 100 → 50 → 25 story at N=9: (N+1)² sectors, halved to
///   the X⊗N order-2 spectral classes (<see cref="XGlobalChargeConjugationPairing"/>, = the banked
///   PrimarySectorCount), quartered to the F1 Π order-4 orbit classes
///   (<see cref="F1PalindromeOrbitPairing"/>, the eig-calls the compute path actually does, since
///   Π² = X⊗N);</item>
///   <item>the full spectrum reconstructed sector-by-sector (<see cref="PerBlockLiouvillianBuilder"/>)
///   for sectors within the live-eig cap, and the F1 palindrome {λ} = {−2σ − λ} checked on it
///   (full at N ≤ 7; at N=9 the central C(9,4)² = 15876² block needs the 3 h SLOW_N9 run, so the
///   live node shows the cap-fitting sub-spectrum, which is itself Π-closed so still pairs);</item>
///   <item>the (0,1) band-edge sector sitting entirely at Re = −2γ (the Absorption floor: every
///   coherence there disagrees in exactly one bit, so L_D = −2γ·I is scalar on the block).</item>
/// </list>
///
/// <para>The full N=9 headline (262144 eigenvalues, the palindrome held bit-exact about −2σ = −9,
/// kernel 10 = N+1, gap 0.0273, 645.95× speedup) is READ live from the committed
/// <c>simulations/results/f1_n8_n9_metrics/chain_N9.json</c> (a [stored] artifact, not recomputed —
/// the run is 3 h), degrading to a "not in this checkout" note if absent. Breadcrumbed from
/// <see cref="F1GeneralTopologyVerifiedClaim"/>; no new claim (this surfaces already-typed results:
/// F1, JointPopcountSectors, the Π-orbit pairing, the Absorption Theorem, the F4 kernel).</para></summary>
public sealed class BlockSpectrumWitness : IInspectable
{
    private static readonly CultureInfo Inv = CultureInfo.InvariantCulture;
    private const double KernelTol = 1e-7;

    /// <summary>Default live-eig cap: a sector block larger than this is skipped in the live
    /// reconstruction (its eigendecomposition is too slow for an interactive inspect). 2048 keeps
    /// N ≤ 7 fully live (max block C(7,3)² = 1225) while the N=9 central blocks (up to 15876²) defer
    /// to the banked headline. The default N is 6 (max block 400²) so a bare render is ~1 s; N=7
    /// is fully live too but its 1225² blocks push an interactive render past ~10 s.</summary>
    public const long DefaultLiveEigCap = 2048;

    public int N { get; }
    public double Gamma { get; }
    public double J { get; }
    public long LiveEigCap { get; }

    public BlockSpectrumWitness(int n = 6, double gamma = 0.5, double j = 1.0, long liveEigCap = DefaultLiveEigCap)
    {
        if (n < 2 || n > 9) throw new ArgumentOutOfRangeException(nameof(n), n, "N in 2..9 for the block-spectrum witness (N=9 is the banked frontier).");
        if (gamma <= 0) throw new ArgumentOutOfRangeException(nameof(gamma), gamma, "gamma must be > 0");
        if (liveEigCap < 1) throw new ArgumentOutOfRangeException(nameof(liveEigCap), liveEigCap, "live-eig cap must be >= 1");
        N = n; Gamma = gamma; J = j; LiveEigCap = liveEigCap;
    }

    // ---- the joint-popcount decomposition (combinatorial, any N) ----

    /// <summary>The sector-decomposition facts: the (N+1)² count, the X⊗N order-2 spectral classes
    /// (= the banked PrimarySectorCount), the F1 Π order-4 orbit classes (the eig-calls), the max
    /// block size + its (p_c, p_r) sector, and the cubic-cost speedup over the dense (4^N)³.</summary>
    public readonly record struct DecompositionFacts(
        int SectorCount, int XnClasses, int PiOrbitClasses, long MaxBlock, int MaxPc, int MaxPr, double CubicSpeedup);

    public static DecompositionFacts Decomposition(int n)
    {
        var decomp = JointPopcountSectorBuilder.Build(n);
        double denseCubic = Math.Pow(Math.Pow(4, n), 3.0);   // ((2^N)²)³ = (4^N)³; double to avoid overflow at N≥9
        double totalCubic = 0.0;
        foreach (var s in decomp.SectorRanges)
            totalCubic += (double)s.Size * s.Size * s.Size;
        return new DecompositionFacts(
            SectorCount: JointPopcountSectors.SectorCount(n),
            XnClasses: XGlobalChargeConjugationPairing.DistinctSpectralClasses(n),
            PiOrbitClasses: F1PalindromeOrbitPairing.DistinctSpectralClasses(n),
            MaxBlock: JointPopcountSectors.MaxSectorSize(n),
            MaxPc: n / 2, MaxPr: n / 2,
            CubicSpeedup: totalCubic > 0 ? denseCubic / totalCubic : 0.0);
    }

    // ---- live spectrum reconstruction, sector by sector ----

    /// <summary>H = (J/4)·Σ_b (X_bX_{b+1}+Y_bY_{b+1}+Z_bZ_{b+1}) on the open chain — the N=9 banked
    /// Heisenberg XXX system. Popcount-conserving, so it lives inside the
    /// <see cref="JointPopcountSectorBuilder"/> block infrastructure.</summary>
    private static ComplexMatrix HeisenbergChain(int n, double j)
    {
        var bonds = Enumerable.Range(0, n - 1).Select(i => new Bond(i, i + 1, 1.0)).ToArray();
        var terms = new (PauliLetter, PauliLetter, Complex)[]
        {
            (PauliLetter.X, PauliLetter.X, j / 4.0),
            (PauliLetter.Y, PauliLetter.Y, j / 4.0),
            (PauliLetter.Z, PauliLetter.Z, j / 4.0),
        };
        return PauliHamiltonian.Bilinear(n, bonds, terms).ToMatrix();
    }

    private static int[] SectorFlat(JointPopcountSectorBuilder.Decomposition decomp, JointPopcountSectorBuilder.SectorRange s)
    {
        var flat = new int[s.Size];
        for (int k = 0; k < s.Size; k++) flat[k] = decomp.Permutation[s.Offset + k];
        return flat;
    }

    /// <summary>Reconstruct the spectrum from the per-sector blocks: for each joint-popcount sector
    /// whose size ≤ <paramref name="cap"/>, build the block (<see cref="PerBlockLiouvillianBuilder.BuildBlockZ"/>)
    /// and concatenate its eigenvalues. With cap ≥ the max block this is the FULL 4^N spectrum
    /// (the blocks are L's exact diagonal blocks). Larger sectors are skipped (too slow to eig live).</summary>
    public static (Complex[] Spectrum, int SectorsUsed, int SectorsSkipped, bool Full) ReconstructSpectrum(
        int n, double gamma, double j, long cap)
    {
        var H = HeisenbergChain(n, j);
        var gammaPerSite = Enumerable.Repeat(gamma, n).ToArray();
        var decomp = JointPopcountSectorBuilder.Build(n);
        var eigs = new List<Complex>();
        int used = 0, skipped = 0;
        foreach (var s in decomp.SectorRanges)
        {
            if (s.Size > cap) { skipped++; continue; }
            var block = PerBlockLiouvillianBuilder.BuildBlockZ(H, gammaPerSite, SectorFlat(decomp, s));
            foreach (var z in block.Evd().EigenValues) eigs.Add(z);
            used++;
        }
        return (eigs.ToArray(), used, skipped, skipped == 0);
    }

    public static double MinReal(Complex[] spectrum)
    {
        double min = double.PositiveInfinity;
        foreach (var z in spectrum) if (z.Real < min) min = z.Real;
        return min;
    }

    public static int KernelDimension(Complex[] spectrum, double tol = KernelTol)
    {
        int k = 0;
        foreach (var z in spectrum) if (z.Magnitude < tol) k++;
        return k;
    }

    /// <summary>The F1 symmetry distance of a spectrum: the (symmetric Hausdorff) distance between
    /// the multiset {λ} and its reflection {−2σ − λ}. ~0 iff {λ} is closed under λ ↦ −2σ − λ (F1).
    /// A nearest-neighbour set distance, NOT an index-pairing — index-pairing on sorted arrays
    /// breaks here because the reflection flips Im while many eigenvalues share a near-equal Re,
    /// so floating-point Re-noise reorders {λ} and {−2σ−λ} inconsistently. Robust to that.</summary>
    public static double PalindromePairingDistance(Complex[] spectrum, double sigma)
    {
        int m = spectrum.Length;
        var reflected = new Complex[m];
        for (int i = 0; i < m; i++) reflected[i] = -2.0 * sigma - spectrum[i];
        return Math.Max(OneSidedHausdorff(spectrum, reflected), OneSidedHausdorff(reflected, spectrum));
    }

    /// <summary>max over <paramref name="from"/> of the nearest-neighbour distance into
    /// <paramref name="to"/>. Sorts <paramref name="to"/> by Re and windows the search (prune once
    /// |ΔRe| ≥ the current best), so once an exact partner is found the window collapses — fast even
    /// with heavy spectral degeneracy.</summary>
    private static double OneSidedHausdorff(Complex[] from, Complex[] to)
    {
        var sorted = (Complex[])to.Clone();
        Array.Sort(sorted, (x, y) => x.Real.CompareTo(y.Real));
        double worst = 0.0;
        foreach (var p in from)
        {
            int lo = 0, hi = sorted.Length;
            while (lo < hi) { int mid = (lo + hi) >> 1; if (sorted[mid].Real < p.Real) lo = mid + 1; else hi = mid; }
            double best = double.PositiveInfinity;
            for (int j = lo; j < sorted.Length; j++)
            {
                if (sorted[j].Real - p.Real >= best) break;
                double d = (p - sorted[j]).Magnitude;
                if (d < best) best = d;
            }
            for (int j = lo - 1; j >= 0; j--)
            {
                if (p.Real - sorted[j].Real >= best) break;
                double d = (p - sorted[j]).Magnitude;
                if (d < best) best = d;
            }
            if (best > worst) worst = best;
        }
        return worst;
    }

    /// <summary>The Re-span of the (p_c=0, p_r=1) band-edge sector — the |1-excitation⟩⟨vacuum⟩
    /// coherences. Every basis element there disagrees in exactly one bit, so with uniform γ the
    /// dissipator is L_D = −2γ·I (scalar) on the block and L_H restricted is anti-Hermitian, so
    /// every eigenvalue has Re = −2γ exactly (the F50 weight-1 Absorption floor).</summary>
    public static (double MinRe, double MaxRe) BandEdgeSectorReSpan(int n, double gamma, double j)
    {
        var H = HeisenbergChain(n, j);
        var gammaPerSite = Enumerable.Repeat(gamma, n).ToArray();
        var decomp = JointPopcountSectorBuilder.Build(n);
        var s = decomp.SectorRanges.First(r => r.PCol == 0 && r.PRow == 1);
        var block = PerBlockLiouvillianBuilder.BuildBlockZ(H, gammaPerSite, SectorFlat(decomp, s));
        double min = double.PositiveInfinity, max = double.NegativeInfinity;
        foreach (var z in block.Evd().EigenValues)
        {
            if (z.Real < min) min = z.Real;
            if (z.Real > max) max = z.Real;
        }
        return (min, max);
    }

    // ---- the banked N=9 headline (live read of the committed artifact) ----

    /// <summary>The full-spectrum invariants from the SLOW_N9 chain run, as banked in
    /// <c>chain_N9.json</c>. A [stored] provenance — the run is ~3 h, not recomputed at inspect.</summary>
    public sealed record BankedHeadline(
        int SpectrumSize, double MinReal, int KernelDimension, double DissipationGap, double MaxImag,
        double MaxPairingDistance, int OutlierPairCount, int SectorCount, int PrimarySectorCount,
        int MaxBlockSize, int MaxPc, int MaxPr, double EffectiveSpeedup, double WallSeconds);

    /// <summary>Locate the committed chain_N9.json by walking up from the runtime base directory
    /// (and the working directory) until the full metrics path exists. Searches for the FILE, not
    /// just <c>simulations/results/</c>, so a shadow copy of that directory in a test bin output
    /// (which lacks the f1_n8_n9_metrics artifact) does not mask the real repo file.</summary>
    private static string? FindBankedN9File()
    {
        foreach (var start in new[] { AppContext.BaseDirectory, Directory.GetCurrentDirectory() })
        {
            var dir = new DirectoryInfo(start);
            while (dir != null)
            {
                var candidate = Path.Combine(dir.FullName, "simulations", "results", "f1_n8_n9_metrics", "chain_N9.json");
                if (File.Exists(candidate)) return candidate;
                dir = dir.Parent;
            }
        }
        return null;
    }

    /// <summary>Read the committed N=9 chain metrics; null if the artifact is not present in this
    /// checkout (e.g. the SLOW_N9 test has never been run here).</summary>
    public static BankedHeadline? ReadBankedN9()
    {
        string? path = FindBankedN9File();
        if (path is null) return null;
        try
        {
            using var doc = JsonDocument.Parse(File.ReadAllText(path));
            var r = doc.RootElement;
            return new BankedHeadline(
                SpectrumSize: r.GetProperty("SpectrumSize").GetInt32(),
                MinReal: r.GetProperty("MinReal").GetDouble(),
                KernelDimension: r.GetProperty("KernelDimension").GetInt32(),
                DissipationGap: r.GetProperty("DissipationGap").GetDouble(),
                MaxImag: r.GetProperty("MaxImag").GetDouble(),
                MaxPairingDistance: r.GetProperty("MaxPairingDistance").GetDouble(),
                OutlierPairCount: r.GetProperty("OutlierPairCount").GetInt32(),
                SectorCount: r.GetProperty("SectorCount").GetInt32(),
                PrimarySectorCount: r.GetProperty("PrimarySectorCount").GetInt32(),
                MaxBlockSize: r.GetProperty("MaxBlockSize").GetInt32(),
                MaxPc: r.GetProperty("MaxBlockSectorPCol").GetInt32(),
                MaxPr: r.GetProperty("MaxBlockSectorPRow").GetInt32(),
                EffectiveSpeedup: r.GetProperty("EffectiveSpeedupOverDense").GetDouble(),
                WallSeconds: r.GetProperty("ComputeSpectrumWallSeconds").GetDouble());
        }
        catch (Exception ex) when (ex is JsonException or IOException) { return null; }
    }

    // ---- the inspect tree ----

    private InspectableNode TheDecompositionNode()
    {
        var d = Decomposition(N);
        var kids = new List<IInspectable>
        {
            InspectableNode.RealScalar("sectors (N+1)²", d.SectorCount),
            InspectableNode.RealScalar("X⊗N order-2 spectral classes (= banked PrimarySectorCount)", d.XnClasses),
            InspectableNode.RealScalar("F1 Π order-4 orbit classes (the eig-calls; Π² = X⊗N)", d.PiOrbitClasses),
            new InspectableNode("max block",
                summary: $"C({N},{N / 2})² = {d.MaxBlock} at sector ({d.MaxPc},{d.MaxPr})"),
            InspectableNode.RealScalar("cubic-cost speedup over dense (4^N)³", d.CubicSpeedup, "0.0"),
        };
        return new InspectableNode("the joint-popcount sector decomposition",
            summary: $"N={N}: (N+1)² = {d.SectorCount} sectors, halved by X⊗N to {d.XnClasses} spectral " +
                     $"classes (Π² = X⊗N, a genuine symmetry → verbatim copy), quartered by the F1 Π " +
                     $"order-4 orbit to {d.PiOrbitClasses} eig-calls; max block C({N},{N / 2})² = {d.MaxBlock}; " +
                     $"{d.CubicSpeedup.ToString("0.0", Inv)}× cubic-cost speedup over the dense (4^N)³.",
            children: kids);
    }

    /// <summary>The navigation hub: each load-bearing joint-popcount sector points to the
    /// sector-specific witness(es) that max-zoom it. blockspectrum is the OVERVIEW; reduction /
    /// ceiling / horizon / survivor / secondclock are the zooms on individual sectors of this same
    /// decomposition. (Sector→witness mappings verified at each witness's source.)</summary>
    private static InspectableNode TheSectorMapNode()
    {
        var kids = new List<IInspectable>
        {
            new InspectableNode("(0,1) — the band edge",
                summary: "inspect --root reduction (SectorReductionWitness): the |1-exc⟩⟨vac| birth-canal " +
                         "boundary mode; the whole sector sits at Re=−2γ (the Absorption floor). Its {0,2} " +
                         "junction at N≥6 crosses into the (2,2) sector."),
            new InspectableNode("(1,1) — the commutant",
                summary: "inspect --root ceiling (StructuralCeilingWitness): the high-Q structural ceiling " +
                         "g2 = strict_gap/2γ from the darkest [H,A]=0 coherence here — the S_N standard-rep " +
                         "sector, g2(K_N)=4/N, g2(star_N)=4/(N−1)."),
            new InspectableNode("single-excitation {0,2} — the EP",
                summary: "inspect --root horizon (CoherenceHorizonWitness): the coherence horizon Q*(N) where " +
                         "the slowest mode stops oscillating — the single-excitation Haken-Strobl √-EP (4^N→N²), " +
                         "= the carbon Frost-Hückel coherent↔incoherent threshold."),
            new InspectableNode("(p,p) half-filling — the second clock",
                summary: "inspect --root survivor (IncompletenessSurvivorWitness): where the longest-lived " +
                         "dissipative mode lives — the interior C=0.5 incompleteness coherence; and inspect " +
                         "--root secondclock (SecondClockRegimeWitness): the {0,2}/half-filling second clock, " +
                         "regime = map(band degeneracy, dispersion). The (2,2) sector is the recurring N=4 " +
                         "anomaly (ceiling's K_4 = 2−2/√3, ring-4 = 1)."),
        };
        return new InspectableNode("the sector map — which live witness zooms each load-bearing sector",
            summary: "the per-sector overview is the index to the sector-specific witnesses: each is a max-zoom " +
                     "on one sector of THIS decomposition. (0,1) → reduction; (1,1) → ceiling; single-excitation " +
                     "{0,2} → horizon; (p,p)/(2,2) half-filling → survivor + secondclock.",
            children: kids);
    }

    private InspectableNode ThePalindromeNode()
    {
        var (spectrum, used, skipped, full) = ReconstructSpectrum(N, Gamma, J, LiveEigCap);
        double sigma = N * Gamma;
        double pairing = PalindromePairingDistance(spectrum, sigma);
        double minRe = MinReal(spectrum);
        int kernel = KernelDimension(spectrum);
        string scope = full
            ? $"all {used} sectors rebuilt → the full {spectrum.Length} = 4^{N} spectrum"
            : $"{used} of {used + skipped} sectors fit the live-eig cap {LiveEigCap} (the central blocks " +
              $"need the SLOW_N9 run) → a {spectrum.Length}-eigenvalue Π-closed sub-spectrum";
        var kids = new List<IInspectable>
        {
            new InspectableNode("F1 max pairing distance |{λ} vs {−2σ−λ}|",
                summary: $"{pairing.ToString("E3", Inv)} (the palindrome about −σ = {(-sigma).ToString("0.###", Inv)})"),
            InspectableNode.RealScalar("MinReal (= −2σ when full: the maximally-dephased coherence)", minRe, "0.000000"),
            new InspectableNode("kernel dimension",
                summary: full
                    ? $"{kernel} = N+1 = {N + 1} (the connected-chain steady states)"
                    : $"{kernel} in this {used}-sector sub-spectrum (the full connected-chain kernel N+1 = {N + 1} is in the banked headline)"),
        };
        return new InspectableNode("the F1 palindrome, reconstructed live",
            summary: $"{scope}; F1 holds: max pairing distance {pairing.ToString("E3", Inv)}, " +
                     $"MinReal {minRe.ToString("0.0000", Inv)}" + (full ? $" = −2σ = {(-2 * sigma).ToString("0.###", Inv)}" : "") +
                     $", kernel {kernel}{(full ? $" = N+1 = {N + 1}" : "")}.",
            children: kids);
    }

    private InspectableNode TheAbsorptionFloorNode()
    {
        var (minRe, maxRe) = BandEdgeSectorReSpan(N, Gamma, J);
        return new InspectableNode("the per-sector Absorption floor (Re = −2γ)",
            summary: $"the (0,1) band-edge sector — the |1-exc⟩⟨vac| coherences, {N}-dim — sits entirely at " +
                     $"Re ∈ [{minRe.ToString("0.0000", Inv)}, {maxRe.ToString("0.0000", Inv)}] = −2γ = " +
                     $"{(-2 * Gamma).ToString("0.###", Inv)} (every coherence disagrees in one bit, so L_D = −2γ·I " +
                     "is scalar there; the F50 weight-1 floor / Absorption Theorem). The decay GAP lives instead " +
                     "in the diagonal (k,k) sectors.");
    }

    private InspectableNode TheBankedN9Node()
    {
        var b = ReadBankedN9();
        if (b is null)
            return new InspectableNode("the N=9 banked headline [chain_N9.json]",
                summary: "not present in this checkout — run the SLOW_N9 test " +
                         "(dotnet test ... --filter \"Category=SLOW_N9\") to bank it (~3 h, writes " +
                         "simulations/results/f1_n8_n9_metrics/chain_N9.json).");
        var kids = new List<IInspectable>
        {
            InspectableNode.RealScalar("spectrum size (= 4^9)", b.SpectrumSize),
            InspectableNode.RealScalar("MinReal (= −2σ = −9)", b.MinReal, "0.000000"),
            InspectableNode.RealScalar("kernel dimension (= N+1)", b.KernelDimension),
            InspectableNode.RealScalar("dissipation gap", b.DissipationGap, "0.00000000"),
            InspectableNode.RealScalar("max |Im|", b.MaxImag, "0.000000"),
            new InspectableNode("F1 pairing (bit-exact)",
                summary: $"max pairing distance {b.MaxPairingDistance.ToString("E3", Inv)}, {b.OutlierPairCount} outliers"),
            new InspectableNode("block structure",
                summary: $"{b.SectorCount} sectors, {b.PrimarySectorCount} X⊗N classes, max block {b.MaxBlockSize}² at ({b.MaxPc},{b.MaxPr})"),
            InspectableNode.RealScalar("effective speedup over dense (4^9)³", b.EffectiveSpeedup, "0.0"),
        };
        return new InspectableNode("the N=9 banked headline [chain_N9.json]",
            summary: $"[stored] the full {b.SpectrumSize}-eigenvalue N=9 chain run (Heisenberg XXX, J=1, γ=0.5): " +
                     $"the F1 palindrome held bit-exact about −2σ = {b.MinReal.ToString("0.#", Inv)} (max pairing " +
                     $"distance {b.MaxPairingDistance.ToString("E2", Inv)}, {b.OutlierPairCount} outliers), kernel " +
                     $"{b.KernelDimension} = N+1, gap {b.DissipationGap.ToString("0.0000", Inv)}, {b.SectorCount} " +
                     $"sectors via {b.EffectiveSpeedup.ToString("0.0", Inv)}× cubic-cost speedup. Wall " +
                     $"{(b.WallSeconds / 3600.0).ToString("0.0", Inv)} h — read from chain_N9.json, not recomputed.",
            children: kids);
    }

    public string DisplayName =>
        $"BlockSpectrumWitness (N={N}, Heisenberg chain, γ={Gamma.ToString("0.###", Inv)}, J={J.ToString("0.###", Inv)})";

    public string Summary =>
        "the joint-popcount block spectrum, live: the (N+1)² sector decomposition (halved by X⊗N, " +
        "quartered by the F1 Π orbit), a sector map indexing the sector-specific witnesses (reduction / " +
        "ceiling / horizon / survivor / secondclock each max-zoom one sector of this decomposition), the " +
        "F1 palindrome {λ} = {−2σ − λ} reconstructed sector-by-sector (full at N ≤ 7), the (0,1) band-edge " +
        "Absorption floor Re = −2γ, and the N=9 banked headline read live from chain_N9.json. The browsable " +
        "overview face of the SLOW_N9 result (arc block_spectrum_n9).";

    public IEnumerable<IInspectable> Children
    {
        get
        {
            yield return TheDecompositionNode();
            yield return TheSectorMapNode();
            yield return ThePalindromeNode();
            yield return TheAbsorptionFloorNode();
            yield return TheBankedN9Node();
        }
    }

    public InspectablePayload Payload => InspectablePayload.Empty;
}
