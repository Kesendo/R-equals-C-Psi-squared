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
///   <item>the (0,1) band-edge sector sitting entirely at Re = −2γ AT UNIFORM γ (the Absorption
///   floor: every coherence there disagrees in exactly one bit, so L_D = −2γ·I is scalar on the
///   block; under a profile it is −2·diag(γ), still normal but no longer scalar);</item>
///   <item>that same sector under a per-site γ PROFILE, where the dissipator is diagonal instead of
///   scalar: the generator is still exactly M = +i·(1/2)·𝓛 − 2·diag(γ) entry-wise, with 𝓛 the
///   WEIGHTED graph Laplacian diag(deg) − A whose degree DIAGONAL is what the ZZ term supplies
///   (XXX only), and the rate law is <see cref="AbsorptionTheoremClaim"/>'s per-channel one read on
///   this block: Re λ_k is the γ-weighted average of the mode's site occupancy, with the Bendixson
///   bracket and the trace (which pins the MEAN, not the bracket's midpoint) as its two corollaries.
///   The arc <c>site_resolved_vacuum_block</c>'s composite: D10 has the degree
///   term at uniform γ, <see cref="VacuumBlockReductionClaim"/> has the profile with an XY H.</item>
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
    private static ComplexMatrix HeisenbergChain(int n, double j) =>
        HeisenbergGraph(n, ChainBonds(n, j));

    /// <summary>The open-chain bond list at uniform coupling J.</summary>
    public static Bond[] ChainBonds(int n, double j) =>
        Enumerable.Range(0, n - 1).Select(i => new Bond(i, i + 1, j)).ToArray();

    /// <summary>H = Σ_b (J_b/4)·(X_aX_b + Y_aY_b + Z_aZ_b) on an arbitrary bond list, the coupling
    /// read PER BOND from <see cref="Bond.Coupling"/>. XXX only: the ZZ coefficient must equal the
    /// XY one, otherwise the (0,1) block's generator is +i·(1/2)·(Δ·diag(deg) − A), which is not a
    /// graph Laplacian. Popcount-conserving, so it lives inside the
    /// <see cref="JointPopcountSectorBuilder"/> block infrastructure.</summary>
    private static ComplexMatrix HeisenbergGraph(int n, IReadOnlyList<Bond> bonds)
    {
        var terms = new List<PauliTerm>(bonds.Count * 3);
        foreach (var b in bonds)
            foreach (var p in new[] { PauliLetter.X, PauliLetter.Y, PauliLetter.Z })
                terms.Add(PauliTerm.TwoSite(n, b.Site1, p, b.Site2, p, b.Coupling / 4.0));
        return new PauliHamiltonian(n, terms).ToMatrix();
    }

    private static int[] SectorFlat(JointPopcountSectorBuilder.Decomposition decomp, JointPopcountSectorBuilder.SectorRange s)
    {
        var flat = new int[s.Size];
        for (int k = 0; k < s.Size; k++) flat[k] = decomp.Permutation[s.Offset + k];
        return flat;
    }

    /// <summary>Reconstruct the spectrum from the per-sector blocks of a prebuilt H: for each
    /// joint-popcount sector whose size ≤ <paramref name="cap"/>, build the block
    /// (<see cref="PerBlockLiouvillianBuilder.BuildBlockZ"/>) and concatenate its eigenvalues. With
    /// cap ≥ the max block this is the FULL 4^N spectrum (the blocks are L's exact diagonal blocks).
    /// Larger sectors are skipped (too slow to eig live). Takes H so a render builds it once and
    /// shares it with the band-edge node (H grows as 2^N — cheap at N≤9, but built once on principle).</summary>
    public static (Complex[] Spectrum, int SectorsUsed, int SectorsSkipped, bool Full) ReconstructSpectrum(
        ComplexMatrix h, int n, double gamma, long cap)
    {
        var gammaPerSite = Enumerable.Repeat(gamma, n).ToArray();
        var decomp = JointPopcountSectorBuilder.Build(n);
        var eigs = new List<Complex>();
        int used = 0, skipped = 0;
        foreach (var s in decomp.SectorRanges)
        {
            if (s.Size > cap) { skipped++; continue; }
            var block = PerBlockLiouvillianBuilder.BuildBlockZ(h, gammaPerSite, SectorFlat(decomp, s));
            foreach (var z in block.Evd().EigenValues) eigs.Add(z);
            used++;
        }
        return (eigs.ToArray(), used, skipped, skipped == 0);
    }

    /// <summary>Convenience overload that builds the Heisenberg chain H (for tests / direct callers).</summary>
    public static (Complex[] Spectrum, int SectorsUsed, int SectorsSkipped, bool Full) ReconstructSpectrum(
        int n, double gamma, double j, long cap) =>
        ReconstructSpectrum(HeisenbergChain(n, j), n, gamma, cap);

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

    /// <summary>The F1 symmetry distance of the spectrum about −σ: ~0 iff {λ} is closed as a
    /// MULTISET under λ ↦ −2σ − λ (the mirror-symmetry F1). Delegates to the canonical
    /// multiplicity-aware greedy matcher <see cref="F1SpectrumStatistics.MaxF1PairingDistance"/>
    /// (the same one the SLOW F1 dogfood tests use). A set/Hausdorff distance is multiplicity-BLIND
    /// — it cannot see a dropped or duplicated eigenvalue that still has a same-valued neighbour,
    /// exactly the failure a reconstruction-wiring bug would produce — so it must not be used here.</summary>
    public static double PalindromePairingDistance(Complex[] spectrum, double sigma) =>
        F1SpectrumStatistics.MaxF1PairingDistance(spectrum, sigma);

    /// <summary>The (p_c=0, p_r=1) band-edge sector block itself, the N×N operator M on the
    /// |1-excitation⟩⟨vacuum| coherences, with a SITE-RESOLVED γ profile
    /// (<paramref name="gammaPerSite"/>[l] is site l, the leftmost Kronecker factor). This is the
    /// operator the <c>site_resolved_vacuum_block</c> arc is about: the engine
    /// (<see cref="PerBlockLiouvillianBuilder.BuildBlockZ"/>) has always returned it site-resolved;
    /// only the callers here were uniform-only.</summary>
    public static ComplexMatrix BandEdgeSectorBlock(ComplexMatrix h, int n, IReadOnlyList<double> gammaPerSite)
    {
        if (gammaPerSite.Count != n)
            throw new ArgumentException($"gamma profile has {gammaPerSite.Count} entries, expected N={n}", nameof(gammaPerSite));
        var decomp = JointPopcountSectorBuilder.Build(n);
        var s = decomp.SectorRanges.First(r => r.PCol == 0 && r.PRow == 1);
        return PerBlockLiouvillianBuilder.BuildBlockZ(h, gammaPerSite, SectorFlat(decomp, s));
    }

    /// <summary>Convenience overload that builds the Heisenberg chain H (for tests / direct callers).</summary>
    public static ComplexMatrix BandEdgeSectorBlock(int n, IReadOnlyList<double> gammaPerSite, double j) =>
        BandEdgeSectorBlock(HeisenbergChain(n, j), n, gammaPerSite);

    /// <summary>Convenience overload on an arbitrary XXX bond graph, the coupling per bond.</summary>
    public static ComplexMatrix BandEdgeSectorBlock(int n, IReadOnlyList<Bond> bonds, IReadOnlyList<double> gammaPerSite) =>
        BandEdgeSectorBlock(HeisenbergGraph(n, bonds), n, gammaPerSite);

    /// <summary>The Re-span of the (p_c=0, p_r=1) band-edge sector under a γ PROFILE. Every basis
    /// element there disagrees in exactly one bit, so the dissipator contributes −2γ_j per element:
    /// diagonal, not scalar. The Hermitian part of M is then exactly −2·diag(γ), which is
    /// <see cref="AbsorptionTheoremClaim"/>'s per-channel law on this block: Re λ_k is the
    /// γ-weighted average of the mode's own site occupancy (<see cref="PerModeAbsorptionResidual"/>).
    /// Two corollaries follow, and they are corollaries: an occupancy is a probability distribution,
    /// so Re λ is a convex combination of the −2γ_l, i.e. the Bendixson bracket
    /// Re λ ∈ [−2·max γ, −2·min γ]; and the trace pins Σ Re λ = −2σ, mean Re = −2·γ̄. At uniform γ
    /// the bracket closes on a point and the whole sector sits at Re = −2γ (the F50 weight-1 floor);
    /// a profile generally opens it into an interval around the same −2·γ̄. GENERALLY, not always:
    /// at N=2 the 2×2 block splits by 2·√(γ̄_diff² − J²/4) and stays closed for J ≥ 2·γ̄_diff, so read
    /// the measured span rather than assuming it opened.</summary>
    public static (double MinRe, double MaxRe) BandEdgeSectorReSpan(ComplexMatrix h, int n, IReadOnlyList<double> gammaPerSite)
    {
        var block = BandEdgeSectorBlock(h, n, gammaPerSite);
        double min = double.PositiveInfinity, max = double.NegativeInfinity;
        foreach (var z in block.Evd().EigenValues)
        {
            if (z.Real < min) min = z.Real;
            if (z.Real > max) max = z.Real;
        }
        return (min, max);
    }

    /// <summary>Uniform-γ overload (the F50 weight-1 floor case: MinRe = MaxRe = −2γ).</summary>
    public static (double MinRe, double MaxRe) BandEdgeSectorReSpan(ComplexMatrix h, int n, double gamma) =>
        BandEdgeSectorReSpan(h, n, Enumerable.Repeat(gamma, n).ToArray());

    /// <summary>Convenience overload that builds the Heisenberg chain H (for tests / direct callers).</summary>
    public static (double MinRe, double MaxRe) BandEdgeSectorReSpan(int n, double gamma, double j) =>
        BandEdgeSectorReSpan(HeisenbergChain(n, j), n, gamma);

    /// <summary>Convenience overload that builds the Heisenberg chain H and takes a γ profile.</summary>
    public static (double MinRe, double MaxRe) BandEdgeSectorReSpan(int n, IReadOnlyList<double> gammaPerSite, double j) =>
        BandEdgeSectorReSpan(HeisenbergChain(n, j), n, gammaPerSite);

    /// <summary>The WEIGHTED graph Laplacian 𝓛 = diag(deg) − A of a bond list, each bond entering
    /// with its own <see cref="Bond.Coupling"/>. At uniform J this is J·(the unweighted Laplacian),
    /// so the generator below reads +i·(J/2)·𝓛_unweighted there.</summary>
    private static ComplexMatrix WeightedLaplacian(int n, IReadOnlyList<Bond> bonds)
    {
        var l = MathNet.Numerics.LinearAlgebra.Matrix<Complex>.Build.Dense(n, n);
        foreach (var b in bonds)
        {
            l[b.Site1, b.Site2] -= b.Coupling; l[b.Site2, b.Site1] -= b.Coupling;
            l[b.Site1, b.Site1] += b.Coupling; l[b.Site2, b.Site2] += b.Coupling;
        }
        return l;
    }

    /// <summary>Which SITE is excited in each row of the (0,1) block, in the block's own basis
    /// order. The order is not site order: the sector's flat indices are row·d + col with
    /// row = 1&lt;&lt;b, sorted ascending, and site l occupies bit n−1−l (l is the leftmost Kronecker
    /// factor), so the block runs from site n−1 down to site 0. Any site-indexed quantity, a γ
    /// profile or a weighted Laplacian, must be permuted through this before it can be compared with
    /// the block entry-wise. What the permutation actually bites on is worth stating: a SYMMETRIC γ
    /// profile cannot see it, and neither can a uniform-J chain or ring Laplacian, which is
    /// reversal-invariant. It takes an asymmetric γ or a non-palindromic per-bond J to make the
    /// difference visible, and the gates use both.
    /// <para>NOT A NEW FINDING. The same site-versus-bit ordering is stated in
    /// <see cref="PerBlockLiouvillianBuilder"/>'s own γ-index-convention paragraph, gated by
    /// <c>PerBlockLiouvillianBuilderGammaOrderTests</c> (which carries a two-sided mutation control
    /// of exactly this shape), and written out in <c>simulations/sacrifice_zone_optics.py</c>. This
    /// method exposes the permutation for entry-wise comparisons against the (0,1) block; it does
    /// not discover the convention.</para></summary>
    public static int[] BandEdgeSectorSiteOrder(int n)
    {
        int d = 1 << n;
        var decomp = JointPopcountSectorBuilder.Build(n);
        var s = decomp.SectorRanges.First(r => r.PCol == 0 && r.PRow == 1);
        var flat = SectorFlat(decomp, s);
        var sites = new int[flat.Length];
        for (int k = 0; k < flat.Length; k++)
        {
            int row = flat[k] / d;
            sites[k] = n - 1 - System.Numerics.BitOperations.TrailingZeroCount((uint)row);
        }
        return sites;
    }

    /// <summary>The entry-wise residual of the site-resolved generator identity on ANY XXX bond
    /// graph, with this witness's convention H = Σ_b (J_b/4)·(XX+YY+ZZ):
    /// <para>M = +i·(1/2)·𝓛 − 2·diag(γ), 𝓛 the WEIGHTED graph Laplacian diag(deg) − A (at uniform
    /// J: +i·(J/2)·𝓛_unweighted).</para>
    /// This is the composite the arc <c>site_resolved_vacuum_block</c> names as owned by neither
    /// side: D10 (<c>D10_W1_DISPERSION.md</c>) has the ZZ degree term but uniform γ;
    /// <see cref="VacuumBlockReductionClaim"/> has the γ profile but an XY H, hence no degree term.
    /// Recomputed live against <see cref="PerBlockLiouvillianBuilder.BuildBlockZ"/>. The residual is
    /// ABSOLUTE and built from cancelling O(J·N) terms, so it is machine zero relative to the entry
    /// scale J·N + 2·max γ, not relative to 1; a gate on it must scale with that, and the gates do.
    /// <para>SCOPE: XXX only. The Laplacian appears because the ZZ coefficient equals the XY one;
    /// for XXZ the generator is +i·(1/2)·(Δ·diag(deg) − A), which is not a Laplacian.</para>
    /// <para>The non-normality this makes visible is not new here: <c>experiments/ANALYTICAL_SPECTRUM.md</c>
    /// (the non-uniform-dephasing paragraph) and <c>docs/ANALYTICAL_FORMULAS.md</c>'s F2 fence both
    /// already state that under a profile the block stays exactly closed while diag(γ) and 𝓛 stop
    /// commuting, so the eigenvalues are not −2γ_j − 2iJ·μ_m. What is new here is the entry-wise
    /// identity under a live witness, gated across graphs, and the per-mode reading below.</para>
    /// <para>BEFORE READING A SIGN CONFLICT, separate the two axes. The i-sign is meaningless on its
    /// own, because 𝓛 = diag(deg) − A carries the opposite off-diagonal sign to the bare hopping h:
    /// <see cref="VacuumBlockReductionClaim"/>'s −iQ·h and this +i·(J/2)·𝓛 AGREE off the diagonal,
    /// and the degree diagonal that ZZ supplies is the whole difference. The other axis is
    /// ORIENTATION, and it is also a label collision. The sector selected here is
    /// (p_c=0, p_r=1) = (popcount of the bra, popcount of the ket) = |1-exc⟩⟨vac|, which carries
    /// <c>+i</c>; the conjugate orientation |vac⟩⟨1-exc| carries <c>−i</c>. D10:139 writes
    /// <c>−2iJ·𝓛</c>, i.e. the conjugate operator, and labels it "(0,1)" too, because the object it
    /// derives is |0⟩⟨j| (D10:112). So the SAME written label names conjugate operators in the two
    /// places. Same Re, mirrored frequencies; never copy a sign between them. The convention here
    /// is <c>JointPopcountSectorBuilder</c>'s (PCol, PRow) = (bra, ket), which agrees with
    /// MirrorWorld's <c>Block</c> (P = bra, Q = ket); D10 uses the label in the other order. Its
    /// normalisation also differs: D10's H = J·Σ(XX+YY+ZZ) makes its J a quarter of this one, so
    /// its −2iJ𝓛 is this −i(J/2)𝓛.</para></summary>
    public static double GeneratorResidual(int n, IReadOnlyList<Bond> bonds, IReadOnlyList<double> gammaPerSite)
    {
        var m = BandEdgeSectorBlock(HeisenbergGraph(n, bonds), n, gammaPerSite);
        var lap = WeightedLaplacian(n, bonds);
        var sites = BandEdgeSectorSiteOrder(n);   // the block's basis order is not site order
        double worst = 0.0;
        for (int r = 0; r < n; r++)
            for (int c = 0; c < n; c++)
            {
                var predicted = Complex.ImaginaryOne * 0.5 * lap[sites[r], sites[c]];
                if (r == c) predicted -= 2.0 * gammaPerSite[sites[r]];
                worst = Math.Max(worst, (m[r, c] - predicted).Magnitude);
            }
        return worst;
    }

    /// <summary>Uniform-J open-chain overload of <see cref="GeneratorResidual"/>.</summary>
    public static double ChainGeneratorResidual(int n, double j, IReadOnlyList<double> gammaPerSite) =>
        GeneratorResidual(n, ChainBonds(n, j), gammaPerSite);

    /// <summary>The entry-wise residual of Herm(M) = (M + M†)/2 = −2·diag(γ): the sharp statement
    /// the Bendixson bracket and the trace identity are both corollaries of, and the one place
    /// where <see cref="AbsorptionTheoremClaim"/>'s per-channel law meets this block. All of the
    /// coupling is anti-Hermitian, whatever the graph and whatever the per-bond J.</summary>
    public static double HermitianPartResidual(int n, IReadOnlyList<Bond> bonds, IReadOnlyList<double> gammaPerSite)
    {
        var m = BandEdgeSectorBlock(HeisenbergGraph(n, bonds), n, gammaPerSite);
        var sites = BandEdgeSectorSiteOrder(n);
        double worst = 0.0;
        for (int r = 0; r < n; r++)
            for (int c = 0; c < n; c++)
            {
                var herm = 0.5 * (m[r, c] + Complex.Conjugate(m[c, r]));
                if (r == c) herm += 2.0 * gammaPerSite[sites[r]];
                worst = Math.Max(worst, herm.Magnitude);
            }
        return worst;
    }

    /// <summary>The residual of <see cref="AbsorptionTheoremClaim"/>'s per-channel absorption law
    /// on this block, mode by mode:
    /// <para>−Re(λ_k) = 2·Σ_l γ_l·⟨Δ_l⟩_k, with ⟨Δ_l⟩_k = ⟨v_k|N_l|v_k⟩/‖v_k‖².</para>
    /// On the (0,1) block ⟨Δ_l⟩_k is just the mode's SITE OCCUPANCY |v_k(l)|²/‖v_k‖², because every
    /// basis element disagrees at exactly its own site. This is the master statement here, and it
    /// is the typed parent's, not a new one: the mode-by-mode equality does not die under a
    /// profile, it becomes a γ-weighted average over the mode's own occupancy. Since the occupancy
    /// is a probability distribution over sites, Re λ_k is a convex combination of the −2γ_l, which
    /// IS the Bendixson bracket, derived sharper. Non-normality of M costs nothing: Re λ is the
    /// Rayleigh quotient of Herm(M) on the RIGHT eigenvector for any eigenvector at all.</summary>
    public static double PerModeAbsorptionResidual(int n, IReadOnlyList<Bond> bonds, IReadOnlyList<double> gammaPerSite)
    {
        var m = BandEdgeSectorBlock(HeisenbergGraph(n, bonds), n, gammaPerSite);
        var sites = BandEdgeSectorSiteOrder(n);
        var evd = m.Evd();
        double worst = 0.0;
        for (int k = 0; k < n; k++)
        {
            var v = evd.EigenVectors.Column(k);
            double norm2 = v.Sum(z => z.Real * z.Real + z.Imaginary * z.Imaginary);
            double predicted = 0.0;
            for (int r = 0; r < n; r++)
                predicted -= 2.0 * gammaPerSite[sites[r]] * (v[r].Real * v[r].Real + v[r].Imaginary * v[r].Imaginary) / norm2;
            worst = Math.Max(worst, Math.Abs(evd.EigenValues[k].Real - predicted));
        }
        return worst;
    }

    /// <summary>A deterministic non-uniform γ profile with the SAME total σ = N·γ as the uniform
    /// one, so the palindrome center −σ is untouched and only the band-edge line splits: the
    /// "deep-edge" shape (both chain ends depressed to γ/4, the bulk flat to make up the sum) that
    /// <see cref="SectorReductionWitness"/> uses as its general-N sibling of the N=5 canal anchor,
    /// here scaled to mean γ. Deep-edge, canal and uniform are three distinct profiles in this repo;
    /// this is the deep-edge one. At γ=1 it is that profile entry for entry, and at N=5 it coincides
    /// with the canal anchor
    /// [0.25, 1.5, 1.5, 1.5, 0.25]. At N=2 there is no bulk to carry the remainder, so the shape
    /// degenerates to [γ/4, 7γ/4]: site 1 is not an end in the sense the shape means, it is
    /// whatever is left. N ≥ 2.</summary>
    public static double[] DeepEdgeProfile(int n, double gamma)
    {
        if (n < 2) throw new ArgumentOutOfRangeException(nameof(n), n, "the deep-edge shape needs at least two sites");
        var p = new double[n];
        double edge = 0.25 * gamma;
        if (n == 2) { p[0] = edge; p[1] = n * gamma - edge; return p; }
        double rest = (n * gamma - 2 * edge) / (n - 2);
        for (int i = 0; i < n; i++) p[i] = (i == 0 || i == n - 1) ? edge : rest;
        return p;
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
                summary: $"C({N},{N / 2})² = {d.MaxBlock} at sector ({d.MaxPc},{d.MaxPr})" +
                         (N % 2 == 0
                            ? " (the unique max at even N)"
                            : ", one of a 4-way tie at odd N — C(N,(N−1)/2)=C(N,(N+1)/2), so (⌊N/2⌋,⌈N/2⌉) etc. tie")),
            InspectableNode.RealScalar("cubic-cost speedup over dense (4^N)³", d.CubicSpeedup, "0.0"),
        };
        return new InspectableNode("the joint-popcount sector decomposition",
            summary: $"N={N}: (N+1)² = {d.SectorCount} sectors, halved by X⊗N to {d.XnClasses} spectral " +
                     $"classes (Π² = X⊗N, a genuine symmetry → verbatim copy), quartered by the F1 Π " +
                     $"order-4 orbit to {d.PiOrbitClasses} eig-calls; max block C({N},{N / 2})² = {d.MaxBlock}; " +
                     $"{d.CubicSpeedup.ToString("0.0", Inv)}× cubic-cost speedup over the dense (4^N)³.",
            children: kids,
            provenance: NodeProvenance.Live);
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
                         "boundary mode; at UNIFORM γ the whole sector sits at Re=−2γ (the Absorption floor), " +
                         "and a per-site γ opens that line into an interval around −2·γ̄ (the site-resolved " +
                         "node below). Its {0,2} junction at N≥6 crosses into the (2,2) sector."),
            new InspectableNode("(1,1) — single-excitation: the {0,2}-coherence, two regimes",
                summary: "the single-excitation populations (n_diff=0) + coherences (n_diff=2) — ONE sector, two " +
                         "regimes: inspect --root ceiling (StructuralCeilingWitness) reads the HIGH-Q dark [H,A]=0 " +
                         "commutant (g2=strict_gap/2γ, the S_N standard rep, g2(K_N)=4/N, star 4/(N−1)); inspect " +
                         "--root horizon (CoherenceHorizonWitness) reads the LOW-Q √-EP Q*(N) where it stops " +
                         "oscillating (the Haken-Strobl reduction 4^N→N², = the carbon Frost-Hückel threshold). " +
                         "inspect --root secondclock stitches the two regimes."),
            new InspectableNode("(2,2)/(p,p) — half-filling: a DISTINCT {0,2}-coherence",
                summary: "NOT the (1,1) mode — same n_diff∈{0,2} histogram, but the TWO-excitation filling sector " +
                         "(the V-Effect seam). inspect --root survivor (IncompletenessSurvivorWitness): the " +
                         "longest-lived interior C=0.5 coherence (ring (2,2)/(N−2,N−2), chain dead-centre). The " +
                         "reduction's {0,2} junction (N≥6) and secondclock's N=4 anomaly (ceiling's K_4 = 2−2/√3, " +
                         "ring-4 = 1) live here, not in (1,1)."),
        };
        return new InspectableNode("the sector map — which live witness zooms each load-bearing sector",
            summary: "the per-sector overview indexes the sector-specific witnesses; each is a max-zoom on one " +
                     "sector of THIS decomposition. (0,1) band edge → reduction; (1,1) single-excitation → ceiling " +
                     "(high-Q) + horizon (low-Q EP), stitched by secondclock; (2,2)/(p,p) half-filling → survivor. " +
                     "The (1,1) and (2,2) {0,2}-coherences are DISTINCT modes (same n_diff, different filling) — do " +
                     "not conflate them.",
            children: kids);
    }

    private InspectableNode ThePalindromeNode(ComplexMatrix h)
    {
        var (spectrum, used, skipped, full) = ReconstructSpectrum(h, N, Gamma, LiveEigCap);
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
            children: kids,
            provenance: NodeProvenance.Live);
    }

    private InspectableNode TheAbsorptionFloorNode(ComplexMatrix h)
    {
        var (minRe, maxRe) = BandEdgeSectorReSpan(h, N, Gamma);
        return new InspectableNode("the per-sector Absorption floor (Re = −2γ)",
            summary: $"at UNIFORM γ, the (0,1) band-edge sector, the |1-exc⟩⟨vac| coherences, {N}-dim, sits entirely at " +
                     $"Re ∈ [{minRe.ToString("0.0000", Inv)}, {maxRe.ToString("0.0000", Inv)}] = −2γ = " +
                     $"{(-2 * Gamma).ToString("0.###", Inv)} (every coherence disagrees in one bit, so L_D = −2γ·I " +
                     "is scalar there; the F50 weight-1 floor / Absorption Theorem). The decay GAP lives instead " +
                     "in the diagonal (k,k) sectors.",
            provenance: NodeProvenance.Live);
    }

    private InspectableNode TheSiteResolvedBandEdgeNode(ComplexMatrix h)
    {
        var bonds = ChainBonds(N, J);
        var g = DeepEdgeProfile(N, Gamma);
        var (minRe, maxRe) = BandEdgeSectorReSpan(h, N, g);
        double residual = GeneratorResidual(N, bonds, g);
        double hermResidual = HermitianPartResidual(N, bonds, g);
        double perMode = PerModeAbsorptionResidual(N, bonds, g);
        double gMin = g.Min(), gMax = g.Max(), gBar = g.Average();
        var m = BandEdgeSectorBlock(h, N, g);
        double traceRe = m.Evd().EigenValues.Sum(z => z.Real);
        // read the trace a second time WITHOUT the eigensolver, so the identity below is not just
        // the eigensolver agreeing with itself
        double traceDirect = Enumerable.Range(0, N).Sum(i => m[i, i].Real);
        double width = maxRe - minRe;
        var kids = new List<IInspectable>
        {
            new InspectableNode("the generator, entry-wise",
                summary: $"M = +i·(J/2)·𝓛 − 2·diag(γ) on the open chain (𝓛 = diag(deg) − A): residual " +
                         $"{residual.ToString("E3", Inv)}. The composite neither owner carries: D10 has the ZZ " +
                         "degree term at uniform γ, VacuumBlockReductionClaim has the γ profile with an XY H " +
                         "(no degree term). XXX only: for XXZ the ZZ coefficient no longer matches the XY one " +
                         "and Δ·diag(deg) − A is not a Laplacian. THE BLOCK GOES BY FOUR NAMES and two of them " +
                         "are the same label on conjugate operators: (0,1) here is (p_c, p_r) = (bra, ket) = " +
                         "|1-exc⟩⟨vac| and carries +i; D10:139's −2iJ𝓛 is the conjugate |vac⟩⟨1-exc|, which " +
                         "D10:112 ALSO calls (0,1); VacuumBlockReductionClaim writes L_(1,0); and " +
                         "simulations/birth_canal_vacuum_block_verifier.py calls it (1,0). Before reading a sign " +
                         "conflict, separate the two axes: ORIENTATION flips i, and 𝓛 = D − A versus the bare " +
                         "hopping h differs on the DEGREE DIAGONAL only. VacuumBlockReduction's −iQ·h and this " +
                         "+i·(J/2)·𝓛 agree off the diagonal; the ZZ degree term is the whole difference, not a " +
                         "sign disagreement."),
            new InspectableNode("the profile",
                summary: $"deep-edge, γ = [{string.Join(", ", g.Select(x => x.ToString("0.###", Inv)))}] by SITE " +
                         "(SectorReductionWitness's general-N deep-edge shape, the sibling of its N=5 canal " +
                         "anchor, scaled to mean γ); " +
                         $"Σγ = σ = {(N * Gamma).ToString("0.###", Inv)} unchanged, so the F1 palindrome center " +
                         "−σ does not move. The block's own basis order is site N−1 down to site 0 " +
                         $"([{string.Join(", ", BandEdgeSectorSiteOrder(N))}]), not site order: any site-indexed " +
                         "quantity must be permuted through it before an entry-wise comparison. What that " +
                         "permutation bites on is the γ diagonal; a uniform-J chain Laplacian is " +
                         "reversal-invariant and cannot see it, which is why the gates also use a " +
                         "non-palindromic per-bond J."),
            new InspectableNode("the per-mode absorption law (the master statement)",
                summary: $"−Re λ_k = 2·Σ_l γ_l·⟨Δ_l⟩_k with ⟨Δ_l⟩_k the mode's site occupancy: residual " +
                         $"{perMode.ToString("E3", Inv)}. This is AbsorptionTheoremClaim's per-channel law " +
                         "(\"the carrier is a vector\"), which is this block's typed parent, not a new statement. " +
                         "The mode-by-mode equality does NOT die under a profile: it becomes a γ-weighted average " +
                         "over the mode's own occupancy. Non-normality costs nothing, because Re λ is the Rayleigh " +
                         $"quotient of Herm(M) = −2·diag(γ) (residual {hermResidual.ToString("E3", Inv)}) on the " +
                         "right eigenvector, eigenvector by eigenvector."),
            new InspectableNode("Bendixson bracket (corollary)",
                summary: $"an occupancy is a probability distribution over sites, so Re λ is a CONVEX COMBINATION " +
                         $"of the −2γ_l and hence lies in [−2·max γ, −2·min γ] = " +
                         $"[{(-2 * gMax).ToString("0.0000", Inv)}, {(-2 * gMin).ToString("0.0000", Inv)}]; " +
                         $"measured span [{minRe.ToString("0.0000", Inv)}, {maxRe.ToString("0.0000", Inv)}]. The " +
                         "same lemma the F89 block lattice uses with n_diff in place of γ."),
            new InspectableNode("the trace (corollary)",
                summary: $"Σ Re λ = {traceRe.ToString("0.000000", Inv)} = Re tr(M) read off the diagonal without " +
                         $"an eigensolver ({traceDirect.ToString("0.000000", Inv)}) = −2σ = " +
                         $"{(-2 * N * Gamma).ToString("0.000", Inv)}, i.e. MEAN Re = −2·γ̄ = " +
                         $"{(-2 * gBar).ToString("0.000", Inv)}, because Re tr(M) = tr(Herm(M)) = −2σ. Note what " +
                         "this is NOT: the per-site occupancies do not close to 1 across modes (measured 1.08 and " +
                         "0.92 at N=6), because a right-eigenvector basis of a non-normal M is not orthonormal. " +
                         "The trace route is the correct one; the occupancy route is not available here. And the " +
                         "object is THIS BLOCK's eigenvalue mean, not the whole spectrum's exactly-Re=−2γ̄ SET, " +
                         "which a profile can empty outright (PROOF_RING_GAP_DOMINANCE's scope section: one site " +
                         "detuned leaves it empty at N=4 and N=5). Do not read the one as rescuing the other."),
        };
        return new InspectableNode("the site-resolved band edge (a γ profile)",
            summary: $"the same (0,1) sector with a per-site γ instead of one number: the generator is still " +
                     $"exactly M = +i·(J/2)·𝓛 − 2·diag(γ) (residual {residual.ToString("E3", Inv)}), but the " +
                     $"dissipator is now diagonal rather than scalar, so the Re = −2γ line " +
                     (width > 1e-9
                        ? $"has opened into the interval [{minRe.ToString("0.0000", Inv)}, {maxRe.ToString("0.0000", Inv)}]"
                        : $"has NOT opened here; the measured span is still a point at {minRe.ToString("0.0000", Inv)} " +
                          "(at N=2 the 2×2 block stays closed for J ≥ 2·γ̄_diff)") +
                     $", bracketed by Bendixson at [{(-2 * gMax).ToString("0.0000", Inv)}, {(-2 * gMin).ToString("0.0000", Inv)}] " +
                     $"and with its MEAN pinned by the trace at −2·γ̄ = {(-2 * gBar).ToString("0.000", Inv)} (the mean " +
                     "of the eigenvalues, NOT the midpoint of the bracket, which is a different number). Both are corollaries " +
                     "of the per-mode law below, which is the Absorption Theorem's own per-channel reading. " +
                     "diag(γ) and 𝓛 do " +
                     "not commute, so M is non-normal and its eigenvalues are NOT −2γ_j − i·(J/2)·μ_m: the slow " +
                     "mode has to be diagonalized for, not read off (arc site_resolved_vacuum_block).",
            children: kids,
            provenance: NodeProvenance.Live);
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
        "Absorption floor Re = −2γ AT UNIFORM γ, and the same sector site-resolved (a γ profile: the generator " +
        "stays exactly +i·(1/2)·𝓛 − 2·diag(γ) on any XXX bond graph, and the rate law is the Absorption Theorem's " +
        "per-channel one, Re λ_k = the γ-weighted average of the mode's site occupancy), and the N=9 " +
        "banked headline read live from chain_N9.json. The browsable " +
        "overview face of the SLOW_N9 result (arc block_spectrum_n9).";

    public IEnumerable<IInspectable> Children
    {
        get
        {
            yield return TheDecompositionNode();
            yield return TheSectorMapNode();
            // Build the Hilbert-space H (2^N × 2^N) once and share it between the two spectrum nodes.
            var h = HeisenbergChain(N, J);
            yield return ThePalindromeNode(h);
            yield return TheAbsorptionFloorNode(h);
            yield return TheSiteResolvedBandEdgeNode(h);
            yield return TheBankedN9Node();
        }
    }

    public InspectablePayload Payload => InspectablePayload.Empty;
}
