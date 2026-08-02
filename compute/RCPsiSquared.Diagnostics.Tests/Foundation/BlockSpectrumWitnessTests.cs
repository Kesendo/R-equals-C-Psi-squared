using System;
using System.Collections.Generic;
using System.Linq;
using RCPsiSquared.Core.ChainSystems;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Diagnostics.Foundation;
using Xunit;

namespace RCPsiSquared.Diagnostics.Tests.Foundation;

/// <summary>Gates for the live block-spectrum witness (<c>inspect --root blockspectrum</c>),
/// which surfaces the N=9 per-joint-popcount-sector Liouvillian spectra the SLOW_N9 test banks.
/// Gate-first: the spectrum reconstruction (Gate 2) and the Absorption floor (Gate 3) are real
/// physics checks that fire if the sector-flat-index wiring or the Hamiltonian is wrong; Gate 4
/// pins the witness's banked-headline read to the committed chain_N9.json so it cannot drift.
///
/// <para>The site-resolved half (Gates 3b..3f) covers the (0,1) block under a per-site γ PROFILE:
/// the generator identity entry-wise on four graph families (3b) and its J-scaling (3b-scale), the
/// basis-order permutation with its mutation control (3b'), Herm(M) = −2·diag(γ) directly (3c), the
/// per-mode Absorption law (3d), its Bendixson and trace corollaries (3e), and the N=2 case where
/// the floor line does NOT open (3f). The graph rows are the point: on a uniform-J chain or ring a
/// global site-versus-bit reversal cancels, so the per-bond-J and star rows are what make the
/// ordering falsifiable.</para></summary>
public class BlockSpectrumWitnessTests
{
    [Fact]
    public void Children_Provenance_LiveReconstructionVsBankedHeadline()
    {
        var w = new BlockSpectrumWitness(n: 5);   // small N -> live reconstruction is cheap
        var kids = ((IInspectable)w).Children.ToList();
        IInspectable Find(string namePart) => kids.Single(c => c.DisplayName.Contains(namePart));

        // Live-computed children (the badge must not contradict their "reconstructed live" prose).
        Assert.Equal(NodeProvenance.Live, Find("decomposition").Provenance);
        Assert.Equal(NodeProvenance.Live, Find("palindrome").Provenance);
        Assert.Equal(NodeProvenance.Live, Find("Absorption floor").Provenance);
        Assert.Equal(NodeProvenance.Live, Find("site-resolved band edge").Provenance);

        // Stored children: the static sector-map prose, and the banked N=9 headline (chain_N9.json).
        Assert.Equal(NodeProvenance.Stored, Find("sector map").Provenance);
        Assert.Equal(NodeProvenance.Stored, Find("banked").Provenance);
    }

    // Gate 1: the joint-popcount decomposition counts (the 100 -> 50 -> 25 story at N=9), pure
    // combinatorial. X(x)N order-2 pairing = the JSON PrimarySectorCount; the F1 Pi order-4 orbit
    // = the eig-calls the compute path actually does (Pi^2 = X(x)N).
    [Theory]
    [InlineData(8, 81, 41, 21, 4900L, 4, 4)]
    [InlineData(9, 100, 50, 25, 15876L, 4, 4)]
    public void Decomposition_CountsMatchTheKnownSectorStructure(
        int n, int sectorCount, int xnClasses, int piClasses, long maxBlock, int maxPc, int maxPr)
    {
        var d = BlockSpectrumWitness.Decomposition(n);
        Assert.Equal(sectorCount, d.SectorCount);
        Assert.Equal(xnClasses, d.XnClasses);
        Assert.Equal(piClasses, d.PiOrbitClasses);
        Assert.Equal(maxBlock, d.MaxBlock);
        Assert.Equal((maxPc, maxPr), (d.MaxPc, d.MaxPr));
    }

    // Gate 1b: the live cubic-cost speedup reproduces the banked N=9 EffectiveSpeedupOverDense
    // (645.95x) -- computed the identical (4^N)^3 / sum(block^3) way, so it must match.
    [Fact]
    public void Decomposition_CubicSpeedup_ReproducesTheBankedN9Value()
    {
        double speedup = BlockSpectrumWitness.Decomposition(9).CubicSpeedup;
        Assert.True(Math.Abs(speedup - 645.9495725858463) < 0.01,
            $"live cubic speedup {speedup} should reproduce the banked 645.95");
    }

    // Gate 2: the witness reconstructs the full spectrum sector-by-sector (every block via
    // BuildBlockZ) and it obeys F1: the multiset is symmetric about the center -sigma. A real
    // gate -- a wrong sector-flat-index extraction or a wrong H gives a spectrum that fails this.
    [Fact]
    public void ReconstructSpectrum_AtN5_IsFull_ObeysTheF1Palindrome_AboutMinusTwoSigma()
    {
        const int n = 5;
        const double gamma = 0.5, j = 1.0;
        var (spectrum, _, skipped, full) = BlockSpectrumWitness.ReconstructSpectrum(n, gamma, j, cap: 2048);
        Assert.True(full);
        Assert.Equal(0, skipped);
        Assert.Equal(1 << (2 * n), spectrum.Length);   // 4^5 = 1024 (sector blocks sum to 4^N)

        double sigma = n * gamma;                       // 2.5
        // F1: {lambda} = {-2sigma - lambda}; the symmetric multiset has ~0 set-distance to its reflection.
        double pairing = BlockSpectrumWitness.PalindromePairingDistance(spectrum, sigma);
        Assert.True(pairing < 1e-9, $"F1 symmetry distance {pairing:E3} should be ~0");
        // the F1 floor: the maximally-dephased coherence (all N bits disagree) sits at Re = -2sigma.
        Assert.Equal(-2.0 * sigma, BlockSpectrumWitness.MinReal(spectrum), 9);
        // F4 kernel: a connected chain has kernel dimension N+1.
        Assert.Equal(n + 1, BlockSpectrumWitness.KernelDimension(spectrum));
    }

    // Gate 2b: the F1 metric must be MULTISET-aware, not set-aware. {+1 x3, -1 x1} is closed under
    // lambda -> -lambda as a SET (sigma=0) but NOT as a multiset (reflected = {-1 x3, +1 x1}). A
    // set-distance metric blindly reports ~0; the real F1 check (greedy NN with removal) must see the
    // gap -- the surplus +1 can only pair to a -1, distance 2. This is the failure mode a reconstruction
    // bug (a dropped/duplicated eigenvalue with a same-valued neighbour) would produce.
    [Fact]
    public void PalindromePairingDistance_DetectsAMultiplicityDefect()
    {
        var defective = new[]
        {
            new System.Numerics.Complex(1, 0), new System.Numerics.Complex(1, 0),
            new System.Numerics.Complex(1, 0), new System.Numerics.Complex(-1, 0),
        };
        double d = BlockSpectrumWitness.PalindromePairingDistance(defective, sigma: 0.0);
        Assert.True(d > 1.0, $"a multiplicity defect must be visible; a set-blind metric returns ~0, got {d:E3}");
    }

    // Gate 2c: the physics-review's STRONGEST counterexample, on a REAL spectrum (not a toy multiset).
    // Take a genuine F1-symmetric N=5 reconstruction, drop one eigenvalue and duplicate a different
    // (non-axis) one -- the exact defect a reconstruction-wiring bug produces. The old set/Hausdorff
    // metric gave ~3e-14 (invisible); the multiplicity-aware metric must flag it large.
    [Fact]
    public void PalindromePairingDistance_FlagsADropAndDuplicateInARealSpectrum()
    {
        const int n = 5;
        const double gamma = 0.5, j = 1.0;
        double sigma = n * gamma;
        var (clean, _, _, _) = BlockSpectrumWitness.ReconstructSpectrum(n, gamma, j, 2048);
        Assert.True(BlockSpectrumWitness.PalindromePairingDistance(clean, sigma) < 1e-9, "clean spectrum: F1 holds");

        // duplicate a non-axis eigenvalue (Re away from -sigma) over a different-valued slot:
        // removes one value, adds a copy of another -> breaks multiset F1 closure.
        int dup = System.Array.FindIndex(clean, z => System.Math.Abs(z.Real - (-sigma)) > 0.5);
        int drop = System.Array.FindIndex(clean, z => z != clean[dup]);
        Assert.True(dup >= 0 && drop >= 0 && clean[dup] != clean[drop], "test setup: distinct non-axis pair exists");
        var defective = (System.Numerics.Complex[])clean.Clone();
        defective[drop] = clean[dup];

        double d = BlockSpectrumWitness.PalindromePairingDistance(defective, sigma);
        Assert.True(d > 1e-6, $"a drop+duplicate defect must be visible; the old set-blind metric gave ~3e-14, got {d:E3}");
    }

    // Gate 3: the (0,1) band-edge sector is entirely at Re = -2gamma. On (p_c=0, p_r=1) every
    // basis coherence disagrees in exactly one bit, so L_D = -2gamma*I (scalar) and L_H restricted
    // is anti-Hermitian -> every eigenvalue has Re = -2gamma exactly (uniform gamma; Absorption).
    [Fact]
    public void BandEdgeSector_01_SitsEntirelyAtMinusTwoGamma()
    {
        const int n = 5;
        const double gamma = 0.5, j = 1.0;
        var (minRe, maxRe) = BlockSpectrumWitness.BandEdgeSectorReSpan(n, gamma, j);
        Assert.Equal(-2.0 * gamma, minRe, 9);
        Assert.Equal(-2.0 * gamma, maxRe, 9);
    }

    // The four graph families the site-resolved gates run on. Chain and ring at uniform J are
    // REVERSAL-INVARIANT, which is what lets a site<->bit mix-up hide; the per-bond J variants
    // below are not, and neither is the star.
    public static IEnumerable<object[]> BandEdgeGraphs()
    {
        foreach (int n in new[] { 3, 4, 5, 6, 7 })
        {
            yield return new object[] { n, "chain-uniform", Chain(n, _ => 1.0) };
            yield return new object[] { n, "chain-perbond", Chain(n, b => 0.4 + 0.9 * b) };   // non-palindromic
            yield return new object[] { n, "ring", Enumerable.Range(0, n).Select(i => new Bond(i, (i + 1) % n, 0.7)).ToArray() };
            yield return new object[] { n, "star", Enumerable.Range(1, n - 1).Select(i => new Bond(0, i, 1.3)).ToArray() };
        }
    }

    private static Bond[] Chain(int n, Func<int, double> j) =>
        Enumerable.Range(0, n - 1).Select(b => new Bond(b, b + 1, j(b))).ToArray();

    // The deliberately asymmetric profile: a reversed gamma cannot pass, and neither can a
    // site-ordered prediction against a block whose basis order is site n-1 down to site 0.
    private static double[] AsymmetricGamma(int n) =>
        Enumerable.Range(0, n).Select(l => 0.13 + 0.31 * l).ToArray();

    // Gate 3b: the SITE-RESOLVED (0,1) block on ANY XXX bond graph. The engine has always returned
    // it with a per-site gamma; the witness caller was uniform-only. The generator identity, in
    // this witness's convention H = sum_b (J_b/4)*(XX+YY+ZZ), is M = +i*(1/2)*WeightedLaplacian
    // - 2*diag(gamma), entry-wise. It fires on a wrong sign, a wrong J factor, a missing ZZ degree
    // term (which would leave the adjacency form, D10's recorded failure mode), a REVERSED gamma
    // profile (the failure PerBlockLiouvillianBuilder's own gamma-index-convention paragraph
    // records and PerBlockLiouvillianBuilderGammaOrderTests gates), and -- only because of the
    // per-bond and star rows -- a site<->bit mix-up in the Laplacian layout.
    //
    // The threshold is a SCALING LAW, not a constant: the residual is an absolute quantity built
    // from cancelling O(J*N) terms, so a fixed 1e-12 would silently become a real tolerance at
    // large J. Scale = J_max*N + 2*max(gamma), and the large-J row below is what proves the scaling
    // is the right one rather than the blind spot being moved one decade.
    [Theory]
    [MemberData(nameof(BandEdgeGraphs))]
    public void SiteResolvedBandEdgeBlock_IsExactlyTheLaplacianPlusDiagGammaGenerator(int n, string label, Bond[] bonds)
    {
        var gamma = AsymmetricGamma(n);
        double scale = bonds.Max(b => Math.Abs(b.Coupling)) * n + 2.0 * gamma.Max();
        double residual = BlockSpectrumWitness.GeneratorResidual(n, bonds, gamma);
        Assert.True(residual < 1e-13 * scale,
            $"[{label} N={n}] entry-wise generator residual {residual:E3} exceeds 1e-13*scale ({1e-13 * scale:E3})");
    }

    // Gate 3b-scale: the same identity at J = 1e5, where an absolute 1e-12 threshold would fail on
    // rounding alone. This is the row that turns the threshold above into a claim about scaling.
    [Theory]
    [InlineData(4, 1e5)]
    [InlineData(6, 1e5)]
    [InlineData(6, 1e-4)]
    public void SiteResolvedBandEdgeBlock_GeneratorIdentity_ScalesWithJ(int n, double j)
    {
        var gamma = AsymmetricGamma(n);
        var bonds = Chain(n, _ => j);
        double scale = j * n + 2.0 * gamma.Max();
        double residual = BlockSpectrumWitness.GeneratorResidual(n, bonds, gamma);
        Assert.True(residual < 1e-13 * scale,
            $"[N={n} J={j}] residual {residual:E3} exceeds 1e-13*scale ({1e-13 * scale:E3})");
    }

    // Gate 3b': the block's basis order is NOT site order. The sector's flat indices row*d + col
    // with row = 1<<b sort ascending, and site l sits at bit n-1-l, so the block runs site n-1 down
    // to site 0. Two halves, because the first alone is arithmetically forced by the implementation
    // and would pass even if the bit<->site semantics were inverted: (a) the permutation is
    // [n-1..0]; (b) a prediction laid out in SITE order must NOT match the block. Half (b) uses a
    // NON-PALINDROMIC per-bond J, so a global site<->bit reversal cannot cancel between the
    // Laplacian and the gamma diagonal the way it does on a uniform chain.
    [Theory]
    [InlineData(3)]
    [InlineData(5)]
    [InlineData(6)]
    public void BandEdgeSectorSiteOrder_RunsFromTheLastSiteDownToTheFirst(int n)
    {
        Assert.Equal(Enumerable.Range(0, n).Reverse().ToArray(), BlockSpectrumWitness.BandEdgeSectorSiteOrder(n));

        var gamma = AsymmetricGamma(n);
        var bonds = Chain(n, b => 0.4 + 0.9 * b);
        var m = BlockSpectrumWitness.BandEdgeSectorBlock(n, bonds, gamma);
        var lapBySite = new double[n, n];
        foreach (var b in bonds)
        {
            lapBySite[b.Site1, b.Site2] -= b.Coupling; lapBySite[b.Site2, b.Site1] -= b.Coupling;
            lapBySite[b.Site1, b.Site1] += b.Coupling; lapBySite[b.Site2, b.Site2] += b.Coupling;
        }
        double worstSiteOrdered = 0.0;
        for (int r = 0; r < n; r++)
            for (int c = 0; c < n; c++)
            {
                var predicted = System.Numerics.Complex.ImaginaryOne * 0.5 * lapBySite[r, c];
                if (r == c) predicted -= 2.0 * gamma[r];
                worstSiteOrdered = Math.Max(worstSiteOrdered, (m[r, c] - predicted).Magnitude);
            }
        Assert.True(worstSiteOrdered > 0.1,
            $"a site-ordered prediction must NOT match the block; got {worstSiteOrdered:E3} -- if this passes, " +
            "the basis-order permutation is a no-op here and Gate 3b is not testing it");
    }

    // Gate 3c: the sharp statement the profile case rests on, gated directly rather than through a
    // consequence with slack -- Herm(M) = -2*diag(gamma) EXACTLY, on every graph. Bendixson and the
    // trace are corollaries of this; asserting only them leaves a theorem standing in for a gate.
    // Threshold scaled like 3b, for the same reason: the entries are O(J*N), so a constant would
    // quietly become a real tolerance at large J even though the cancellation here is exact today.
    [Theory]
    [MemberData(nameof(BandEdgeGraphs))]
    public void SiteResolvedBandEdge_HermitianPartIsExactlyMinusTwoDiagGamma(int n, string label, Bond[] bonds)
    {
        var gamma = AsymmetricGamma(n);
        double scale = bonds.Max(b => Math.Abs(b.Coupling)) * n + 2.0 * gamma.Max();
        double residual = BlockSpectrumWitness.HermitianPartResidual(n, bonds, gamma);
        Assert.True(residual < 1e-14 * scale,
            $"[{label} N={n}] Herm(M) + 2*diag(gamma) residual {residual:E3} exceeds 1e-14*scale ({1e-14 * scale:E3})");
    }

    // Gate 3c-scale: the row that earns 3c's scaled threshold. Without it, scaling by J*N + 2*max
    // gamma is a 4.5x to 38x LOOSENING on every case actually run, justified by a large-J argument
    // no row exercises. Here J = 1e5, where a flat 1e-14 would fail on rounding alone.
    [Theory]
    [InlineData(4, 1e5)]
    [InlineData(6, 1e5)]
    public void SiteResolvedBandEdge_HermitianPart_ScalesWithJ(int n, double j)
    {
        var gamma = AsymmetricGamma(n);
        var bonds = Chain(n, _ => j);
        double scale = j * n + 2.0 * gamma.Max();
        double residual = BlockSpectrumWitness.HermitianPartResidual(n, bonds, gamma);
        Assert.True(residual < 1e-14 * scale,
            $"[N={n} J={j}] Herm residual {residual:E3} exceeds 1e-14*scale ({1e-14 * scale:E3})");
    }

    // Gate 3d: the MASTER statement, and it belongs to AbsorptionTheoremClaim, not to this witness:
    // -Re(lambda_k) = 2*sum_l gamma_l*<Delta_l>_k with <Delta_l>_k the mode's site occupancy, on the
    // RIGHT eigenvectors of a non-normal M.
    //
    // WHAT THIS CAN AND CANNOT CATCH, stated because it would otherwise read stronger than it is:
    // given Gate 3c, this is a linear-algebra identity (Re lambda is the Rayleigh quotient of
    // Herm(M) on an eigenvector), so it CANNOT fail independently -- the only things that fire it
    // are a broken eigensolver or the same gamma-index error 3c already catches. It is kept as the
    // cross-check that ties this block to its typed parent's per-channel law, and as the executable
    // form of "non-normality costs nothing here", not as independent evidence.
    //
    // Its threshold is scaled like 3b/3c, but note the scale is only a PROXY here: this residual is
    // an eigensolver quantity on a deliberately non-normal M, so the governing size is
    // eps*||M||*cond(V), not ||M||. J*N + 2*max gamma tracks ||M|| and leaves the conditioning out,
    // which is why this gate is given a decade more room than 3c rather than the same.
    [Theory]
    [MemberData(nameof(BandEdgeGraphs))]
    public void SiteResolvedBandEdge_ObeysThePerModeAbsorptionLaw(int n, string label, Bond[] bonds)
    {
        var gamma = AsymmetricGamma(n);
        double scale = bonds.Max(b => Math.Abs(b.Coupling)) * n + 2.0 * gamma.Max();
        double residual = BlockSpectrumWitness.PerModeAbsorptionResidual(n, bonds, gamma);
        Assert.True(residual < 1e-13 * scale,
            $"[{label} N={n}] per-mode absorption residual {residual:E3} exceeds 1e-13*scale ({1e-13 * scale:E3})");
    }

    // Gate 3e: what the profile does to the floor, as the two corollaries. Bendixson brackets every
    // Re in [-2*max gamma, -2*min gamma] and the trace pins the MEAN at -2*gamma_bar; the deep-edge
    // profile keeps sigma fixed, so the F1 center -sigma does not move either. N >= 3: at N=2 the
    // block is 2x2 and the line stays CLOSED for J >= 2*(gamma difference), so the "has opened"
    // half below is false there -- which is why the witness renders that sentence conditionally.
    [Theory]
    [InlineData(3)]
    [InlineData(4)]
    [InlineData(5)]
    [InlineData(6)]
    public void SiteResolvedBandEdge_BendixsonBracketAndTraceHold_WhileTheLineOpens(int n)
    {
        const double gamma = 0.5, j = 1.0;
        var g = BlockSpectrumWitness.DeepEdgeProfile(n, gamma);
        Assert.Equal(n * gamma, g.Sum(), 9);                       // sigma unchanged: the F1 center holds

        var (minRe, maxRe) = BlockSpectrumWitness.BandEdgeSectorReSpan(n, g, j);

        // Bendixson: the Hermitian part is exactly -2*diag(gamma), so no Re escapes its spectrum.
        Assert.True(minRe >= -2.0 * g.Max() - 1e-9, $"minRe {minRe} below the Bendixson floor {-2 * g.Max()}");
        Assert.True(maxRe <= -2.0 * g.Min() + 1e-9, $"maxRe {maxRe} above the Bendixson ceiling {-2 * g.Min()}");

        // the trace: sum of Re = -2*sigma exactly (the Absorption floor, in the mean).
        var m = BlockSpectrumWitness.BandEdgeSectorBlock(n, g, j);
        double traceRe = m.Evd().EigenValues.Sum(z => z.Real);
        Assert.Equal(-2.0 * n * gamma, traceRe, 9);

        // and the line really has opened -- otherwise the two assertions above are vacuous.
        Assert.True(maxRe - minRe > 1e-3, $"the profile must split the floor line; span width {maxRe - minRe:E3}");

        // the uniform case, same sigma, is the collapsed limit: span width 0.
        var (uMin, uMax) = BlockSpectrumWitness.BandEdgeSectorReSpan(n, gamma, j);
        Assert.Equal(-2.0 * gamma, uMin, 9);
        Assert.Equal(uMin, uMax, 9);
    }

    // Gate 3f: the N=2 exception the prose must not paper over. The 2x2 block splits by
    // 2*sqrt(d^2 - J^2/4) with d the half-difference of the two rates, so for J >= 2d the two
    // eigenvalues have EQUAL Re and the line does not open at all. So "a profile opens the line" is
    // false in general, and the witness renders that sentence conditionally on the measured width.
    [Theory]
    [InlineData(1.0)]            // below 2d = 1.5: the line opens
    [InlineData(1.4)]
    [InlineData(1.5)]            // exactly at the coalescence
    [InlineData(2.5)]            // above: closed
    public void SiteResolvedBandEdge_AtNEqualsTwo_TheLineNeedNotOpen(double j)
    {
        const int n = 2;
        const double gamma = 0.5;
        var g = BlockSpectrumWitness.DeepEdgeProfile(n, gamma);
        Assert.Equal(n * gamma, g.Sum(), 9);
        double d = g[1] - g[0];                                    // 0.75 at gamma = 0.5
        var (minRe, maxRe) = BlockSpectrumWitness.BandEdgeSectorReSpan(n, g, j);
        Assert.Equal(-2.0 * n * gamma, minRe + maxRe, 9);          // the trace holds either way
        double expectedWidth = j < 2.0 * d ? 2.0 * Math.Sqrt(d * d - j * j / 4.0) : 0.0;
        // J = 2d is an exact coalescence (a defective 2x2), where the numerically resolved split
        // goes as sqrt(eps), not eps -- so the tolerance there is sqrt-scaled, not loosened by
        // taste. Away from it the identity is exact to 1e-9.
        double tol = Math.Abs(j - 2.0 * d) < 1e-12 ? 1e-6 : 1e-9;
        Assert.True(Math.Abs((maxRe - minRe) - expectedWidth) < tol,
            $"J={j}: width {maxRe - minRe:E6} vs predicted {expectedWidth:E6}");
    }

    // The deep-edge profile is the sibling witness's canal shape, scaled to mean gamma -- an
    // attribution, so it is gated: at gamma=1 it must reproduce SectorReductionWitness's profile
    // entry for entry, including the N=5 canal anchor [0.25, 1.5, 1.5, 1.5, 0.25].
    [Fact]
    public void DeepEdgeProfile_ReproducesTheSiblingCanalShapeAtGammaOne()
    {
        Assert.Equal(new[] { 0.25, 1.5, 1.5, 1.5, 0.25 }, BlockSpectrumWitness.DeepEdgeProfile(5, 1.0));
        for (int n = 3; n <= 8; n++)
        {
            double edge = 0.25, rest = (n - 2 * edge) / (n - 2);
            var expected = Enumerable.Range(0, n).Select(i => (i == 0 || i == n - 1) ? edge : rest).ToArray();
            Assert.Equal(expected, BlockSpectrumWitness.DeepEdgeProfile(n, 1.0));
        }
        // and the invariant the profile exists for: sum = N*gamma at every N >= 2.
        for (int n = 2; n <= 9; n++)
            Assert.Equal(n * 0.37, BlockSpectrumWitness.DeepEdgeProfile(n, 0.37).Sum(), 9);
        Assert.Throws<ArgumentOutOfRangeException>(() => BlockSpectrumWitness.DeepEdgeProfile(1, 0.5));
    }

    // Gate 4: the banked N=9 headline node reads the committed chain_N9.json and matches it. This
    // pins the witness's "stored" provenance to the artifact so it cannot silently drift.
    [Fact]
    public void ReadBankedN9_MatchesTheCommittedArtifact()
    {
        var b = BlockSpectrumWitness.ReadBankedN9();
        Assert.NotNull(b);
        Assert.Equal(262144, b!.SpectrumSize);
        Assert.Equal(-9.0, b.MinReal, 6);                 // = -2sigma, the F1 floor at N=9 (sigma=4.5)
        Assert.Equal(10, b.KernelDimension);              // = N+1, the F4 connected-chain kernel
        Assert.Equal(0.02727562511208863, b.DissipationGap, 9);
        Assert.Equal(5.736321706379341, b.MaxImag, 6);
        Assert.True(b.MaxPairingDistance < 1e-12);        // 3.48e-13: the F1 pairing held bit-exact
        Assert.Equal(0, b.OutlierPairCount);
        Assert.Equal(100, b.SectorCount);
        Assert.Equal(50, b.PrimarySectorCount);           // X(x)N order-2 classes
        Assert.Equal(15876, b.MaxBlockSize);
        Assert.True(Math.Abs(b.EffectiveSpeedup - 645.9495725858463) < 0.01);
    }

    // Render smoke: the default (N=6) witness renders six nodes without crashing, and the live
    // reconstruction node at N=6 is full (max block C(6,3)^2 = 400 < the 2048 cap).
    [Fact]
    public void Witness_RendersSixNodes_WithoutCrash()
    {
        var w = new BlockSpectrumWitness();   // default N=6, gamma=0.5, J=1
        Assert.Contains("N=6", w.DisplayName);
        Assert.False(string.IsNullOrWhiteSpace(w.Summary));

        var children = w.Children.ToList();
        Assert.Equal(6, children.Count);
        Assert.All(children, c => Assert.False(string.IsNullOrWhiteSpace(c.DisplayName)));
        Assert.All(children, c => Assert.False(string.IsNullOrWhiteSpace(c.Summary)));
    }

    // The sector map turns blockspectrum into the navigation hub: each load-bearing sector points
    // to the sector-specific witness(es) that zoom it. The five roots must all be reachable from it.
    [Fact]
    public void SectorMap_PointsToTheSectorSpecificWitnesses()
    {
        var w = new BlockSpectrumWitness();
        var map = w.Children.Single(c => c.DisplayName.Contains("sector map"));
        var text = string.Join(" | ", map.Children.Select(c => $"{c.DisplayName} {c.Summary}"));
        foreach (var root in new[] { "reduction", "ceiling", "horizon", "survivor", "secondclock" })
            Assert.Contains(root, text);
    }

    // The reverse leg of the hub: each sector-specific witness points back to blockspectrum (the
    // per-sector overview), so the navigation is bidirectional.
    [Fact]
    public void TheFiveSectorWitnesses_PointBackToBlockSpectrum()
    {
        var summaries = new[]
        {
            new SectorReductionWitness(n: 5).Summary,
            new StructuralCeilingWitness().Summary,
            new CoherenceHorizonWitness().Summary,
            new IncompletenessSurvivorWitness(6, 1.5).Summary,
            new SecondClockRegimeWitness().Summary,
        };
        Assert.All(summaries, s => Assert.Contains("blockspectrum", s));
    }
}
