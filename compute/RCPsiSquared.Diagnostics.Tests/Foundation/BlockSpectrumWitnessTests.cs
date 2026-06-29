using System;
using System.Linq;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Diagnostics.Foundation;
using Xunit;

namespace RCPsiSquared.Diagnostics.Tests.Foundation;

/// <summary>Gates for the live block-spectrum witness (<c>inspect --root blockspectrum</c>),
/// which surfaces the N=9 per-joint-popcount-sector Liouvillian spectra the SLOW_N9 test banks.
/// Gate-first: the spectrum reconstruction (Gate 2) and the Absorption floor (Gate 3) are real
/// physics checks that fire if the sector-flat-index wiring or the Hamiltonian is wrong; Gate 4
/// pins the witness's banked-headline read to the committed chain_N9.json so it cannot drift.</summary>
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

    // Render smoke: the default (N=6) witness renders five nodes without crashing, and the live
    // reconstruction node at N=6 is full (max block C(6,3)^2 = 400 < the 2048 cap).
    [Fact]
    public void Witness_RendersFiveNodes_WithoutCrash()
    {
        var w = new BlockSpectrumWitness();   // default N=6, gamma=0.5, J=1
        Assert.Contains("N=6", w.DisplayName);
        Assert.False(string.IsNullOrWhiteSpace(w.Summary));

        var children = w.Children.ToList();
        Assert.Equal(5, children.Count);
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
