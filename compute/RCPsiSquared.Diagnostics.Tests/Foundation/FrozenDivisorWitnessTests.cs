using RCPsiSquared.Diagnostics.Foundation;
using Xunit;

namespace RCPsiSquared.Diagnostics.Tests.Foundation;

/// <summary>Gate for <see cref="FrozenDivisorWitness"/>, the live lab of the R90 frozen divisor
/// (F140, proof <c>docs/proofs/PROOF_R90_FROZEN_DIVISOR.md</c>): on the anti-palindromic dephasing
/// locus of F91, the single-excitation corner block holds λ = −4γ̄ with multiplicity at least
/// ⌊N/2⌋ at every coupling, one frozen mode per balanced pair, and exactly ⌊N/2⌋ away from the
/// finitely many exceptional couplings (of which J = 0, where it doubles, is one). Every read here
/// is on the taxed stratum γ̄ &gt; 0, which the witness enforces; at γ̄ = 0 the count is N and the
/// stratum is pinned by the gate's G16 and MirrorWorld's <c>Divisor</c> instead.
///
/// <para>Two-sided throughout: every ⌊N/2⌋ read sits beside a zero that a broken rank routine, a
/// mislabelled root or a lost locus condition could not produce, plus the J = 0 corollary
/// 2⌊N/2⌋, which is neither.</para></summary>
public class FrozenDivisorWitnessTests
{
    [Theory]
    [InlineData(3, 1)]
    [InlineData(4, 2)]
    [InlineData(5, 2)]
    [InlineData(6, 3)]
    [InlineData(7, 3)]
    public void CornerBlock_FreezesFloorNOverTwo_AtEveryCoupling(int n, int expected)
    {
        var w = new FrozenDivisorWitness(n);

        Assert.Equal(expected, w.PredictedFrozen);
        Assert.Equal(2 * expected, w.FixedCells);          // τQ's anti-diagonal fixed cells
        Assert.Equal(expected, w.Tax);                     // dim D₋, even under τQ, eats half

        Assert.Equal(expected, w.KernelDimension);
        Assert.Equal(expected, w.SecondCouplingKernelDimension);
        Assert.Equal(expected, w.SecondProfileKernelDimension);
        Assert.True(w.Holds, $"N={n}: {w.Summary}");
    }

    [Theory]
    [InlineData(3)]
    [InlineData(4)]
    [InlineData(5)]
    [InlineData(6)]
    [InlineData(7)]
    public void TheFalsifiers_AreZero(int n)
    {
        var w = new FrozenDivisorWitness(n);

        // The value is the value: one grid step off the root carries nothing.
        Assert.Equal(0, w.WrongRootKernelDimension);

        // Off the locus at UNCHANGED σ and γ̄: partial balance yields nothing.
        Assert.Equal(0, w.OffLocusKernelDimension);
        Assert.Equal(w.Gamma.Sum(), w.OffLocusGamma.Sum(), 12);

        // A non-corner block carries neither root, on the Heisenberg chain, which is the
        // default and the chain the proof document's census was verified on.
        if (w.ControlBlock is { } cb)
        {
            Assert.Equal(0, cb.AtRoot);
            Assert.Equal(0, cb.AtFoldRoot);
        }
    }

    /// <summary>The confinement to the four corners is a HEISENBERG statement. Drop the ZZ
    /// diagonal and the same root is carried, at the same multiplicity, by blocks that are not
    /// corners: the bound and the fold parity survive the change of chain, "only the corners" does
    /// not. Found by this witness, 2026-07-25; the proof document's census ran on Heisenberg
    /// alone.</summary>
    [Theory]
    [InlineData(4)]
    [InlineData(5)]
    public void TheCornerConfinement_IsHeisenbergOnly(int n)
    {
        var heis = new FrozenDivisorWitness(n, chain: DivisorChain.Heisenberg);
        var xy = new FrozenDivisorWitness(n, chain: DivisorChain.Xy);

        Assert.True(heis.CensusRan);
        Assert.True(xy.CensusRan);

        // Both chains freeze the same number of modes in the corner, at the same root.
        Assert.Equal(n / 2, heis.KernelDimension);
        Assert.Equal(n / 2, xy.KernelDimension);
        Assert.Equal(heis.Root, xy.Root, 12);

        // Heisenberg: exactly the four corners. XY: strictly more.
        Assert.True(heis.CornersAreTheOnlyCarriers, $"Heisenberg N={n}: {heis.Summary}");
        Assert.Equal(4, heis.Census.Count);
        Assert.False(xy.CornersAreTheOnlyCarriers, $"XY N={n}: {xy.Summary}");
        Assert.True(xy.Census.Count > 4, $"XY N={n} census has only {xy.Census.Count} carriers");

        // Every XY carrier carries the full ⌊N/2⌋, never a partial count: what spreads is the
        // whole divisor, not a remnant of it. And every block carrying the UNFOLDED root −4γ̄ has
        // p + q even (measured over the full census at these N, not derived; at odd N a one-sided
        // fold flips that parity, which is why the fold root's carriers sit on the odd blocks).
        foreach (var c in xy.Census)
        {
            Assert.Equal(n / 2, Math.Max(c.AtRoot, c.AtFoldRoot));
            if (c.AtRoot > 0)
                Assert.True((c.P + c.Q) % 2 == 0,
                    $"XY block ({c.P},{c.Q}) carries −4γ̄ but has odd p+q");
        }

        // Both chains still obey the theorem as this witness scopes it per chain.
        Assert.True(heis.Holds, heis.Summary);
        Assert.True(xy.Holds, xy.Summary);
    }

    [Theory]
    [InlineData(3)]
    [InlineData(4)]
    [InlineData(5)]
    [InlineData(6)]
    [InlineData(7)]
    public void ZeroCoupling_DoublesTheMultiplicity(int n)
    {
        var w = new FrozenDivisorWitness(n);
        Assert.Equal(2 * w.PredictedFrozen, w.ZeroCouplingKernelDimension);
    }

    /// <summary>The four corners p, q ∈ {1, N−1}, each carrying ⌊N/2⌋ at the root its fold parity
    /// selects. At N = 4 the fold sends the root to itself ((4 − 2N)γ̄ = −4γ̄), so the parity is
    /// invisible and all four corners carry the single root.</summary>
    [Theory]
    [InlineData(3)]
    [InlineData(4)]
    [InlineData(5)]
    [InlineData(6)]
    public void FourCorners_SplitByFoldParity(int n)
    {
        var w = new FrozenDivisorWitness(n);
        Assert.Equal(4, w.Corners.Count);

        foreach (var c in w.Corners)
        {
            if (w.RootsCoincide)
            {
                Assert.Equal(w.PredictedFrozen, c.AtRoot);
                Assert.Equal(w.PredictedFrozen, c.AtFoldRoot);
            }
            else if (c.Folds % 2 == 0)
            {
                Assert.Equal(w.PredictedFrozen, c.AtRoot);
                Assert.Equal(0, c.AtFoldRoot);
            }
            else
            {
                Assert.Equal(0, c.AtRoot);
                Assert.Equal(w.PredictedFrozen, c.AtFoldRoot);
            }
        }

        Assert.Equal(n == 4, w.RootsCoincide);
    }

    /// <summary>The root value itself: −4γ̄ = −4σ/N, the mean dephasing rate and nothing else.</summary>
    [Fact]
    public void Root_IsMinusFourTimesTheMeanWatching()
    {
        var w = new FrozenDivisorWitness(6);
        Assert.Equal(-4.0 * w.Sigma / w.N, w.Root, 12);
        Assert.Equal(4.0 * w.GammaBar - 2.0 * w.Sigma, w.FoldRoot, 12);

        // The locus, cell by cell: every reflection pair carries the same total.
        for (int l = 0; l < w.N; l++)
            Assert.Equal(2.0 * w.GammaBar, w.Gamma[l] + w.Gamma[w.N - 1 - l], 12);

        // ... and the profile is genuinely non-uniform, or the locus would be trivial.
        Assert.True(w.Gamma.Distinct().Count() > 1);
    }

    /// <summary>The ladder, closed form: distances N + 1 − 2c, departure orders twice that, total
    /// 2⌊N²/4⌋.</summary>
    [Theory]
    [InlineData(5, 12)]
    [InlineData(6, 18)]
    [InlineData(7, 24)]
    public void Ladder_MatchesTheClosedForm(int n, int totalValuation)
    {
        var w = new FrozenDivisorWitness(n);
        Assert.Equal(n / 2, w.PairDistances.Count);
        for (int c = 1; c <= n / 2; c++)
        {
            Assert.Equal(n + 1 - 2 * c, w.PairDistances[c - 1]);
            Assert.Equal(2 * (n + 1 - 2 * c), w.DepartureOrders[c - 1]);
        }
        Assert.Equal(totalValuation, w.TotalValuation);
        Assert.Equal(2 * w.PairDistances.Sum(), w.TotalValuation);
    }

    [Theory]
    [InlineData(FrozenDivisorWitness.MinN - 1)]
    [InlineData(FrozenDivisorWitness.MaxN + 1)]
    public void OutOfRangeN_Throws(int n) =>
        Assert.Throws<ArgumentOutOfRangeException>(() => new FrozenDivisorWitness(n));

    [Fact]
    public void GammaBarBelowTheHalfProfileSteps_Throws() =>
        Assert.Throws<ArgumentOutOfRangeException>(() => new FrozenDivisorWitness(6, gammaBar: 0.001));
}
