using System.Numerics;
using RCPsiSquared.Core.BlockSpectrum;
using RCPsiSquared.Core.F89PathK;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Diagnostics.Foundation;
using RCPsiSquared.Diagnostics.Knowledge;
using Xunit;

namespace RCPsiSquared.Diagnostics.Tests.Foundation;

/// <summary>F153, gated from below. The entry-wise route has no eigensolver in it, so every residual
/// here is compared to 0.0 EXACTLY; the one place a tolerance appears is the ZZ frequency, where the
/// assertion is that a spread is strictly positive rather than that it equals anything.</summary>
public class PinnedBlockFloorClaimTests
{
    [Theory]
    [InlineData(2)]
    [InlineData(3)]
    [InlineData(4)]
    [InlineData(5)]
    [InlineData(6)]
    public void Criterion_Counts_Exactly_Four_N_Blocks(int n)
    {
        int counted = 0;
        for (int p = 0; p <= n; p++)
            for (int q = 0; q <= n; q++)
                if (PinnedBlockFloorClaim.IsPinned(n, p, q)) counted++;

        Assert.Equal(4 * n, counted);
        Assert.Equal(4 * n, PinnedBlockFloorClaim.PinnedBlockCount(n));
        Assert.Equal((n + 1) * (n + 1), JointPopcountSectors.SectorCount(n));
    }

    /// <summary>The criterion IS the window's zero-width case, which is the half PROOF_CODIM1 §6 owns.
    /// Both sides are computed here, and the window itself is checked against the enumerated cells so
    /// the closed form is not merely asserted.</summary>
    [Theory]
    [InlineData(3)]
    [InlineData(4)]
    [InlineData(5)]
    public void Window_Is_The_Enumerated_Range_And_Zero_Width_Is_The_Criterion(int n)
    {
        for (int p = 0; p <= n; p++)
            for (int q = 0; q <= n; q++)
            {
                var kets = WeightCoherenceBlock.Configs(n, p);
                var bras = WeightCoherenceBlock.Configs(n, q);
                int lo = int.MaxValue, hi = int.MinValue;
                foreach (var a in kets)
                    foreach (var b in bras)
                    {
                        int nd = System.Numerics.BitOperations.PopCount((uint)(a ^ b));
                        lo = Math.Min(lo, nd); hi = Math.Max(hi, nd);
                    }

                var (predLo, predHi) = PinnedBlockFloorClaim.DisagreementWindow(n, p, q);
                Assert.Equal(lo, predLo);
                Assert.Equal(hi, predHi);
                Assert.Equal(lo == hi, PinnedBlockFloorClaim.IsPinned(n, p, q));
            }
    }

    /// <summary>The entry-wise certificate, both directions, bit-exact. On a pinned block: constant REAL
    /// diagonal, no real part off the diagonal, off-diagonal symmetric, i.e. L = c·Id + i·S with S real
    /// symmetric. On the other (N−1)² blocks the real diagonal spread is STRICTLY positive, which is the
    /// half that earns the criterion: a one-sided gate would pass on "every block is pinned".</summary>
    [Theory]
    [InlineData(4, 1.0, 1.0)]
    [InlineData(4, 0.25, 1.0)]
    [InlineData(4, 3.0, 1.0)]
    [InlineData(5, 1.0, 1.0)]
    [InlineData(5, 1.0, 0.5)]
    [InlineData(5, 1.0, 2.0)]
    [InlineData(5, 1.0, -1.0)]
    [InlineData(5, 1.0, 0.0)]
    [InlineData(6, 1.0, 1.0)]
    public void Entry_Wise_Certificate_Holds_In_Both_Directions(int n, double j, double delta)
    {
        int whole = 0;
        foreach (var r in PinnedBlockFloorWitness.Sweep(n, j, delta))
        {
            Assert.True(r.OffDiagonalRealMass == 0.0,
                $"({r.P},{r.Q}) has a real part off the diagonal: {r.OffDiagonalRealMass}");
            Assert.True(r.OffDiagonalAsymmetry == 0.0,
                $"({r.P},{r.Q}) off-diagonal is not symmetric: {r.OffDiagonalAsymmetry}");

            if (r.Pinned)
            {
                Assert.True(r.RealDiagonalSpread == 0.0,
                    $"pinned block ({r.P},{r.Q}) has a real-diagonal spread: {r.RealDiagonalSpread}");
                Assert.True(PinnedBlockFloorWitness.ObservedFloor(n, r.P, r.Q, j, delta) == r.Floor,
                    $"pinned block ({r.P},{r.Q}) does not sit on −2|p−q| = {r.Floor}");
                whole++;
            }
            else
            {
                Assert.True(r.RealDiagonalSpread > 0.0,
                    $"free block ({r.P},{r.Q}) has a CONSTANT real diagonal, which would break the criterion");
            }
        }

        Assert.Equal(4 * n, whole);
        Assert.Equal(4 * n, PinnedBlockFloorWitness.WholeCount(n, j, delta));
    }

    /// <summary>The correction that must travel with a Heisenberg carrier. MirrorWorld's gate certifies
    /// the XY-only form c·Id + i·J·A and passes bit-exactly because its rate atom has no ZZ. With Δ ≠ 0
    /// the pinned diagonal is NOT c·Id: it is constant on the real axis and configuration-dependent on the
    /// imaginary one. At Δ = 0 the same spread collapses to exactly 0.0, which is the XY statement, so
    /// this test also shows the two forms are distinguishable rather than asserting it.</summary>
    [Fact]
    public void The_ZZ_Diagonal_Is_Imaginary_And_Not_Constant()
    {
        var heis = PinnedBlockFloorWitness.Read(4, 0, 1, j: 1.0, delta: 1.0);
        Assert.True(heis.RealDiagonalSpread == 0.0, "the ZZ term must not move the real part");
        Assert.True(heis.ImaginaryDiagonalSpread > 0.0,
            "at Δ=1 the pinned diagonal must NOT be c·Id: the ZZ term is a configuration-dependent frequency");

        var xy = PinnedBlockFloorWitness.Read(4, 0, 1, j: 1.0, delta: 0.0);
        Assert.True(xy.ImaginaryDiagonalSpread == 0.0,
            "at Δ=0 the pinned diagonal IS c·Id, the form the MirrorWorld gate certifies");
    }

    /// <summary>The two frequencies the arc recorded for the N=4 (0,1) block, reproduced from the builder
    /// rather than restated: the cells read −2−2i and −2−4i at γ = 1 (the arc's −1−2i and −1−4i are the
    /// same frequencies at γ = ½). Every value exact, so compared exactly.</summary>
    [Fact]
    public void The_N4_Vacuum_Block_Reads_Two_Frequencies_On_One_Floor()
    {
        var l = WeightCoherenceBlock.Build(4, 0, 1, new Complex(1.0, 0.0), 1.0);
        var seen = new SortedSet<double>();
        for (int r = 0; r < l.GetLength(0); r++)
        {
            Assert.True(l[r, r].Real == -2.0, $"cell {r} left the floor: {l[r, r].Real}");
            seen.Add(l[r, r].Imaginary);
        }

        Assert.Equal(new[] { -4.0, -2.0 }, seen);
    }

    /// <summary>The converse's premise, checked where it is cheap: the real trace equals −2·Σ_cells n_diff,
    /// which is the identity the ⟸ argument reads the block through. If this ever stopped holding the
    /// trace proof would be about a different matrix than the one the witness builds.</summary>
    [Theory]
    [InlineData(4, 1.0)]
    [InlineData(4, 0.0)]
    [InlineData(4, 2.0)]
    [InlineData(5, 1.0)]
    [InlineData(5, -1.0)]
    public void Real_Trace_Counts_The_Disagreements(int n, double delta)
    {
        for (int p = 0; p <= n; p++)
            for (int q = 0; q <= n; q++)
            {
                var l = WeightCoherenceBlock.Build(n, p, q, new Complex(1.0, 0.0), delta);
                double reTrace = 0.0;
                for (int r = 0; r < l.GetLength(0); r++) reTrace += l[r, r].Real;

                double counted = 0.0;
                foreach (var a in WeightCoherenceBlock.Configs(n, p))
                    foreach (var b in WeightCoherenceBlock.Configs(n, q))
                        counted += -2.0 * System.Numerics.BitOperations.PopCount((uint)(a ^ b));

                Assert.True(reTrace == counted,
                    $"({p},{q}): real trace {reTrace} is not the disagreement count {counted}");
            }
    }

    /// <summary>The REAL-q fence, as a negative control. §6 states its Edge lemma "at real q" and the
    /// criterion inherits that: off the real axis the hops put a real part off the diagonal and the ZZ term
    /// puts one on it, so a pinned block leaves its floor. It can even gain: the (0,1) block at N=4, Δ=1,
    /// q = 1+0.3i reaches Re λ ≈ +0.05 against a floor of −2. Read entry-wise, no eigensolver: the two
    /// residuals that are exactly 0.0 at real q are both nonzero here.</summary>
    [Fact]
    public void Complex_Coupling_Breaks_The_Certificate()
    {
        var l = WeightCoherenceBlock.Build(4, 0, 1, new Complex(1.0, 0.3), 1.0);
        int d = l.GetLength(0);

        double reLo = double.PositiveInfinity, reHi = double.NegativeInfinity, realMass = 0.0;
        for (int r = 0; r < d; r++)
        {
            reLo = Math.Min(reLo, l[r, r].Real); reHi = Math.Max(reHi, l[r, r].Real);
            for (int c = 0; c < d; c++)
                if (r != c) realMass = Math.Max(realMass, Math.Abs(l[r, c].Real));
        }

        Assert.True(reHi - reLo > 0.0, "at complex q the ZZ term must put a spread on the REAL diagonal");
        Assert.True(realMass > 0.0, "at complex q the hops must put a real part OFF the diagonal");

        // the same block at real q, for contrast: both exactly 0.0 (asserted elsewhere, restated here so
        // the control is a comparison and not an isolated inequality)
        var real = PinnedBlockFloorWitness.Read(4, 0, 1, j: 1.0, delta: 1.0);
        Assert.True(real.RealDiagonalSpread == 0.0 && real.OffDiagonalRealMass == 0.0);
    }

    /// <summary>The UNIFORM-γ fence, gated entry-wise. Pinning fixes the disagreement COUNT, not the
    /// disagreeing SITES, so under a profile a pinned block's cell rates differ and the mechanism fails.
    /// This needs no Hamiltonian at all (at real q neither the hops nor the ZZ term touch the real
    /// diagonal), which is why it can be gated exactly while the SPECTRAL consequence cannot: the
    /// committed Python verifier takes a scalar γ, so the spectral numbers in the claim are measured and
    /// not gated. The blocks chosen are pinned ones with more than one cell; a corner block (dim 1, or one
    /// whose cells all disagree everywhere) has nothing to spread.</summary>
    [Theory]
    [InlineData(4, 0, 1)]
    [InlineData(4, 0, 2)]
    [InlineData(4, 1, 4)]
    [InlineData(5, 0, 1)]
    [InlineData(5, 5, 2)]
    public void A_Gamma_Profile_Breaks_The_Pinning_Mechanism(int n, int p, int q)
    {
        var uniform = Enumerable.Repeat(1.0, n).ToArray();
        Assert.True(PinnedBlockFloorWitness.ProfileRealDiagonalSpread(n, p, q, uniform) == 0.0,
            "at uniform γ a pinned block's cell rates are all equal, exactly");

        var profile = Enumerable.Range(0, n).Select(k => 0.1 + 0.9 * k / (n - 1.0)).ToArray();
        Assert.True(PinnedBlockFloorWitness.ProfileRealDiagonalSpread(n, p, q, profile) > 0.0,
            "under a profile the cell rates must differ: pinning fixes the COUNT, not the SITES");

        Assert.Throws<ArgumentException>(
            () => PinnedBlockFloorWitness.ProfileRealDiagonalSpread(n, p, q, new double[n - 1]));
    }

    [Fact]
    public void Claim_Is_Registered_With_Both_Typed_Parents()
    {
        var registry = KnowledgeRegistryFactory.BuildDefault();
        Assert.True(registry.Contains<PinnedBlockFloorClaim>());
        Assert.Equal(Tier.Tier1Derived, registry.Get<PinnedBlockFloorClaim>().Tier);

        var ancestors = registry.AncestorsOf<PinnedBlockFloorClaim>()
            .Select(c => c.GetType()).ToHashSet();
        Assert.Contains(typeof(AbsorptionTheoremClaim), ancestors);
        Assert.Contains(typeof(JointPopcountSectors), ancestors);
    }

    /// <summary>The claim's prose must carry the two fences a reader could otherwise drop: uniform γ, and
    /// that the S includes the ZZ diagonal. Both were review findings before the carrier existed.</summary>
    [Fact]
    public void Statement_Carries_The_Uniform_Gamma_And_ZZ_Fences()
    {
        var claim = PinnedBlockFloorClaim.Shared;
        Assert.Contains("UNIFORM γ", claim.Name);
        Assert.Contains("ZZ diagonal", claim.Name);
        Assert.Contains("PROOF_CODIM1_BY_ADDITIVITY", claim.Name);
        Assert.Contains("inspect --root pinned", claim.Anchor);
        Assert.Contains("F153", claim.Anchor);
    }

    [Fact]
    public void Guards_Reject_Blocks_Outside_The_Grading()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => PinnedBlockFloorClaim.IsPinned(4, -1, 0));
        Assert.Throws<ArgumentOutOfRangeException>(() => PinnedBlockFloorClaim.IsPinned(4, 0, 5));
        Assert.Throws<ArgumentOutOfRangeException>(() => new PinnedBlockFloorWitness(n: 1));
        Assert.Throws<ArgumentOutOfRangeException>(() => new PinnedBlockFloorWitness(n: 8));
    }
}
