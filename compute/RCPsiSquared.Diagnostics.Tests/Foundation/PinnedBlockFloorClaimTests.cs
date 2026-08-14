using System.Numerics;
using RCPsiSquared.Core.BlockSpectrum;
using RCPsiSquared.Core.F89PathK;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Pauli;
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
    /// diagonal), which is why it can be gated exactly while the SPECTRAL consequence was long thought not to be (it is, since 2026-08-14, further down this file): the
    /// committed Python verifier takes a scalar γ, so the spectral numbers in the claim are measured and
    /// not gated. The blocks chosen are pinned ones with more than one cell; a one-cell block (dim 1, or one
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

    // ---- the SPECTRAL half of the γ fence ----
    //
    // Every profile used by an EXACT gate below is DYADIC on purpose (the off-locus and linspace profiles
    // further down are not; they feed only eigensolver readings, where the question does not arise). The fence turns on the cell rates being CONSTANT, and
    // subset sums that are equal in exact arithmetic need not be bit-equal in doubles, so a comparison
    // to 0.0 on a non-dyadic profile would be reading the summation order. Dyadic values make every
    // subset sum exact, which is the only honest way to keep the comparison at 0.0.

    private static double[] DyadicRamp(int n) =>
        Enumerable.Range(0, n).Select(k => 0.25 + 0.25 * k).ToArray();

    private static double[] Uniform(int n, double g) => Enumerable.Repeat(g, n).ToArray();

    /// <summary>THE MASTER INVARIANT, from below and by TWO ROUTES: the block comes from the general
    /// Liouvillian builder acting on a Hilbert-space H, and the prediction −2·diag(rate) is assembled
    /// from the cells alone, so this is not one code path agreeing with itself. Herm(L) = −2·diag(rate)
    /// exactly, for every one of the (N+1)² blocks, at every Δ including the anisotropic and negative
    /// ones, and under a profile. It is what makes the whole fence exact: Herm is DIAGONAL, so its
    /// spectrum is read off and no eigensolver is needed to know the real parts.</summary>
    [Theory]
    [InlineData(3, 1.0, 1.0)]
    [InlineData(3, 0.25, 0.0)]
    [InlineData(4, 1.0, 1.0)]
    [InlineData(4, 3.0, -0.7)]
    [InlineData(4, 0.25, 2.0)]
    [InlineData(5, 1.0, 1.0)]
    public void ProfileHermitianPart_IsExactlyMinusTwoDiagRate_OnEveryBlock(int n, double j, double delta)
    {
        var gamma = DyadicRamp(n);
        for (int p = 0; p <= n; p++)
            for (int q = 0; q <= n; q++)
            {
                double residual = PinnedBlockFloorWitness.ProfileHermitianPartResidual(n, p, q, j, delta, gamma);
                Assert.True(residual == 0.0,
                    $"[N={n} J={j} Δ={delta} block ({p},{q})] Herm(L) + 2·diag(rate) is not bit-exact 0: {residual:E3}");
            }
    }

    /// <summary>The control the invariant above needs, because a residual that cannot move is not a
    /// gate. Perturbing ONE site's rate in the PREDICTION by one dyadic step must make it fire, on a
    /// block whose cells actually see that site. It fires by exactly twice the perturbation, since the
    /// prediction enters as −2·rate.</summary>
    [Fact]
    public void ProfileHermitianPart_Fires_WhenThePredictionIsPerturbed()
    {
        const int n = 4;
        var gamma = DyadicRamp(n);
        var wrong = (double[])gamma.Clone();
        wrong[2] += 0.5;

        var l = PinnedBlockFloorWitness.ProfileBlock(n, 0, 2, 1.0, 1.0, gamma);
        var cells = PinnedBlockFloorWitness.Cells(n, 0, 2);
        double worst = 0.0;
        for (int r = 0; r < cells.Count; r++)
            worst = Math.Max(worst, Math.Abs(l[r, r].Real + 2.0 * PinnedBlockFloorWitness.CellRate(n, cells[r].Ket, cells[r].Bra, wrong)));

        // Exactly 1.0, not to twelve places: every input is dyadic, so the perturbation enters the
        // prediction as −2·0.5 with no rounding anywhere. A tolerance here would be an error code.
        Assert.True(worst == 1.0, $"the perturbed prediction must fire by exactly 2·0.5: {worst:R}");
    }

    /// <summary>The SECOND control the invariant needs, and it is the one that keeps the class's
    /// convention paragraph honest. <see cref="PinnedBlockFloorWitness.CellRate"/> claims site l is bit
    /// n−1−l; if that were merely decorative, feeding the REVERSED profile would pass just as well. It
    /// does not: on a non-palindromic profile the reversed reading fires. Contrast
    /// <see cref="PinnedBlockFloorWitness.ProfileRealDiagonalSpread"/>, which is a max−min over a cell
    /// set closed under bit permutations and therefore returns the SAME number either way, asserted here
    /// so the difference between the two members is a gated fact and not a claim in a doc-comment.</summary>
    [Fact]
    public void TheSiteConvention_IsObservable_EntryWise_AndInvisible_InTheSpread()
    {
        const int n = 4;
        var gamma = DyadicRamp(n);
        var reversed = gamma.Reverse().ToArray();

        Assert.True(PinnedBlockFloorWitness.ProfileHermitianPartResidual(n, 0, 1, 1.0, 1.0, gamma) == 0.0);

        var l = PinnedBlockFloorWitness.ProfileBlock(n, 0, 1, 1.0, 1.0, gamma);
        var cells = PinnedBlockFloorWitness.Cells(n, 0, 1);
        double worstReversed = 0.0;
        for (int r = 0; r < cells.Count; r++)
            worstReversed = Math.Max(worstReversed,
                Math.Abs(l[r, r].Real + 2.0 * PinnedBlockFloorWitness.CellRate(n, cells[r].Ket, cells[r].Bra, reversed)));
        Assert.True(worstReversed > 0.0,
            "the reversed profile must break the entry-wise reading: the convention is load-bearing here");

        Assert.Equal(PinnedBlockFloorWitness.ProfileRealDiagonalSpread(n, 0, 1, gamma),
                     PinnedBlockFloorWitness.ProfileRealDiagonalSpread(n, 0, 1, reversed));
        Assert.Equal(PinnedBlockFloorWitness.CellRateSpread(n, 0, 1, gamma),
                     PinnedBlockFloorWitness.CellRateSpread(n, 0, 1, reversed));
    }

    /// <summary>The trace, the second exact route and the one that carries the CONVERSE. Σ Re λ = Re tr L
    /// holds for defective blocks too, the trace being the characteristic polynomial's, so this is the
    /// direction "entirely on the floor ⇒ the rates are constant" with no eigensolver in it.</summary>
    [Theory]
    [InlineData(3, 1.0, 1.0)]
    [InlineData(4, 1.0, 1.0)]
    [InlineData(4, 0.25, -0.7)]
    [InlineData(5, 3.0, 0.0)]
    public void ProfileTrace_IsTheExactRateSum_OnEveryBlock(int n, double j, double delta)
    {
        var gamma = DyadicRamp(n);
        for (int p = 0; p <= n; p++)
            for (int q = 0; q <= n; q++)
            {
                double residual = PinnedBlockFloorWitness.ProfileTraceResidual(n, p, q, j, delta, gamma);
                Assert.True(residual == 0.0,
                    $"[N={n} J={j} Δ={delta} block ({p},{q})] Re tr L + 2·Σ rate is not bit-exact 0: {residual:E3}");
            }
    }

    /// <summary>THE FENCE ITSELF, and the sharp form is not the one F153's scope paragraph states. The
    /// criterion turns on the cell rates being CONSTANT, which on a pinned block of dimension greater
    /// than one happens exactly when γ is uniform; the four ONE-CELL blocks (0,0), (0,N), (N,0), (N,N) are one-dimensional and
    /// therefore stay on their floor under ANY profile, which is the fence's only exception and it is a
    /// dimensional one. Both halves are exact comparisons to 0.0.</summary>
    [Theory]
    [InlineData(3)]
    [InlineData(4)]
    [InlineData(5)]
    public void TheFence_IsRateConstancy_AndTheFourOneCellBlocksAreItsOnlyException(int n)
    {
        var profile = DyadicRamp(n);
        var uniform = Uniform(n, 0.5);
        int oneCell = 0, spread = 0;

        for (int p = 0; p <= n; p++)
            for (int q = 0; q <= n; q++)
            {
                if (!PinnedBlockFloorClaim.IsPinned(n, p, q)) continue;

                Assert.True(PinnedBlockFloorWitness.CellRateSpread(n, p, q, uniform) == 0.0,
                    $"block ({p},{q}) at uniform γ must have constant cell rates, exactly");

                double underProfile = PinnedBlockFloorWitness.CellRateSpread(n, p, q, profile);
                if (PinnedBlockFloorWitness.Cells(n, p, q).Count == 1)
                {
                    oneCell++;
                    // NOT "the spread of a one-element set is zero", which no input could break. The
                    // content is that the BUILT block, Hamiltonian and all, still sits on its floor under
                    // a profile: its one diagonal entry is the floor exactly, and its imaginary part is
                    // whatever the ZZ term makes it.
                    var oneCellBlock = PinnedBlockFloorWitness.ProfileBlock(n, p, q, 1.0, 1.0, profile);
                    var cell = PinnedBlockFloorWitness.Cells(n, p, q)[0];
                    double floor = -2.0 * PinnedBlockFloorWitness.CellRate(n, cell.Ket, cell.Bra, profile);
                    Assert.True(oneCellBlock[0, 0].Real == floor,
                        $"one-cell block ({p},{q}) must sit on its floor under ANY profile: {oneCellBlock[0, 0].Real:R} vs {floor:R}");
                }
                else
                {
                    spread++;
                    Assert.True(underProfile > 0.0,
                        $"pinned block ({p},{q}) with more than one cell must leave its floor under a profile");
                }
            }

        Assert.Equal(4, oneCell);
        Assert.Equal(4 * n - 4, spread);
    }

    /// <summary>The remaining direction, and the ONE place an eigensolver appears: at uniform γ every
    /// real part of a pinned block IS the floor. There is no exact route to an eigenvalue, so this is
    /// gated as a LAW rather than a number. The model is the eigensolver's backward error, ε·‖L‖_F, and
    /// the gate is that the RATIO to it stays bounded as J crosses four decades; a threshold that merely
    /// passed at one coupling would say nothing. The ratios are printed on failure.</summary>
    [Theory]
    [InlineData(4, 0, 1)]
    [InlineData(4, 0, 2)]
    [InlineData(4, 4, 1)]
    [InlineData(5, 0, 2)]
    public void AtUniformGamma_EveryRealPartIsTheFloor_AtTheBackwardErrorLaw(int n, int p, int q)
    {
        const double gammaValue = 0.5;
        double floor = -2.0 * gammaValue * Math.Abs(p - q);
        var gamma = Uniform(n, gammaValue);
        var ratios = new List<double>();

        foreach (double j in new[] { 0.125, 1.0, 8.0, 64.0, 512.0 })
        {
            var l = PinnedBlockFloorWitness.ProfileBlock(n, p, q, j, 1.0, gamma);
            double frob = l.FrobeniusNorm();
            double worst = l.Evd().EigenValues.Max(z => Math.Abs(z.Real - floor));
            ratios.Add(worst / (Eps * frob));
        }

        Assert.True(ratios.Max() < 100.0,
            $"block ({p},{q}) at N={n}: max|Re λ − floor| must track ε·‖L‖_F across decades; " +
            $"ratios {string.Join(", ", ratios.Select(r => r.ToString("0.0")))}");
    }

    /// <summary>WHAT THE GENERIC FLAT SET ACTUALLY IS, and TWO earlier readings of it were wrong: first
    /// that it tracks the coupling, then that it is the EVEN |p−q| blocks. The second reading is an N=4
    /// ACCIDENT and this gate is built to kill it, which is why N=6 is a row: there the flat blocks carry
    /// |p−q| = 3, ODD, and every even one is spread.
    ///
    /// <para>The invariant is the block's FREE INDEX, the one of p, q that is not 0 or N. Flattening needs
    /// an antilinear involution of the block about the trace constant, i.e. a pairing of cells sending
    /// rate ↦ 2·mean(rate) − rate. When the free index is N/2 the bitwise COMPLEMENT of the free config
    /// stays inside the block and does exactly that, since rate + rate(complement) = Σγ = 2·mean(rate),
    /// for ANY profile. So: four such blocks at even N, none at odd N, and no parity anywhere in it.</para>
    ///
    /// <para>Flatness is an eigensolver reading with no exact route, but it needs no chosen threshold:
    /// "flat" means AT the eigensolver's backward error ε·‖L‖_F, never below a fixed number, while every
    /// spread one is above 1e-2, and the gate asserts BOTH sides of that gap. An absolute cut would not do,
    /// since ‖L‖_F grows with J and flat spans really do reach 1.7e-11 by J = 1000. The rows where the flat
    /// set is EMPTY (odd N) are the discriminating ones.</para></summary>
    [Theory]
    [InlineData(3)]
    [InlineData(4)]
    [InlineData(5)]
    [InlineData(6)]
    public void TheGenericFlatSet_IsTheFreeIndexAtHalfFilling_AndNotAParity(int n)
    {
        var gamma = OffLocusProfile(n);
        AssertOffTheR90Locus(n, gamma);

        var (flat, spread) = FlatAndSpreadPinnedBlocks(n, PauliJ, 1.0, gamma);

        var expected = AllPinnedAboveDimensionOne(n)
            .Where(b => FreeIndex(n, b.P, b.Q) * 2 == n)
            .ToHashSet();
        Assert.Equal(expected, flat.Select(f => f.Block).ToHashSet());
        Assert.Equal(n % 2 == 0 ? 4 : 0, expected.Count);

        Assert.All(flat, f => Assert.True(f.AtBackwardError,
            $"block {f.Block} is called flat at {f.Ratio:F1}× the eigensolver's error scale"));
        Assert.All(spread, s => Assert.True(s.Span > 1e-2,
            $"block {s.Block} is called spread at span {s.Span:E3}; the verdict must not rest on the cut"));
    }

    /// <summary>The N=4 accident, named as such. At N=4 the free-index-N/2 blocks happen to carry
    /// |p−q| = 2 and at N=6 they carry 3, so a parity read off N=4 alone inverts one step later. This
    /// gate holds the coincidence AND its failure, so neither can be mistaken for the law again.</summary>
    [Fact]
    public void TheEvenParityReading_HoldsAtN4_AndIsInvertedAtN6()
    {
        Assert.All(AllPinnedAboveDimensionOne(4).Where(b => FreeIndex(4, b.P, b.Q) == 2),
            b => Assert.Equal(0, Math.Abs(b.P - b.Q) % 2));
        Assert.All(AllPinnedAboveDimensionOne(6).Where(b => FreeIndex(6, b.P, b.Q) == 3),
            b => Assert.Equal(1, Math.Abs(b.P - b.Q) % 2));
    }

    /// <summary>The SECOND source, and the only one that reaches all 4N: γ ANTI-PALINDROMIC,
    /// γ_l + γ_{N−1−l} constant, the R₉₀ locus of the parameter Klein group (F91). There the chain
    /// reflection supplies the involution for EVERY pinned block, so the whole 4N coalesce, at odd N as
    /// well, where the free-index source gives nothing at all. Both linspace(0.1, 1, N) and the gates'
    /// dyadic ramp lie on this locus, which is the source that CAN reach all 4N; whether it does at a given
    /// coupling is the onset's business, and F153's quoted numbers are taken below it.
    ///
    /// <para>NOT F140, which lives on the same locus and is a different object: F140 pins −4γ̄ on the
    /// (1,1)-type corner blocks, which are not pinned in F153's sense, and its proof says the spectrum
    /// there is NOT palindromic about the root. What the two share is the locus, and that is F91's.</para></summary>
    [Theory]
    [InlineData(4)]
    [InlineData(5)]
    [InlineData(6)]
    public void OnTheR90Locus_EveryPinnedBlockFlattens_IncludingAtOddN(int n)
    {
        var gamma = R90Profile(n);
        var pairSums = Enumerable.Range(0, n).Select(l => gamma[l] + gamma[n - 1 - l]).ToArray();
        Assert.True(pairSums.Max() == pairSums.Min(), "the profile is dyadic, so locus membership is EXACT");

        var (flat, spread) = FlatAndSpreadPinnedBlocks(n, PauliJ, 1.0, gamma);
        Assert.Empty(spread);
        Assert.Equal(4 * n - 4, flat.Count);            // the four one-cell blocks are excluded by dimension
        Assert.All(flat, f => Assert.True(f.AtBackwardError,
            $"block {f.Block} span {f.Span:E3} is {f.Ratio:F1}× the eigensolver's error scale"));
    }

    /// <summary>WHAT IS Δ-DEPENDENT IS THE COLLAPSE, NOT THE INVOLUTION, and the difference matters
    /// because an earlier wording said the source "fails at Δ = 0", which reads as the symmetry failing.
    /// It does not: the palindrome is machine-exact at every Δ here (gated next door). What Δ decides is
    /// whether the collapse has an ONSET at all. At Δ = 1 and Δ = 2 it does. At Δ = 0 and Δ = 0.5 the span
    /// SATURATES instead: it stops moving with J entirely, so there is no threshold anywhere, and calling
    /// it "below threshold" would be wrong.</summary>
    [Theory]
    [InlineData(4, 0, 2, 1.0, true)]
    [InlineData(4, 0, 2, 2.0, true)]
    [InlineData(4, 0, 2, 0.0, false)]
    [InlineData(6, 0, 3, 0.5, false)]
    public void TheCollapse_HasAnOnset_ForSomeDeltaOnly(int n, int p, int q, double delta, bool hasOnset)
    {
        var gamma = OffLocusProfile(n);
        double Span(double j)
        {
            var (lo, hi) = PinnedBlockFloorWitness.ProfileReSpan(n, p, q, j, delta, gamma);
            return hi - lo;
        }

        double atThousand = Span(1000.0);
        if (hasOnset)
        {
            // J = 0.2, not 1: at Δ = 1 this profile's onset is already below J = 1, which is exactly the
            // point that the onset moves with the profile and with Δ and is therefore not quoted anywhere.
            Assert.True(Span(0.2) > 1e-2, $"below the onset the block must be wide: {Span(0.2):E3}");
            Assert.True(atThousand < 1e-9, $"above it, collapsed: {atThousand:E3}");
        }
        else
        {
            // Saturation, not a distant threshold: the span at J = 1000 and at J = 100000 agree to four
            // figures, so raising the coupling further cannot collapse it.
            Assert.True(atThousand > 1e-2, $"at Δ={delta} the block must NOT collapse: {atThousand:E3}");
            Assert.Equal(atThousand, Span(100_000.0), 4);
        }
    }

    /// <summary>FLAT AND ON THE FLOOR ARE TWO PREDICATES, and this keeps them apart. A flattened pinned
    /// block sits at the TRACE's constant −2·|p−q|·γ̄, which is strictly BELOW its floor −2·min(rate)
    /// whenever the profile is non-uniform, so flattening is never the criterion being met. The gap is
    /// asserted against its own closed form rather than a chosen slack.</summary>
    [Fact]
    public void AFlattenedBlock_SitsAtTheTraceConstant_StrictlyBelowItsFloor()
    {
        // |p−q| = 2, NOT 1, and that is the whole point: at |p−q| = 1 the trace constant −2·|p−q|·γ̄
        // coincides with the withdrawn −2·mean(γ), so a row there would pass under the mislabel this gate
        // exists to catch. Here they are −2.5 and −1.25.
        const int n = 4, p = 0, q = 2;
        var gamma = R90Profile(n);
        var cells = PinnedBlockFloorWitness.Cells(n, p, q);
        var l = PinnedBlockFloorWitness.ProfileBlock(n, p, q, PauliJ, 1.0, gamma);
        var re = l.Evd().EigenValues.Select(z => z.Real).ToArray();
        double floor = -2.0 * cells.Min(c => PinnedBlockFloorWitness.CellRate(n, c.Ket, c.Bra, gamma));
        double traceConstant = -2.0 * Math.Abs(p - q) * gamma.Average();

        // FIRST that the block has actually FLATTENED. Without this the test would assert only the trace,
        // which equals the trace constant on every pinned block, flat or spread, and is gated next door;
        // the whole point here is a block that IS flat and still not on its floor.
        Assert.True(re.Max() - re.Min() < 100.0 * Eps * l.FrobeniusNorm(),
            $"this row needs a flattened block: span {re.Max() - re.Min():E3}");
        Assert.Equal(traceConstant, re.Average(), 9);
        Assert.NotEqual(-2.0 * gamma.Average(), re.Average(), 3);   // the withdrawn formula must FAIL here
        Assert.True(traceConstant < floor,
            $"the trace constant {traceConstant:F6} must sit strictly below the floor {floor:F6}");
    }

    /// <summary>WHAT THE INVOLUTION ACTUALLY BUYS, which is less than flatness and is why a J threshold
    /// exists at all. The pairing sends rate ↦ 2·mean(rate) − rate, so it makes the block's Re-spectrum
    /// PALINDROMIC about the trace constant, and it does so at EVERY coupling: the defect below is machine
    /// zero at J = 0.05 exactly as at J = 10. FLATNESS is the UNBROKEN case of that symmetry, all
    /// eigenvalues sitting on its axis, and that sets in only above a coupling. Below it the block is
    /// symmetric and wide; above it, symmetric and collapsed. Stating "the involution makes it flatten"
    /// would skip precisely this step.
    ///
    /// <para>The palindrome has no exact route (it is an eigensolver reading), so it is gated against the
    /// backward-error model ε·‖L‖_F rather than a chosen number, and the gate is that the ratio stays
    /// bounded while the SPAN moves over four orders across the same couplings. The span assertions on
    /// either side of the threshold are what keep the two predicates apart.</para></summary>
    [Theory]
    [InlineData(4, 0, 2, 1.0, true)]
    [InlineData(6, 0, 3, 1.0, true)]
    [InlineData(4, 0, 2, 0.0, false)]    // no onset at this Δ, and the palindrome must hold anyway
    [InlineData(6, 0, 3, 0.5, false)]
    public void TheInvolution_GivesAPalindromeAtEveryCoupling_FlatnessIsItsUnbrokenCase(
        int n, int p, int q, double delta, bool hasOnset)
    {
        var gamma = OffLocusProfile(n);
        AssertOffTheR90Locus(n, gamma);
        double constant = -2.0 * Math.Abs(p - q) * gamma.Average();

        double widestBelow = 0.0, widestAbove = 0.0, worstRatio = 0.0;
        foreach (double j in new[] { 0.05, 0.2, 0.5, 10.0 })
        {
            var l = PinnedBlockFloorWitness.ProfileBlock(n, p, q, j, delta, gamma);
            var re = l.Evd().EigenValues.Select(z => z.Real).OrderBy(x => x).ToArray();
            var mirrored = re.Select(x => 2.0 * constant - x).OrderBy(x => x).ToArray();

            double defect = re.Zip(mirrored, (a, b) => Math.Abs(a - b)).Max();
            worstRatio = Math.Max(worstRatio, defect / (Eps * l.FrobeniusNorm()));

            double span = re[^1] - re[0];
            // Above the onset the span is compared to the eigensolver's own error scale, not to a fixed
            // number: ‖L‖_F grows with J, so an absolute cut would quietly become a coupling-dependent one.
            if (j < 1.0) widestBelow = Math.Max(widestBelow, span);
            else widestAbove = Math.Max(widestAbove, span / (Eps * l.FrobeniusNorm()));
        }

        Assert.True(worstRatio < 100.0,
            $"the Re-palindrome about {constant:F4} must hold at EVERY coupling, tracking ε·‖L‖_F: worst ratio {worstRatio:F1}");
        Assert.True(widestBelow > 0.1, $"below the onset the block must be symmetric and WIDE: {widestBelow:E3}");
        if (hasOnset)
            Assert.True(widestAbove < 100.0,
                $"above the onset, symmetric and COLLAPSED: span is {widestAbove:F1}× the backward error");
        else
            Assert.True(widestAbove > 1e6,
                $"at Δ={delta} there is no onset: the block stays wide, yet the palindrome above still holds");
    }

    /// <summary>THE TWO SOURCES ARE NOT THE SAME KIND OF THING, and a uniform-J chain hides it. Source (2)
    /// is the chain REFLECTION, so it needs H to be reflection-invariant too: give the bonds a
    /// non-palindromic J profile and the palindrome breaks outright, six orders above machine. Source (1)
    /// is the global complement, i.e. X^N, which commutes with XX+YY+ZZ bond by bond, so it is blind to
    /// the bond profile and its palindrome survives untouched. Every other gate here runs uniform J, which
    /// is exactly why this scope clause needs its own row.</summary>
    [Fact]
    public void SourceTwoNeedsAReflectionSymmetricH_SourceOneDoesNot()
    {
        const int n = 6;
        double[] palindromicBonds = [5.0, 11.5, 2.0, 11.5, 5.0];
        double[] skewBonds = [5.0, 11.5, 2.0, 8.5, 4.5];

        // Source (2): the R₉₀ locus, on a block with no free-index help (|p−q| = 1, free index 1 ≠ N/2).
        Assert.True(PalindromeDefectRatio(n, 0, 1, palindromicBonds, R90Profile(n)) < 100.0);
        Assert.True(PalindromeDefectRatio(n, 0, 1, skewBonds, R90Profile(n)) > 1e6,
            "a non-palindromic bond profile must break the reflection source outright");

        // Source (1): the free index, off the locus, under the SAME skew bonds.
        Assert.True(PalindromeDefectRatio(n, 0, 3, skewBonds, OffLocusProfile(n)) < 100.0,
            "the complement source is blind to the bond profile");
    }

    /// <summary>The Re-palindrome defect about the trace constant, in units of the eigensolver's backward
    /// error, for a chain with per-bond couplings.</summary>
    private static double PalindromeDefectRatio(int n, int p, int q, double[] bondJ, double[] gamma)
    {
        var terms = new List<PauliTerm>();
        for (int b = 0; b < n - 1; b++)
        {
            terms.Add(PauliTerm.TwoSite(n, b, PauliLetter.X, b + 1, PauliLetter.X, bondJ[b]));
            terms.Add(PauliTerm.TwoSite(n, b, PauliLetter.Y, b + 1, PauliLetter.Y, bondJ[b]));
            terms.Add(PauliTerm.TwoSite(n, b, PauliLetter.Z, b + 1, PauliLetter.Z, bondJ[b]));
        }
        var l = PinnedBlockFloorWitness.ProfileBlock(n, p, q, new PauliHamiltonian(n, terms).ToMatrix(), gamma);
        var re = l.Evd().EigenValues.Select(z => z.Real).OrderBy(x => x).ToArray();
        double constant = -2.0 * Math.Abs(p - q) * gamma.Average();
        var mirrored = re.Select(x => 2.0 * constant - x).OrderBy(x => x).ToArray();
        return re.Zip(mirrored, (a, b) => Math.Abs(a - b)).Max() / (Eps * l.FrobeniusNorm());
    }

    private const double PauliJ = 10.0;

    /// <summary>The one of p, q that is not pinned to 0 or N. Defined only above dimension one, where
    /// exactly one of the two is free; the four one-cell blocks have neither.</summary>
    private static int FreeIndex(int n, int p, int q) => p == 0 || p == n ? q : p;

    private static IEnumerable<(int P, int Q)> AllPinnedAboveDimensionOne(int n) =>
        from p in Enumerable.Range(0, n + 1)
        from q in Enumerable.Range(0, n + 1)
        where PinnedBlockFloorClaim.IsPinned(n, p, q) && PinnedBlockFloorWitness.Cells(n, p, q).Count > 1
        select (p, q);

    /// <summary>Profiles chosen once and asserted off / on the locus in the gates that use them, so the
    /// classification is never a comment.</summary>
    private static double[] OffLocusProfile(int n) => n switch
    {
        3 => [0.31, 0.87, 0.13],
        4 => [0.31, 0.87, 0.13, 0.52],
        5 => [0.31, 0.87, 0.13, 0.52, 0.19],
        6 => [0.31, 0.87, 0.13, 0.52, 0.19, 0.71],
        _ => throw new ArgumentOutOfRangeException(nameof(n)),
    };

    /// <summary>DYADIC anti-palindromic profiles, so the locus membership γ_l + γ_{N−1−l} = const is an
    /// EXACT equality of doubles and not a tolerance: the class's own floating-point paragraph applies
    /// here too, and being ON the locus is a premise the gates lean on.</summary>
    private static double[] R90Profile(int n) => n switch
    {
        4 => [0.25, 0.75, 0.5, 1.0],
        5 => [0.25, 0.75, 0.625, 0.5, 1.0],
        6 => [0.25, 0.75, 0.5, 0.75, 0.5, 1.0],
        _ => throw new ArgumentOutOfRangeException(nameof(n)),
    };

    private static void AssertOffTheR90Locus(int n, double[] gamma)
    {
        var pairSums = Enumerable.Range(0, n).Select(l => gamma[l] + gamma[n - 1 - l]).ToArray();
        // A separation guard, not a derived bound: off-locus-ness has no exact route on a non-dyadic
        // profile, and this only has to exclude a profile that is accidentally near the locus.
        Assert.True(pairSums.Max() - pairSums.Min() > 1e-3,
            $"this profile must be OFF the R₉₀ locus, else the row tests the other source: pair sums {string.Join(", ", pairSums)}");
    }

    /// <summary>Pinned blocks of dimension above one, split by whether the Re-span coalesces, with the
    /// span AND the eigensolver's error scale ε·‖L‖_F carried along. The scale is the point: a span is an
    /// eigensolver output, so "flat" can only mean "at the backward error", and an ABSOLUTE cut silently
    /// becomes a coupling-dependent one, ‖L‖_F growing with J. Flat spans really do reach 1.7e-11 at
    /// J = 1000 while remaining the same phenomenon.</summary>
    private static (List<Reading> Flat, List<Reading> Spread) FlatAndSpreadPinnedBlocks(
        int n, double j, double delta, double[] gamma)
    {
        var flat = new List<Reading>();
        var spread = new List<Reading>();
        var h = PinnedBlockFloorWitness.ChainHamiltonian(n, j, delta);
        foreach (var (p, q) in AllPinnedAboveDimensionOne(n))
        {
            var l = PinnedBlockFloorWitness.ProfileBlock(n, p, q, h, gamma);
            var re = l.Evd().EigenValues.Select(z => z.Real).ToArray();
            double span = re.Max() - re.Min();
            var reading = new Reading((p, q), span, Eps * l.FrobeniusNorm());
            (reading.AtBackwardError ? flat : spread).Add(reading);
        }
        return (flat, spread);
    }

    private const double Eps = 2.220446049250313e-16;

    private readonly record struct Reading((int P, int Q) Block, double Span, double ErrorScale)
    {
        /// <summary>Flat means "at the eigensolver's backward error", never "below a fixed number".</summary>
        public bool AtBackwardError => Span < 100.0 * ErrorScale;
        public double Ratio => Span / ErrorScale;
    }

    /// <summary>THE SCOPE SEAM, gated so it cannot be lost again: the entry-wise identity holds for ANY
    /// Hermitian H, but the SPECTRAL reading needs a NUMBER-CONSERVING one, because only then is the
    /// (p, q) index set an invariant block whose submatrix eigenvalues are the Liouvillian's. With a
    /// transverse field added, Herm(L) = −2·diag(rate) still holds bit-exactly while H itself connects
    /// the block to configurations outside it, so reading a floor off that submatrix would be reading a
    /// number the dynamics does not have.</summary>
    [Fact]
    public void WithATransverseField_TheEntryWiseIdentityHolds_ButTheBlockIsNoLongerInvariant()
    {
        const int n = 4;
        var gamma = DyadicRamp(n);
        var terms = new List<PauliTerm>();
        for (int b = 0; b < n - 1; b++)
        {
            terms.Add(PauliTerm.TwoSite(n, b, PauliLetter.X, b + 1, PauliLetter.X, 1.0));
            terms.Add(PauliTerm.TwoSite(n, b, PauliLetter.Y, b + 1, PauliLetter.Y, 1.0));
            terms.Add(PauliTerm.TwoSite(n, b, PauliLetter.Z, b + 1, PauliLetter.Z, 1.0));
        }
        for (int s = 0; s < n; s++)
            terms.Add(PauliTerm.SingleSite(n, s, PauliLetter.X, 0.7));
        var h = new PauliHamiltonian(n, terms).ToMatrix();

        var block = PerBlockLiouvillianBuilder.BuildBlockZ(h, gamma, PinnedBlockFloorWitness.FlatIndices(n, 0, 1));
        var cells = PinnedBlockFloorWitness.Cells(n, 0, 1);
        double worst = 0.0;
        for (int r = 0; r < cells.Count; r++)
            for (int c = 0; c < cells.Count; c++)
            {
                var herm = 0.5 * (block[r, c] + Complex.Conjugate(block[c, r]));
                if (r == c) herm += 2.0 * PinnedBlockFloorWitness.CellRate(n, cells[r].Ket, cells[r].Bra, gamma);
                worst = Math.Max(worst, herm.Magnitude);
            }
        Assert.True(worst == 0.0, $"the entry-wise identity must survive a non-number-conserving H: {worst:E3}");

        // And the consequence, which is the half that was missing: the submatrix's eigenvalues are then
        // NOT the Liouvillian's. Building the FULL L on all 4^N indices and asking how far each submatrix
        // eigenvalue sits from the nearest true one is the only way to show that reading it is reading a
        // number the dynamics does not have. (Asserting merely that H has popcount-changing entries would
        // re-state the construction three lines above, which no input could break.)
        var full = PerBlockLiouvillianBuilder.BuildBlockZ(h, gamma, Enumerable.Range(0, 1 << (2 * n)).ToArray());
        var trueSpectrum = full.Evd().EigenValues.ToArray();
        double worstDistance = block.Evd().EigenValues
            .Max(z => trueSpectrum.Min(t => (t - z).Magnitude));
        Assert.True(worstDistance > 0.1,
            $"with a transverse field the submatrix spectrum must LEAVE the Liouvillian's: worst distance {worstDistance:E3}");
    }

    /// <summary>The four spectral numbers F153 carries in prose, given a generator at last. Until this
    /// test they existed only as three copies of the same sentence, in the entry, the claim and the
    /// witness, with no code in the repository that produced them. They reproduce only in the PAULI book
    /// (H = J·Σ(XX+YY+Δ·ZZ)), which is WeightCoherenceBlock's book and which the entry never names; in
    /// the spin book next door the same J is four times smaller and none of them lands.
    ///
    /// <para>The last assertion is the correction this gate was built for. The entry reads the flat
    /// constant −2.200 as −2·mean(γ), which is −1.100 for this profile. It is −2·|p−q|·γ̄, which is what
    /// the trace pins, and at |p−q| = 1 the two expressions coincide, which is where the slip came
    /// from.</para></summary>
    [Fact]
    public void TheMeasuredSpectralNumbers_Reproduce_InThePauliBook()
    {
        double[] Ramp(int n) => Enumerable.Range(0, n).Select(k => 0.1 + 0.9 * k / (n - 1.0)).ToArray();

        var (lo4, hi4) = PinnedBlockFloorWitness.ProfileReSpan(4, 0, 1, 1.0, 1.0, Ramp(4));
        Assert.Equal(0.4903, hi4 - lo4, 4);

        var (lo5, hi5) = PinnedBlockFloorWitness.ProfileReSpan(5, 0, 1, 1.0, 1.0, Ramp(5));
        Assert.Equal(0.9477, hi5 - lo5, 4);

        // The Δ=0 re-flattening is not the coupling's doing alone: this ramp is anti-palindromic, i.e. on
        // the R₉₀ locus, which is what lets the whole block coalesce. The span is machine zero and not
        // "exactly 0.0": it is an eigensolver output, so it is read against ε·‖L‖_F rather than a chosen
        // number, and the off-locus control below shows what the locus is buying.
        var (loXY, hiXY) = PinnedBlockFloorWitness.ProfileReSpan(4, 0, 1, 1.0, 0.0, Ramp(4));
        double frob = PinnedBlockFloorWitness.ProfileBlock(4, 0, 1, 1.0, 0.0, Ramp(4)).FrobeniusNorm();
        Assert.True(hiXY - loXY < 100.0 * Eps * frob,
            $"at Δ=0, J=1 on the R₉₀ locus the block coalesces: span {hiXY - loXY:E3} against ε·‖L‖_F = {Eps * frob:E3}");

        var (loOff, hiOff) = PinnedBlockFloorWitness.ProfileReSpan(4, 0, 1, 1.0, 0.0, new[] { 0.13, 0.62, 0.29, 0.94 });
        Assert.True(hiOff - loOff > 0.01,
            $"off the R₉₀ locus the same Δ and J do NOT coalesce: span {hiOff - loOff:E3}");

        var (loSlow, hiSlow) = PinnedBlockFloorWitness.ProfileReSpan(4, 0, 1, 0.05, 0.0, Ramp(4));
        Assert.Equal(1.77, hiSlow - loSlow, 2);

        // the flat-at-the-wrong-constant block, and what that constant actually is
        var gamma = Ramp(4);
        var mean = PinnedBlockFloorWitness.ProfileBlock(4, 0, 2, 1.0, 1.0, gamma)
                                          .Evd().EigenValues.Average(z => z.Real);
        Assert.Equal(-2.2, mean, 9);
        Assert.Equal(-2.0 * 2 * gamma.Average(), mean, 9);
        Assert.NotEqual(-2.0 * gamma.Average(), mean, 3);
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
