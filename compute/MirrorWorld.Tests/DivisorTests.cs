using MirrorWorld;

namespace MirrorWorldTests;

// From-below guard for the frozen divisor (adopted 2026-07-25 from
// docs/proofs/PROOF_R90_FROZEN_DIVISOR.md, registry F140). Everything here is a count, an entry-wise
// residual, or a rank by elimination: no eigensolver anywhere, the Mirror genre.
public class DivisorTests
{
    const long JNum = 90, Den = 100;                 // J = 0.9, rates in hundredths: exact integers
    static readonly World W = new();

    static Divisor Make(int n, long jNum = JNum)
        => new(new Mirror(W, n, 1.0, 0.5), n, jNum, Divisor.Locus(n, 9, 2, -3, 5), Den);

    // ---- the rooms: the shortage IS the mirror's fixed-cell count ----

    [Theory]
    [InlineData(3)]
    [InlineData(4)]
    [InlineData(5)]
    [InlineData(6)]
    [InlineData(7)]
    [InlineData(8)]
    [InlineData(9)]
    [InlineData(12)]
    public void Surplus_Equals_The_Mirrors_Fixed_Rooms(int n)
    {
        var (fixedCells, dimOPlus, dimOMinus, surplus, _, _) = Make(n).Rooms();
        Assert.Equal(fixedCells, surplus);                 // the shortage is the trace
        Assert.Equal(dimOPlus - dimOMinus, surplus);
        Assert.Equal(2 * (n / 2), fixedCells);             // one per site, minus the odd-N centre
    }

    // the frozen count is what survives the diagonal tax: one mode per balanced pair.
    [Theory]
    [InlineData(3, 1)]
    [InlineData(4, 2)]
    [InlineData(5, 2)]
    [InlineData(6, 3)]
    [InlineData(7, 3)]
    [InlineData(10, 5)]
    [InlineData(15, 7)]
    public void Frozen_Is_Surplus_Minus_Tax(int n, int expected)
    {
        var (_, _, _, surplus, tax, frozen) = Make(n).Rooms();
        Assert.Equal(expected, frozen);
        Assert.Equal(n / 2, frozen);
        Assert.Equal(surplus - tax, frozen);
    }

    // the odd-N centre site is invisible to the divisor: its self-coherence is a POPULATION, not a
    // coherence, so it never enters O. Hence N and N-1 (N odd) hold identical room counts.
    [Theory]
    [InlineData(5)]
    [InlineData(7)]
    [InlineData(9)]
    [InlineData(11)]
    public void Odd_N_Counts_The_Same_As_The_Even_N_Below(int n)
    {
        var odd = Make(n).Rooms();
        var even = Make(n - 1).Rooms();
        Assert.Equal(even.Fixed, odd.Fixed);
        Assert.Equal(even.Frozen, odd.Frozen);
    }

    // ---- the hypothesis, entry by entry ----

    [Theory]
    [InlineData(3)]
    [InlineData(4)]
    [InlineData(5)]
    [InlineData(6)]
    [InlineData(7)]
    public void On_The_Locus_Both_Parts_Are_Odd_Under_The_Cell_Mirror(int n)
    {
        var (hop, rate) = Make(n).OddnessResidual();
        Assert.True(hop < 1e-12, $"hop residual {hop:E2}");
        Assert.True(rate < 1e-12, $"rate residual {rate:E2}");
    }

    // the control that matters: step OFF the locus and the rate part stops being odd. The hop part
    // never cared (it is odd for any R-symmetric h), so this is exactly what the locus buys.
    [Theory]
    [InlineData(4)]
    [InlineData(5)]
    [InlineData(6)]
    public void Off_The_Locus_The_Rate_Part_Stops_Being_Odd(int n)
    {
        var g = Divisor.Locus(n, 9, 2, -3, 5);
        g[0] += 17;                                        // break one reflection pair only
        var off = new Divisor(new Mirror(W, n, 1.0, 0.5), n, JNum, g, Den);
        var (hop, rate) = off.OddnessResidual();
        Assert.True(hop < 1e-12, "the hop part is odd regardless of the profile");
        Assert.True(rate > 1e-3, $"the rate part must break off the locus, got {rate:E2}");
    }

    // ---- the multiplicity, by rank, never by a spectrum ----

    [Theory]
    [InlineData(3)]
    [InlineData(4)]
    [InlineData(5)]
    [InlineData(6)]
    [InlineData(7)]
    [InlineData(8)]
    public void Kernel_Dimension_Matches_The_Room_Count(int n)
    {
        var d = Make(n);
        Assert.Equal(d.Rooms().Frozen, d.KernelDimension());
    }

    // the coupling cannot move it: the same count at wildly different J.
    // The frozen root is exact at EVERY coupling, so the count must not move -- including where a
    // floating-point rank would collapse: small J and a long chain, the corner the other eigenvalues
    // crowd at spacing J^(2d). A float rank returned 2 instead of 10 at N=20, J=1/1000 (review catch);
    // the exact GF(p) rank does not.
    [Theory]
    [InlineData(1)]          // J = 1/100
    [InlineData(100)]        // J = 1
    [InlineData(750)]        // J = 7.5
    public void The_Count_Does_Not_Depend_On_The_Coupling(long jNum)
    {
        for (int n = 3; n <= 10; n++)
            Assert.Equal(n / 2, Make(n, jNum).KernelDimension());
    }

    [Theory]
    [InlineData(12, 1)]      // the exact corner a floating-point rank got wrong
    [InlineData(16, 1)]
    [InlineData(20, 1)]
    [InlineData(20, 10)]
    public void Small_Coupling_On_A_Long_Chain_Still_Counts_Right(int n, long jNum)
        => Assert.Equal(n / 2, Make(n, jNum).KernelDimension());

    // ---- the zero-mean stratum: where the tax is not there to be paid ----
    // The locus never asks the site rates to be positive, so gbar = 0 is a point of it, and there the
    // even defect 4*gbar*P_D is the zero matrix: the whole block is odd, the diagonal cells charge
    // nothing, and the count is N (one per SITE) rather than floor(N/2) (one per balanced pair). The
    // room count has to say so on its own, without the rank being consulted.
    static Divisor Zero(int n, long jNum = JNum)
        => new(new Mirror(W, n, 1.0, 0.5), n, jNum, Divisor.Locus(n, 0, 2, -3, 5), Den);

    [Theory]
    [InlineData(3)]
    [InlineData(4)]
    [InlineData(5)]
    [InlineData(6)]
    [InlineData(7)]
    [InlineData(8)]
    public void Zero_Mean_Is_Still_The_Locus(int n)
    {
        var d = Zero(n);
        Assert.Equal(0.0, d.Sigma);
        Assert.Equal(0.0, d.GBar);
        Assert.Equal(0.0, d.Root);
        Assert.Equal(d.Root, d.FoldRoot);                  // the fold partner collapses onto it
        Assert.True(d.OddnessResidual().Rate < 1e-12);     // still on the locus, entry by entry
    }

    // the mechanism, and the read that can actually SEE it. OddnessResidual skips the diagonal
    // cells, so it reads machine zero on both strata and discriminates nothing; the whole-block
    // residual is zero only where the even defect is absent, and O(gbar) where it is not.
    [Theory]
    [InlineData(3)]
    [InlineData(4)]
    [InlineData(5)]
    [InlineData(6)]
    [InlineData(7)]
    [InlineData(8)]
    public void Only_The_Zero_Mean_Block_Is_Odd_All_The_Way_Through(int n)
    {
        Assert.True(Zero(n).WholeBlockOddnessResidual() < 1e-12);
        // not merely nonzero: the defect is exactly 8*gbar, the proof's own G1 constant.
        Assert.Equal(8 * 0.09, Make(n).WholeBlockOddnessResidual(), 12);
        Assert.True(Make(n).OddnessResidual().Rate < 1e-12);   // the blind read, both strata
        Assert.True(Zero(n).OddnessResidual().Rate < 1e-12);
    }

    [Theory]
    [InlineData(3)]
    [InlineData(4)]
    [InlineData(5)]
    [InlineData(6)]
    [InlineData(7)]
    [InlineData(8)]
    public void Zero_Mean_Charges_No_Tax_And_Freezes_N(int n)
    {
        var (_, _, _, surplus, tax, frozen) = Zero(n).Rooms();
        Assert.Equal(2 * (n / 2), surplus);                 // the surplus is the mirror's, unchanged
        Assert.Equal(-(n % 2), tax);                        // the diagonal cells PAY here, they do
        Assert.Equal(surplus - tax, frozen);                // not charge; the identity still holds
        Assert.Equal(n, frozen);                            // one per site, and NOT 2*floor(n/2)
        var taxed = Make(n).Rooms();
        Assert.Equal(n / 2, taxed.Tax);                     // the taxed stratum is untouched
        Assert.Equal(taxed.Surplus - taxed.Tax, taxed.Frozen);
        Assert.Equal(n / 2, taxed.Frozen);
    }

    [Theory]
    [InlineData(3)]
    [InlineData(4)]
    [InlineData(5)]
    [InlineData(6)]
    [InlineData(7)]
    [InlineData(8)]
    public void Zero_Mean_Kernel_Matches_The_Untaxed_Room_Count(int n)
    {
        var d = Zero(n);
        Assert.Equal(n, d.KernelDimension());
        Assert.Equal(d.Rooms().Frozen, d.KernelDimension());
    }

    // and the coupling still cannot move it, on this stratum either.
    [Theory]
    [InlineData(1)]          // J = 1/100
    [InlineData(100)]        // J = 1
    [InlineData(750)]        // J = 7.5
    public void Zero_Mean_Count_Does_Not_Depend_On_The_Coupling(long jNum)
    {
        for (int n = 3; n <= 9; n++)
            Assert.Equal(n, Zero(n, jNum).KernelDimension());
    }

    // past the wall: the spectrum died at N=8, the corner block is N^2, the divisor does not notice.
    [Theory]
    [InlineData(10)]
    [InlineData(12)]
    [InlineData(16)]
    public void It_Walks_Past_The_Wall(int n)
    {
        var d = Make(n);
        Assert.Equal(n / 2, d.KernelDimension());
        Assert.Equal(n / 2, d.Rooms().Frozen);
    }

    // ---- the four corners, by fold parity ----

    [Theory]
    [InlineData(5)]
    [InlineData(6)]
    [InlineData(7)]
    public void Four_Corners_With_The_Root_Chosen_By_Fold_Parity(int n)
    {
        var d = Make(n);
        var corners = d.Corners();
        Assert.Equal(4, corners.Length);
        foreach (var (p, q, folds, root) in corners)
        {
            Assert.True(p == 1 || p == n - 1);
            Assert.True(q == 1 || q == n - 1);
            Assert.Equal(folds % 2 == 0 ? d.Root : d.FoldRoot, root, 12);
        }
        // each one-sided fold sends r to -r - 2*sigma, so two of them return to where they started
        Assert.Equal(d.Root, -d.FoldRoot - 2 * d.Sigma, 12);
    }

    // ---- the ladder ----

    [Theory]
    [InlineData(4)]
    [InlineData(5)]
    [InlineData(6)]
    [InlineData(7)]
    [InlineData(8)]
    public void The_Ladder_Is_Twice_The_Site_Distance_And_Sums_To_The_Valuation(int n)
    {
        var (dist, per, total) = Make(n).Ladder();
        Assert.Equal(n / 2, dist.Length);
        for (int c = 1; c <= dist.Length; c++)
        {
            Assert.Equal(n + 1 - 2 * c, dist[c - 1]);      // the pair's site distance
            Assert.Equal(2 * dist[c - 1], per[c - 1]);     // twice it, one power per hop each way
        }
        Assert.Equal(2 * (n * n / 4), total);
        Assert.Equal(per.Sum(), total);
        // the outermost pair is the last to let go
        Assert.Equal(per.Max(), per[0]);
    }

    [Theory]
    [InlineData(6)]
    [InlineData(9)]
    public void The_Boundary_Clock_Is_N_Plus_1_Without_The_ZZ_Diagonal_And_N_With_It(int n)
    {
        var mir = new Mirror(W, n, 1.0, 0.5);
        var g = Divisor.Locus(n, 9, 2, -3, 5);
        Assert.Equal(n + 1, new Divisor(mir, n, JNum, g, Den).ClockModulus);
        Assert.Equal(n, new Divisor(mir, n, JNum, g, Den, zz: true).ClockModulus);
    }

    // ---- the ontology: the divisor cannot hold its own reason ----

    [Fact]
    public void The_Rooms_Are_Inherited_From_The_Mirror_Not_Owned()
    {
        var d = Make(6);
        Assert.Equal(new[] { "frozen", "root", "ladder" }, d.Own);
        // neither count that goes into the subtraction is the divisor's: both are tauQ fixed-counts
        Assert.DoesNotContain("rooms", d.Own);
        Assert.DoesNotContain("surplus", d.Own);
        Assert.DoesNotContain("tax", d.Own);
        // the parent is the mirror, not the frame (Marginal reached a non-frame parent first)
        Assert.IsType<Mirror>(d.Parent);
        foreach (var leg in new[] { "legs", "price", "orbit" })
            Assert.Contains(leg, d.Inherited);
        foreach (var frame in new[] { "x", "y", "z" })
            Assert.Contains(frame, d.Inherited);           // the frame still arrives, one step further up
    }
}
