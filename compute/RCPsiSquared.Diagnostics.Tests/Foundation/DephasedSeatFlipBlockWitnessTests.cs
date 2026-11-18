using RCPsiSquared.Diagnostics.Foundation;

namespace RCPsiSquared.Diagnostics.Tests.Foundation;

/// <summary>From-below pins for <see cref="DephasedSeatFlipBlockWitness"/> (the bullets "The block from one excitation to
/// N − 1, at every N" and "At N = 4, a seat that is not blind" of <c>docs/proofs/PROOF_ONE_DEPHASED_SEAT.md</c>): rows the
/// proof and its gates <c>simulations/blind_seat_block_lemma_gate.py</c> and <c>simulations/blind_seat_descent_gate.py</c>
/// read, the grids, a planted mutation (the jump's equation dropped) that the block reading must catch, and readings
/// that the N = 4 rule and condition (II) must call disagreements.</summary>
public class DephasedSeatFlipBlockWitnessTests
{
    static DephasedSeatFlipBlockWitness W(int n, string bonds, string fields = "") => new(n, bonds, fields);

    [Theory]
    [InlineData(3, "0-1:1,0-2:2,1-2:3", 0, 1)]                  // generic triangle: not blind, X^3's block alone
    [InlineData(3, "0-1:3,0-2:6,1-2:-2", 1, 2)]                 // C's couplings, e2 = 0: blind
    [InlineData(3, "0-1:2,0-2:2,1-2:3", 1, 2)]                  // B's couplings, J_du = J_dv: blind
    [InlineData(3, "0-1:3,1-2:-2", 0, 1)]                       // a path with d at an end
    [InlineData(4, "0-1:1,0-2:-1,1-3:-1,2-3:1,1-2:2", 1, 2)]    // the N = 4 family's couplings (a = 1, e = 2)
    [InlineData(4, "0-1:1,0-2:2,0-3:-1,1-2:3,2-3:1", 0, 1)]
    [InlineData(5, "0-1:2,1-2:-1,2-3:3,3-4:1,0-4:2,1-3:1", 1, 2)]   // blind at these couplings: the commutant count
    public void Pinned_Blocks(int n, string bonds, int blind, int block)
    {
        var r = W(n, bonds).Block;
        Assert.Equal((blind, block), (r.Blind, r.BlockDimension));
        Assert.Equal(r.BlockDimension, r.CommutantDimension);
        Assert.True(r.XiSolves);
        Assert.True(r.Agrees);
    }

    /// <summary>The mutation: without the jump's equation the bond condition alone admits more blocks at a seat that
    /// is not blind, so the jump is what pins the block to X^N.</summary>
    [Fact]
    public void Dropping_The_Jump_Breaks_The_Lemma()
    {
        var bonds = DephasedSeatFlipBlockWitness.ParseBonds("0-1:1,0-2:2,1-2:3", 3);
        var r = DephasedSeatFlipBlockWitness.ReadBlock(3, bonds, withoutJump: true);
        Assert.True(r.BlockDimension > 1);
        Assert.False(r.Agrees);
        Assert.True(DephasedSeatFlipBlockWitness.ReadBlockGrid(4, withoutJump: true).Disagreements > 0);
    }

    [Theory]
    [InlineData(3)]
    [InlineData(4)]
    [InlineData(5)]
    public void Block_Grid_Agrees_On_Both_Kinds(int n)
    {
        var g = DephasedSeatFlipBlockWitness.ReadBlockGrid(n);
        Assert.True(g.NotBlind > 0 && g.Blind > 0, $"N = {n}: the grid must hold rows of both kinds");
        Assert.Equal(0, g.Disagreements);
    }

    [Theory]
    [InlineData("0-1:1,0-2:2,1-2:3", "1:3,4;2:5,0", 0, 0, false)]     // not blind, fields off one line: none
    [InlineData("0-1:1,0-2:2,1-2:3", "1:3,4;2:-3,-4", 1, 0, true)]    // antiparallel: class A's colouring, odd
    [InlineData("0-1:2,0-2:-1,1-2:1", "1:1,1;2:2,2", 1, 0, true)]     // parallel
    [InlineData("0-1:1,0-2:2,1-2:3", "", 4, 0, true)]                 // no field: (m·σ)^3 for every planar m spans four, odd; even none from bonds alone
    public void Far_At_N3_Not_Blind(string bonds, string fields, int odd, int even, bool line)
    {
        var f = W(3, bonds, fields).Far;
        Assert.Equal((odd, even, line), (f.Odd, f.Even, f.OnLine));
        Assert.True(f.Agrees);
    }

    /// <summary>The gate's reading beside the proof: on e2 = 0 with J_du != J_dv four even far elements without
    /// fields, none with f_u = (3, 4), f_v = (5, 0); the seat is blind there and nothing is predicted.</summary>
    [Fact]
    public void Reciprocal_Row_Even_Height_Reading()
    {
        var bare = W(3, "0-1:3,0-2:6,1-2:-2").Far;
        var lit = W(3, "0-1:3,0-2:6,1-2:-2", "1:3,4;2:5,0").Far;
        Assert.Equal(4, bare.Even);
        Assert.Equal(0, lit.Even);
        Assert.Null(bare.Agrees);
        Assert.Equal(1, lit.Odd);                                    // C's far element, of odd height
    }

    [Fact]
    public void Far_Grids()
    {
        var g3 = DephasedSeatFlipBlockWitness.ReadFarGrid(3, 30);
        Assert.True(g3.NotBlindOffLine > 0 && g3.NotBlindOnLine > 0);
        Assert.Equal(0, g3.Disagreements);
        Assert.Equal(0, g3.OffLineWithElement);
        Assert.Equal(g3.NotBlindOnLine, g3.OnLineWithElement);
        var g4 = DephasedSeatFlipBlockWitness.ReadFarGrid(4, 20);     // predicted by the bullet at N = 4
        Assert.True(g4.NotBlindOffLine > 0);
        Assert.Equal(0, g4.Disagreements);
        Assert.Equal(0, g4.OffLineWithElement);
        Assert.Equal(g4.NotBlindOnLine, g4.OnLineWithElement);
    }

    /// <summary>The bullet "At N = 4, a seat that is not blind": off one line none of either height, with all three
    /// fields nonzero and with one zero; on one line, not all zero, exactly class A's colouring, of even height, at a row
    /// off the antiparallel locus. Both conditions of the descent hold at each row.</summary>
    [Theory]
    [InlineData("1:3,4;2:5,0;3:1,1", 0, 0, true)]
    [InlineData("1:3,4;2:5,0", 0, 0, true)]                         // site 3 without a field, 1 and 2 off one line
    [InlineData("1:1,0;2:-2,0;3:3,0", 1, 0, true)]                  // class A, all three nonzero
    [InlineData("1:1,0;2:2,0", 1, 0, true)]                         // on one line with a field zero
    [InlineData("1:1,0;2:-1,0", 1, 0, true)]                        // opposite fields, the bonds off the locus
    [InlineData("2:3,4", 1, 0, true)]                               // one field alone
    public void Far_At_N4_Not_Blind(string fields, int even, int odd, bool predicted)
    {
        var w = W(4, "0-1:1,0-2:2,0-3:-1,1-2:3,2-3:1", fields);
        Assert.Equal((even, odd), (w.Far.Even, w.Far.Odd));
        Assert.Equal(predicted, w.Far.Predicted);
        Assert.Equal(predicted ? true : null, w.Far.Agrees);
        Assert.Equal(16, w.Descent.Cyclic);
        Assert.Equal(new[] { 0, 0, 0 }, w.Descent.Fence);
        Assert.True(w.Descent.FenceAgrees);
    }

    /// <summary>With no field at all nothing is predicted at N = 4 (X^⊗4 and Y^⊗4 are both far elements there).</summary>
    [Fact]
    public void No_Field_Is_Not_Predicted_At_N4()
    {
        var w = W(4, "0-1:1,0-2:2,0-3:-1,1-2:3,2-3:1");
        Assert.False(w.Far.AnyField);
        Assert.False(w.Far.Predicted);
        Assert.Null(w.Far.Agrees);
        Assert.True(w.Far.Even >= 2);
    }

    /// <summary>The antiparallel locus (the bullet's steps (vi) and (vii)): one site without a field, opposite fields on
    /// the other two, and the field-free site coupled to them exactly as 0 is, oppositely to each, J03 ≠ 0. Seat 0 is not
    /// blind, the cyclic space has dimension 14, and the far elements are two of even height, class A's colouring and
    /// the hop on the triplet state; at J03 = 0 the seat turns blind and nothing is predicted. The labels are free.</summary>
    [Theory]
    [InlineData("0-1:-1,0-2:1,1-3:-1,2-3:1,0-3:2,1-2:3", "1:1,0;2:-1,0", true)]
    [InlineData("0-1:2,0-2:-2,1-3:2,2-3:-2,0-3:-1", "1:3,4;2:-3,-4", true)]
    [InlineData("0-2:1,1-2:1,0-3:-1,1-3:-1,0-1:2,2-3:3", "2:1,1;3:-1,-1", true)]   // site 1 without a field
    [InlineData("0-1:-1,0-2:1,1-3:-1,2-3:1,1-2:3", "1:1,0;2:-1,0", false)]          // J03 = 0: blind
    public void Antiparallel_Locus_Carries_Two(string bonds, string fields, bool notBlind)
    {
        var w = W(4, bonds, fields);
        Assert.True(w.Far.OnLocus);
        Assert.Equal(notBlind, w.Far.Blind == 0);
        if (!notBlind) { Assert.Null(w.Far.Agrees); return; }
        Assert.Equal((2, 0), (w.Far.Even, w.Far.Odd));
        Assert.Equal(true, w.Far.Agrees);
        Assert.Equal(14, w.Descent.Cyclic);
        Assert.Equal(new[] { 0, 0, 0 }, w.Descent.Fence);
    }

    /// <summary>The gate's control: a blind seat bonded to two sites, where |1⟩_0 ⊗ (J01·e1 − J02·e2) is an eigenvector of
    /// H_b killed by F⁻, so the cyclic space misses exactly that line (15 = 16 − 1, the vector bounding it from above).</summary>
    [Fact]
    public void Blind_Control_Cyclic_Space_Misses_A_Line()
    {
        var w = W(4, "0-1:1,0-2:2,1-2:-2,1-3:4,2-3:2", "1:2,0;2:1,0;3:0,1");
        Assert.Equal(1, w.Far.Blind);
        Assert.Equal(15, w.Descent.Cyclic);
        Assert.Null(w.Far.Agrees);
    }

    /// <summary>The bullet's N = 5 row: seat 0 not blind, the second condition fails at three excitations (Ξ applied to
    /// |0⟩_0 ⊗ 𝓈12 ⊗ 𝓈34), the first holds, and still no homogeneous far element, at fields antiparallel within the pairs.</summary>
    [Fact]
    public void N5_Row_Where_The_Fence_Fails_Carries_None()
    {
        var w = W(5, "0-1:1,0-2:1,0-3:2,0-4:2,1-3:1,2-4:2,2-3:3,1-2:-1,3-4:-2", "1:1,0;2:-1,0;3:0,1;4:0,-1");
        Assert.Equal(0, w.Far.Blind);
        Assert.Equal(new[] { 0, 0, 1, 0 }, w.Descent.Fence);
        Assert.True(w.Descent.FenceAgrees);                         // sectors 1, 2 and N - 1 empty, as proven
        Assert.Equal(32, w.Descent.Cyclic);
        Assert.Equal((0, 0), (w.Far.Odd, w.Far.Even));
        Assert.Null(w.Far.Agrees);
    }

    /// <summary>Anti-vacuity: readings that break the N = 4 rule are called disagreements, and at a blind seat with a
    /// kernel vector of 𝓛 vanishing at 0 (the triangle 1, 2, 3 on e2 = 0, 0 a leaf on 1) the sector N − 1 holds a
    /// vector with 0 excited at the vacuum's energy, which FenceAgrees would reject at a seat that is not blind.</summary>
    [Fact]
    public void Disagreements_Are_Seen()
    {
        Assert.Equal(false, new DephasedSeatFlipBlockWitness.FarReading(4, 0, false, true, false, 0, 1).Agrees);
        Assert.Equal(false, new DephasedSeatFlipBlockWitness.FarReading(4, 0, true, true, false, 1, 1).Agrees);
        Assert.Equal(false, new DephasedSeatFlipBlockWitness.FarReading(4, 0, true, true, false, 2, 0).Agrees);
        Assert.Equal(false, new DephasedSeatFlipBlockWitness.FarReading(4, 0, true, true, true, 0, 1).Agrees);
        Assert.Equal(false, new DephasedSeatFlipBlockWitness.FarReading(4, 0, true, true, true, 1, 0).Agrees);
        Assert.Null(new DephasedSeatFlipBlockWitness.FarReading(4, 0, true, false, false, 0, 2).Agrees);
        var d = W(4, "0-1:1,1-2:3,1-3:6,2-3:-2").Descent;
        Assert.Equal(new[] { 0, 0, 1 }, d.Fence);
        Assert.False(new DephasedSeatFlipBlockWitness.DescentReading(0, d.Cyclic, d.Full, d.Fence).FenceAgrees);
    }

    [Theory]
    [InlineData(4)]
    [InlineData(5)]
    public void Fence_Empty_At_Sectors_One_Two_And_N_Minus_One_On_The_Grid(int n)
    {
        int notBlind = 0;
        foreach (var b in DephasedSeatFlipBlockWitness.BondGrid(n, 12, 99UL))
        {
            var r = DephasedSeatFlipBlockWitness.ReadDescent(n, b, new Dictionary<int, (RCPsiSquared.Core.Numerics.BigRational, RCPsiSquared.Core.Numerics.BigRational)>());
            if (r.Blind == 0) notBlind++;
            Assert.True(r.FenceAgrees);
        }
        Assert.True(notBlind > 0);
    }

    [Fact]
    public void Arguments_Are_Validated()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => new DephasedSeatFlipBlockWitness(6));
        Assert.Throws<ArgumentException>(() => W(3, "0-1:1,0-2:2", "0:1,0"));
        Assert.Throws<ArgumentException>(() => W(3, "0-1:1"));
        Assert.Throws<ArgumentException>(() => W(3, "0-1:1,0-1:2,1-2:1"));
    }

    [Fact]
    public void Default_Row_Reads()
    {
        var w = new DephasedSeatFlipBlockWitness();
        Assert.Contains("holds", w.Summary);
        Assert.Equal(6, w.Children.Count());
    }
}
