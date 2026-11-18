using RCPsiSquared.Diagnostics.Foundation;

namespace RCPsiSquared.Diagnostics.Tests.Foundation;

/// <summary>From-below pins for <see cref="DephasedSeatFlipBlockWitness"/> (the bullet "The block from one excitation to
/// N − 1, at every N" of <c>docs/proofs/PROOF_ONE_DEPHASED_SEAT.md</c>): rows the proof and its gate
/// <c>simulations/blind_seat_block_lemma_gate.py</c> read, the grids, and a planted mutation (the jump's equation
/// dropped) that the block reading must catch.</summary>
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
        var g4 = DephasedSeatFlipBlockWitness.ReadFarGrid(4, 20);     // measured, as the proof measures it
        Assert.True(g4.NotBlindOffLine > 0);
        Assert.Equal(0, g4.OffLineWithElement);
        Assert.Equal(g4.NotBlindOnLine, g4.OnLineWithElement);
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
        Assert.Equal(5, w.Children.Count());
    }
}
