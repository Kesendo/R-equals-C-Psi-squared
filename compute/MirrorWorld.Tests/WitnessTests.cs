using MirrorWorld;

namespace MirrorWorldTests;

// The witness reading (F135 + F136) pinned from below: every family, letter, sign and price the
// adopted laws force, at the numbers the proof's gates measured (docs/proofs/PROOF_RECORD_LETTER_LAW.md,
// 87/87). The class computes pure cosine-parity arithmetic; these tests pin it against the proof's
// exact values, including the two review-round regions (even dressers in the pointer family; the
// pendant role-swap behind the hinge) and the fan-out sighting (anti-pointer redundancy is not
// bounded by deg(S)).
public class WitnessTests
{
    static readonly World W = new();
    const double G = 0.05;                                   // canonical gamma
    const double PointerPriced = 0.7680396679759238;         // 1 - h2((1+e^{-2g t*})/2), the F135 number
    const double BellPriced = 0.6241464301142645;            // 1 - h2((1+e^{-4g t*})/2), both sites pay

    static readonly (int, int, double)[] Chain3 = { (0, 1, 1.0), (1, 2, 1.0) };
    static readonly (int, int, double)[] TriangleOdd = { (0, 1, 1.0), (0, 2, 1.0), (1, 2, 1.0) };
    static readonly (int, int, double)[] TriangleEven = { (0, 1, 1.0), (0, 2, 1.0), (1, 2, 2.0) };
    static readonly (int, int, double)[] Square = { (0, 1, 1.0), (1, 2, 1.0), (2, 3, 1.0), (3, 0, 1.0) };

    [Fact]
    public void The_Leaf_Records_The_Pointer_And_The_Pendant_Reads_Backwards()
    {
        var leaf = new Witness(W, Chain3, s: 1, j: 0);       // S interior, j the unwatched leaf
        Assert.Equal(Witness.Kind.Pointer, leaf.Family);
        Assert.Equal(+1, leaf.Sign);
        Assert.Equal(1.0, leaf.Bits, 12);
        var swap = new Witness(W, Chain3, s: 0, j: 1);       // S the pendant: the role swap
        Assert.Equal(Witness.Kind.RoleSwap, swap.Family);
        Assert.Equal('Z', swap.Letter);                      // j's pointer, recorded in S's equator (YZ)
        Assert.Equal(1.0, swap.Bits, 12);
    }

    [Fact]
    public void The_Triangle_Splits_By_Its_Far_Bonds_Parity()
    {
        var odd = new Witness(W, TriangleOdd, 0, 1);         // far bond 1: the m = 1 Bell row
        Assert.Equal(Witness.Kind.Bell, odd.Family);
        Assert.Equal('Y', odd.Letter);
        Assert.Equal(+1, odd.Sign);
        var even = new Witness(W, TriangleEven, 0, 1);       // far bond 2: a pointer record, rotated by pi
        Assert.Equal(Witness.Kind.Pointer, even.Family);
        Assert.Equal(-1, even.Sign);                         // Prod_{D u P} (-1)^{r/2} at r = 2
        Assert.Equal(1.0, even.Bits, 12);
    }

    [Fact]
    public void The_Letter_Is_The_Dresser_Parity_With_Closed_Form_Signs()
    {
        var plaquette = new Witness(W, Square, 0, 2);        // m = 2 -> XX
        Assert.Equal(Witness.Kind.Bell, plaquette.Family);
        Assert.Equal('X', plaquette.Letter);
        Assert.Equal(+1, plaquette.Sign);
        var k21r3 = new Witness(W, new[] { (0, 2, 1.0), (1, 2, 3.0) }, 0, 1);   // m = 1, r = 3
        Assert.Equal('Y', k21r3.Letter);
        Assert.Equal(-1, k21r3.Sign);                        // sigma_2 = (-1)^{(1-3)/2}
        var k22mixed = new Witness(W, new[] { (0, 2, 1.0), (0, 3, 1.0), (1, 2, 1.0), (1, 3, 3.0) }, 0, 1);
        Assert.Equal('X', k22mixed.Letter);                  // m = 2 -> XX, sign carries the odd-3 dresser
        Assert.Equal(-1, k22mixed.Sign);
        var dressedPrivate = new Witness(W, new[] { (0, 1, 1.0), (0, 2, 1.0), (1, 2, 1.0), (1, 3, 2.0) }, 0, 1);
        Assert.Equal(Witness.Kind.Bell, dressedPrivate.Family);   // private r = 2 forgives with a pi rotation
        Assert.Equal(-1, dressedPrivate.Sign);
    }

    [Fact]
    public void Integer_Violations_Are_Exact_Kills()
    {
        var pentagon = new[] { (0, 1, 1.0), (1, 2, 1.0), (2, 3, 1.0), (3, 4, 1.0), (4, 0, 1.0) };
        Assert.Equal(Witness.Kind.Dark, new Witness(W, pentagon, 0, 1).Family);  // odd private watcher
        Assert.Equal(Witness.Kind.Dark, new Witness(W, pentagon, 0, 2).Family);  // nonempty Q
        Assert.Equal(0.0, new Witness(W, pentagon, 0, 1).Bits, 12);
        var hexagon = new[] { (0, 1, 1.0), (1, 2, 1.0), (2, 3, 1.0), (3, 4, 1.0), (4, 5, 1.0), (5, 0, 1.0) };
        Assert.Equal(Witness.Kind.Dark, new Witness(W, hexagon, 0, 3).Family);
        var mixedD = new Witness(W, new[] { (0, 2, 1.0), (0, 3, 1.0), (1, 2, 1.0), (1, 3, 2.0) }, 0, 1);
        Assert.Equal(Witness.Kind.Dark, mixedD.Family);      // mixed-parity dressers kill both families
        var noWriter = new Witness(W, new[] { (0, 2, 1.0), (0, 4, 1.0), (1, 3, 2.0) }, 0, 1);
        Assert.Equal(Witness.Kind.Dark, noWriter.Family);    // D = empty, P even, no write bond: no writer
    }

    [Fact]
    public void The_Hinge_Decides_And_The_Pendant_Splits_By_Watcher_Parity()
    {
        var entangled = new Witness(W, new[] { (0, 1, 1.0), (1, 2, 2.0) }, 0, 1);
        Assert.Equal(Witness.Kind.Entangled, entangled.Family);   // even watcher keeps I = 2 exactly
        Assert.Equal(2.0, entangled.Bits, 12);
        var restored = new Witness(W, new[] { (0, 1, 1.0), (1, 2, 1.0), (0, 3, 1.0) }, 0, 1);
        Assert.Equal(Witness.Kind.Dark, restored.Family);         // any second S-neighbor closes the channel
        var mixed = new Witness(W, new[] { (0, 1, 1.0), (1, 2, 1.0), (1, 3, 2.0) }, 0, 1);
        Assert.Equal(Witness.Kind.RoleSwap, mixed.Family);        // EXISTENTIAL: one odd watcher suffices
        Assert.Equal(1.0, mixed.Bits, 12);
        var mixedNonInt = new Witness(W, new[] { (0, 1, 1.0), (1, 2, 1.0), (1, 3, 1.5) }, 0, 1);
        Assert.Equal(Witness.Kind.RoleSwap, mixedNonInt.Family);  // non-integer bystanders do not spoil it
        var offContract = new Witness(W, new[] { (0, 1, 2.0), (1, 2, 3.0) }, 0, 1);
        Assert.Equal(Witness.Kind.Generic, offContract.Family);   // write bond != Delta_S: not classified here
        var entangledPriced = new Witness(W, new[] { (0, 1, 1.0), (1, 2, 2.0) }, 0, 1, gammaS: 0.2, gammaJ: 0.2);
        Assert.True(double.IsNaN(entangledPriced.Bits));          // the entangled face is pinned at gamma = 0 only
    }

    [Fact]
    public void The_Prices_Are_Exact_And_Asymmetric()
    {
        var bell = new Witness(W, Square, 0, 2, gammaS: G, gammaJ: G);
        Assert.Equal(BellPriced, bell.Bits, 10);             // the Bell record pays BOTH sites
        var pointer = new Witness(W, Chain3, 1, 0, gammaS: 0.0, gammaJ: G);
        Assert.Equal(PointerPriced, pointer.Bits, 10);       // the pointer record pays the witness
        var pointerFreeS = new Witness(W, Chain3, 1, 0, gammaS: 0.3, gammaJ: 0.0);
        Assert.Equal(1.0, pointerFreeS.Bits, 12);            // gamma_S exactly invisible to it (F135 Law C face)
        var swapFreeJ = new Witness(W, Chain3, 0, 1, gammaS: 0.0, gammaJ: 0.3);
        Assert.Equal(1.0, swapFreeJ.Bits, 12);               // and the swap ignores the swapped side's gamma
        var swapPriced = new Witness(W, Chain3, 0, 1, gammaS: G, gammaJ: 0.0);
        Assert.Equal(PointerPriced, swapPriced.Bits, 10);    // the swap pays its witness's rate: same formula
        var bellOneSided = new Witness(W, Square, 0, 2, gammaS: G, gammaJ: 0.0);
        Assert.True(bellOneSided.Bits < 1.0);                // even one-sided watching prices the Bell record
    }

    [Fact]
    public void The_Anti_Pointer_Fans_Out_Beyond_Law_Bs_Bound()
    {
        // K_{4,2}: S = 0 and corners 1, 2, 3 all share the two dressers 4, 5. Law B bounds pointer
        // redundancy by deg(S) = 2; the Bell family ignores that bound: three perfect XX bits, and
        // the corners Bell-record each other (the clique). Gated from below in the proof's substrate
        // before this class printed it.
        var k42 = Enumerable.Range(0, 4).SelectMany(i => new[] { (i, 4, 1.0), (i, 5, 1.0) }).ToArray();
        foreach (int j in new[] { 1, 2, 3 })
        {
            var w = new Witness(W, k42, 0, j);
            Assert.Equal(Witness.Kind.Bell, w.Family);
            Assert.Equal('X', w.Letter);
            Assert.Equal(+1, w.Sign);
            Assert.Equal(1.0, w.Bits, 12);
        }
        Assert.Equal(Witness.Kind.Bell, new Witness(W, k42, 1, 2).Family);
        var alignedChain = new[] { (2, 1, 1.0), (2, 3, 1.0), (1, 0, 2.0), (3, 4, 2.0) };
        var aligned = new Witness(W, alignedChain, 2, 1);
        Assert.Equal(Witness.Kind.Pointer, aligned.Family);  // Law B's alignment: interior chain 0 -> 2
        Assert.Equal(-1, aligned.Sign);                      // the r = 2 watcher rotates by pi
    }
}
