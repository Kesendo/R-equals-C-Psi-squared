using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;

namespace RCPsiSquared.Core.Tests.Symmetry;

/// <summary>Welle 7 (2026-05-26): structural tests for the 9 typed inheritance
/// edges wiring KleinFourCellClaim → KleinEightCellClaim → 8 YParity Claims.
///
/// <para>Before Welle 7: KleinEightCellClaim was structurally orphan (no typed
/// parents, no typed children). After Welle 7: a single typed inheritance path
/// from Klein2 (4-cell decomposition) → Cubic3 (8-cell decomposition) → all 8
/// YParity-axis Claims (F102/F103/F105/F106/F107/F109/F110/F111).</para></summary>
public class Welle7Klein2Cubic3InheritanceTests
{
    private static KleinEightCellClaim BuildKlein8() =>
        new KleinEightCellClaim(new KleinFourCellClaim());

    // ------------------------------------------------------------------
    // B1. KleinFourCellClaim → KleinEightCellClaim
    // ------------------------------------------------------------------

    [Fact]
    public void KleinEightCellClaim_Klein4Parent_IsNotNull()
    {
        var klein8 = BuildKlein8();
        Assert.NotNull(klein8.Klein4Parent);
        Assert.IsType<KleinFourCellClaim>(klein8.Klein4Parent);
    }

    [Fact]
    public void KleinEightCellClaim_Klein4Parent_HasKlein2Axis()
    {
        var klein8 = BuildKlein8();
        Assert.Equal(Z2Axis.Klein2, klein8.Klein4Parent.Z2Axis);
    }

    [Fact]
    public void KleinEightCellClaim_NullKlein4Parent_Throws()
    {
        Assert.Throws<ArgumentNullException>(() =>
            new KleinEightCellClaim(null!));
    }

    [Fact]
    public void KleinEightCellClaim_Z2Axis_IsCubic3()
    {
        Assert.Equal(Z2Axis.Cubic3, BuildKlein8().Z2Axis);
    }

    // ------------------------------------------------------------------
    // B2. YParityIndependenceAtK3 takes KleinEightCellClaim
    // ------------------------------------------------------------------

    [Fact]
    public void F102_KleinEightParent_IsNotNull()
    {
        var klein8 = BuildKlein8();
        var f102 = new YParityIndependenceAtK3(klein8);
        Assert.Same(klein8, f102.KleinEightParent);
    }

    [Fact]
    public void F102_NullKleinEightParent_Throws()
    {
        Assert.Throws<ArgumentNullException>(() =>
            new YParityIndependenceAtK3(null!));
    }

    // ------------------------------------------------------------------
    // B3. F87Z2CubedRefinementN4K3 (F103) takes KleinEightCellClaim
    // ------------------------------------------------------------------

    [Fact]
    public void F103_KleinEightParent_IsNotNull()
    {
        var klein8 = BuildKlein8();
        var f103 = new F87Z2CubedRefinementN4K3(klein8);
        Assert.Same(klein8, f103.KleinEightParent);
    }

    [Fact]
    public void F103_NullKleinEightParent_Throws()
    {
        Assert.Throws<ArgumentNullException>(() =>
            new F87Z2CubedRefinementN4K3(null!));
    }

    // ------------------------------------------------------------------
    // B4. F87Z2CubedRefinementN5K3 (F105) takes KleinEightCellClaim
    // ------------------------------------------------------------------

    [Fact]
    public void F105_KleinEightParent_IsNotNull()
    {
        var klein8 = BuildKlein8();
        var f105 = new F87Z2CubedRefinementN5K3(klein8);
        Assert.Same(klein8, f105.KleinEightParent);
    }

    [Fact]
    public void F105_NullKleinEightParent_Throws()
    {
        Assert.Throws<ArgumentNullException>(() =>
            new F87Z2CubedRefinementN5K3(null!));
    }

    // ------------------------------------------------------------------
    // B5. F87Z2CubedRefinementN4K4 (F106) takes KleinEightCellClaim
    // ------------------------------------------------------------------

    [Fact]
    public void F106_KleinEightParent_IsNotNull()
    {
        var klein8 = BuildKlein8();
        var f106 = new F87Z2CubedRefinementN4K4(klein8);
        Assert.Same(klein8, f106.KleinEightParent);
    }

    [Fact]
    public void F106_NullKleinEightParent_Throws()
    {
        Assert.Throws<ArgumentNullException>(() =>
            new F87Z2CubedRefinementN4K4(null!));
    }

    // ------------------------------------------------------------------
    // B6. TrulyYParityZeroPurity (F107) takes KleinEightCellClaim
    // ------------------------------------------------------------------

    [Fact]
    public void F107_KleinEightParent_IsNotNull()
    {
        var klein8 = BuildKlein8();
        var f107 = new TrulyYParityZeroPurity(klein8);
        Assert.Same(klein8, f107.KleinEightParent);
    }

    [Fact]
    public void F107_NullKleinEightParent_Throws()
    {
        Assert.Throws<ArgumentNullException>(() =>
            new TrulyYParityZeroPurity(null!));
    }

    // ------------------------------------------------------------------
    // B7. MotherSoftYParityOnePurity (F109) takes KleinEightCellClaim
    // ------------------------------------------------------------------

    [Fact]
    public void F109_KleinEightParent_IsNotNull()
    {
        var klein8 = BuildKlein8();
        var f109 = new MotherSoftYParityOnePurity(klein8);
        Assert.Same(klein8, f109.KleinEightParent);
    }

    [Fact]
    public void F109_NullKleinEightParent_Throws()
    {
        Assert.Throws<ArgumentNullException>(() =>
            new MotherSoftYParityOnePurity(null!));
    }

    // ------------------------------------------------------------------
    // B8. HardCellYInversionPattern (F110) takes KleinEightCellClaim
    // ------------------------------------------------------------------

    [Fact]
    public void F110_KleinEightParent_IsNotNull()
    {
        var klein8 = BuildKlein8();
        var f110 = new HardCellYInversionPattern(klein8);
        Assert.Same(klein8, f110.KleinEightParent);
    }

    [Fact]
    public void F110_NullKleinEightParent_Throws()
    {
        Assert.Throws<ArgumentNullException>(() =>
            new HardCellYInversionPattern(null!));
    }

    // ------------------------------------------------------------------
    // B9. HardCellPureDTemplate (F111) takes KleinEightCellClaim
    // ------------------------------------------------------------------

    [Fact]
    public void F111_KleinEightParent_IsNotNull()
    {
        var klein8 = BuildKlein8();
        var f111 = new HardCellPureDTemplate(klein8);
        Assert.Same(klein8, f111.KleinEightParent);
    }

    [Fact]
    public void F111_NullKleinEightParent_Throws()
    {
        Assert.Throws<ArgumentNullException>(() =>
            new HardCellPureDTemplate(null!));
    }

    // ------------------------------------------------------------------
    // Inheritance chain: all 8 YParity Claims see the same KleinEightCellClaim
    // ------------------------------------------------------------------

    [Fact]
    public void Klein2_Cubic3_YParity_Chain_TraversableEndToEnd()
    {
        var klein4 = new KleinFourCellClaim();
        var klein8 = new KleinEightCellClaim(klein4);

        var f102 = new YParityIndependenceAtK3(klein8);
        var f103 = new F87Z2CubedRefinementN4K3(klein8);
        var f105 = new F87Z2CubedRefinementN5K3(klein8);
        var f106 = new F87Z2CubedRefinementN4K4(klein8);
        var f107 = new TrulyYParityZeroPurity(klein8);
        var f109 = new MotherSoftYParityOnePurity(klein8);
        var f110 = new HardCellYInversionPattern(klein8);
        var f111 = new HardCellPureDTemplate(klein8);

        // All 8 YParity Claims share the same KleinEightCellClaim parent instance
        var yParityClaims = new IZ2AxisClaim[] { f102, f103, f105, f106, f107, f109, f110, f111 };
        foreach (var c in yParityClaims)
            Assert.Equal(Z2Axis.YParity, c.Z2Axis);

        Assert.Same(klein8, f102.KleinEightParent);
        Assert.Same(klein8, f103.KleinEightParent);
        Assert.Same(klein8, f105.KleinEightParent);
        Assert.Same(klein8, f106.KleinEightParent);
        Assert.Same(klein8, f107.KleinEightParent);
        Assert.Same(klein8, f109.KleinEightParent);
        Assert.Same(klein8, f110.KleinEightParent);
        Assert.Same(klein8, f111.KleinEightParent);

        // KleinEightCellClaim sees KleinFourCellClaim as its single Klein2 parent
        Assert.Same(klein4, klein8.Klein4Parent);
        Assert.Equal(Z2Axis.Klein2, klein4.Z2Axis);
        Assert.Equal(Z2Axis.Cubic3, klein8.Z2Axis);
    }
}
