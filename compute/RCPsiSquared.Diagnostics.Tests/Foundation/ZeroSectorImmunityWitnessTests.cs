using System.Linq;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Diagnostics.Foundation;
using Xunit;

namespace RCPsiSquared.Diagnostics.Tests.Foundation;

/// <summary>The live Zero-Sector Immunity witness: for a RANDOM parity-violating 2-body
/// Hamiltonian under uniform Z-dephasing, the centered palindromic residual M vanishes on the
/// (w=0,w=0) block ({I,Z}^⊗N) and the (w=N,w=N) block ({X,Y}^⊗N), while the full M is NON-zero
/// (the non-triviality gate: the H genuinely breaks the global F1 palindrome).</summary>
public class ZeroSectorImmunityWitnessTests
{
    private static System.Collections.Generic.List<IInspectable> Children(ZeroSectorImmunityWitness w) =>
        ((IInspectable)w).Children.ToList();

    [Theory]
    [InlineData(2)]
    [InlineData(3)]
    [InlineData(4)]
    public void ZeroSector_Vanishes_ForParityViolatingH(int n)
    {
        var w = new ZeroSectorImmunityWitness(n);
        Assert.True(w.ZeroSectorNorm < 1e-9,
            $"‖M|(w=0)‖_F = {w.ZeroSectorNorm:e3} must vanish at N={n} (Zero-Sector Immunity)");
    }

    [Theory]
    [InlineData(2)]
    [InlineData(3)]
    [InlineData(4)]
    public void MirrorSector_Vanishes_ByPiSymmetry(int n)
    {
        var w = new ZeroSectorImmunityWitness(n);
        Assert.True(w.MirrorSectorNorm < 1e-9,
            $"‖M|(w=N)‖_F = {w.MirrorSectorNorm:e3} must vanish at N={n} (Π-mirror of w=0)");
    }

    [Theory]
    [InlineData(2)]
    [InlineData(3)]
    [InlineData(4)]
    public void FullResidual_IsNonZero_TheNonTrivialGate(int n)
    {
        var w = new ZeroSectorImmunityWitness(n);
        Assert.True(w.IsParityViolating,
            $"the witness H must be genuinely parity-violating at N={n}, else the test is trivial");
        Assert.True(w.FullResidualNorm > 1e-3,
            $"‖M‖_F = {w.FullResidualNorm:e3} must be non-zero at N={n} (the H breaks the global palindrome)");
    }

    [Fact]
    public void SectorDimension_Is2PowerN()
    {
        Assert.Equal(1 << 3, new ZeroSectorImmunityWitness(n: 3).SectorDim);   // 8 = |{I,Z}^⊗3|
    }

    [Fact]
    public void Constructor_RejectsBadArgs()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => new ZeroSectorImmunityWitness(n: 0));
        Assert.Throws<ArgumentOutOfRangeException>(() => new ZeroSectorImmunityWitness(gamma: 0.0));
    }

    [Fact]
    public void Witness_SurfacesGateAndSectorChildren()
    {
        var labels = Children(new ZeroSectorImmunityWitness(n: 3)).Select(c => c.DisplayName).ToList();
        Assert.Contains(labels, l => l.Contains("immunity") || l.Contains("w=0"));
        Assert.Contains(labels, l => l.Contains("non-trivial") || l.Contains("full M"));
    }

    [Fact]
    public void Witness_RendersToJson()
    {
        var json = InspectionJsonExporter.ToJson(new ZeroSectorImmunityWitness(n: 3));
        Assert.Contains("w=0", json);
    }
}
