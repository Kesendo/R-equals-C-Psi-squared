using System.Linq;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Diagnostics.Foundation;
using Xunit;

namespace RCPsiSquared.Diagnostics.Tests.Foundation;

/// <summary>The F89 octic branch locus is a palindrome: its EP/diabolic collisions are mirror-symmetric
/// about Re λ = −4, forced by the F1 palindrome carried on the block as an ANTIUNITARY symmetry. The
/// two-sided gate: the octic roots close under the antilinear mirror λ → −λ̄ − 2σ (reflect Re about −σ,
/// keep Im), but NOT under the bare linear palindrome λ → −λ − 2σ (which would flip Im) — so a sign or
/// construction error fails it. And every branch point is on the line or in a mirror pair, no orphan.</summary>
public class BranchLocusPalindromeWitnessTests
{
    [Fact]
    public void OcticRoots_CloseUnderTheAntiunitaryMirror_NotTheLinearPalindrome()
    {
        var (antiunitary, linear) = BranchLocusPalindromeWitness.MirrorClosureResiduals(q0: 2.0);
        Assert.True(antiunitary < 1e-9, $"octic closes under λ→−λ̄−2σ (residual {antiunitary:E2})");
        Assert.True(linear > 1.0, $"octic does NOT close under the linear λ→−λ−2σ (residual {linear:E2}); two-sided gate");
    }

    [Fact]
    public void EveryBranchPoint_IsOnTheLineOrInAMirrorPair_NoOrphan()
    {
        var (total, onLine, paired, orphans) = BranchLocusPalindromeWitness.BranchLocusStructure();
        Assert.True(total >= 6, $"found enough branch points to test ({total})");
        Assert.Equal(0, orphans);
        Assert.Equal(total, onLine + paired);
    }

    // The line is not the silence. Under a sampled XXZ anisotropy the q_EP crossing splits into two Jordan EP2s,
    // and the mirror keeps both on Re λ = −4 at real q. Location model: an EP2 is fixed to (g/s)² ≈ 1e-15 by its
    // residual gap g ≈ 1e-8 (s ≈ 0.67 at Δ = 0.02, 1.5 at 0.10), the pair mean to u·‖M‖ ≈ 1e-15; 1e-12 leaves 1000.
    [Theory]
    [InlineData(0.02)]
    [InlineData(0.10)]
    public void DeltaControl_StaysOnTheLine_YetDefects(double delta)
    {
        var (first, second) = BranchLocusPalindromeWitness.DeltaControl(delta);
        foreach (var r in new[] { first, second })
        {
            Assert.Equal(XxzCoherenceBlock.DeltaFlipVerdict.Defective, r.Verdict);
            Assert.Equal(2, r.Algebraic);
            Assert.Equal(1, r.Geometric);
            Assert.Equal(1, r.DiscriminantZeroOrder);
            Assert.True(System.Math.Abs(r.QCandidate.Imaginary) <= 1e-12, $"Im q {r.QCandidate.Imaginary:E2}");
            Assert.True(System.Math.Abs(r.LambdaCandidate.Real + 4) <= 1e-12, $"Re lambda + 4 = {r.LambdaCandidate.Real + 4:E2}");
        }
        Assert.True(System.Math.Abs(first.QCandidate.Real - second.QCandidate.Real) > 1e-3, "two distinct EP2s");
    }

    [Fact]
    public void Witness_Renders_ThePalindromeVerdict()
    {
        var children = ((IInspectable)new BranchLocusPalindromeWitness()).Children.ToList();
        Assert.Contains(children, c =>
            c.Summary.Contains("palindrome", System.StringComparison.OrdinalIgnoreCase) ||
            c.Summary.Contains("mirror", System.StringComparison.OrdinalIgnoreCase));
    }
}
