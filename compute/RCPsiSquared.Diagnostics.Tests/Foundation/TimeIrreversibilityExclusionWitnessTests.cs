using System;
using System.Linq;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Diagnostics.Foundation;
using Xunit;

namespace RCPsiSquared.Diagnostics.Tests.Foundation;

/// <summary>The live witness for TIME_IRREVERSIBILITY_EXCLUSION. The gate is two-sided so it
/// cannot pass trivially: the anticommutator {L_H, L_Dc} must VANISH at N=2 (the orthogonality)
/// AND be nonzero for N&gt;2 (the break), and the commutator [L_H, L_Dc] must be nonzero already at
/// N=2 (the caveat: the vanishing is Frobenius-orthogonality, not reversibility). The live matrix
/// norms are cross-checked against the F49 closed form.</summary>
public class TimeIrreversibilityExclusionWitnessTests
{
    [Fact]
    public void Constructor_RejectsBadArgs()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => new TimeIrreversibilityExclusionWitness(nMax: 1));
        Assert.Throws<ArgumentOutOfRangeException>(() =>
            new TimeIrreversibilityExclusionWitness(nMax: TimeIrreversibilityExclusionWitness.MaxN + 1));
        Assert.Throws<ArgumentOutOfRangeException>(() => new TimeIrreversibilityExclusionWitness(gamma: 0.0));
    }

    [Fact]
    public void Anticommutator_VanishesAtN2_ButBreaksForN3()   // the two-sided orthogonality gate
    {
        var w = new TimeIrreversibilityExclusionWitness(nMax: 3, gamma: 0.05, j: 1.0);
        var n2 = w.Rows.Single(r => r.N == 2);
        var n3 = w.Rows.Single(r => r.N == 3);
        Assert.True(n2.AntiNorm < 1e-8, $"‖{{L_H,L_Dc}}‖ should vanish at N=2, got {n2.AntiNorm:E3}");
        Assert.True(n3.AntiNorm > 1e-3, $"‖{{L_H,L_Dc}}‖ should break (be nonzero) at N=3, got {n3.AntiNorm:E3}");
    }

    [Fact]
    public void Commutator_IsNonzeroAtN2_TheReversibilityCaveat()   // vanishing anti is orthogonality, NOT reversibility
    {
        var w = new TimeIrreversibilityExclusionWitness(nMax: 2, gamma: 0.05, j: 1.0);
        var n2 = w.Rows.Single(r => r.N == 2);
        Assert.True(n2.CommNorm > 1e-4,
            $"[L_H,L_Dc] must be nonzero at N=2 (the caveat: separability is the commutator, not the anticommutator), got {n2.CommNorm:E3}");
        // and it must dwarf the (vanishing) anticommutator: they are genuinely different objects.
        Assert.True(n2.CommNorm > 1e6 * n2.AntiNorm,
            $"the commutator {n2.CommNorm:E3} must dwarf the vanishing anticommutator {n2.AntiNorm:E3} at N=2");
    }

    [Fact]
    public void LiveAnticommutator_MatchesF49ClosedForm_EveryN()   // two independent computations meet
    {
        var w = new TimeIrreversibilityExclusionWitness(nMax: 4, gamma: 0.05, j: 1.0);
        foreach (var r in w.Rows)
            Assert.True(Math.Abs(r.AntiNorm - r.ClosedFormAntiNorm) < 1e-6,
                $"N={r.N}: live ‖{{L_H,L_Dc}}‖={r.AntiNorm:E6} vs F49 closed form {r.ClosedFormAntiNorm:E6}");
    }

    [Fact]
    public void AntiNorm_GrowsMonotonicallyAfterN2()
    {
        var w = new TimeIrreversibilityExclusionWitness(nMax: 4, gamma: 0.05, j: 1.0);
        var byN = w.Rows.OrderBy(r => r.N).ToList();
        for (int i = 1; i < byN.Count; i++)
            Assert.True(byN[i].AntiNorm >= byN[i - 1].AntiNorm,
                $"anticommutator norm should not decrease with N: N={byN[i - 1].N}->{byN[i].N}");
    }

    [Fact]
    public void Witness_SurfacesPerNAndCaveatChildren()
    {
        var labels = ((IInspectable)new TimeIrreversibilityExclusionWitness(nMax: 3))
            .Children.Select(c => c.DisplayName).ToList();
        Assert.Contains(labels, l => l.Contains("caveat") && l.Contains("commutator"));
        Assert.Contains(labels, l => l.Contains("N=2"));
        Assert.Contains(labels, l => l.Contains("N=3"));
    }
}
