using System.Linq;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Diagnostics.Foundation;
using Xunit;

namespace RCPsiSquared.Diagnostics.Tests.Foundation;

/// <summary>The live filling-threshold witness (inspect --root fillcsr), the persistent evidence for the F89
/// Door-C decisive follow-up. Gates that the live recomputation reproduces the qualitative result: the dilute
/// (1,2) block stays Poisson while the dense (3,4) block develops GinUE angular repulsion that grows with N.</summary>
public class FillingThresholdWitnessTests
{
    /// <summary>The witness materializes its children without error and the verdict node reads CONFIRMED — the
    /// live read reproduces the dilute-flat / dense-rising-with-N filling-threshold pattern.</summary>
    [Fact]
    public void Witness_VerdictReadsConfirmed()
    {
        var children = new FillingThresholdWitness().Children.ToList();
        Assert.NotEmpty(children);
        var verdict = children.Single(c => c.DisplayName.Contains("the verdict"));
        Assert.Contains("CONFIRMED", verdict.Summary);
    }
}
