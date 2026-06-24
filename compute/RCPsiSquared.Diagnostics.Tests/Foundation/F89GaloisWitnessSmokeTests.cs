using System.Collections.Generic;
using System.Linq;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Diagnostics.Foundation;
using Xunit;

namespace RCPsiSquared.Diagnostics.Tests.Foundation;

/// <summary>Smoke tests for the live H_B-mixed Galois witnesses. They exercise the full inspect-time
/// orchestration (block → Berkowitz over Z[i] → reconstruct the AT factor from the rate-confined
/// invariant subspace → isolate F_d via the validation triple → Frobenius cycle types → JordanVerdict)
/// and assert the rendered verdict node. The Core-level tests assert the MATH (live F_d == oracle,
/// JordanVerdict gates); these assert the WITNESS WIRING the Core tests do not touch. Kept to path-3 +
/// path-4 to stay cheap (~2 s); path-5/6 live isolation is asserted by F89FullDReconstructionTests.</summary>
public class F89GaloisWitnessSmokeTests
{
    private static List<IInspectable> Children(IInspectable w) => w.Children.ToList();

    [Fact]
    public void Path3_OcticWitness_RecomputesLive_AndRendersS8NonSolvable()
    {
        var children = Children(new F89OcticGaloisWitness());
        Assert.Contains(children, c => c.Summary.Contains("S_8") && c.Summary.Contains("non-solvable"));
    }

    [Fact]
    public void Path4_LiveWitness_RecomputesLive_AndRendersS18NonSolvable()
    {
        var children = Children(new F89PathKLiveGaloisWitness(4));
        Assert.Contains(children, c => c.Summary.Contains("S_18") && c.Summary.Contains("non-solvable"));
    }
}
