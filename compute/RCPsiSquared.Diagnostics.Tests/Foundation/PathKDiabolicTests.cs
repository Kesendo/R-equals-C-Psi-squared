using System;
using System.Linq;
using System.Numerics;
using RCPsiSquared.Diagnostics.Foundation;
using Xunit;

namespace RCPsiSquared.Diagnostics.Tests.Foundation;

/// <summary>The path-k diabolic investigation (the N=4→N=5 forward edge of the zeros_connecting_structure
/// arc). The N=4 diabolic is the +2.349 σ_T twin pair merging onto the fold at λ_EP=−4+1.318i; these tests
/// build and validate the generalized path-k diabolic tooling against that KNOWN path-3 case before it is
/// trusted at path-4. See docs/superpowers/plans/2026-06-27-f89-path-k-diabolic.md.</summary>
public class PathKDiabolicTests
{
    // Task 1 — ResidualRootsAt(k,q) is the path-k analogue of GaloisMonodromyWitness.OcticRootsAt, valid at
    // q≈2 (R-7 scope). The trusted cross-check: at q=2 the path-3 residual MUST be the path-3 octic exactly.
    [Fact]
    public void ResidualRootsAt_Path3_MatchesOcticWitness()
    {
        var pk = PathKMonodromyScout.ResidualRootsAt(3, new Complex(2, 0))
            .OrderBy(z => z.Real).ThenBy(z => z.Imaginary).ToArray();
        var oc = GaloisMonodromyWitness.OcticRootsAt(new Complex(2, 0))
            .OrderBy(z => z.Real).ThenBy(z => z.Imaginary).ToArray();
        Assert.Equal(8, pk.Length);
        for (int i = 0; i < 8; i++)
            Assert.True((pk[i] - oc[i]).Magnitude < 1e-9, $"strand {i}: {pk[i]} vs {oc[i]}");
    }
}
