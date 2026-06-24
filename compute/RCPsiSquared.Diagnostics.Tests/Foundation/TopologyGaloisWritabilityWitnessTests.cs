using System.Linq;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Diagnostics.Foundation;
using Xunit;

namespace RCPsiSquared.Diagnostics.Tests.Foundation;

/// <summary>The topology-controls-Galois-writability witness recomputes, from the (SE,DE) block + its
/// symmetry, the N-independent factor-degree caps that classify the relaxation: complete K_N → 4
/// (writable), star → 9 (a fixed S_9 scramble); and confirms the chain has no such large-group symmetry.</summary>
public class TopologyGaloisWritabilityWitnessTests
{
    [Theory]
    [InlineData(5)]
    [InlineData(6)]
    public void CompleteGraph_StandardMultiplicity_Is4_NIndependent(int n)
    {
        Assert.Equal(4, TopologyGaloisWritabilityWitness.CompleteCap(n));
    }

    [Theory]
    [InlineData(5)]
    [InlineData(6)]
    [InlineData(7)]
    public void Star_StandardMultiplicity_Is9_NIndependent(int n)
    {
        Assert.Equal(9, TopologyGaloisWritabilityWitness.StarCap(n));
    }

    [Theory]
    [InlineData(5, "complete")]
    [InlineData(6, "complete")]
    [InlineData(5, "star")]
    [InlineData(6, "star")]
    public void Symmetry_CommutesWithLiouvillian_ForLargeGroupTopologies(int n, string topo)
    {
        Assert.True(TopologyGaloisWritabilityWitness.SymmetryHolds(n, topo));
    }

    [Fact]
    public void Chain_HasNoAdjacentTranspositionSymmetry()
    {
        // the chain's only nontrivial automorphism is the full reversal, not an adjacent transposition,
        // so it has no large-group cap (its multiplicities grow => S_8/S_18/S_32/S_53).
        Assert.False(TopologyGaloisWritabilityWitness.SymmetryHolds(5, "chain"));
    }

    [Fact]
    public void Witness_Renders_TheWritableCompleteAndBoundedStarVerdicts()
    {
        var children = ((IInspectable)new TopologyGaloisWritabilityWitness()).Children.ToList();
        Assert.Contains(children, c => c.DisplayName.Contains("cap = 4") && c.DisplayName.Contains("writable"));
        Assert.Contains(children, c => c.DisplayName.Contains("S_9 scramble"));
        Assert.Contains(children, c => c.Summary.Contains("for ALL N"));
    }
}
