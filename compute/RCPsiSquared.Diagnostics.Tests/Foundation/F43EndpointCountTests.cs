using RCPsiSquared.Core.F1;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Diagnostics.Foundation;

namespace RCPsiSquared.Diagnostics.Tests.Foundation;

/// <summary>The F43 endpoint multiplicity, counted instead of returned.
///
/// <para><see cref="F43BandSffPairingPi2Inheritance.EndpointMultiplicity"/> states that a connected
/// uniform chain at γ &gt; 0 carries N+1 modes at each palindrome end, λ = 0 and λ = −2Nγ. The
/// count is F158's two-end count, and <see cref="PalindromeTwoEndCountWitness"/> reads it from the
/// actual Liouvillian by exact ranks over GF(p) (Heisenberg bonds, uniform Z-dephasing at
/// γ = 1/20, so 2σ = 2Nγ): dim ker L and dim ker(L + 2Nγ). Nothing compared here is a float.
/// A formula that drifted from the spectrum, or a spectrum that lost a mode at either end, would
/// fail the equality; the witness's own suite carries the row where a single Z field empties the
/// far end, so the instrument is known to report counts other than N+1.</para></summary>
public class F43EndpointCountTests
{
    private static F43BandSffPairingPi2Inheritance BuildClaim()
    {
        var ladder = new Pi2DyadicLadderClaim();
        var f1 = new F1Pi2Inheritance(new F1PalindromeIdentity(), ladder, new Pi2I4MemoryLoopClaim());
        return new F43BandSffPairingPi2Inheritance(ladder, f1);
    }

    [Theory]
    [InlineData(2, "chain")]
    [InlineData(3, "chain")]
    [InlineData(4, "chain")]
    [InlineData(3, "ring")]
    [InlineData(4, "ring")]
    public void EndpointMultiplicity_EqualsTheExactCountAtBothEnds(int n, string topology)
    {
        var claim = BuildClaim();
        var reading = new PalindromeTwoEndCountWitness(n, topology: topology).Read();

        int stated = claim.EndpointMultiplicity(n, gamma: 0.05);
        Assert.Equal(stated, reading.FarCount);
        Assert.Equal(stated, reading.NearCount);
        // The far end read a second way, from the operator conditions (anticommutant with the
        // jump sign flipped), which never forms L.
        Assert.Equal(stated, reading.FarFromOperatorConditions);
        // The unnormalized endpoint SFF is the square of that count, (N+1)².
        Assert.Equal((double)reading.FarCount * reading.FarCount,
            claim.EndpointUnnormalizedFrequencySff(n, gamma: 0.05, time: 1.0));
    }
}
