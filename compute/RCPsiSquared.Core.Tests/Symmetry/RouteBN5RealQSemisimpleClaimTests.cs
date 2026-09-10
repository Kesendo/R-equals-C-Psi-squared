using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;

namespace RCPsiSquared.Core.Tests.Symmetry;

public sealed class RouteBN5RealQSemisimpleClaimTests
{
    [Fact]
    public void ParameterlessClaimHasNoEmbeddedClaimParentAndCarriesTheCurrentEvidenceSplit()
    {
        var claim = new RouteBN5RealQSemisimpleClaim();

        Assert.DoesNotContain(claim.Children, child => child is Claim);
        Assert.Contains("only real-q points in the N=5/N=6 A2 inventories", claim.Summary);
        Assert.Contains("q is the settable parameter", claim.Summary);
        Assert.Contains("lambda is the resulting spectral eigenvalue", claim.Summary);
        Assert.DoesNotContain("observed", claim.Summary, StringComparison.OrdinalIgnoreCase);
        Assert.Contains("24 imaginary-q directly Hermitian-certified", claim.Summary);
        Assert.Contains("F164 certifies all 26 R-odd loci", claim.Summary);
        Assert.Contains("38 distinct exact-certified loci in total", claim.Summary);
        Assert.Contains("20 nonreal-q R-even loci numerical-only", claim.Summary);
    }
}
