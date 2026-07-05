using RCPsiSquared.Diagnostics.Foundation;

namespace RCPsiSquared.Diagnostics.Tests.Foundation;

/// <summary>From-below gates for <see cref="WatchedLetterRoutingWitness"/> (the label-layer
/// witness, <c>inspect --root label</c>): the 4^N Pauli strings are ONE shared eigenbasis of all
/// three letter dissipators with three different price lists (rate −2γ·n_anti(S, P), the
/// disagreement with the held letter alone); the letter swap relocates which cells pay; only the
/// identity rides free under every watcher. Two-sided throughout: the dense-vs-closed-form
/// residuals must vanish AND the repriced/max-rate controls must not.</summary>
public class WatchedLetterRoutingWitnessTests
{
    private const double Gamma = 0.05;

    private static WatchedLetterRoutingWitness Build() => new(n: 3, gamma: Gamma);

    [Fact]
    public void SharedEigenbasis_DenseMatchesClosedForm()
    {
        var w = Build();
        // Every Pauli string is an exact eigenvector of every letter dissipator (dense side),
        // and the dense eigenvalue equals −2γ·n_anti (closed-form side), all 3·4^N cases.
        Assert.True(w.MaxEigenResidual < 1e-12, $"eigen residual {w.MaxEigenResidual}");
        Assert.True(w.MaxClosedFormDeviation < 1e-12, $"closed-form dev {w.MaxClosedFormDeviation}");
    }

    [Fact]
    public void Routing_ThePriceListFollowsTheHeldLetter()
    {
        var w = Build();
        // Z^⊗N: free under the Z-watcher, maximal under the X-watcher; X^⊗N mirrored.
        Assert.Equal(0.0, w.RateOfAllZUnderZ, 12);
        Assert.Equal(-2.0 * Gamma * 3, w.RateOfAllZUnderX, 12);
        Assert.Equal(-2.0 * Gamma * 3, w.RateOfAllXUnderZ, 12);
        Assert.Equal(0.0, w.RateOfAllXUnderX, 12);
    }

    [Fact]
    public void Routing_RepricedCount_MatchesTheCombinatorics()
    {
        var w = Build();
        // rate_Z(S) ≠ rate_X(S) iff #X-letters ≠ #Z-letters in S. Equal-count strings at N=3 are
        // the constant term of (2 + x + 1/x)^3 = 8 + 12 = 20, so 4^3 − 20 = 44 strings repriced.
        Assert.Equal(44, w.RepricedCountZtoX);
    }

    [Fact]
    public void OnlyTheIdentityRidesFreeUnderEveryWatcher()
    {
        var w = Build();
        Assert.Equal(8, w.FreeCountZ);          // {I,Z}^⊗3
        Assert.Equal(8, w.FreeCountX);          // {I,X}^⊗3
        Assert.Equal(8, w.FreeCountY);          // {I,Y}^⊗3
        Assert.Equal(1, w.UniversalFreeCount);  // I^⊗3 alone
    }

    [Fact]
    public void LetterSwap_IsAnExactTransport()
    {
        var w = Build();
        // Ad(h_zx)·L_Z·Ad(h_zx)† = L_X and Ad(h_yz)·L_Z·Ad(h_yz)† = L_Y, entry-exact.
        Assert.True(w.TransportDevZtoX < 1e-12, $"h_zx dev {w.TransportDevZtoX}");
        Assert.True(w.TransportDevZtoY < 1e-12, $"h_yz dev {w.TransportDevZtoY}");
    }

    [Fact]
    public void TwoSided_TheNonzeroControlHolds()
    {
        var w = Build();
        Assert.Equal(2.0 * Gamma * 3, w.MaxRateMagnitude, 12);  // 2γN > 0: someone always pays
        Assert.True(w.RepricedCountZtoX > 0);                    // the routing is real, not vacuous
    }

    [Fact]
    public void Guards_RejectOversizeNAndNonpositiveGamma()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => new WatchedLetterRoutingWitness(n: 6));
        Assert.Throws<ArgumentOutOfRangeException>(() => new WatchedLetterRoutingWitness(n: 3, gamma: 0.0));
    }
}
