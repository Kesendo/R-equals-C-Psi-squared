using RCPsiSquared.Core.ChainSystems;
using RCPsiSquared.Core.Receivers;
using RCPsiSquared.Core.States;

namespace RCPsiSquared.Core.Tests.Receivers;

public class ReceiverTests
{
    [Fact]
    public void Constructor_RejectsUnnormalizedState()
    {
        int N = 3;
        var psi = MathNet.Numerics.LinearAlgebra.Vector<System.Numerics.Complex>.Build.Dense(1 << N);
        psi[0] = 2.0;
        Assert.Throws<ArgumentException>(() => new Receiver(psi));
    }

    [Fact]
    public void FromUnnormalized_NormalizesAndAccepts()
    {
        var psi = MathNet.Numerics.LinearAlgebra.Vector<System.Numerics.Complex>.Build.Dense(8);
        psi[0] = 3.0;
        psi[7] = 4.0;
        var r = Receiver.FromUnnormalized(psi);
        Assert.Equal(3, r.N);
        Assert.Equal(1.0, r.Psi.ConjugateDotProduct(r.Psi).Real, 10);
    }

    [Fact]
    public void BondingMode_IsF71Eigenstate_AtBalancedN()
    {
        // F65 mode k=1 at N=5: ψ_1(j) = sin(πj/6) for j=1..5 → symmetric pattern → F71 +1.
        var chain = new ChainSystem(N: 5, J: 1.0, GammaZero: 0.05);
        var r = Receiver.BondingMode(chain, k: 1);
        Assert.Equal(+1, r.F71Eigenvalue);
    }

    [Fact]
    public void Signature_AtN5_BalancedSplit_ReportsCapacityOptimal()
    {
        // F71-eigenstate on balanced odd-N split → capacity-optimal (per J_BLIND_RECEIVER_CLASSES.md).
        var chain = new ChainSystem(N: 5, J: 1.0, GammaZero: 0.05);
        var r = Receiver.BondingMode(chain, k: 2);
        var sig = r.Signature();
        Assert.Equal(5, sig.N);
        Assert.True(sig.BondBlockBalanced);
        Assert.Contains("capacity-optimal", sig.Prediction);
    }

    [Fact]
    public void Signature_AtN6_UnbalancedSplit_ReportsCapacitySuboptimal()
    {
        var chain = new ChainSystem(N: 6, J: 1.0, GammaZero: 0.05);
        var r = Receiver.BondingMode(chain, k: 2);
        var sig = r.Signature();
        Assert.False(sig.BondBlockBalanced);
        Assert.Contains("capacity-suboptimal", sig.Prediction);
    }
}
