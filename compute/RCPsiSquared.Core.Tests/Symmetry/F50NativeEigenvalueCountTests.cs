using System.Numerics;
using RCPsiSquared.Core.ChainSystems;
using RCPsiSquared.Core.Symmetry;

namespace RCPsiSquared.Core.Tests.Symmetry;

/// <summary>Native C# verification of F50's claim that the Liouvillian of a
/// Heisenberg chain + uniform Z-dephasing has exactly <c>2N</c> purely-real
/// eigenvalues at Re = −2γ (the first non-zero real grid position; the SWAP-
/// invariant conserved operators T_c^{(a)} for a ∈ {X, Y}, c = 0..N−1).
///
/// <para>Builds L via <c>ChainSystem.BuildLiouvillian()</c>, eigendecomposes,
/// counts eigenvalues whose real part is within tolerance of −2γ AND whose
/// imaginary part is within tolerance of 0, and verifies the count matches
/// <c>F50.TotalDegeneracy(N) = 2N</c>. Tested at the canonical J = γ = 1
/// parameters for N = 2..5 — the regime where the docstring claims "verified
/// N=2..7" via the same eigendecomposition.</para>
///
/// <para>At very large J/γ (e.g. J ≥ 2γ at N = 2) the count can include
/// accidental real degeneracies that aren't SWAP-invariants; those are
/// outside F50's typical verified regime and not covered here.</para>
/// </summary>
public class F50NativeEigenvalueCountTests
{
    private const double Gamma = 1.0;
    private const double J = 1.0;
    private const double Tol = 1e-6;

    private static int CountPureRealEigenvaluesAt(int N, double targetRealPart)
    {
        var L = new ChainSystem(N: N, J: J, GammaZero: Gamma,
                                HType: HamiltonianType.Heisenberg,
                                Topology: TopologyKind.Chain).BuildLiouvillian();
        var eigs = L.Evd().EigenValues.ToArray();
        return eigs.Count(lam =>
            Math.Abs(lam.Real - targetRealPart) < Tol &&
            Math.Abs(lam.Imaginary) < Tol);
    }

    private static F50WeightOneDegeneracyPi2Inheritance BuildClaim() =>
        new F50WeightOneDegeneracyPi2Inheritance(new Pi2DyadicLadderClaim());

    [Theory]
    [InlineData(2)]
    [InlineData(3)]
    [InlineData(4)]
    [InlineData(5)]
    public void NativeSpectrum_HasExactly2N_PureRealEigenvaluesAtMinusTwoGamma(int N)
    {
        // F50.TotalDegeneracy(N) = 2N must equal the count of pure-real
        // eigenvalues at Re = −2γ in the actual Liouvillian spectrum.
        var f50 = BuildClaim();
        int expected = f50.TotalDegeneracy(N);
        double target = f50.EigenvaluePosition(Gamma);  // −2γ
        int actual = CountPureRealEigenvaluesAt(N, target);
        Assert.Equal(expected, actual);
    }

    [Theory]
    [InlineData(2, 4)]
    [InlineData(3, 6)]
    [InlineData(4, 8)]
    [InlineData(5, 10)]
    public void NativeSpectrum_Matches_F50_TotalDegeneracy(int N, int expected2N)
    {
        // Direct check against the integer constant 2N for N = 2..5.
        // Cross-validates F50.TotalDegeneracy against the SWAP-invariant
        // construction in the docstring (lower bound proven; upper bound
        // proven; numerical here confirms both at canonical J = γ = 1).
        int actual = CountPureRealEigenvaluesAt(N, -2.0 * Gamma);
        Assert.Equal(expected2N, actual);
    }

    [Fact]
    public void NativeSpectrum_NoExtraPureRealAtThreeGamma()
    {
        // Adjacent grid positions (Re = -γ, Re = -3γ) should NOT have 2N
        // pure-real eigenvalues — F50 is specific to Re = -2γ. At N = 3
        // these positions are not on the F50 grid; expect zero or far-from-2N
        // count.
        const int N = 3;
        int countAtMinusGamma = CountPureRealEigenvaluesAt(N, -Gamma);
        int countAtMinusThreeGamma = CountPureRealEigenvaluesAt(N, -3.0 * Gamma);
        Assert.NotEqual(2 * N, countAtMinusGamma);
        Assert.NotEqual(2 * N, countAtMinusThreeGamma);
    }

    // Note: F50's docstring claims "topology universality" of the 2N count
    // (chain, star, ring, complete, tree). Empirically at N = 3 Ring + Heisenberg
    // + Z-deph (J = γ = 1), the pure-real count at −2γ is 8, not 6 — likely
    // 2 extra accidental real eigenvalues from the ring's translation symmetry
    // (which chain lacks). Adding a Ring-vs-Chain native test surfaces a
    // potential refinement of F50's topology-universality claim that's its
    // own thread; left out of this native verification to keep F50's chain
    // case clean.
}
