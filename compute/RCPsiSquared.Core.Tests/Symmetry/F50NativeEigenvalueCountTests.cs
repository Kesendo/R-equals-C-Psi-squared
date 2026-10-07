using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.ChainSystems;
using RCPsiSquared.Core.Symmetry;
using ComplexMatrix = MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex>;

namespace RCPsiSquared.Core.Tests.Symmetry;

/// <summary>Native C# verification of F50's claim that the Liouvillian of a
/// Heisenberg chain + uniform Z-dephasing has exactly <c>2N</c> purely-real
/// eigenvalues at Re = −2γ (the weight-1 pure-sector line, not necessarily the
/// full spectrum's first nonzero decay rate; the SWAP-
/// invariant conserved operators T_c^{(a)} for a ∈ {X, Y}, c = 0..N−1).
///
/// <para>Builds L via <c>ChainSystem.BuildLiouvillian()</c>, eigendecomposes,
/// counts eigenvalues whose real part is within tolerance of −2γ AND whose
/// imaginary part is within tolerance of 0, and verifies the count matches
/// <c>F50.TotalDegeneracy(N) = 2N</c>. Tested at the canonical J = γ = 1
/// parameters for N = 2..5, the regime where the docstring claims "verified
/// N=2..7" via the same eigendecomposition.</para>
///
/// <para>The count 2N holds at every γ/J outside a finite exceptional set,
/// where a diagonal joint-popcount block (p, p) puts real modes that are not
/// SWAP-invariants exactly on −2γ. In this C# spin convention J_s = 4J (Pauli),
/// so J_s = 2γ is γ = 2J, the N = 2 exceptional point, where the count is 6. Two tests
/// (ExceptionalCoupling_*) pin that point and the N = 3 chain's γ/J = √3. J = γ = 1 here is γ/J = 4 (Pauli), outside every
/// exceptional set computed for N ≤ 4 and the N = 5 chain (PROOF_WEIGHT1_DEGENERACY § The count at
/// exceptional couplings).</para>
/// </summary>
public class F50NativeEigenvalueCountTests
{
    private const double Gamma = 1.0;
    private const double J = 1.0;
    private const double Tol = 1e-6;

    private static int CountPureRealEigenvaluesAt(int N, double targetRealPart, double j = J)
    {
        var L = new ChainSystem(N: N, J: j, GammaZero: Gamma,
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
        // measured; numerical here confirms both at canonical J = γ = 1).
        int actual = CountPureRealEigenvaluesAt(N, -2.0 * Gamma);
        Assert.Equal(expected2N, actual);
    }

    [Fact]
    public void ExceptionalCoupling_N2_AtJsEqualTwoGamma_CountIsSix()
    {
        // N = 2 at J_s = 2γ (γ = 2J Pauli): the (1,1) pair λ² + 4γλ + 16J² has its double root on
        // −2γ (a defective pair; in floating point it splits by ~1e-8, inside the window), so the
        // line holds 2N + 2 = 6 real eigenvalues.
        Assert.Equal(6, CountPureRealEigenvaluesAt(2, -2.0 * Gamma, j: 2.0 * Gamma));
    }

    [Fact]
    public void ExceptionalCoupling_N3Chain_AtGammaOverJEqualSqrt3_CountIsEight()
    {
        // N = 3 chain at γ/J = √3 (Pauli), J_s = 4γ/√3: one real mode from the (1,1) block and one
        // from (2,2) join the line, semisimple, so the count is 2N + 2 = 8.
        Assert.Equal(8, CountPureRealEigenvaluesAt(3, -2.0 * Gamma, j: 4.0 * Gamma / Math.Sqrt(3.0)));
    }

    [Fact]
    public void NativeSpectrum_NoExtraPureRealAtThreeGamma()
    {
        // Adjacent grid positions (Re = -γ, Re = -3γ) should NOT have 2N
        // pure-real eigenvalues; F50 is specific to Re = -2γ. At N = 3
        // these positions are not on the F50 grid; expect zero or far-from-2N
        // count.
        const int N = 3;
        int countAtMinusGamma = CountPureRealEigenvaluesAt(N, -Gamma);
        int countAtMinusThreeGamma = CountPureRealEigenvaluesAt(N, -3.0 * Gamma);
        Assert.NotEqual(2 * N, countAtMinusGamma);
        Assert.NotEqual(2 * N, countAtMinusThreeGamma);
    }

    // ────────────────────────────────────────────────────────────────────
    // K_3 N=3 anomaly (2026-05-17 finding)
    // ────────────────────────────────────────────────────────────────────
    //
    // F50's "topology universality" claim (any connected graph → count = 2N)
    // is empirically violated at N=3 K_3 (= ring = triangle = complete on 3
    // vertices): the count is 8, not 6. The 2 extras are weight-1 operators
    // that commute with H_K_3 globally (the three bond commutators cancel
    // pairwise) but NOT with H_chain alone. They correspond to the 2-dim
    // standard irreducible representation of S_3 = Aut(K_3) acting on the
    // weight-1 c=1 sector. Adding any external bond to K_3 (paw, bowtie, book)
    // breaks the S_3 symmetry and restores the count to 2N.
    //
    // The original proof's Step 5 derivation conflates matrix-commutator
    // [SWAP_b, v] with conjugation action SWAP_b·v·SWAP_b^dagger − v; these
    // are not equivalent, and the proof's triangle-inequality argument
    // applies only to the conjugation form. See PROOF_WEIGHT1_DEGENERACY §
    // Appendix (2026-05-17) for the full structural reading.

    [Fact]
    public void NativeSpectrum_K3N3_Anomaly_Has8PureRealAtMinusTwoGamma()
    {
        // K_3 N=3 (= ring at N=3 = triangle = complete on 3 vertices) has
        // 8 pure-real Liouvillian eigenvalues at Re = -2γ, NOT the 2N = 6
        // that F50 predicts. The 2 extras are S_3-standard-rep operators
        // that exist because K_3 has full S_3 symmetry.
        var L = new ChainSystem(N: 3, J: J, GammaZero: Gamma,
                                HType: HamiltonianType.Heisenberg,
                                Topology: TopologyKind.Ring).BuildLiouvillian();
        int count = L.Evd().EigenValues.Count(lam =>
            Math.Abs(lam.Real - (-2.0 * Gamma)) < Tol &&
            Math.Abs(lam.Imaginary) < Tol);
        Assert.Equal(8, count);  // 2N + 2 = 6 + 2; F50 predicts only 6 here
    }

    [Fact]
    public void NativeSpectrum_K3N3_AllExtrasAreInKerOfHCommutator()
    {
        // The 2 K_3-extras commute with H_K3 globally (||[H_K3, A]|| ≤ 1e-12),
        // which means F50's "lower bound dim(ker) ≥ 2N" is satisfied (it's
        // actually dim ≥ 2N + 2 here) but the "upper bound = 2N" is empirically
        // violated. This test confirms the extras are real H-commutants, not
        // numerical noise.
        var H = BuildHeisenbergRingN3();
        var L = BuildLiouvillianN3FromH(H);
        var evd = L.Evd();
        var eigs = evd.EigenValues.ToArray();
        var vecs = evd.EigenVectors;

        // Build commutator superop [H, ·] in vec form.
        var I = Matrix<Complex>.Build.DenseIdentity(8);
        var C = I.KroneckerProduct(H) - H.Transpose().KroneckerProduct(I);

        int countInKer = 0;
        for (int k = 0; k < eigs.Length; k++)
        {
            var lam = eigs[k];
            if (Math.Abs(lam.Real - (-2.0 * Gamma)) > Tol ||
                Math.Abs(lam.Imaginary) > Tol) continue;
            var v = vecs.Column(k);
            var Hv = C * v;
            double norm = Hv.L2Norm();
            if (norm < 1e-10) countInKer++;
        }
        // All 8 should be in ker[H_K3, ·], not just the 6 F50-predicted ones.
        Assert.Equal(8, countInKer);
    }

    [Fact]
    public void NativeSpectrum_LowerBoundAt2N_HoldsForAllTestedTopologies()
    {
        // F50's lower bound `dim(ker) ≥ 2N` holds for every connected graph
        // tested. The K_3 N=3 case exceeds the lower bound (8 > 6) but never
        // falls below it. This test verifies the lower bound on 3 topologies
        // at N = 3 (chain, ring/K_3, star at N=3 = chain by relabeling).
        var topologies = new[]
        {
            (Label: "Chain N=3", Topo: TopologyKind.Chain, ExpectedAtLeast: 6),
            (Label: "Ring/K_3 N=3", Topo: TopologyKind.Ring, ExpectedAtLeast: 6),
            (Label: "Star N=3", Topo: TopologyKind.Star, ExpectedAtLeast: 6),
        };
        foreach (var t in topologies)
        {
            var L = new ChainSystem(N: 3, J: J, GammaZero: Gamma,
                                    HType: HamiltonianType.Heisenberg,
                                    Topology: t.Topo).BuildLiouvillian();
            int count = L.Evd().EigenValues.Count(lam =>
                Math.Abs(lam.Real - (-2.0 * Gamma)) < Tol &&
                Math.Abs(lam.Imaginary) < Tol);
            Assert.True(count >= t.ExpectedAtLeast,
                $"{t.Label}: F50 lower bound 2N = {t.ExpectedAtLeast} violated; got {count}");
        }
    }

    private static ComplexMatrix BuildHeisenbergRingN3() =>
        new ChainSystem(N: 3, J: J, GammaZero: 0.0,
                        HType: HamiltonianType.Heisenberg,
                        Topology: TopologyKind.Ring).BuildHamiltonian();

    private static ComplexMatrix BuildLiouvillianN3FromH(ComplexMatrix H) =>
        new ChainSystem(N: 3, J: J, GammaZero: Gamma,
                        HType: HamiltonianType.Heisenberg,
                        Topology: TopologyKind.Ring).BuildLiouvillian();
}
