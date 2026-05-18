using RCPsiSquared.Core.F1;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.Tests.F1;

public class F49NonUniformCrossTermClaimTests
{
    // Phase 1 verified bondNormSquared values (J = 1, full N-qubit operator-space convention).
    // ‖L_H^bond‖²_F scales as 4 × 4^(N−3) above the N=3 anchor (one extra spectator factor tr(I_4) = 4).
    private const double HeisenbergBondNormSqN3 = 384.0;
    private const double HeisenbergBondNormSqN4 = 1536.0;
    private const double HeisenbergBondNormSqN5 = 6144.0;
    private const double IsingBondNormSqN3 = 128.0;
    private const double XyBondNormSqN3 = 256.0;
    private const double XyBondNormSqN4 = 1024.0;

    // Canonical N=3 γ triple used across the Phase-1 anchor row tests.
    private static double[] GammaN3() => new[] { 0.1, 0.2, 0.3 };

    // Per-bond (bondNormSq, gFraction) arrays for a graph with <paramref name="bondCount"/>
    // Heisenberg bonds sharing the same per-bond superoperator norm.
    private static (double[] bondNorms, double[] gFractions) UniformHeisenberg(int bondCount, double bondNormSq) =>
        (Enumerable.Repeat(bondNormSq, bondCount).ToArray(),
         Enumerable.Repeat(F49NonUniformCrossTermClaim.GHeisenbergFraction, bondCount).ToArray());

    [Fact]
    public void Tier_IsTier1Derived()
    {
        var claim = new F49NonUniformCrossTermClaim();
        Assert.Equal(Tier.Tier1Derived, claim.Tier);
    }

    [Fact]
    public void GHeisenbergFraction_IsExactly4Over3()
    {
        // ZZ is 1/3 of the Heisenberg bond norm ⟹ G/‖L_H^bond‖² = 4·(1/3).
        Assert.Equal(4.0 / 3.0, F49NonUniformCrossTermClaim.GHeisenbergFraction);
    }

    [Fact]
    public void GIsingFraction_IsExactly4()
    {
        // ZZ is 100% of the Ising bond norm ⟹ G/‖L_H^bond‖² = 4·1.
        Assert.Equal(4.0, F49NonUniformCrossTermClaim.GIsingFraction);
    }

    [Fact]
    public void GXyFraction_IsZero()
    {
        // No ZZ in XX+YY ⟹ G/‖L_H^bond‖² = 0.
        Assert.Equal(0.0, F49NonUniformCrossTermClaim.GXyFraction);
    }

    [Fact]
    public void GSoftXyYxFraction_IsZero()
    {
        // No ZZ in XY+YX ⟹ G/‖L_H^bond‖² = 0.
        Assert.Equal(0.0, F49NonUniformCrossTermClaim.GSoftXyYxFraction);
    }

    [Fact]
    public void SpectatorPrefactor_IsExactly4()
    {
        // (d_α^c + d_β^c)² = 4·(spectator-sum)² + A²·(γ_i−γ_j)² after the bond-sum rule.
        Assert.Equal(4.0, F49NonUniformCrossTermClaim.SpectatorPrefactor);
    }

    [Fact]
    public void PredictHeisenbergChain_MatchesPhase1Anchor()
    {
        // Phase 1 anchor at N=3 γ=[0.1, 0.2, 0.3], bondNormSquared=384:
        //   spectator: 4·384·(γ_2²+γ_0²) over bonds (0,1),(1,2) = 4·384·(0.09+0.01) = 153.60
        //   asymmetry: (4/3)·384·((γ_0−γ_1)²+(γ_1−γ_2)²) = (4/3)·384·(0.01+0.01) = 10.24
        //   total: 163.84 (matches truth at F1-centered L_Dc, Section 1 of verify script).
        double predicted = F49NonUniformCrossTermClaim.PredictHeisenbergChain(
            N: 3, gamma: GammaN3(), bondNormSquared: HeisenbergBondNormSqN3);
        Assert.Equal(163.84, predicted, 12);
    }

    [Theory]
    // N-scan with γ_l = 0.05·(l+1), Heisenberg J = 1; matches Section 3 of the verify script bit-exact.
    [InlineData(3, HeisenbergBondNormSqN3, 40.96)]
    [InlineData(4, HeisenbergBondNormSqN4, 737.28)]
    [InlineData(5, HeisenbergBondNormSqN5, 8437.76)]
    public void PredictHeisenbergChain_MatchesPhase1NScan(int N, double bondNormSquared, double expected)
    {
        var gamma = new double[N];
        for (int l = 0; l < N; l++) gamma[l] = 0.05 * (l + 1);
        double predicted = F49NonUniformCrossTermClaim.PredictHeisenbergChain(N, gamma, bondNormSquared);
        Assert.Equal(expected, predicted, 10);
    }

    [Fact]
    public void PredictIsingChain_MatchesPhase1Anchor()
    {
        // N=3 γ=[0.1, 0.2, 0.3] Ising J=1, bondNormSquared=128:
        //   spectator: 4·128·(0.09+0.01) = 51.20
        //   asymmetry: 4·128·((0.01)+(0.01)) = 10.24
        //   total: 61.44 (independently verified in section 4 of the verify script at N=4 cross-class).
        double predicted = F49NonUniformCrossTermClaim.PredictIsingChain(
            N: 3, gamma: GammaN3(), bondNormSquared: IsingBondNormSqN3);
        Assert.Equal(61.44, predicted, 12);
    }

    [Fact]
    public void PredictXyChain_HasZeroAsymmetryPart()
    {
        // XY G-fraction = 0 ⟹ formula reduces to spectator-only. N=3 γ=[0.1, 0.2, 0.3] XY (XX+YY) J=1,
        // bondNormSquared=256: spectator = 4·256·(0.09+0.01) = 102.40; asymmetry = 0.
        double predicted = F49NonUniformCrossTermClaim.PredictXyChain(
            N: 3, gamma: GammaN3(), bondNormSquared: XyBondNormSqN3);
        Assert.Equal(102.40, predicted, 12);
    }

    [Fact]
    public void PredictSoftXyYxChain_HasZeroAsymmetryPart()
    {
        // Soft XY+YX has the same G-fraction = 0 and the same per-bond intrinsic norm structure
        // as XY at N=3 (256 for J=1). Predictions therefore coincide at γ=[0.1, 0.2, 0.3].
        double predictedSoft = F49NonUniformCrossTermClaim.PredictSoftXyYxChain(
            N: 3, gamma: GammaN3(), bondNormSquared: XyBondNormSqN3);
        double predictedXy = F49NonUniformCrossTermClaim.PredictXyChain(
            N: 3, gamma: GammaN3(), bondNormSquared: XyBondNormSqN3);
        Assert.Equal(predictedXy, predictedSoft, 12);
        Assert.Equal(102.40, predictedSoft, 12);
    }

    [Fact]
    public void Predict_GeneralFormula_MatchesChainSpecialization()
    {
        // A Heisenberg chain at N=4 encoded as bondEdges (0,1)(1,2)(2,3) with uniform G-fraction 4/3
        // must give the same value as PredictHeisenbergChain.
        const int N = 4;
        var gamma = new[] { 0.05, 0.10, 0.15, 0.20 };
        var bondEdges = new (int i, int j)[] { (0, 1), (1, 2), (2, 3) };
        var (bondNorms, gFractions) = UniformHeisenberg(bondEdges.Length, HeisenbergBondNormSqN4);
        double general = F49NonUniformCrossTermClaim.Predict(N, gamma, bondEdges, bondNorms, gFractions);
        double specialised = F49NonUniformCrossTermClaim.PredictHeisenbergChain(N, gamma, HeisenbergBondNormSqN4);
        Assert.Equal(specialised, general, 10);
        Assert.Equal(737.28, general, 10);
    }

    [Fact]
    public void Predict_GeneralFormula_HandlesRingTopology()
    {
        // Heisenberg ring at N=4 with 4 bonds (0,1),(1,2),(2,3),(0,3), uniform Heisenberg fractions.
        // Bond (0,3) spectator pair is {1, 2}; (γ_0 − γ_3)² = (0.05 − 0.20)² = 0.0225 carries the
        // largest asymmetry weight.
        // Hand calc with γ = [0.05, 0.10, 0.15, 0.20] and bondNormSquared = 1536:
        //   spec contributions: 4·1536·(0.0625 + 0.0425 + 0.0125 + 0.0325) = 4·1536·0.15 = 921.60
        //   asym contributions: (4/3)·1536·(0.0025 + 0.0025 + 0.0025 + 0.0225) = (4/3)·1536·0.03 = 61.44
        //   total: 983.04
        const int N = 4;
        var gamma = new[] { 0.05, 0.10, 0.15, 0.20 };
        var bondEdges = new (int i, int j)[] { (0, 1), (1, 2), (2, 3), (0, 3) };
        var (bondNorms, gFractions) = UniformHeisenberg(bondEdges.Length, HeisenbergBondNormSqN4);
        double predicted = F49NonUniformCrossTermClaim.Predict(N, gamma, bondEdges, bondNorms, gFractions);
        Assert.Equal(983.04, predicted, 10);
    }

    [Fact]
    public void PredictXyChain_MatchesPhase1CrossClassAnchor()
    {
        // Section 4 of the verify script: N=4 γ=[0.05, 0.10, 0.15, 0.20] XY chain ‖{L_H, L_Dc}‖² = 481.28.
        double predicted = F49NonUniformCrossTermClaim.PredictXyChain(
            N: 4, gamma: new[] { 0.05, 0.10, 0.15, 0.20 }, bondNormSquared: XyBondNormSqN4);
        Assert.Equal(481.28, predicted, 10);
    }

    [Fact]
    public void Predict_RejectsTooSmallN()
    {
        var exChain = Assert.Throws<ArgumentOutOfRangeException>(
            () => F49NonUniformCrossTermClaim.PredictHeisenbergChain(N: 1, new[] { 0.1 }, HeisenbergBondNormSqN3));
        Assert.Contains("N must be ≥ 2", exChain.Message);

        var exGeneral = Assert.Throws<ArgumentOutOfRangeException>(
            () => F49NonUniformCrossTermClaim.Predict(
                N: 1, new[] { 0.1 },
                bondEdges: Array.Empty<(int, int)>(),
                bondNormSquaredPerBond: Array.Empty<double>(),
                gFractionPerBond: Array.Empty<double>()));
        Assert.Contains("N must be ≥ 2", exGeneral.Message);
    }

    [Fact]
    public void Predict_RejectsMismatchedArrayLengths()
    {
        // Chain: gamma.Count != N.
        var exChain = Assert.Throws<ArgumentException>(
            () => F49NonUniformCrossTermClaim.PredictHeisenbergChain(
                N: 3, gamma: new[] { 0.1, 0.2 }, bondNormSquared: HeisenbergBondNormSqN3));
        Assert.Contains("must equal N", exChain.Message);

        // General: bondNormSquaredPerBond.Count != bondEdges.Count.
        var heisenbergG = F49NonUniformCrossTermClaim.GHeisenbergFraction;
        var exBondNorm = Assert.Throws<ArgumentException>(
            () => F49NonUniformCrossTermClaim.Predict(
                N: 3, gamma: GammaN3(),
                bondEdges: new (int, int)[] { (0, 1), (1, 2) },
                bondNormSquaredPerBond: new[] { HeisenbergBondNormSqN3 },
                gFractionPerBond: new[] { heisenbergG, heisenbergG }));
        Assert.Contains("must equal bondEdges length", exBondNorm.Message);

        // General: gFractionPerBond.Count != bondEdges.Count.
        var exGFraction = Assert.Throws<ArgumentException>(
            () => F49NonUniformCrossTermClaim.Predict(
                N: 3, gamma: GammaN3(),
                bondEdges: new (int, int)[] { (0, 1), (1, 2) },
                bondNormSquaredPerBond: new[] { HeisenbergBondNormSqN3, HeisenbergBondNormSqN3 },
                gFractionPerBond: new[] { heisenbergG }));
        Assert.Contains("must equal bondEdges length", exGFraction.Message);

        // General: gamma.Count != N.
        var exGamma = Assert.Throws<ArgumentException>(
            () => F49NonUniformCrossTermClaim.Predict(
                N: 3, gamma: new[] { 0.1, 0.2 },
                bondEdges: new (int, int)[] { (0, 1), (1, 2) },
                bondNormSquaredPerBond: new[] { HeisenbergBondNormSqN3, HeisenbergBondNormSqN3 },
                gFractionPerBond: new[] { heisenbergG, heisenbergG }));
        Assert.Contains("must equal N", exGamma.Message);
    }

    [Fact]
    public void Predict_RejectsInvalidBondIndices()
    {
        var gamma = GammaN3();
        var bondNorms = new[] { HeisenbergBondNormSqN3 };
        var gFractions = new[] { F49NonUniformCrossTermClaim.GHeisenbergFraction };

        // Out of range (j >= N).
        var exOob = Assert.Throws<ArgumentException>(
            () => F49NonUniformCrossTermClaim.Predict(
                N: 3, gamma,
                bondEdges: new (int, int)[] { (0, 3) },
                bondNormSquaredPerBond: bondNorms, gFractionPerBond: gFractions));
        Assert.Contains("out of range", exOob.Message);

        // Out of range (negative i).
        var exNeg = Assert.Throws<ArgumentException>(
            () => F49NonUniformCrossTermClaim.Predict(
                N: 3, gamma,
                bondEdges: new (int, int)[] { (-1, 1) },
                bondNormSquaredPerBond: bondNorms, gFractionPerBond: gFractions));
        Assert.Contains("out of range", exNeg.Message);

        // Violates i < j ordering (i == j).
        var exEq = Assert.Throws<ArgumentException>(
            () => F49NonUniformCrossTermClaim.Predict(
                N: 3, gamma,
                bondEdges: new (int, int)[] { (1, 1) },
                bondNormSquaredPerBond: bondNorms, gFractionPerBond: gFractions));
        Assert.Contains("violates i < j ordering", exEq.Message);

        // Violates i < j ordering (i > j).
        var exSwap = Assert.Throws<ArgumentException>(
            () => F49NonUniformCrossTermClaim.Predict(
                N: 3, gamma,
                bondEdges: new (int, int)[] { (2, 1) },
                bondNormSquaredPerBond: bondNorms, gFractionPerBond: gFractions));
        Assert.Contains("violates i < j ordering", exSwap.Message);
    }

    [Fact]
    public void PredictHeisenbergChain_RejectsNegativeBondNormSquared()
    {
        var ex = Assert.Throws<ArgumentOutOfRangeException>(
            () => F49NonUniformCrossTermClaim.PredictHeisenbergChain(
                N: 3, gamma: GammaN3(), bondNormSquared: -1.0));
        Assert.Contains("bondNormSquared must be ≥ 0", ex.Message);
    }

    [Fact]
    public void Anchor_References_AllRequiredSources()
    {
        var claim = new F49NonUniformCrossTermClaim();
        Assert.Contains("PROOF_F49_NONUNIFORM_GAMMA_EXTENSION", claim.Anchor);
        Assert.Contains("PROOF_CROSS_TERM_FORMULA", claim.Anchor);
        Assert.Contains("PROOF_F1_NONUNIFORM_GAMMA", claim.Anchor);
        Assert.Contains("F1T1ResidualClosedForm", claim.Anchor);
        Assert.Contains("F1DepolResidualClosedForm", claim.Anchor);
    }

    [Fact]
    public void ExtraChildren_ExposeKeyStructuralNodes()
    {
        var claim = new F49NonUniformCrossTermClaim();
        var names = claim.Children.Select(c => c.DisplayName).ToList();
        Assert.Contains("tier", names);
        Assert.Contains("anchor", names);
        Assert.Contains("statement", names);
        Assert.Contains("spectator prefactor", names);
        Assert.Contains("G(bond, Heisenberg) / ‖L_H^bond‖² = 4/3 (ZZ is 1/3 of bond norm)", names);
        Assert.Contains("G(bond, Ising) / ‖L_H^bond‖² = 4 (ZZ is 100% of bond norm)", names);
        Assert.Contains("G(bond, XY) / ‖L_H^bond‖² = 0 (no ZZ content)", names);
        Assert.Contains("G(bond, soft XY+YX) / ‖L_H^bond‖² = 0 (no ZZ content)", names);
        Assert.Contains("spectator-vs-asymmetry split", names);
        Assert.Contains("A-classification per bond Pauli class", names);
        Assert.Contains("bondNormSquared convention", names);
        Assert.Contains("uniform-γ recovery", names);
        Assert.Contains("verification", names);
        Assert.Contains("closes follow-up", names);
    }
}
