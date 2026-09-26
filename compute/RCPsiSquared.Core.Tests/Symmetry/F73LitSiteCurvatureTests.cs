using RCPsiSquared.Core.Numerics;
using RCPsiSquared.Core.Symmetry;
using Bond = RCPsiSquared.Core.Symmetry.F73LitSiteCurvature.Bond;

namespace RCPsiSquared.Core.Tests.Symmetry;

/// <summary>The F73 lit-site corollary, gated exactly over ℚ(i): at t = 0 the F73 term and the
/// site purity of one lit site both have first derivative −2γ₀ (no H in it) and second derivative
/// 8γ₀² − V_l with V_l = 4·Σ_{b∋l} J_b², whatever Δ_b and the longitudinal fields are, and the
/// F73 term is ½e^(−4γ₀t) times a γ₀-free function (gated to sixth order), and the purity
/// agrees with it through the third derivative only. Every comparison is ==,
/// since every quantity is an exact rational. The pairing term G(XX − YY) on the lit bond is the
/// mutation (it breaks the law at G = 1/2, 3, −1 and restores it at G = 2J); a transverse field is
/// recorded as leaving the second derivatives standing.</summary>
public class F73LitSiteCurvatureTests
{
    private static BigRational Q(long p, long q = 1) => new(p, q);

    private static Bond B(int a, int b, BigRational j, BigRational delta) => new(a, b, j, delta);

    private static IReadOnlyList<Bond> Chain(int n, BigRational delta) =>
        Enumerable.Range(0, n - 1).Select(i => B(i, i + 1, 1, delta)).ToList();

    private static IReadOnlyList<Bond> Star(int n, BigRational delta) =>
        Enumerable.Range(1, n - 1).Select(i => B(0, i, 1, delta)).ToList();

    private static IReadOnlyList<Bond> Complete(int n, BigRational delta)
    {
        var list = new List<Bond>();
        for (int a = 0; a < n; a++) for (int b = a + 1; b < n; b++) list.Add(B(a, b, 1, delta));
        return list;
    }

    // Ring N = 4, inhomogeneous J and a different Δ on every bond.
    private static IReadOnlyList<Bond> InhomogeneousRing() => new List<Bond>
    {
        B(0, 1, 1, Q(0)), B(1, 2, Q(3, 2), Q(2)), B(2, 3, Q(1, 2), Q(-7, 3)), B(3, 0, 2, 1),
    };

    public static IEnumerable<object[]> Cases()
    {
        // (graph, n, Δ numerator, Δ denominator, lit site, degree of the lit site)
        foreach (var (dn, dd) in new[] { (0L, 1L), (1L, 1L), (2L, 1L), (-7L, 3L) })
        {
            yield return new object[] { "chain", 4, dn, dd, 0, 1 };
            yield return new object[] { "chain", 4, dn, dd, 1, 2 };
            yield return new object[] { "star", 5, dn, dd, 0, 4 };
            yield return new object[] { "star", 5, dn, dd, 3, 1 };
            yield return new object[] { "complete", 4, dn, dd, 2, 3 };
        }
    }

    private static IReadOnlyList<Bond> Graph(string label, int n, BigRational delta) => label switch
    {
        "chain" => Chain(n, delta),
        "star" => Star(n, delta),
        "complete" => Complete(n, delta),
        _ => throw new ArgumentException(label),
    };

    private static void AssertCorollary(int n, IReadOnlyList<Bond> bonds, int lit, BigRational gamma,
        IReadOnlyList<BigRational>? hz = null)
    {
        var h = F73LitSiteCurvature.Hamiltonian(n, bonds, hz);
        var v = F73LitSiteCurvature.PredictedVariance(bonds, lit);
        Assert.True(F73LitSiteCurvature.SingleExcitationVariance(h, lit) == v, "V_l read off H ≠ 4·Σ_{b∋l} J_b²");

        var t = F73LitSiteCurvature.TaylorAtZero(n, h, gamma, lit);
        var second = F73LitSiteCurvature.PredictedSecondDerivative(gamma, v);
        Assert.True(t.P0Term == Q(1, 2), $"p_l(0) = {t.P0Term}");
        Assert.True(t.P0Purity == 1, $"P_l(0) = {t.P0Purity}");
        Assert.True(t.P1Term == -2 * gamma, $"p_l'(0) = {t.P1Term}, expected {-2 * gamma}");
        Assert.True(t.P1Purity == -2 * gamma, $"P_l'(0) = {t.P1Purity}, expected {-2 * gamma}");
        Assert.True(t.P2Term == second, $"p_l''(0) = {t.P2Term}, expected {second}");
        Assert.True(t.P2Purity == second, $"P_l''(0) = {t.P2Purity}, expected {second}");
    }

    [Theory]
    [MemberData(nameof(Cases))]
    public void LitSite_AtGammaZero_LosesTwoDegJSquaredTSquared(string label, int n, long dn, long dd, int lit, int degree)
    {
        var bonds = Graph(label, n, Q(dn, dd));
        AssertCorollary(n, bonds, lit, BigRational.Zero);
        // Uniform J = 1: loss = (V_l/2)·t² = 2·deg(l)·t², i.e. P'' = −4·deg(l).
        var t = F73LitSiteCurvature.TaylorAtZero(n, F73LitSiteCurvature.Hamiltonian(n, bonds), BigRational.Zero, lit);
        Assert.True(t.P2Purity == -4 * degree, $"{label} deg {degree}: P'' = {t.P2Purity}");
        Assert.True(t.P1Purity.IsZero, "at γ₀ = 0 the first derivative must vanish: no rate");
    }

    [Theory]
    [MemberData(nameof(Cases))]
    public void LitSite_AtGammaPositive_RateIsTwoGammaAndCurvatureIsEightGammaSquaredMinusV(string label, int n, long dn, long dd, int lit, int degree)
    {
        var gamma = Q(1, 20);
        var bonds = Graph(label, n, Q(dn, dd));
        AssertCorollary(n, bonds, lit, gamma);
        var t = F73LitSiteCurvature.TaylorAtZero(n, F73LitSiteCurvature.Hamiltonian(n, bonds), gamma, lit);
        Assert.True(t.P2Purity == 8 * gamma * gamma - 4 * degree, $"{label} deg {degree}: P'' = {t.P2Purity}");
    }

    [Fact]
    public void InhomogeneousRing_WithLongitudinalFields_HoldsAtEverySite()
    {
        var bonds = InhomogeneousRing();
        var hz = new[] { Q(3, 5), Q(-1), Q(0), Q(9, 4) };
        for (int lit = 0; lit < 4; lit++)
        {
            AssertCorollary(4, bonds, lit, BigRational.Zero, hz);
            AssertCorollary(4, bonds, lit, Q(1, 7), hz);
        }
        // Site 1 carries J = 1 and J = 3/2: V = 4·(1 + 9/4) = 13.
        Assert.True(F73LitSiteCurvature.PredictedVariance(bonds, 1) == 13);
    }

    [Fact]
    public void DeltaIndependence_SecondDerivativeIdenticalAcrossDelta()
    {
        var reference = F73LitSiteCurvature.TaylorAtZero(4,
            F73LitSiteCurvature.Hamiltonian(4, Chain(4, 0)), BigRational.Zero, 1).P2Purity;
        foreach (var delta in new[] { Q(1), Q(5), Q(-13, 6) })
        {
            var p2 = F73LitSiteCurvature.TaylorAtZero(4,
                F73LitSiteCurvature.Hamiltonian(4, Chain(4, delta)), BigRational.Zero, 1).P2Purity;
            Assert.True(p2 == reference, $"Δ = {delta}: P'' = {p2}, Δ = 0 gave {reference}");
        }
    }

    [Theory]
    [InlineData(0L, 1L, 1L, 2L)]
    [InlineData(1L, 20L, 1L, 2L)]
    [InlineData(0L, 1L, 3L, 1L)]
    [InlineData(1L, 20L, -1L, 1L)]
    public void Mutation_PairingTermOutsideU1_BreaksTheSecondDerivativeLaw(long gn, long gd, long pn, long pd)
    {
        // J_x ≠ J_y on the lit site's bond: G·(XX − YY) creates |1_l 1_j⟩ out of |vac⟩, so
        // [H, N_total] ≠ 0. On this chain the lit bond's J² becomes (J − G)², so the in-class law
        // 8γ₀² − 4·Σ J_b² fails on both readings for these G (none is 0 or 2J).
        var gamma = Q(gn, gd);
        var g = Q(pn, pd);
        var inClass = Chain(3, 1);
        var mutated = inClass.Select((b, i) => i == 0 ? b with { G = g } : b).ToList();
        var predicted = F73LitSiteCurvature.PredictedSecondDerivative(gamma, F73LitSiteCurvature.PredictedVariance(inClass, 0));
        var t = F73LitSiteCurvature.TaylorAtZero(3, F73LitSiteCurvature.Hamiltonian(3, mutated), gamma, 0);
        Assert.True(t.P2Term != predicted, $"p'' = {t.P2Term} equals the in-class prediction {predicted} outside the class");
        Assert.True(t.P2Purity != predicted, $"P'' = {t.P2Purity} equals the in-class prediction {predicted} outside the class");
        // What the mutation does instead, read exactly: J² → (J − G)² on the lit bond.
        var shifted = 8 * gamma * gamma - 4 * (1 - g) * (1 - g);
        Assert.True(t.P2Term == shifted && t.P2Purity == shifted, $"p''={t.P2Term}, P''={t.P2Purity}, expected {shifted}");
        // The same bond with G switched off is back inside the class and the law holds again.
        var t0 = F73LitSiteCurvature.TaylorAtZero(3, F73LitSiteCurvature.Hamiltonian(3, inClass), gamma, 0);
        Assert.True(t0.P2Term == predicted && t0.P2Purity == predicted);
    }

    [Theory]
    [InlineData(0L, 1L)]
    [InlineData(1L, 20L)]
    public void PairingTermAtGEqualsTwoJ_RestoresTheSecondDerivativeLaw(long gn, long gd)
    {
        // The break is not general: at G = 2J (J_x = 3J, J_y = −J) the lit bond's (J − G)² = J²,
        // and the in-class value returns although [H, N_total] ≠ 0.
        var gamma = Q(gn, gd);
        var inClass = Chain(3, 1);
        var mutated = inClass.Select((b, i) => i == 0 ? b with { G = 2 } : b).ToList();
        var predicted = F73LitSiteCurvature.PredictedSecondDerivative(gamma, F73LitSiteCurvature.PredictedVariance(inClass, 0));
        var t = F73LitSiteCurvature.TaylorAtZero(3, F73LitSiteCurvature.Hamiltonian(3, mutated), gamma, 0);
        Assert.True(t.P2Term == predicted && t.P2Purity == predicted, $"p''={t.P2Term}, P''={t.P2Purity}, predicted {predicted}");
    }

    [Theory]
    [InlineData("chain", 4, 1, 1L, 20L)]
    [InlineData("star", 4, 0, 3L, 7L)]
    [InlineData("complete", 4, 2, 1L, 3L)]
    public void TermIsHalfExpMinusFourGammaTTimesAGammaFreeFunction_ToSixthOrder(string label, int n, int lit, long gn, long gd)
    {
        // p_l(t) = ½e^(−4γ₀t)·S(t) with S γ₀-free, so every derivative of p_l at γ₀ follows from the
        // γ₀ = 0 derivatives by Leibniz. Gated exactly through k = 6 at Δ = 2.
        var gamma = Q(gn, gd);
        var h = F73LitSiteCurvature.Hamiltonian(n, Graph(label, n, 2));
        var free = F73LitSiteCurvature.Derivatives(n, h, BigRational.Zero, lit, 6);
        var damped = F73LitSiteCurvature.Derivatives(n, h, gamma, lit, 6);
        for (int k = 0; k <= 6; k++)
        {
            var expected = F73LitSiteCurvature.DampedDerivative(free.Term, gamma, k);
            Assert.True(damped.Term[k] == expected, $"{label}: p^({k})(0) = {damped.Term[k]}, factorized {expected}");
        }
        // The purity does not factor this way: its excitation weight has n_l'''(0) = 4γ₀·V_l, and the
        // purity agrees with the term through k = 3, departing at k = 4 by 12·V_l².
        var v = F73LitSiteCurvature.PredictedVariance(Graph(label, n, 2), lit);
        foreach (var (series, g) in new[] { (free, BigRational.Zero), (damped, gamma) })
        {
            Assert.True(series.Weight[3] == 4 * g * v, $"n_l'''(0) = {series.Weight[3]}, expected {4 * g * v}");
            for (int k = 1; k <= 3; k++)
                Assert.True(series.Purity[k] == series.Term[k], $"k={k}: P={series.Purity[k]}, p={series.Term[k]}");
            Assert.True(series.Purity[4] - series.Term[4] == 12 * v * v, $"P''''−p'''' = {series.Purity[4] - series.Term[4]}, expected {12 * v * v}");
        }
    }

    [Theory]
    [InlineData(0L, 1L)]
    [InlineData(1L, 20L)]
    public void TransverseField_LeavesTheLawStanding_SoU1IsSufficientNotNecessary(long gn, long gd)
    {
        // Recorded, not assumed: a transverse field on every site, the lit one included, is outside
        // F73's class and still leaves both second derivatives at 8γ₀² − V_l on these graphs.
        var gamma = Q(gn, gd);
        foreach (var (n, bonds) in new[] { (3, Chain(3, 1)), (4, Star(4, 2)) })
        {
            var hx = Enumerable.Range(0, n).Select(i => Q(i + 1, 3)).ToArray();
            var predicted = F73LitSiteCurvature.PredictedSecondDerivative(gamma, F73LitSiteCurvature.PredictedVariance(bonds, 0));
            var t = F73LitSiteCurvature.TaylorAtZero(n, F73LitSiteCurvature.Hamiltonian(n, bonds, hx: hx), gamma, 0);
            Assert.True(t.P2Term == predicted && t.P2Purity == predicted, $"n={n}: p''={t.P2Term}, P''={t.P2Purity}, predicted {predicted}");
        }
    }
}
