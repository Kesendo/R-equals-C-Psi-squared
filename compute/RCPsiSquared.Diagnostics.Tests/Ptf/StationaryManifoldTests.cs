using RCPsiSquared.Core.ChainSystems;
using RCPsiSquared.Core.Lindblad;
using RCPsiSquared.Core.Pauli;
using RCPsiSquared.Diagnostics.Ptf;
using ComplexMatrix = MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex>;

namespace RCPsiSquared.Diagnostics.Tests.Ptf;

/// <summary>The closure guardian's six pinned cases (Edge 5, 2026-06-10): the PTF closure
/// holds iff the Liouvillian keeps its complete N+1 stationary sector manifold, counted in
/// any frame.
///
/// <para>All six integer counts below were pinned by a Python convention check at N = 4,
/// γ = 0.05, J = 1 against the Edge-5 probe's builders
/// (<c>simulations/ptf_carrier_seam_retrodiction.py</c>: H_A = (J/2)·Σ_b (XX+YY),
/// perturbations ε·V with V at unit coefficient per term on every bond, replacement chain
/// H = V alone, stationary tolerance |λ| &lt; 10⁻¹⁰). The measured gap is wide on both sides
/// of the tolerance: surviving projectors sit at |λ| ≤ 2·10⁻¹⁵, lifted ones at
/// |λ| ≥ 3.3·10⁻³ (X-field) / 1.41·10⁻¹ (XY+YX perturbation). Pinned counts:
/// baseline 5, J-defect 5, Z-field 5, X-field 1, +ε(XY+YX) 2 (the parity remnant),
/// replacement XY+YX 5 (the hidden sublattice-flip U(1)).</para>
///
/// <para>Sources: <c>experiments/PTF_PALINDROME_BREAKING_PERTURBATIONS.md</c> ("Edge 5
/// result"), <c>reflections/ON_HOW_THE_CARRIER_SHOWS_ITSELF.md</c> ("Seen again
/// 2026-06-10").</para></summary>
public class StationaryManifoldTests
{
    private const int N = 4;
    private const double Gamma = 0.05;
    private const double Eps = 0.10;

    private static readonly double[] Gammas = Enumerable.Repeat(Gamma, N).ToArray();
    private static readonly IReadOnlyList<Bond> ChainBonds =
        Enumerable.Range(0, N - 1).Select(i => new Bond(i, i + 1, 1.0)).ToArray();

    /// <summary>Baseline H_A = (J/2)·Σ_b (XX+YY) at J = 1, the canonical chain.</summary>
    private static ComplexMatrix BaselineH() => PauliHamiltonian.XYChain(N, 1.0).ToMatrix();

    /// <summary>Edge-5 perturbation convention: V = Σ_bonds Σ_terms σσ at unit coefficient,
    /// scaled here by ε directly in the term coefficients.</summary>
    private static ComplexMatrix Bilinear(PauliLetter a, PauliLetter b, double coeff) =>
        PauliHamiltonian.Bilinear(N, ChainBonds,
            new[] { (a, b, (System.Numerics.Complex)coeff), (b, a, (System.Numerics.Complex)coeff) }).ToMatrix();

    private static void AssertPinned(ComplexMatrix liouvillian, int expectedCount, bool expectedHolds)
    {
        Assert.Equal(expectedCount, StationaryManifold.Count(liouvillian));
        var forecast = StationaryManifold.Forecast(N, liouvillian);
        Assert.Equal(expectedCount, forecast.StationaryCount);
        Assert.Equal(N + 1, forecast.FullManifold);
        Assert.Equal(expectedHolds, forecast.PredictsClosureHolds);
    }

    [Fact]
    public void Baseline_UniformXYChain_KeepsFullManifold()
    {
        // Pinned: 5 stationary (= N+1, the F4 excitation-sector projectors); next |λ| = 0.160.
        var l = new ChainSystem(N, J: 1.0, GammaZero: Gamma).BuildLiouvillian();
        AssertPinned(l, expectedCount: 5, expectedHolds: true);
    }

    [Fact]
    public void JDefect_OnBond01_KeepsFullManifold()
    {
        // δJ = 0.1 on bond (0, 1): the hopping bilinear conserves U(1), all 5 projectors
        // survive. Pinned: 5 stationary; next |λ| = 0.157. The canonical hold control of the
        // Edge-5 probe (S = +0.097 vs +8.3..+11.3 for the breaking bilinears).
        var h = PauliHamiltonian.XYChain(N, new[] { 1.0 + Eps, 1.0, 1.0 }).ToMatrix();
        AssertPinned(PauliDephasingDissipator.BuildZ(h, Gammas), expectedCount: 5, expectedHolds: true);
    }

    [Fact]
    public void ZField_OnSite0_KeepsFullManifold()
    {
        // ε·Z₀ commutes with N̂ (and with the dephasers): full manifold intact.
        // Pinned: 5 stationary; next |λ| = 0.161.
        var h = BaselineH() + Eps * PauliString.SiteOp(N, 0, PauliLetter.Z);
        AssertPinned(PauliDephasingDissipator.BuildZ(h, Gammas), expectedCount: 5, expectedHolds: true);
    }

    [Fact]
    public void XField_OnSite0_CollapsesToSingleSteadyState()
    {
        // ε·X₀ breaks U(1) with ΔN = ±1 content: only the trace steady state survives.
        // Pinned: 1 stationary; lifted projectors start at |λ| = 3.3·10⁻³ (still far above
        // the 10⁻¹⁰ tolerance: the verdict is robust).
        var h = BaselineH() + Eps * PauliString.SiteOp(N, 0, PauliLetter.X);
        AssertPinned(PauliDephasingDissipator.BuildZ(h, Gammas), expectedCount: 1, expectedHolds: false);
    }

    [Fact]
    public void XYplusYXPerturbation_KeepsOnlyParityRemnant_AndThatBuysNothing()
    {
        // H_A + ε·Σ_b (XY+YX): ΔN = ±2, parity-conserving, so exactly the 2 parity projectors
        // survive (pinned: 2 stationary; lifted at |λ| = 0.141). This is the refutation case
        // for the weaker guardian: despite keeping the parity remnant it breaks the closure
        // FASTEST of the five tour bilinears (canonical protocol S = +11.3). The forecast is
        // the equality test, so 2 < 5 correctly predicts the break.
        var h = BaselineH() + Bilinear(PauliLetter.X, PauliLetter.Y, Eps);
        AssertPinned(PauliDephasingDissipator.BuildZ(h, Gammas), expectedCount: 2, expectedHolds: false);
    }

    [Fact]
    public void ReplacementPureXYplusYXChain_HiddenU1_KeepsFullManifold()
    {
        // The tour's replacement reading: H = Σ_b (XY+YX) at unit coefficient INSTEAD of
        // XX+YY. The bare [H, N̂] ≠ 0 lens misfires here; the sublattice X-flip frame
        // (mask 0101) maps this chain to a U(1)-conserving one while commuting with the
        // Z-dephasers, and the frame-independent count sees it without any frame search.
        // Pinned: 5 stationary (full manifold); next |λ| = 0.160.
        var h = Bilinear(PauliLetter.X, PauliLetter.Y, 1.0);
        AssertPinned(PauliDephasingDissipator.BuildZ(h, Gammas), expectedCount: 5, expectedHolds: true);
    }
}
