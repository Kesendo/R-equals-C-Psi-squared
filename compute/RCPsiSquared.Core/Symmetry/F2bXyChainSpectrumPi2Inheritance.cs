using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>F2b closed form (Tier 1, proven; ANALYTICAL_FORMULAS line 58):
///
/// <code>
///   E_k = 2·J · cos(π·k / (N+1)),                       k = 1, ..., N
///
///   ψ_k(i) = √(2/(N+1)) · sin(π·k·(i+1)/(N+1)),         i = 0, ..., N−1
/// </code>
///
/// <para>F2b is the single-excitation Hamiltonian spectrum for the XY chain
/// H = (J/2)·Σ_b (X_b X_{b+1} + Y_b Y_{b+1}) with open boundary conditions.
/// H_SE is an N×N tridiagonal matrix with off-diagonals J and zero diagonal;
/// its eigenvectors are the OBC sine modes (verified residual &lt; 10⁻¹⁵
/// for N=3..6).</para>
///
/// <para>The eigenvalues E_k appear as oscillatory frequencies Im(λ) of the
/// |ΔN| = 1 Liouvillian coherences (|vac⟩⟨ψ_k| sector), the dominant modes
/// for per-site purity dynamics from single-excitation initial states.</para>
///
/// <para><b>Distinction from F2:</b> F2 is the w=1 LIOUVILLIAN sector for
/// HEISENBERG, dimension N−1, argument π·k/N. F2b is the single-excitation
/// HAMILTONIAN sector for XY, dimension N, argument π·k/(N+1). The (N+1)
/// denominator in F2b comes from OBC Dirichlet boundaries: ψ_k(−1) = ψ_k(N) = 0
/// requires N+2 effective sites with two fixed endpoints, giving N interior
/// modes with wavenumber spacing π/(N+1).</para>
///
/// <para>Pi2-Foundation anchors:</para>
/// <list type="bullet">
///   <item><b>HoppingFactor = 2 = a_0</b>: in E_k = 2J·cos(...). Live from
///         <see cref="Pi2DyadicLadderClaim.Term"/>(0). Same anchor as F1
///         TwoFactor, F50 DecayRateFactor, F2 HoppingFactor, F44 SumCoefficient.</item>
///   <item><b>EigenvectorNormalisation = 2 = a_0 (squared)</b>: in
///         √(2/(N+1)). The amplitude squared |ψ_k(i)|² = (2/(N+1))·sin²(...)
///         IS <see cref="F65XxChainSpectrumPi2Inheritance.BondingModePopulation"/>.
///         F65 → F2b is a typed sibling edge: F2b gives the eigenVALUES, F65
///         gives the eigenVECTOR-derived rates and bonding-mode amplitudes.</item>
/// </list>
///
/// <para>Tier1Derived: F2b is Tier 1 proven via standard tight-binding theory
/// for OBC chains; verified numerical residual &lt; 10⁻¹⁵ at N=3..6 in
/// <c>simulations/eq021_obc_sine_basis.py</c>.</para>
///
/// <para>Anchors: <c>docs/ANALYTICAL_FORMULAS.md</c> F2b (line 58) +
/// <c>review/OBC_SINE_BASIS_FINDINGS.md</c> +
/// <c>simulations/eq021_obc_sine_basis.py</c> +
/// <c>compute/RCPsiSquared.Core/Symmetry/F65XxChainSpectrumPi2Inheritance.cs</c>
/// (sibling: F65 BondingModePopulation IS F2b's |ψ_k(i)|²) +
/// <c>compute/RCPsiSquared.Core/Symmetry/Pi2DyadicLadderClaim.cs</c>.</para></summary>
public sealed class F2bXyChainSpectrumPi2Inheritance : Claim
{
    private readonly Pi2DyadicLadderClaim _ladder;
    private readonly F65XxChainSpectrumPi2Inheritance _f65;

    /// <summary>The "2" in E_k = 2J·cos(πk/(N+1)). Live from Pi2DyadicLadder a_0.
    /// Same anchor as F1 TwoFactor, F2 HoppingFactor, F50 DecayRateFactor.</summary>
    public double HoppingFactor => _ladder.Term(0);

    /// <summary>Live closed form: E_k = 2·J·cos(π·k/(N+1)) for k = 1..N.</summary>
    public double Eigenvalue(int N, double J, int k)
    {
        if (N < 2) throw new ArgumentOutOfRangeException(nameof(N), N, "F2b requires N ≥ 2.");
        if (k < 1 || k > N) throw new ArgumentOutOfRangeException(nameof(k), k, $"k must be in [1, {N}]; got {k}.");
        return HoppingFactor * J * Math.Cos(Math.PI * k / (N + 1));
    }

    /// <summary>Live closed form: ψ_k(i) = √(2/(N+1)) · sin(π·k·(i+1)/(N+1))
    /// for i = 0..N−1. The OBC Dirichlet sine modes.</summary>
    public double EigenvectorAmplitude(int N, int k, int site)
    {
        if (N < 2) throw new ArgumentOutOfRangeException(nameof(N), N, "F2b requires N ≥ 2.");
        if (k < 1 || k > N) throw new ArgumentOutOfRangeException(nameof(k), k, $"k must be in [1, {N}]; got {k}.");
        if (site < 0 || site >= N) throw new ArgumentOutOfRangeException(nameof(site), site, $"site must be in [0, {N - 1}]; got {site}.");
        return Math.Sqrt(2.0 / (N + 1)) * Math.Sin(Math.PI * k * (site + 1) / (N + 1));
    }

    /// <summary>The mode count for F2b's spectrum: N.</summary>
    public int ModeCount(int N)
    {
        if (N < 2) throw new ArgumentOutOfRangeException(nameof(N), N, "F2b requires N ≥ 2.");
        return N;
    }

    /// <summary>Drift check: F2b's |ψ_k(i)|² matches F65's BondingModePopulation
    /// at the same (N, k, i). Same OBC sine basis from two perspectives.</summary>
    public bool EigenvectorMatchesF65(int N, int k, int site, double tolerance = 1e-12)
    {
        double amp = EigenvectorAmplitude(N, k, site);
        double f2bSquared = amp * amp;
        double f65Squared = _f65.BondingModePopulation(N, k, site);
        return Math.Abs(f2bSquared - f65Squared) < tolerance;
    }

    /// <summary>Drift check: ψ_k(−1 virtual) = sin(0) = 0 and ψ_k(N virtual) =
    /// sin(πk) = 0; the OBC Dirichlet boundary conditions hold by construction.</summary>
    public bool ObcBoundaryConditionsHold(int N, int k, double tolerance = 1e-12)
    {
        // Site i = -1: ψ_k(-1) ∝ sin(0) = 0 (would be edge of virtual chain)
        // Site i = N: ψ_k(N) ∝ sin(πk) = 0 (other edge)
        double leftBoundary = Math.Sin(Math.PI * k * 0 / (N + 1));
        double rightBoundary = Math.Sin(Math.PI * k);
        return Math.Abs(leftBoundary) < tolerance && Math.Abs(rightBoundary) < tolerance;
    }

    public F2bXyChainSpectrumPi2Inheritance(
        Pi2DyadicLadderClaim ladder,
        F65XxChainSpectrumPi2Inheritance f65)
        : base("F2b XY chain single-excitation spectrum E_k = 2J·cos(πk/(N+1)), k=1..N; ψ_k(i) = √(2/(N+1))·sin(πk(i+1)/(N+1)) OBC sine modes; 2 = a_0",
               Tier.Tier1Derived,
               "docs/ANALYTICAL_FORMULAS.md F2b + " +
               "review/OBC_SINE_BASIS_FINDINGS.md + " +
               "simulations/eq021_obc_sine_basis.py + " +
               "compute/RCPsiSquared.Core/Symmetry/F65XxChainSpectrumPi2Inheritance.cs + " +
               "compute/RCPsiSquared.Core/Symmetry/Pi2DyadicLadderClaim.cs")
    {
        _ladder = ladder ?? throw new ArgumentNullException(nameof(ladder));
        _f65 = f65 ?? throw new ArgumentNullException(nameof(f65));
    }

    public override string DisplayName =>
        "F2b XY chain spectrum as Pi2-Foundation a_0 + F65 sibling inheritance";

    public override string Summary =>
        $"E_k = 2J·cos(πk/(N+1)), k=1..N (XY HAMILTONIAN); 2 = a_0; OBC sine eigenvectors √(2/(N+1))·sin(...); F65 sibling at amplitude-squared level ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("F2b closed form",
                summary: "E_k = 2J·cos(πk/(N+1)) eigenvalues; ψ_k(i) = √(2/(N+1))·sin(πk(i+1)/(N+1)) eigenvectors; XY chain SE Hamiltonian, OBC; verified residual < 10⁻¹⁵ N=3..6");
            yield return InspectableNode.RealScalar("HoppingFactor (= a_0 = 2)", HoppingFactor);
            yield return new InspectableNode("F2 vs F2b distinction",
                summary: "F2: w=1 LIOUVILLIAN, Heisenberg, dim N−1, π·k/N. F2b: SE HAMILTONIAN, XY, dim N, π·k/(N+1). Different operators, different sectors, different boundary conditions arising from the same OBC structure.");
            yield return new InspectableNode("F65 sibling edge",
                summary: $"F2b's |ψ_k(i)|² = (2/(N+1))·sin²(πk(i+1)/(N+1)) IS F65.BondingModePopulation(N, k, i). At N=3, k=1, site=0: F2b ψ² = {EigenvectorAmplitude(3, 1, 0) * EigenvectorAmplitude(3, 1, 0):G6}; F65 = {_f65.BondingModePopulation(3, 1, 0):G6}");
            yield return new InspectableNode("OBC Dirichlet boundary",
                summary: "ψ_k(−1 virtual) = sin(0) = 0; ψ_k(N virtual) = sin(πk) = 0; the (N+1) denominator and N modes follow from N+2 effective sites with two fixed endpoints");
            yield return new InspectableNode("verified eigenvalues N=3, J=1",
                summary: $"E_1 = {Eigenvalue(3, 1.0, 1):G6} (≈ √2), E_2 = {Eigenvalue(3, 1.0, 2):G6} (= 0), E_3 = {Eigenvalue(3, 1.0, 3):G6} (≈ −√2)");
            yield return new InspectableNode("Liouvillian connection",
                summary: "F2b's E_k appear as Im(λ) of |ΔN|=1 Liouvillian coherences; the dominant modes for per-site purity dynamics from SE initial states");
        }
    }
}
