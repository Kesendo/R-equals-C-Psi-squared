using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Spectrum;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>F2 closed form (Tier 1, proven D10; ANALYTICAL_FORMULAS line 32):
///
/// <code>
///   ω_k = 4·J · (1 − cos(π·k / N)),      k = 1, ..., N−1
///
///   N−1 distinct frequencies for the Heisenberg chain w=1 Liouvillian sector
///   Tight-binding hopping 2·J
/// </code>
///
/// <para>F2 is the registry-formal entry for the w=1 Liouvillian dispersion.
/// The IDENTICAL closed form is exposed through
/// <see cref="W1Dispersion"/> as a typed Spectrum primitive (Tier1Derived in
/// D10). F2 inherits the dispersion from W1Dispersion and adds the Pi2-Foundation
/// anchoring readings.</para>
///
/// <para><b>Distinction from F2b:</b> F2 is the w=1 LIOUVILLIAN sector (Pauli
/// strings with exactly one X or Y factor) for the HEISENBERG Hamiltonian,
/// dimension N−1, argument π·k/N. F2b is the single-excitation HAMILTONIAN
/// sector for the XY Hamiltonian, dimension N, argument π·k/(N+1). Different
/// mathematical objects in different Hamiltonians.</para>
///
/// <para>Three independent validations of F2: (1) eigenvalue match &lt; 10⁻¹²,
/// (2) Poisson spacing in w=1 sector (RMT signature), (3) SFF modulation peak
/// at ω_1 matches to &lt;1% for N=2-4, 6.</para>
///
/// <para>Pi2-Foundation anchors:</para>
/// <list type="bullet">
///   <item><b>BandwidthPrefactor = 4 = a_{−1}</b>: in ω_k = 4·J·(1−cos(...)).
///         Live from <see cref="Pi2DyadicLadderClaim.Term"/>(−1). Same anchor
///         as F25 decay rate (e^{−4γt}), F65/F76 numerators, F73 spatial-sum
///         closure, F23 4^N denominator base.</item>
///   <item><b>HoppingFactor = 2 = a_0</b>: in tight-binding hopping 2·J.
///         Live from <see cref="Pi2DyadicLadderClaim.Term"/>(0). Same anchor
///         as F1 TwoFactor, F50 DecayRateFactor, F44 SumCoefficient.</item>
/// </list>
///
/// <para>The (N−1) mode count and the π/N argument come from the (N−1)-
/// dimensional reduction of the w=1 Liouvillian sector under the Heisenberg
/// commutator action (D10 step 4): the zero-frequency stationary mode is
/// excluded, leaving N−1 oscillation modes.</para>
///
/// <para>Tier1Derived: F2 is Tier 1 proven in D10 via reduction to nearest-
/// neighbour tight-binding on (N−1) frequency sites. Pi2-Foundation anchoring
/// is composition through W1Dispersion + Pi2DyadicLadder.</para>
///
/// <para>Anchors: <c>docs/ANALYTICAL_FORMULAS.md</c> F2 (line 32) +
/// <c>experiments/ANALYTICAL_SPECTRUM.md</c> +
/// <c>docs/proofs/derivations/D10_W1_DISPERSION.md</c> +
/// <c>compute/RCPsiSquared.Core/Spectrum/W1Dispersion.cs</c> (typed primitive
/// with same closed form) +
/// <c>compute/RCPsiSquared.Core/Symmetry/Pi2DyadicLadderClaim.cs</c>.</para></summary>
public sealed class F2W1DispersionPi2Inheritance : Claim
{
    private readonly Pi2DyadicLadderClaim _ladder;
    private readonly W1Dispersion _w1;

    /// <summary>The "4" prefactor in ω_k = 4J·(1−cos(πk/N)). Live from
    /// Pi2DyadicLadder a_{−1}. Same anchor as F25/F65/F73/F76.</summary>
    public double BandwidthPrefactor => _ladder.Term(-1);

    /// <summary>The "2" in tight-binding hopping 2J. Live from Pi2DyadicLadder
    /// a_0. Same anchor as F1 TwoFactor, F50 DecayRateFactor.</summary>
    public double HoppingFactor => _ladder.Term(0);

    /// <summary>Live closed form: ω_k = 4·J·(1 − cos(π·k/N)) for k = 1..N−1.</summary>
    public double Frequency(int N, double J, int k)
    {
        if (N < 2) throw new ArgumentOutOfRangeException(nameof(N), N, "F2 requires N ≥ 2.");
        if (k < 1 || k > N - 1) throw new ArgumentOutOfRangeException(nameof(k), k, $"k must be in [1, {N - 1}]; got {k}.");
        if (J < 0) throw new ArgumentOutOfRangeException(nameof(J), J, "J must be ≥ 0.");
        return BandwidthPrefactor * J * (1.0 - Math.Cos(Math.PI * k / N));
    }

    /// <summary>The mode count for F2's dispersion: N − 1.</summary>
    public int ModeCount(int N)
    {
        if (N < 2) throw new ArgumentOutOfRangeException(nameof(N), N, "F2 requires N ≥ 2.");
        return N - 1;
    }

    /// <summary>Drift check: F2's Frequency(N, J, k) matches W1Dispersion's
    /// stored Frequencies[k − 1] for the parent W1Dispersion's (N, J).</summary>
    public bool MatchesW1DispersionParent(double tolerance = 1e-12)
    {
        for (int k = 1; k <= _w1.N - 1; k++)
        {
            double f = Frequency(_w1.N, _w1.J, k);
            double w1f = _w1.Frequencies[k - 1];
            if (Math.Abs(f - w1f) > tolerance) return false;
        }
        return true;
    }

    public F2W1DispersionPi2Inheritance(Pi2DyadicLadderClaim ladder, W1Dispersion w1)
        : base("F2 w=1 Liouvillian dispersion ω_k = 4J·(1−cos(πk/N)), k=1..N−1; 4 = a_{-1}, 2J hopping = a_0; N−1 distinct frequencies",
               Tier.Tier1Derived,
               "docs/ANALYTICAL_FORMULAS.md F2 + " +
               "experiments/ANALYTICAL_SPECTRUM.md + " +
               "docs/proofs/derivations/D10_W1_DISPERSION.md + " +
               "compute/RCPsiSquared.Core/Spectrum/W1Dispersion.cs + " +
               "compute/RCPsiSquared.Core/Symmetry/Pi2DyadicLadderClaim.cs")
    {
        _ladder = ladder ?? throw new ArgumentNullException(nameof(ladder));
        _w1 = w1 ?? throw new ArgumentNullException(nameof(w1));
    }

    public override string DisplayName =>
        "F2 w=1 dispersion as Pi2-Foundation a_{-1} + a_0 + W1Dispersion inheritance";

    public override string Summary =>
        $"ω_k = 4J·(1−cos(πk/N)), k=1..N−1; 4 = a_{{-1}} (= {BandwidthPrefactor}); 2J hopping = a_0 (= {HoppingFactor}); N−1 distinct frequencies; Heisenberg w=1 LIOUVILLIAN sector ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("F2 closed form",
                summary: "ω_k = 4J·(1−cos(πk/N)) for k=1..N−1; Heisenberg w=1 Liouvillian sector; tight-binding hopping 2J on (N−1) frequency sites; D10 derivation");
            yield return InspectableNode.RealScalar("BandwidthPrefactor (= a_{-1} = 4)", BandwidthPrefactor);
            yield return InspectableNode.RealScalar("HoppingFactor (= a_0 = 2)", HoppingFactor);
            yield return new InspectableNode("W1Dispersion source-claim edge",
                summary: $"F2 IS W1Dispersion's content at the F-formula registry level. W1Dispersion at (N={_w1.N}, J={_w1.J}, γ={_w1.GammaZero}) gives ω_1 = {_w1.Frequencies[0]:G6}; F2.Frequency(N, J, 1) = {Frequency(_w1.N, _w1.J, 1):G6}; drift check: MatchesW1DispersionParent = {MatchesW1DispersionParent()}");
            yield return new InspectableNode("F2 vs F2b distinction",
                summary: "F2: w=1 LIOUVILLIAN, Heisenberg, dim N−1, argument πk/N. F2b: single-excitation HAMILTONIAN, XY, dim N, argument πk/(N+1). Different operators, different sectors, different boundary conditions.");
            yield return new InspectableNode("three independent validations",
                summary: "(1) eigenvalue match < 10⁻¹², (2) Poisson spacing in w=1 (RMT), (3) SFF modulation peak at ω_1 < 1% N=2..4, 6");
            yield return new InspectableNode("verified frequencies",
                summary: $"N=4: ω_1 = {Frequency(4, 1.0, 1):G6}, ω_2 = {Frequency(4, 1.0, 2):G6}, ω_3 = {Frequency(4, 1.0, 3):G6}; mode count = {ModeCount(4)}");
        }
    }
}
