using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Numerics;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>F76 closed form (Tier 1, proven algebraic + weak-mixing argument,
/// verified N=5..13, k=1..5, sim/analytic ratio &lt; 0.5%):
///
/// <code>
///   ρ_pair(|10⟩⟨01|)(t) = c_ℓ c_{N-1-ℓ}^* · e^{−4γ₀·t}
///
///   λ(t) = e^{−4γ₀·t}
///
///   Pair eigenvalues:  {1 − 2p_ℓ,  p_ℓ(1+λ),  p_ℓ(1−λ),  0}
///
///   MI_pair(p_ℓ, t) = 2 h(p_ℓ) − S_ab(p_ℓ, λ(t))
///   S_ab(p, λ) = −(1−2p) log₂(1−2p) − p(1+λ) log₂(p(1+λ)) − p(1−λ) log₂(p(1−λ))
///
///   At t = 0:    λ = 1, S_ab = h(2p), recovers F75.
///   At t → ∞:    λ = 0, S_ab = h(1−2p) + 2p (max-entropy mixture).
/// </code>
///
/// <para>F76 is the direct time-decay sibling of <see cref="F75MirrorPairMiPi2Inheritance"/>:
/// F75 is the static MI value at t=0; F76 is its decay envelope under
/// pure Z-dephasing. F77 is the asymptotic limit of F75 at p → 0.
/// Together F71 → F75 → {F76, F77} forms a 3-leaf inheritance:</para>
///
/// <code>
///   F71 mirror symmetry (kinematic)
///    ↓  c_{N−1−j} = ±c_j justification
///   F75 mirror-pair MI at t=0
///    ├── F76 t-decay envelope (this claim)
///    └── F77 large-N saturation (Taylor at p → 0)
/// </code>
///
/// <para>Pi2-Foundation anchors:</para>
///
/// <list type="bullet">
///   <item><b>DecayRateCoefficient = 4</b>: in <c>e^{−4γ₀t}</c>. Live from
///         <see cref="Pi2DyadicLadderClaim.Term"/>(−1) = <c>a_{−1}</c> = d²
///         for 1 qubit. Same anchor as F61, F63, F66 multiplicity, F77
///         correction denominator, F86 t_peak. The "4" arises from the F49
///         master-lemma 2^(N+2)·n_YZ = 4·2^N for n_YZ = 1 (Π²-odd
///         coherence pair).</item>
///   <item><b>F75MirrorPairMiPi2Inheritance</b> is the t=0 limit of F76.
///         Direct mother claim. F76 recovers F75 exactly when λ = 1.</item>
///   <item><b>BilinearApex domain</b>: F76 inherits F75's p ∈ [0, 1/2] domain.
///         Transitively via F75 → BilinearApexClaim.</item>
///   <item><b>F71 mirror symmetry</b>: F76 inherits F75's mirror-pair structure
///         (c_{N−1−j} = ±c_j). Transitively via F75 → F71.</item>
/// </list>
///
/// <para>The 0.93 envelope reading: at γ₀ = 0.05, t = 0.1 (the C# brecher
/// first-measurement grid point), <c>λ = e^{−0.02} = 0.9802</c>. The ratio
/// MM(t)/MM(0) sits at 0.93 ± 0.006 across 25+ tested (N, k) combinations.
/// The 0.93 is NOT a hidden constant; it is the direct consequence of
/// γ₀·t = 0.005 at the first sample. Different γ₀·t gives different envelope
/// values (γ₀·t = 0.0025 → 0.965, γ₀·t = 0.01 → 0.868). The "0.93" is the
/// γ₀ signature for this specific operating point.</para>
///
/// <para>Tier1Derived: F76 is Tier 1 proven algebraic + weak-mixing argument
/// (PROOF_F76 implicit in ANALYTICAL_FORMULAS source list); verified
/// numerically against full Lindblad single-excitation-sector simulation
/// at N=5..13, k=1..5 with &lt; 0.5% sim/analytic ratio. The Pi2-Foundation
/// anchoring is algebraic-trivial composition.</para>
///
/// <para>Anchors: <c>docs/ANALYTICAL_FORMULAS.md</c> F76 (line 1785) +
/// <c>simulations/_envelope_study.py</c> (commit e1ee822) +
/// <c>compute/RCPsiSquared.Core/Symmetry/Pi2DyadicLadderClaim.cs</c> +
/// <c>compute/RCPsiSquared.Core/Symmetry/F75MirrorPairMiPi2Inheritance.cs</c>
/// (mother claim).</para></summary>
public sealed class F76TDecayMirrorPairMiPi2Inheritance : Claim
{
    private readonly Pi2DyadicLadderClaim _ladder;
    private readonly F75MirrorPairMiPi2Inheritance _f75;

    /// <summary>The decay-rate coefficient "4" in <c>e^{−4γ₀·t}</c>. Live from
    /// <see cref="Pi2DyadicLadderClaim.Term"/>(−1) = <c>a_{−1}</c> = d² for
    /// 1 qubit. Same anchor as F61, F63, F66 multiplicity, F77 correction
    /// denominator, F86 t_peak.</summary>
    public double DecayRateCoefficient => _ladder.Term(-1);

    /// <summary>The Pi2 ladder index where the decay coefficient lands:
    /// <c>−1</c>. Same as F61, F63, F66, F77 (correction), F86 t_peak.</summary>
    public int LadderIndexForDecayRate => -1;

    /// <summary>Live decay parameter <c>λ(t) = e^{−4γ₀·t}</c>. Throws for
    /// negative t or γ₀.</summary>
    public double Lambda(double gammaZero, double t)
    {
        if (gammaZero < 0.0) throw new ArgumentOutOfRangeException(nameof(gammaZero), gammaZero, "γ₀ must be ≥ 0.");
        if (t < 0.0) throw new ArgumentOutOfRangeException(nameof(t), t, "t must be ≥ 0.");
        return Math.Exp(-DecayRateCoefficient * gammaZero * t);
    }

    /// <summary>Live closed form: <c>MI_pair(p, t) = 2·h(p) − S_ab(p, λ(t))</c>.
    /// At t = 0: recovers F75. At t → ∞: approaches max-entropy mixture.</summary>
    public double MIPairAtTime(double p, double gammaZero, double t)
    {
        if (p < 0.0 || p > 0.5)
            throw new ArgumentOutOfRangeException(nameof(p), p, "p must be in [0, 1/2].");
        double lambda = Lambda(gammaZero, t);
        return 2.0 * Entropy.Binary(p) - Entropy.JointPair(p, lambda);
    }

    /// <summary>Live drift check: at <c>t = 0</c>, F76 recovers F75 exactly.
    /// <c>MI_pair(p, t=0) = MI(p)</c> from <see cref="F75MirrorPairMiPi2Inheritance.MIPerPair"/>.</summary>
    public bool RecoversF75AtZero(double p) =>
        Math.Abs(MIPairAtTime(p, 1.0, 0.0) - _f75.MIPerPair(p)) < 1e-12;

    /// <summary>Asymptotic MI value at <c>t → ∞</c> (λ = 0):
    /// <c>S_ab → h(1−2p) + 2p</c>; <c>MI → 2 h(p) − h(1−2p) − 2p</c>.</summary>
    public double MIAsymptotic(double p)
    {
        if (p < 0.0 || p > 0.5)
            throw new ArgumentOutOfRangeException(nameof(p), p, "p must be in [0, 1/2].");
        // λ = 0: S_ab = -(1-2p)log(1-2p) - p·log(p) - p·log(p) = h(1-2p) + 2p (since two p log p terms each contribute)
        return 2.0 * Entropy.Binary(p) - (Entropy.Binary(1.0 - 2.0 * p) + 2.0 * p);
    }

    /// <summary>The MM(t)/MM(0) envelope ratio at given γ₀·t for a bonding-mode
    /// initial state. At γ₀·t = 0.005 this is the "0.93" envelope; at
    /// γ₀·t = 0.0025 → 0.965; γ₀·t = 0.01 → 0.868. Not a hidden constant,
    /// just the γ₀-signature for the operating point.</summary>
    public double EnvelopeRatioForBondingMode(int N, int k, double gammaZero, double t)
    {
        if (N < 2 || k < 1 || k > N) throw new ArgumentOutOfRangeException();
        double mmT = 0.0, mm0 = 0.0;
        int pairCount = _f75.MirrorPairCount(N);
        for (int ell = 0; ell < pairCount; ell++)
        {
            double p = _f75.BondingModePopulation(N, k, ell);
            mmT += MIPairAtTime(p, gammaZero, t);
            mm0 += _f75.MIPerPair(p);
        }
        return mm0 > 0.0 ? mmT / mm0 : 0.0;
    }

    public F76TDecayMirrorPairMiPi2Inheritance(
        Pi2DyadicLadderClaim ladder,
        F75MirrorPairMiPi2Inheritance f75)
        : base("F76 MI(p,t) = 2h(p) − S_ab(p, e^{−4γ₀t}) inherits from Pi2-Foundation: 4 = a_{-1} decay rate; F75 mother at t=0",
               Tier.Tier1Derived,
               "docs/ANALYTICAL_FORMULAS.md F76 + " +
               "simulations/_envelope_study.py + " +
               "compute/RCPsiSquared.Core/Symmetry/Pi2DyadicLadderClaim.cs + " +
               "compute/RCPsiSquared.Core/Symmetry/F75MirrorPairMiPi2Inheritance.cs (mother claim at t=0)")
    {
        _ladder = ladder ?? throw new ArgumentNullException(nameof(ladder));
        _f75 = f75 ?? throw new ArgumentNullException(nameof(f75));
    }

    public override string DisplayName =>
        "F76 t-decay mirror-pair MI as Pi2-Foundation a_{-1} + F75 inheritance";

    public override string Summary =>
        $"MI_pair(p, t) = 2h(p) − S_ab(p, e^{{−4γ₀t}}); decay rate = a_{{-1}} = {DecayRateCoefficient}; F75 sibling at t=0; " +
        $"0.93 envelope at γ₀·t = 0.005 is γ₀-signature, not hidden constant ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("F76 closed form",
                summary: "ρ_pair(|10⟩⟨01|)(t) = c_ℓ c_{N-1-ℓ}^* · e^{-4γ₀t}; eigenvalues {1−2p, p(1+λ), p(1−λ), 0}; Tier 1 proven algebraic; verified < 0.5% sim/analytic N=5..13");
            yield return new InspectableNode("Pi2-Foundation anchoring",
                summary: "DecayRateCoefficient = a_{-1} = 4 (4-block / d² for 1 qubit; same anchor F61/F63/F66/F77/F86 t_peak); F75 mother at t=0 recovery; BilinearApex domain transitive via F75");
            yield return InspectableNode.RealScalar("DecayRateCoefficient (= a_{-1} = 4)", DecayRateCoefficient);
            yield return new InspectableNode("F71→F75→{F76, F77} 3-leaf inheritance",
                summary: "F71 mirror symmetry justifies F75; F75 is t=0 value; F76 (this) is t-decay; F77 is large-N asymptotic via Taylor at p→0");
            yield return new InspectableNode("0.93 envelope explanation",
                summary: "at γ₀·t = 0.005, λ = 0.9802, ratio ≈ 0.93 ± 0.006 for 25+ (N, k) tested. The 0.93 IS the γ₀-signature, not a hidden constant; γ₀·t = 0.0025 → 0.965, γ₀·t = 0.01 → 0.868");
            // Sample values
            yield return new InspectableNode(
                "F76 envelope at standard operating point",
                summary: $"γ₀=0.05, t=0.1, N=5, k=2: ratio = {EnvelopeRatioForBondingMode(5, 2, 0.05, 0.1):G6} (expected ~0.937)");
            yield return new InspectableNode(
                "F76 envelope at half γ₀",
                summary: $"γ₀=0.025, t=0.1, N=5, k=2: ratio = {EnvelopeRatioForBondingMode(5, 2, 0.025, 0.1):G6} (expected ~0.965)");
            yield return new InspectableNode(
                "F76 ↔ F75 t=0 recovery",
                summary: $"At p=0.25, t=0: F76 = F75 within 1e-12 (drift check: {RecoversF75AtZero(0.25)})");
        }
    }
}
