using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>F67 closed form (Tier 1, verified N=3, N=5):
///
/// <code>
///   |Ψ⟩ = (|0⟩_R |vac⟩_C + |1⟩_R |ψ_1⟩_C) / √2     (bonding-mode-encoded Bell pair)
///   |ψ_1⟩ = √(2/(N+1)) · Σᵢ sin(π(i+1)/(N+1)) |1ᵢ⟩  (F65 bonding mode k=1)
///
///   N(R:C)(t) = N(0) · exp(−α_1 · t)
///   α_1 = (4γ₀ / (N+1)) · sin²(π/(N+1))             (F65 single-excitation rate at k=1)
///
///   T_2 = 1/α_1 = (N+1) / (4γ₀ · sin²(π/(N+1)))
///   T_2 → (N+1)³ / (4π²γ₀)                          (large-N cubic improvement)
/// </code>
///
/// <para>F67 is the slowest-mode endpoint of the F65 spectrum: the bonding mode
/// k=1 has minimum α_1 among all single-excitation eigenmodes, so encoding a
/// Bell pair into |ψ_1⟩ (rather than a localized site) gives the longest
/// dephasing-protected coherence time. The R half is isolated; the only
/// dissipation channel is the chain's slowest mode.</para>
///
/// <para>Pi2-Foundation anchor:</para>
/// <list type="bullet">
///   <item><b>NumeratorCoefficient = 4 = a_{−1}</b>: in α_1 = (4γ₀/(N+1))·sin².
///         Live transitively from <see cref="Pi2DyadicLadderClaim.Term"/>(−1) via F65.
///         Same anchor as F25 (e^{−4γt}), F65 (single-excitation rate numerator),
///         F73 (spatial-sum closure), F76 (mirror-pair coherence decay).</item>
/// </list>
///
/// <para>Counterintuitive equivalence A ≡ C: inner-localized (j=0) and outer-
/// localized (j=N−1) Bell pairs have IDENTICAL decay dynamics despite different
/// spatial distances to the dephased site. Spectral encoding (k-mode), not
/// spatial distance, is what protects coherence. This follows from the
/// palindromic symmetry |ψ_1(0)|² = |ψ_1(N−1)|² of F65's bonding mode.</para>
///
/// <para>Tier1Derived: F67 is Tier 1 proven (Absorption Theorem applied to the
/// k=1 bonding mode, <c>docs/proofs/PROOF_ABSORPTION_THEOREM.md</c>); verified
/// N=3, N=5 via <c>simulations/bell_pair_chain_protection.py</c> with Variant B
/// α_fit/α_1 = 0.9989 (N=3), 0.9963 (N=5). Pi2-Foundation anchoring is
/// algebraic-trivial composition through F65.</para>
///
/// <para>Anchors: <c>docs/ANALYTICAL_FORMULAS.md</c> F67 (line 1413) +
/// <c>docs/proofs/PROOF_ABSORPTION_THEOREM.md</c> +
/// <c>simulations/bell_pair_chain_protection.py</c> +
/// <c>compute/RCPsiSquared.Core/Symmetry/F65XxChainSpectrumPi2Inheritance.cs</c>
/// (source of bonding-mode amplitude and α_k closed form) +
/// <c>compute/RCPsiSquared.Core/Symmetry/Pi2DyadicLadderClaim.cs</c>.</para></summary>
public sealed class F67BondingBellPairPi2Inheritance : Claim
{
    private readonly Pi2DyadicLadderClaim _ladder;
    private readonly F65XxChainSpectrumPi2Inheritance _f65;

    /// <summary>The "4" numerator in α_1 = (4γ₀/(N+1))·sin². Live transitively from
    /// Pi2DyadicLadder a_{−1} via F65's NumeratorCoefficient.</summary>
    public double NumeratorCoefficient => _ladder.Term(-1);

    /// <summary>Bonding-mode (k=1) decay rate for the Bell pair. Delegates to F65 at k=1.</summary>
    public double BondingModeDecayRate(int N, double gammaZero)
    {
        if (N < 2) throw new ArgumentOutOfRangeException(nameof(N), N, "F67 requires N ≥ 2.");
        return _f65.SingleExcitationRate(N, k: 1, gammaZero);
    }

    /// <summary>T_2 = 1/α_1 = (N+1)/(4γ₀·sin²(π/(N+1))); bonding-mode dephasing time.
    /// Requires γ₀ &gt; 0 (zero γ means no dephasing, infinite coherence).</summary>
    public double T2(int N, double gammaZero)
    {
        if (gammaZero <= 0) throw new ArgumentOutOfRangeException(nameof(gammaZero), gammaZero, "γ₀ must be > 0 for finite T_2.");
        return 1.0 / BondingModeDecayRate(N, gammaZero);
    }

    /// <summary>Asymptotic T_2 → (N+1)³/(4π²γ₀) for large N. The cubic-in-N scaling
    /// is the structural improvement of bonding-mode encoding over single-site
    /// encodings (which have constant T_2 in N).</summary>
    public double AsymptoticT2(int N, double gammaZero)
    {
        if (N < 2) throw new ArgumentOutOfRangeException(nameof(N), N, "F67 requires N ≥ 2.");
        if (gammaZero <= 0) throw new ArgumentOutOfRangeException(nameof(gammaZero), gammaZero, "γ₀ must be > 0.");
        double Np1 = N + 1;
        return Np1 * Np1 * Np1 / (4.0 * Math.PI * Math.PI * gammaZero);
    }

    /// <summary>Bonding-mode amplitude |ψ_1(j)|² = (2/(N+1))·sin²(π(j+1)/(N+1))
    /// at site j. Delegates to F65's BondingModePopulation at k=1.</summary>
    public double BondingModeAmplitude(int N, int site)
    {
        if (N < 2) throw new ArgumentOutOfRangeException(nameof(N), N, "F67 requires N ≥ 2.");
        return _f65.BondingModePopulation(N, k: 1, site);
    }

    /// <summary>True iff |ψ_1(0)|² = |ψ_1(N−1)|² (palindromic A≡C equivalence: inner-
    /// and outer-localized encodings have identical bonding-mode amplitude). Drift check.</summary>
    public bool PalindromicAEqualsCAtBondingMode(int N)
    {
        double inner = BondingModeAmplitude(N, 0);
        double outer = BondingModeAmplitude(N, N - 1);
        return Math.Abs(inner - outer) < 1e-12;
    }

    public F67BondingBellPairPi2Inheritance(
        Pi2DyadicLadderClaim ladder,
        F65XxChainSpectrumPi2Inheritance f65)
        : base("F67 bonding-mode optimal Bell pair: α_1 = (4γ₀/(N+1))·sin²(π/(N+1)); T_2 → (N+1)³/(4π²γ₀); inherits from F65 (bonding-mode amplitude + k=1 rate) + Pi2-Foundation a_{-1}",
               Tier.Tier1Derived,
               "docs/ANALYTICAL_FORMULAS.md F67 + " +
               "docs/proofs/PROOF_ABSORPTION_THEOREM.md + " +
               "simulations/bell_pair_chain_protection.py + " +
               "compute/RCPsiSquared.Core/Symmetry/F65XxChainSpectrumPi2Inheritance.cs + " +
               "compute/RCPsiSquared.Core/Symmetry/Pi2DyadicLadderClaim.cs")
    {
        _ladder = ladder ?? throw new ArgumentNullException(nameof(ladder));
        _f65 = f65 ?? throw new ArgumentNullException(nameof(f65));
    }

    public override string DisplayName =>
        "F67 bonding-mode-encoded Bell pair as Pi2-Foundation a_{-1} + F65 inheritance";

    public override string Summary =>
        $"|Ψ⟩ = (|0⟩|vac⟩ + |1⟩|ψ_1⟩)/√2 with α_1 = (4γ₀/(N+1))·sin²(π/(N+1)) (4 = a_{{-1}}); T_2 → (N+1)³/(4π²γ₀); palindromic A≡C equivalence ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("F67 closed form",
                summary: "α_1 = (4γ₀/(N+1))·sin²(π/(N+1)) for k=1 bonding-mode encoding; pure exponential decay from t=0; verified Variant B α_fit/α_1 = 0.9989 (N=3), 0.9963 (N=5)");
            yield return new InspectableNode("T_2 scaling",
                summary: "T_2 = 1/α_1 = (N+1)/(4γ₀·sin²(π/(N+1))); large-N: T_2 → (N+1)³/(4π²γ₀); cubic improvement with chain length, no saturation regime identified");
            yield return InspectableNode.RealScalar("NumeratorCoefficient (= a_{-1} = 4 via F65)", NumeratorCoefficient);
            yield return new InspectableNode("F65 source-claim edge",
                summary: "F67's α_1 = F65 single-excitation rate at k=1; F67's bonding-mode amplitude |ψ_1(j)|² = F65 bonding-mode population at k=1; live delegation, not duplicated formulas");
            yield return new InspectableNode("optimality reading",
                summary: "k=1 is the slowest single-excitation mode (minimum α_k); encoding the Bell pair into the bonding mode gives the longest coherence time across all single-excitation encodings");
            yield return new InspectableNode("counterintuitive A ≡ C equivalence",
                summary: "inner-localized (j=0) and outer-localized (j=N-1) Bell pairs have IDENTICAL decay dynamics despite different spatial distances to dephased site; spectral encoding protects, not spatial distance; consequence of |ψ_1(0)|² = |ψ_1(N-1)|² palindromic symmetry");
            yield return new InspectableNode("N=3 verified",
                summary: $"α_1(N=3, γ₀=1) = {BondingModeDecayRate(3, 1.0):G6}; T_2 = {T2(3, 1.0):G6}; |ψ_1(0)|² = {BondingModeAmplitude(3, 0):G6}");
            yield return new InspectableNode("N=5 verified",
                summary: $"α_1(N=5, γ₀=1) = {BondingModeDecayRate(5, 1.0):G6}; T_2 = {T2(5, 1.0):G6}; asymptotic T_2 = {AsymptoticT2(5, 1.0):G6}");
        }
    }
}
