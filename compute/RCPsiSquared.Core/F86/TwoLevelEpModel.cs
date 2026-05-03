using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Resonance;

namespace RCPsiSquared.Core.F86;

/// <summary>F86 Statement 1: the 2-level effective Liouvillian for the channel pair
/// (HD = 2k−1, HD = 2k+1), parametrised by the level index <see cref="K"/>. k=1 is the
/// slowest pair (the one F86 Statement 1's t_peak refers to); higher k decay faster.
///
/// <para>Eigenvalues: <c>λ_±(k) = −4γ₀·k ± √(4γ₀² − J²·g_eff²)</c>. The discriminant is
/// k-independent — the EP location <c>J·g_eff = 2γ₀</c> is the same for every level pair —
/// but the trace centre <c>−4γ₀·k</c> shifts to faster decay with k.</para>
///
/// <para>The decay time at the EP for level pair k is <c>1/(4γ₀·k)</c>: 1/(4γ₀) for the
/// slowest k=1, 1/(8γ₀) for k=2, 1/(12γ₀) for k=3, etc. The slowest k=1 dominates the
/// long-time response — this is the universality basis of <see cref="TPeakLaw"/>.</para>
///
/// <para>Algebraic note: the same-sign-imaginary off-diagonal form gives an EP at finite
/// J·g_eff = 2γ₀. This is "PT-like phenomenology" but algebraically inside class AIII chiral
/// (Π linear; classical PT anti-linear) — see <c>experiments/PT_SYMMETRY_ANALYSIS.md</c>
/// and <c>hypotheses/FRAGILE_BRIDGE.md</c>.</para>
/// </summary>
public sealed class TwoLevelEpModel : F86Claim
{
    public int K { get; }
    public double GammaZero { get; }
    public double J { get; }
    public double GEff { get; }

    public double Discriminant => 4.0 * GammaZero * GammaZero - J * J * GEff * GEff;
    public double TraceCentre => -4.0 * GammaZero * K;
    public (double LamPlus, double LamMinus) Eigenvalues { get; }
    public EpRegime Regime { get; }

    /// <summary>Decay time at the EP for this level pair: <c>1/(4γ₀·k)</c>.</summary>
    public double DecayTimeAtEp => 1.0 / (4.0 * GammaZero * K);

    /// <summary>|J·g_eff − 2γ₀| / (2γ₀) — relative distance to the EP, zero at coalescence.</summary>
    public double RelativeDistanceToEp { get; }

    /// <summary>Eigenvector rotation parameter <c>τ² = (Q − Q_EP)/(Q + Q_EP)</c> (real pre-EP,
    /// switches to phase parameterisation post-EP). Probe overlap depends only on this ratio.</summary>
    public double TauSquared
    {
        get
        {
            double sum = J * GEff + 2 * GammaZero;
            return sum == 0 ? 0 : (J * GEff - 2 * GammaZero) / sum;
        }
    }

    public TwoLevelEpModel(double gammaZero, double j, double gEff, int k = 1)
        : base($"2-level EP model (k={k})", Tier.Tier1Derived, "docs/proofs/PROOF_F86_QPEAK.md Statement 1")
    {
        if (k < 1) throw new ArgumentOutOfRangeException(nameof(k), $"k must be ≥ 1; got {k}");
        K = k;
        GammaZero = gammaZero;
        J = j;
        GEff = gEff;
        // Slowest-pair eigenvalues are around centre -4γ₀; level-k eigenvalues are around -4γ₀·k.
        var slowest = EpAlgebra.SlowestPairEigenvalues(gammaZero, j, gEff);
        // Shift by -(k-1)·4γ₀ to land at level k's trace centre.
        double shift = -4.0 * gammaZero * (k - 1);
        Eigenvalues = (slowest.LamPlus + shift, slowest.LamMinus + shift);
        const double epTol = 1e-10;
        if (Math.Abs(Discriminant) < epTol) Regime = EpRegime.AtEp;
        else if (Discriminant > 0) Regime = EpRegime.PreEp;
        else Regime = EpRegime.PostEp;
        RelativeDistanceToEp = Math.Abs(j * gEff - 2 * gammaZero) / (2 * gammaZero);
    }

    /// <summary>Build at a chosen Q (uniform-J): <c>J = Q·γ₀</c>.</summary>
    public static TwoLevelEpModel AtQ(double gammaZero, double q, double gEff, int k = 1) =>
        new(gammaZero, q * gammaZero, gEff, k);

    public override string DisplayName => $"2-level EP k={K} @ J·g_eff/(2γ₀) = {(J * GEff) / (2 * GammaZero):F4}";
    public override string Summary =>
        $"regime = {Regime}, λ₊ = {Eigenvalues.LamPlus:F4}, λ₋ = {Eigenvalues.LamMinus:F4}, " +
        $"t_decay@EP = {DecayTimeAtEp:G4}, discriminant = {Discriminant:E2}";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return InspectableNode.RealScalar("k (level pair index)", K);
            yield return InspectableNode.RealScalar("γ₀", GammaZero, "G6");
            yield return InspectableNode.RealScalar("J", J, "G6");
            yield return InspectableNode.RealScalar("g_eff", GEff, "G6");
            yield return new InspectableNode("regime", summary: Regime.ToString());
            yield return InspectableNode.RealScalar("discriminant", Discriminant, "E3");
            yield return InspectableNode.RealScalar("trace centre (-4γ₀·k)", TraceCentre, "G6");
            yield return InspectableNode.RealScalar("decay time at EP (1/(4γ₀·k))", DecayTimeAtEp, "G6");
            yield return InspectableNode.RealScalar("λ₊", Eigenvalues.LamPlus, "G6");
            yield return InspectableNode.RealScalar("λ₋", Eigenvalues.LamMinus, "G6");
            yield return InspectableNode.RealScalar("relative distance to EP", RelativeDistanceToEp, "G6");
            yield return InspectableNode.RealScalar("τ² (eigenvector rotation)", TauSquared, "G6");
        }
    }
}

public enum EpRegime
{
    PreEp,
    AtEp,
    PostEp,
}
