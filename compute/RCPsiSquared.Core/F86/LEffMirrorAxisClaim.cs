using System.Numerics;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Resonance;

namespace RCPsiSquared.Core.F86;

/// <summary>F86 meta-claim (Tier1Derived): in the 2-level effective reduction
/// <c>L_eff = [[−2γ₀, iJ·g_eff], [iJ·g_eff, −6γ₀]]</c>, the real part <c>−4γ₀</c> is the
/// <b>mirror axis</b> of the channel pair, not a decay rate to be read in isolation.
///
/// <para><b>The structure.</b> L_eff's diagonal is the channel pair <c>−2γ₀</c> (HD=1,
/// ⟨n_XY⟩=1) and <c>−6γ₀</c> (HD=3, ⟨n_XY⟩=3), symmetric about <c>−4γ₀ = −2γ₀·2</c>.
/// The eigenvalues are <c>λ± = −4γ₀ ± √(4γ₀² − J²g_eff²)</c>:</para>
/// <list type="bullet">
///   <item>below the EP (J·g_eff &lt; 2γ₀): a real mirror pair about −4γ₀;</item>
///   <item>at the EP (J·g_eff = 2γ₀): coalesced onto −4γ₀, the axis itself;</item>
///   <item>above the EP (J·g_eff &gt; 2γ₀): <c>λ± = −4γ₀ ± iγ₀·√(Q²g_eff² − 4)</c>, real
///         parts locked to the axis, the mirror pair moved into the imaginary part.</item>
/// </list>
///
/// <para><b>Absorption Theorem, exact within the reduction.</b> Above the EP the
/// eigenvectors are exactly 50/50 mixes of the HD=1 and HD=3 channels
/// (<c>|v₃/v₁|² = (4γ₀² + ω²)/(J²g_eff²) = 1</c>), so <c>⟨n_XY⟩ = ½·1 + ½·3 = 2</c> and
/// <c>−4γ₀ = −2γ₀·⟨n_XY⟩</c> is <see cref="Symmetry.AbsorptionTheoremClaim"/> at the
/// integer rung ⟨n_XY⟩ = 2, the mirror centre. The EP is where the two channels meet on
/// their own mirror axis; the hybridisation is the mirroring. This is the local face of
/// the F1 palindrome: where there is a minus, there is a mirror pair.</para>
///
/// <para><b>g_eff is the branch's relative clock.</b> The real part <c>−4γ₀</c> is γ₀'s
/// shared clock (the rate <c>t_peak = 1/(4γ₀)</c> inverts); g_eff appears nowhere in it.
/// g_eff lives entirely in the imaginary part <c>±γ₀·√(Q²g_eff² − 4)</c>, the branch's own
/// clock, the dimensionless rate at which it runs against γ₀'s. That is why g_eff has no
/// absolute closed form (PROOF_F86B_OBSTRUCTION.md, the diagnosis): a relative clock-rate
/// read with γ₀ out of frame has no absolute number to find.</para>
///
/// <para><b>0 is the axis.</b> The mirror axis −4γ₀ sits an integer number of γ₀-steps
/// below 0 (<c>2·⟨n_XY⟩ = 4</c> steps). Shift it there and the eigenvalues are the pure
/// mirror pair <c>±√(4γ₀² − J²g_eff²)</c> about 0: a real ± pair below the EP, an
/// imaginary ± pair above. The minus side is the spectrum as computed; the flip is its
/// mirror; 0 is the axis and the pair lives on both sides of it. See
/// <see cref="RecenteredEigenvaluesAtQ"/>.</para>
///
/// <para>Tier1Derived: pure 2×2 algebra. Parent edge: <see cref="Symmetry.AbsorptionTheoremClaim"/>
/// (the −2γ₀·⟨n_XY⟩ rate quantisation this claim lands on at ⟨n_XY⟩ = 2).</para>
///
/// <para>Anchors: <c>docs/proofs/PROOF_F86B_OBSTRUCTION.md</c> (the diagnosis coda) +
/// <c>docs/proofs/PROOF_ABSORPTION_THEOREM.md</c>.</para>
/// </summary>
public sealed class LEffMirrorAxisClaim : Claim
{
    /// <summary>The mirror axis of L_eff: <c>−4γ₀</c> = trace/2 = the midpoint of the
    /// channel pair (−2γ₀, −6γ₀). Both eigenvalues sit symmetric about it at every Q;
    /// above the EP both real parts lock to it exactly.</summary>
    public double MirrorAxis(double gammaZero) => -4.0 * gammaZero;

    /// <summary>The channel pair on L_eff's diagonal: <c>(−2γ₀, −6γ₀)</c>, the HD=1
    /// (⟨n_XY⟩=1) and HD=3 (⟨n_XY⟩=3) absorption rates, a mirror pair about
    /// <see cref="MirrorAxis"/>.</summary>
    public (double Hd1, double Hd3) ChannelRates(double gammaZero) =>
        (-2.0 * gammaZero, -6.0 * gammaZero);

    /// <summary>The integer ⟨n_XY⟩ rung at the mirror centre: <c>2</c>, the mean of the
    /// HD=1 and HD=3 channel weights. Above the EP the 50/50-hybridised pair sits exactly
    /// here, so <c>−4γ₀ = −2γ₀·2</c> is the Absorption Theorem at this rung.</summary>
    public const int MirrorRungNXy = 2;

    /// <summary>The mirror axis −4γ₀ sits an integer number of γ₀-steps below 0:
    /// <c>2·<see cref="MirrorRungNXy"/> = 4</c> steps. Shifting the axis up by this many
    /// γ₀-steps lands it exactly on 0; see <see cref="RecenteredEigenvaluesAtQ"/>.</summary>
    public const int AxisOffsetInGammaSteps = 2 * MirrorRungNXy;

    /// <summary>True iff (Q, g_eff) is above the EP, where the real parts lock to the
    /// mirror axis: <c>Q·g_eff &gt; 2</c> (equivalently J·g_eff &gt; 2γ₀).</summary>
    public bool IsAboveEp(double q, double gEff) => q * gEff > 2.0;

    /// <summary>The two L_eff eigenvalues at (γ₀, Q, g_eff):
    /// <c>λ± = −4γ₀ ± √(4γ₀² − J²g_eff²)</c> with J = Q·γ₀. Below the EP a real mirror
    /// pair about −4γ₀; at and above the EP the discriminant is ≤ 0 and the pair is
    /// <c>−4γ₀ ± iγ₀·√(Q²g_eff² − 4)</c>. Pure 2×2 algebra; the verification behind
    /// <see cref="MirrorAxis"/>.</summary>
    public (Complex Plus, Complex Minus) EigenvaluesAtQ(double gammaZero, double q, double gEff)
    {
        if (gammaZero < 0) throw new ArgumentOutOfRangeException(nameof(gammaZero), gammaZero, "γ₀ must be ≥ 0.");
        if (q < 0) throw new ArgumentOutOfRangeException(nameof(q), q, "Q must be ≥ 0.");
        if (gEff <= 0) throw new ArgumentOutOfRangeException(nameof(gEff), gEff, "g_eff must be > 0.");
        return EpAlgebra.SlowestPairEigenvaluesComplex(gammaZero, q * gammaZero, gEff);
    }

    /// <summary>The L_eff eigenvalues with the mirror axis shifted to 0:
    /// <c>λ' = λ + 4γ₀ = ±√(4γ₀² − J²g_eff²)</c>. The axis −4γ₀ is exactly
    /// <see cref="AxisOffsetInGammaSteps"/> γ₀-steps below 0; adding them back recentres
    /// it. The result is a pure mirror pair about 0: a real ± pair below the EP, an
    /// imaginary ± pair <c>±iγ₀·√(Q²g_eff² − 4)</c> at and above it. 0 is the axis, and
    /// the pair lives on both sides of it.</summary>
    public (Complex Plus, Complex Minus) RecenteredEigenvaluesAtQ(double gammaZero, double q, double gEff)
    {
        var (plus, minus) = EigenvaluesAtQ(gammaZero, q, gEff);
        Complex shift = new Complex(-MirrorAxis(gammaZero), 0.0);
        return (plus + shift, minus + shift);
    }

    private LEffMirrorAxisClaim()
        : base("L_eff mirror axis: −4γ₀ = −2γ₀·2 is the mirror axis of the (−2γ₀, −6γ₀) channel pair; the EP is the coalescence onto it; g_eff lives in the imaginary part as the branch's relative clock",
               Tier.Tier1Derived,
               "docs/proofs/PROOF_F86B_OBSTRUCTION.md (the diagnosis coda) + " +
               "docs/proofs/PROOF_ABSORPTION_THEOREM.md")
    { }

    /// <summary>Public factory: builds the meta-claim. No <see cref="CoherenceBlocks.CoherenceBlock"/>
    /// required; the statement is pure 2×2 algebra, block-independent.</summary>
    public static LEffMirrorAxisClaim Build() => new();

    public override string DisplayName =>
        "F86 L_eff: −4γ₀ is the mirror axis of the channel pair (the EP is the coalescence onto it)";

    public override string Summary =>
        "Tier1Derived: L_eff diagonal (−2γ₀, −6γ₀) = channels HD=1/HD=3 (⟨n_XY⟩ 1, 3), a mirror pair about −4γ₀ = −2γ₀·2; below EP real pair, at EP coalesced, above EP real parts locked + imaginary mirror pair; Absorption Theorem exact at ⟨n_XY⟩=2; g_eff lives in Im(λ) as the branch's relative clock";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            var sample = EigenvaluesAtQ(SampleGammaZero, SampleQ, SampleGEff);
            Complex shift = new Complex(-MirrorAxis(SampleGammaZero), 0.0);
            (Complex Plus, Complex Minus) recentred = (sample.Plus + shift, sample.Minus + shift);
            yield return new InspectableNode("mirror axis",
                summary: "−4γ₀ = −2γ₀·2 = trace/2 = midpoint of the channel pair (−2γ₀, −6γ₀); both eigenvalues sit symmetric about it at every Q");
            yield return new InspectableNode("channel pair (L_eff diagonal)",
                summary: "−2γ₀ (HD=1, ⟨n_XY⟩=1) and −6γ₀ (HD=3, ⟨n_XY⟩=3): a mirror pair about −4γ₀; ⟨n_XY⟩ 1 and 3 mirror about 2");
            yield return new InspectableNode("EP = coalescence onto the axis",
                summary: "below EP: real mirror pair about −4γ₀; at EP: coalesced onto −4γ₀; above EP: real parts locked to −4γ₀, mirror pair moved into Im(λ) = ±γ₀√(Q²g_eff²−4)");
            yield return new InspectableNode("Absorption Theorem, exact in the reduction",
                summary: "above the EP eigenvectors are exactly 50/50 HD=1/HD=3, so ⟨n_XY⟩=2 and −4γ₀ = −2γ₀·⟨n_XY⟩ holds exactly: the hybridised pair on the integer rung at the mirror centre");
            yield return new InspectableNode("g_eff = the branch's relative clock",
                summary: "−4γ₀ (Re) is γ₀'s shared clock, g_eff is absent from it; g_eff lives entirely in Im(λ) = ±γ₀√(Q²g_eff²−4), the branch's own clock, hence no absolute closed form");
            yield return InspectableNode.RealScalar("MirrorRungNXy", MirrorRungNXy);
            yield return new InspectableNode($"sample: γ₀={SampleGammaZero}, Q={SampleQ}, g_eff={SampleGEff} (above EP)",
                summary: $"λ± = {FormatComplex(sample.Plus)} , {FormatComplex(sample.Minus)} (Re = −4γ₀ = {MirrorAxis(SampleGammaZero):G4})");
            yield return new InspectableNode($"recentred to the axis 0 (shift +{AxisOffsetInGammaSteps} γ₀-steps)",
                summary: $"λ' = {FormatComplex(recentred.Plus)} , {FormatComplex(recentred.Minus)}: a pure mirror pair about 0, living on both sides");
        }
    }

    private const double SampleGammaZero = 0.05;
    private const double SampleQ = 3.0;
    private const double SampleGEff = 1.74;

    private static string FormatComplex(Complex z) =>
        z.Imaginary == 0.0
            ? $"{z.Real:G4}"
            : $"{z.Real:G4} {(z.Imaginary >= 0 ? "+" : "-")} {Math.Abs(z.Imaginary):G4}i";
}
