using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>F43 closed form (Tier 1, proven D09):
///
/// <code>
///   K_freq(w, t) = K_freq(N − w, t)        identical SFF for paired sectors
///
///   Special case: XOR sector (w = N): K = 1.000
///   (all eigenvalues degenerate at rate 2·N·γ)
/// </code>
///
/// <para>F43 is the sector-level reading of F1's palindromic symmetry. The
/// Π conjugation maps Pauli-string sector-weight w to N − w; spectral
/// statistics (the sector spectral form factor K_freq) of paired sectors
/// are identical. F1 lifts to four readings of the palindromic structure:</para>
///
/// <list type="bullet">
///   <item><see cref="F44CrooksLikeRateIdentityPi2Inheritance"/>: pair RATIO
///         at the eigenvalue level (artanh closed form).</item>
///   <item><see cref="F68PalindromicPartnerPi2Inheritance"/>: pair SUM at the
///         eigenvalue level (α_b + α_p = 2γ₀).</item>
///   <item><see cref="F41PalindromicTimePi2Inheritance"/>: pair DIFFERENCE in
///         the time domain (period t_Pi).</item>
///   <item><b>F43 (this claim)</b>: SECTOR-LEVEL identity (K_freq(w) =
///         K_freq(N−w)): the palindromic mirror at the Pauli-string-weight
///         partition.</item>
/// </list>
///
/// <para><b>Mirror axis structure:</b> the pairing w ↔ N−w has its mirror axis
/// at w = N/2. For odd N the axis is at a half-integer (no Pauli string sits
/// on it); for even N the axis is at an integer (Pauli strings ON the axis
/// exist). This matches <see cref="HalfIntegerMirrorClaim"/>'s regime
/// classification.</para>
///
/// <para><b>XOR sector specialization (w = N):</b> the partner is w = 0
/// (identity sector). All Pauli strings in the XOR sector have all sites
/// carrying X or Y (full bit_a = 1 strings). Each gets the same dissipator
/// hit 2·N·γ, so all eigenvalues in the sector are degenerate at this rate;
/// hence K_freq = 1 (delta-spike at zero frequency, normalised to 1).</para>
///
/// <para>Pi2-Foundation anchor:</para>
/// <list type="bullet">
///   <item><b>XorRateCoefficient = 2 = a_0</b>: in the XOR-sector decay rate
///         2·N·γ. Live from <see cref="Pi2DyadicLadderClaim.Term"/>(0). Same
///         anchor as F1's TwoFactor, F50's DecayRateFactor, F66's
///         UpperPoleCoefficient.</item>
/// </list>
///
/// <para>F43 itself is sector-pairing-symmetric (K_freq(w) = K_freq(N−w))
/// without a numerical "2" anchor; that lives in the XOR-sector specialization
/// above. The sector-pairing identity is purely structural (Π palindrome at
/// sector level), inherited from F1.</para>
///
/// <para>Tier1Derived: F43 is Tier 1 proven in D09 (one-line consequence of
/// Π's commutation with the sector projector at the spectral level); valid
/// for Heisenberg chain, Z-dephasing, all N. Pi2-Foundation anchoring is
/// algebraic-trivial composition through F1 + Pi2DyadicLadder.</para>
///
/// <para>Anchors: <c>docs/ANALYTICAL_FORMULAS.md</c> F43 (line 867) +
/// <c>experiments/SPECTRAL_FORM_FACTOR.md</c> +
/// <c>docs/proofs/derivations/D09</c> +
/// <c>compute/RCPsiSquared.Core/Symmetry/F1Pi2Inheritance.cs</c> (palindrome
/// identity at sector-decomposition level) +
/// <c>compute/RCPsiSquared.Core/Symmetry/Pi2DyadicLadderClaim.cs</c>.</para></summary>
public sealed class F43SectorSffPairingPi2Inheritance : Claim
{
    private readonly Pi2DyadicLadderClaim _ladder;
    private readonly F1Pi2Inheritance _f1;

    /// <summary>The "2" coefficient in the XOR-sector decay rate 2·N·γ. Live from
    /// Pi2DyadicLadder a_0. Same anchor as F1's TwoFactor.</summary>
    public double XorRateCoefficient => _ladder.Term(0);

    /// <summary>The mirror partner of sector weight w under Π: <c>N − w</c>.</summary>
    public int PartnerSector(int w, int N)
    {
        if (N < 1) throw new ArgumentOutOfRangeException(nameof(N), N, "F43 requires N ≥ 1.");
        if (w < 0 || w > N) throw new ArgumentOutOfRangeException(nameof(w), w, $"w must be in [0, {N}]; got {w}.");
        return N - w;
    }

    /// <summary>True iff sector w is its own mirror partner (w = N − w ⟺ w = N/2).
    /// Self-paired at the mirror axis only when N is even (then w = N/2 is integer);
    /// for odd N no sector is self-paired.</summary>
    public bool IsSelfPaired(int w, int N)
    {
        return w == PartnerSector(w, N);
    }

    /// <summary>The mirror axis position w = N/2. For even N this is an integer
    /// (some Pauli strings sit ON the axis); for odd N this is a half-integer
    /// (no Pauli string is at the axis).</summary>
    public double MirrorAxis(int N)
    {
        if (N < 1) throw new ArgumentOutOfRangeException(nameof(N), N, "F43 requires N ≥ 1.");
        return N / 2.0;
    }

    /// <summary>True iff the mirror axis at w = N/2 is half-integer (odd N).
    /// Matches <see cref="HalfIntegerMirrorClaim.IsHalfIntegerRegime"/>.</summary>
    public bool IsHalfIntegerMirrorRegime(int N)
    {
        if (N < 1) throw new ArgumentOutOfRangeException(nameof(N), N, "F43 requires N ≥ 1.");
        return N % 2 == 1;
    }

    /// <summary>The XOR-sector decay rate <c>2·N·γ</c>: all Pauli strings in the
    /// w = N sector get the same dissipator hit (all sites carry X or Y).
    /// Per-site rate 2γ from Z-dephasing (per F50) times N sites.</summary>
    public double XorSectorRate(int N, double gammaZero)
    {
        if (N < 1) throw new ArgumentOutOfRangeException(nameof(N), N, "F43 requires N ≥ 1.");
        if (gammaZero < 0) throw new ArgumentOutOfRangeException(nameof(gammaZero), gammaZero, "γ₀ must be ≥ 0.");
        return XorRateCoefficient * N * gammaZero;
    }

    /// <summary>The XOR sector SFF value: K_freq(w = N, t) = 1 exactly (delta-spike
    /// from full eigenvalue degeneracy at the per-sector decay rate). Time-independent.</summary>
    public const double XorSectorSffValue = 1.0;

    /// <summary>Drift check: w + PartnerSector(w, N) = N for any valid w.</summary>
    public bool PartnerSumEqualsN(int w, int N)
    {
        return w + PartnerSector(w, N) == N;
    }

    public F43SectorSffPairingPi2Inheritance(
        Pi2DyadicLadderClaim ladder,
        F1Pi2Inheritance f1)
        : base("F43 sector SFF pairing K_freq(w) = K_freq(N−w); palindromic mirror at sector-weight level; XOR sector (w=N) at rate 2·N·γ with K=1; both '2's = a_0",
               Tier.Tier1Derived,
               "docs/ANALYTICAL_FORMULAS.md F43 + " +
               "experiments/SPECTRAL_FORM_FACTOR.md + " +
               "docs/proofs/derivations/D09 + " +
               "compute/RCPsiSquared.Core/Symmetry/F1Pi2Inheritance.cs + " +
               "compute/RCPsiSquared.Core/Symmetry/Pi2DyadicLadderClaim.cs")
    {
        _ladder = ladder ?? throw new ArgumentNullException(nameof(ladder));
        _f1 = f1 ?? throw new ArgumentNullException(nameof(f1));
    }

    public override string DisplayName =>
        "F43 sector SFF pairing as Pi2-Foundation a_0 + F1 (sector-level palindrome) inheritance";

    public override string Summary =>
        $"K_freq(w) = K_freq(N−w) sector pairing under Π; XOR sector (w=N) has K=1, decay rate 2·N·γ (= a_0·N·γ); F1 palindrome at sector-decomposition level ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("F43 closed form",
                summary: "K_freq(w, t) = K_freq(N−w, t); palindromic Π maps sector w → N−w; identical spectral statistics for paired sectors; valid all N, Heisenberg chain, Z-dephasing");
            yield return InspectableNode.RealScalar("XorRateCoefficient (= a_0 = 2)", XorRateCoefficient);
            yield return new InspectableNode("F1 palindrome at sector level",
                summary: $"F43's pairing w ↔ N−w is the F1 Π conjugation acting at the Pauli-string-sector level. F1's TwoFactor (= {_f1.TwoFactor}) is the same '2' as the XOR-sector decay coefficient.");
            yield return new InspectableNode("XOR sector specialization",
                summary: $"w = N (all sites X or Y): K_freq = {XorSectorSffValue} (delta-spike); all eigenvalues degenerate at rate 2·N·γ; partner sector w = 0 (identity sector) trivially has K = 1");
            yield return new InspectableNode("mirror axis at N/2",
                summary: "axis at w = N/2; even N: integer mirror (Pauli strings ON axis exist); odd N: half-integer mirror (no Pauli string on axis); matches HalfIntegerMirrorClaim regime");
            yield return new InspectableNode("F1 family of palindromic readings",
                summary: "F44 (eigenvalue ratio, artanh) + F68 (eigenvalue sum, 2γ₀) + F41 (eigenvalue difference period, t_Pi) + F43 (sector-weight pairing): four typed readings of one F1 palindromic identity");
            yield return new InspectableNode("N=3 verified",
                summary: $"Pairs: (0,3), (1,2). MirrorAxis = {MirrorAxis(3)} (half-integer, odd N). XOR rate 2·3·γ; PartnerSector(1, 3) = {PartnerSector(1, 3)}");
            yield return new InspectableNode("N=4 verified",
                summary: $"Pairs: (0,4), (1,3); self-paired w=2 = N/2 (integer mirror). MirrorAxis = {MirrorAxis(4)}; PartnerSector(2, 4) = {PartnerSector(2, 4)} (= 2, self-paired)");
        }
    }
}
