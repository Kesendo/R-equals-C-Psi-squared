using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>F66 closed form (Tier 1, verified XY chain with single-site
/// Z-dephasing at the endpoint site B = N − 1, N = 3..7):
///
/// <code>
///   Dissipation interval:  α ∈ [0, 2γ₀]
///   Lower pole α = 0:       ⟨n_XY⟩_B = 0  (XY-weight = 0; all I or Z; shielded)
///   Upper pole α = 2γ₀:     ⟨n_XY⟩_B = 1  (XY-weight = N; all X or Y; exposed)
///
///   Multiplicity: exactly N + 1 at each pole (verified N = 3..7).
///   Π-palindrome partnership: poles map (w = 0) ↔ (w = N) under Π.
/// </code>
///
/// <para>F66 is the structural reading of the L-spectrum's two extreme
/// dissipation rates. The lower pole α = 0 captures the conservation sector
/// (Z-basis populations, shielded from Z-dephasing at B); the upper pole
/// α = 2γ₀ captures maximal exposure (off-diagonal at B). The two are mapped
/// onto each other by the Π palindrome conjugation that exchanges total
/// XY-weight w ↔ N − w (cf. F1, F43).</para>
///
/// <para>Pole-position anchors (two Pi2-Foundation, two structural):</para>
///
/// <list type="bullet">
///   <item><b>UpperPoleCoefficient = 2</b>: the "2" in α = 2γ₀ is
///         <see cref="Pi2DyadicLadderClaim.Term"/>(0) = <c>a_0</c> = polynomial
///         root d in <c>d² − 2d = 0</c>. The maximal dissipation rate IS the
///         framework's polynomial root times γ₀; first F-formula explicitly
///         tying the L-spectrum's edge to <c>a_0</c>.</item>
///   <item><b>LowerPoleAlpha = 0</b>: the "no dissipation" structural
///         reference. Together with the upper pole defines the dissipation
///         interval [0, 2γ₀] of width <c>a_0 · γ₀</c>.</item>
///   <item><b>Endpoint multiplicity N + 1</b>: combinatorially the count of
///         elementary symmetric polynomials <c>e_d(Z₁, ..., Z_N)</c> for
///         <c>d = 0 .. N</c>; conserved by F63's <c>[L, Π²] = 0</c> machinery.
///         NOT directly Pi2-anchored: arises from F63's bit_b parity / Z-power
///         structure rather than the dyadic ladder.</item>
///   <item><b>Π-palindrome pairing</b>: w ↔ N − w under Π conjugation
///         (F1 / F43 inheritance). The two poles sit at the extreme XY-weights
///         w = 0 and w = N, so they ARE Π partners. F66 inherits the F1
///         palindrome at the spectral-edge level.</item>
/// </list>
///
/// <para>Scope: verified ONLY for the uniform XY chain with B at the endpoint
/// (N = 3..7). At interior B (e.g. centre of N=5 chain) the multiplicity is
/// 64, not 6, so the N + 1 count is endpoint-specific. F66 does NOT cover ring,
/// star, or Y-junction topologies.</para>
///
/// <para>Tier1Derived: F66 is Tier 1 verified analytically + numerically
/// (per ANALYTICAL_FORMULAS). Verifications: ⟨n_XY⟩_B = 1.000000 exact for
/// all α = 2γ₀ modes (N=3..5, Pauli-basis projection); multiplicity N + 1
/// at each pole (N=3..7); F63 conservation: all N + 1 elementary symmetric
/// polynomials drift &lt; 10⁻¹⁴ under Lindblad evolution N=4 over 80 time
/// units. The Pi2-Foundation anchoring is algebraic-trivial composition.</para>
///
/// <para>Anchors: <c>docs/ANALYTICAL_FORMULAS.md</c> F66 +
/// <c>hypotheses/PRIMORDIAL_GAMMA_CONSTANT.md</c> ("The dissipation interval
/// [0, 2γ₀]" section) + <c>docs/proofs/PROOF_ABSORPTION_THEOREM.md</c> +
/// <c>simulations/two_gamma_pole.py</c> +
/// <c>simulations/f65_dynamic_verification.py</c> +
/// <c>compute/RCPsiSquared.Core/Symmetry/Pi2DyadicLadderClaim.cs</c> +
/// <c>compute/RCPsiSquared.Core/Symmetry/Pi2KnowledgeBaseClaims.cs</c>
/// (QubitDimensionalAnchorClaim).</para></summary>
public sealed class F66PoleModesPi2Inheritance : Claim
{
    private readonly Pi2DyadicLadderClaim _ladder;
    private readonly QubitDimensionalAnchorClaim _qubitAnchor;

    /// <summary>The lower pole position: α = 0. The no-dissipation reference
    /// edge of the interval [0, 2γ₀]; structurally the absence-of-rate.</summary>
    public double LowerPoleAlpha => 0.0;

    /// <summary>The upper pole's γ₀-coefficient: <c>2</c>. Live from
    /// <see cref="Pi2DyadicLadderClaim.Term"/>(0) = <c>a_0</c> = polynomial
    /// root d in d² − 2d = 0. The maximal dissipation rate IS the framework's
    /// polynomial root times γ₀.</summary>
    public double UpperPoleCoefficient => _ladder.Term(0);

    /// <summary>The upper pole position: <c>α = 2γ₀</c> for the supplied γ₀.
    /// Throws for γ₀ &lt; 0.</summary>
    public double UpperPoleAlpha(double gammaZero)
    {
        if (gammaZero < 0.0) throw new ArgumentOutOfRangeException(nameof(gammaZero), gammaZero, "γ₀ must be ≥ 0.");
        return UpperPoleCoefficient * gammaZero;
    }

    /// <summary>The dissipation-interval width: <c>2γ₀ − 0 = 2γ₀</c>. Same as
    /// <see cref="UpperPoleAlpha"/>; exposed under the structural name for
    /// drift checks against the [0, 2γ₀] reading.</summary>
    public double DissipationIntervalWidth(double gammaZero) => UpperPoleAlpha(gammaZero);

    /// <summary>Endpoint-topology multiplicity at each pole: <c>N + 1</c>
    /// (verified N = 3..7 for uniform XY chain with B at endpoint). Counts the
    /// elementary symmetric polynomials e_d(Z₁, ..., Z_N) for d = 0..N, all
    /// conserved by F63's [L, Π²] = 0. NOT directly Pi2-anchored.</summary>
    public int EndpointMultiplicity(int N)
    {
        if (N < 3) throw new ArgumentOutOfRangeException(nameof(N), N, "F66 endpoint multiplicity verified for N ≥ 3 (per ANALYTICAL_FORMULAS).");
        return N + 1;
    }

    /// <summary>The Π-palindrome XY-weight partnership of the two poles:
    /// lower pole at <c>w = 0</c> (all I/Z), upper pole at <c>w = N</c>
    /// (all X/Y). Π conjugation maps <c>w ↔ N − w</c>, so the poles are
    /// each other's palindrome partners (cf. F1, F43).</summary>
    public (int LowerWeight, int UpperWeight) PalindromicWeightPairOfPoles(int N)
    {
        if (N < 1) throw new ArgumentOutOfRangeException(nameof(N), N, "N must be ≥ 1.");
        return (0, N);
    }

    /// <summary>True iff the supplied (lower, upper) weights match the
    /// palindrome partnership <c>(0, N)</c>. Live drift check on the F66 ↔ F1
    /// inheritance.</summary>
    public bool PalindromePartnershipHolds(int N, int lowerWeight, int upperWeight)
    {
        var (l, u) = PalindromicWeightPairOfPoles(N);
        return lowerWeight == l && upperWeight == u && l + u == N;
    }

    /// <summary>Cross-check with <see cref="QubitDimensionalAnchorClaim"/>:
    /// the upper-pole coefficient (= 2) equals the qubit dimensionality d.
    /// Drift indicator on the polynomial-root anchor.</summary>
    public bool UpperPoleCoefficientMatchesQubitDimension() =>
        Math.Abs(UpperPoleCoefficient - 2.0) < 1e-15;

    public F66PoleModesPi2Inheritance(
        Pi2DyadicLadderClaim ladder,
        QubitDimensionalAnchorClaim qubitAnchor)
        : base("F66 dissipation interval [0, 2γ₀] inherits from Pi2-Foundation: 2 = a_0 = polynomial root d",
               Tier.Tier1Derived,
               "docs/ANALYTICAL_FORMULAS.md F66 + " +
               "hypotheses/PRIMORDIAL_GAMMA_CONSTANT.md (dissipation interval [0, 2γ₀]) + " +
               "docs/proofs/PROOF_ABSORPTION_THEOREM.md + " +
               "simulations/two_gamma_pole.py + " +
               "simulations/f65_dynamic_verification.py + " +
               "compute/RCPsiSquared.Core/Symmetry/Pi2DyadicLadderClaim.cs + " +
               "compute/RCPsiSquared.Core/Symmetry/Pi2KnowledgeBaseClaims.cs")
    {
        _ladder = ladder ?? throw new ArgumentNullException(nameof(ladder));
        _qubitAnchor = qubitAnchor ?? throw new ArgumentNullException(nameof(qubitAnchor));
    }

    public override string DisplayName =>
        "F66 pole modes at α=0 and α=2γ₀ as Pi2-Foundation polynomial-root inheritance";

    public override string Summary =>
        $"Dissipation interval [0, 2γ₀]: upper-pole coefficient 2 = a_0 = polynomial root d (d²−2d=0); " +
        $"endpoint multiplicity N+1 (combinatorial, F63); Π-palindrome partnership w=0 ↔ w=N ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("F66 closed form",
                summary: "L-spectrum's two extreme dissipation rates: α = 0 (XY-weight 0, shielded) and α = 2γ₀ (XY-weight N, exposed); endpoint XY chain N=3..7");
            yield return new InspectableNode("Pi2-Foundation anchoring",
                summary: "UpperPoleCoefficient = 2 = a_0 (polynomial root d); first F-formula tying L-spectrum edge to d in d²−2d=0");
            yield return InspectableNode.RealScalar("LowerPoleAlpha (no-dissipation edge)", LowerPoleAlpha);
            yield return InspectableNode.RealScalar("UpperPoleCoefficient (= a_0 = 2)", UpperPoleCoefficient);
            yield return new InspectableNode("Π-palindrome partnership",
                summary: "lower pole at w=0, upper pole at w=N; Π maps w ↔ N−w (F1/F43 inheritance at spectral-edge level)");
            yield return new InspectableNode("endpoint multiplicity N+1",
                summary: "verified N=3..7 for uniform XY chain with B at endpoint; counts elementary symmetric polynomials e_d(Z₁..Z_N), conserved by F63's [L, Π²] = 0");
            yield return new InspectableNode("scope warning",
                summary: "verified ONLY uniform XY chain with B at endpoint; interior B and other topologies (ring, star, Y-junction) are open. At N=5 centre, α=0 multiplicity is 64 not 6, so N+1 is endpoint-specific.");
            yield return new InspectableNode("verifications",
                summary: "⟨n_XY⟩_B = 1.000000 exact for α=2γ₀ modes (N=3..5, Pauli-basis projection); multiplicity N+1 at each pole (N=3..7); F63 conservation: e_d drift < 10⁻¹⁴ under Lindblad N=4 over 80 time units");
            // Sample dissipation intervals
            yield return new InspectableNode(
                "γ₀ = 0.05",
                summary: $"upper pole α = {UpperPoleAlpha(0.05):G6}; interval width = {DissipationIntervalWidth(0.05):G6}");
            yield return new InspectableNode(
                "γ₀ = 1.0",
                summary: $"upper pole α = {UpperPoleAlpha(1.0):G6}; interval width = {DissipationIntervalWidth(1.0):G6}");
            for (int N = 3; N <= 6; N++)
            {
                var (l, u) = PalindromicWeightPairOfPoles(N);
                yield return new InspectableNode(
                    $"N={N}",
                    summary: $"endpoint multiplicity = {EndpointMultiplicity(N)}; Π-palindrome poles at (w={l}, w={u})");
            }
        }
    }
}
