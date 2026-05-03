using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>The framework's primitive half-anchor 1/2 = 1/d (Tier 1 derived; root claim).
/// At the rawest layer, before any of the higher abstractions camouflage it: d = 2 is the
/// qubit dimension (the only non-zero solution of <c>d² − 2d = 0</c>), so 1/d = 1/2.
///
/// <para>Every recurring 0.5 in the framework descends from this primitive: F81 50/50 split,
/// F83 anti-fraction limit r=0, bilinear apex p·(1−p) maximum, half-integer mirror w_XY = N/2,
/// qubit purity floor, slow-mode Π²_Z balance. They are camouflages at successive abstraction
/// layers; the dimensional anchor is where the trail starts.</para>
///
/// <para>This claim exists so future readings of "0.5 in the framework" can locate the
/// primitive without rediscovering it. The lineage is explicit; the anchor is unmasked.</para>
/// </summary>
public sealed class QubitDimensionalAnchorClaim : Claim
{
    public QubitDimensionalAnchorClaim()
        : base("1/2 = 1/d (qubit dimensional anchor; root of all framework half-anchors)",
               Tier.Tier1Derived,
               "docs/EXCLUSIONS.md:251 (d²−2d=0 with C=1/2) + experiments/ORTHOGONALITY_SELECTION_FAMILY.md:357 (apex synthesis)")
    { }

    public override string DisplayName => "1/2 = 1/d (root anchor; qubit dimensionality)";

    public override string Summary =>
        "d=2 is the only non-zero solution of d²−2d=0; 1/d=1/2 is the primitive half-anchor; every framework 0.5 descends from it";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("layer −1 (root)",
                summary: "d²−2d=0 has solutions {0, 2}; d=2 is qubit dim; 1/d=1/2 is primitive");
            yield return new InspectableNode("layer 0: Pauli-basis normalisation",
                summary: "a_α = Tr(σ_α·ρ)/d = Tr(σ_α·ρ)/2; the 1/2 enters every Pauli-vector decomposition");
            yield return new InspectableNode("layer 0: maximally mixed state",
                summary: "ρ_mm = I/d = I/2; diagonal entries of the most-uncertain qubit state are exactly 1/2");
            yield return new InspectableNode("layer 0: Bloch form",
                summary: "ρ = (I + r·σ)/2 for any single qubit; the 2 in the denominator is d");
            yield return new InspectableNode("layer 1: F1 palindrome shift",
                summary: "Π·L·Π⁻¹ = −L − 2σ·I; the factor 2 carries d=2 dimensional structure");
            yield return new InspectableNode("layer 1: qubit purity floor",
                summary: "C ≥ Ψ²/2 + 1/2 for qubits (CORE_ALGEBRA); concurrence cannot fall below 1/d");
            yield return new InspectableNode("layer 2: F81 50/50 split",
                summary: "‖M_sym‖² = ‖M_anti‖² = ‖M‖²/2 for pure Π²-odd; F-chain's flagship 1/2");
            yield return new InspectableNode("layer 2: F83 anti-fraction r=0 limit",
                summary: "anti-fraction = 1/(2+4r) → 1/2 at r=0; F83 inherits F81's split as a special case");
            yield return new InspectableNode("layer 3: bilinear apex",
                summary: "p·(1−p) maximised at p=1/2 universally (ORTHOGONALITY_SELECTION_FAMILY:357)");
            yield return new InspectableNode("layer 4: half-integer mirror",
                summary: "w_XY = N/2 at odd N: 1/2, 3/2, 5/2, … all share fractional 0.5; no Pauli string on mirror axis");
            yield return new InspectableNode("layer 5: slow-mode Klein apex (Schicht 3)",
                summary: "within active Π²_X axis, Pp/(Pp+Mp) ≈ 1/2 (and Pm/(Pm+Mm) ≈ 1/2); the 1/2 surfaces in dynamics too");
            yield return new InspectableNode("ontological anchor",
                summary: "EXCLUSIONS.md:251 — setting Ψ=0 and C=1/2 in R=C(Ψ+R)² collapses to R(R−2)=0; the qubit IS the framework polynomial at C=1/2");
        }
    }
}

/// <summary>F1 → Π² commutes with L (Tier 1 derived). Squaring the F1 palindrome
/// Π·L·Π⁻¹ = −L − 2σ·I gives Π²·L·Π⁻² = L. Therefore each L-eigenspace lives in one
/// Π²-eigenspace, and L is block-diagonal in the Π²-eigenbasis.</summary>
public sealed class Pi2InvolutionClaim : Claim
{
    public Pi2InvolutionClaim()
        : base("Π² commutes with L (consequence of F1 palindrome)",
               Tier.Tier1Derived,
               "docs/ANALYTICAL_FORMULAS.md F1 + F88")
    { }

    public override string DisplayName => "Π² commutes with L";

    public override string Summary =>
        "Π²·L·Π⁻² = L; L is block-diagonal in the Π²-eigenbasis; each L-eigenmode lives in one Π² eigenspace";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("derivation",
                summary: "F1: Π·L·Π⁻¹ = −L − 2σ·I; squaring gives Π²·L·Π⁻² = L");
            yield return new InspectableNode("operator-level statement",
                summary: "L block-diagonalises in the ±1 eigenspaces of Π²");
            yield return new InspectableNode("eigenvector caveat",
                summary: "for degenerate L-eigenvalues spanning both Π² subspaces, generic eigenvector basis can mix; the eigenSPACE always respects Π²");
        }
    }
}

/// <summary>F88 two-axis Π² decomposition (Tier 1 derived). The 4^N Pauli operator space
/// decomposes into 4 Klein cells under (Π²_Z, Π²_X). The Z-axis counts bit_b parity
/// (Z- and Y-dephasing share this), the X-axis counts bit_a parity (X-dephasing). For
/// 2-body bilinears the cells are: Pp = {XX, YY, ZZ} (truly), Pm = {YZ, ZY} (Π²-even
/// non-truly), Mp = {XY, YX} (Π²-odd subgroup A), Mm = {XZ, ZX} (Π²-odd subgroup B).</summary>
public sealed class KleinFourCellClaim : Claim
{
    public KleinFourCellClaim()
        : base("F88 two-axis Π² decomposition: (Π²_Z, Π²_X) → 4 Klein cells",
               Tier.Tier1Derived,
               "docs/ANALYTICAL_FORMULAS.md F88")
    { }

    public override string DisplayName => "F88 Klein 4-cell decomposition";

    public override string Summary =>
        "Pauli operator space splits into Pp / Pm / Mp / Mm via (Π²_Z, Π²_X). 2-body bilinears: 3 truly in Pp, 2 non-truly Π²-even in Pm, 4 Π²-odd split across Mp + Mm";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("Pp = (+, +)", summary: "truly bilinears: XX, YY, ZZ");
            yield return new InspectableNode("Pm = (+, −)", summary: "Π²-even non-truly: YZ, ZY");
            yield return new InspectableNode("Mp = (−, +)", summary: "Π²-odd subgroup A: XY, YX");
            yield return new InspectableNode("Mm = (−, −)", summary: "Π²-odd subgroup B: XZ, ZX");
            yield return new InspectableNode("F80 universality scope",
                summary: "F80's 4-Π²-odd-case universality spans Mp + Mm (two Klein cells, not one)");
        }
    }
}

/// <summary>Bilinear apex 1/2 (Tier 1 derived). The framework's recurring 0.5 anchor: any
/// bilinear form p·(1−p) in a probability variable is maximised at p = 1/2. Manifestations:
/// F81 50/50 split for pure Π²-odd; F83 anti-fraction at r=0; balanced Π² partition
/// (each subspace = 4^N / 2 strings); slow-mode Π²_Z balance within active Π²_X axis
/// (Schicht 3 finding).</summary>
public sealed class BilinearApexClaim : Claim
{
    public BilinearApexClaim()
        : base("Bilinear apex p = 1/2 (recurring framework half-anchor)",
               Tier.Tier1Derived,
               "experiments/ORTHOGONALITY_SELECTION_FAMILY.md:357 + F81 50/50 + F83 anti-fraction r=0")
    { }

    public override string DisplayName => "Apex 1/2 of any bilinear form p·(1−p)";

    public override string Summary =>
        "1/2 is the universal maximum of any bilinear form in a probability variable; appears as F81 50/50, F83 anti-fraction at r=0, balanced Π² partition, and slow-mode apex within active Π²_X axis";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("F81 50/50",
                summary: "pure Π²-odd: ‖M_sym‖² = ‖M_anti‖² = ‖M‖²/2");
            yield return new InspectableNode("F83 anti-fraction r=0",
                summary: "anti-fraction = 1/(2 + 4r) → 1/2 at r = 0 (pure Π²-odd)");
            yield return new InspectableNode("Π² partition",
                summary: "each Π²_Z eigenspace contains 4^N / 2 Pauli strings; balanced 50/50 by combinatorics");
            yield return new InspectableNode("slow-mode apex (Schicht 3)",
                summary: "within active Π²_X axis, slow-mode Pp/(Pp+Mp) ≈ 1/2 (and Pm/(Pm+Mm) ≈ 1/2)");
        }
    }
}

/// <summary>Half-integer mirror w_XY = N/2 (Tier 1 derived). For odd N the mirror axis is
/// at a half-integer (1/2, 3/2, 5/2, …) and no Pauli string sits on the axis. For even N
/// the mirror is at an integer and 4^N · C(N, N/2) / 2^N strings sit on the axis. The
/// half-integer-mirror family shares the fractional part 1/2: Tom's "alte Bekannte 0.5" as
/// a sequence rather than a single number.</summary>
public sealed class HalfIntegerMirrorClaim : Claim
{
    public int N { get; }
    public double WXY => N / 2.0;
    public bool IsHalfIntegerRegime => N % 2 == 1;

    public HalfIntegerMirrorClaim(int N)
        : base($"Half-integer mirror w_XY = N/2 (regime classifier; N = {N})",
               Tier.Tier1Derived,
               "experiments/OPERATOR_RIGIDITY_ACROSS_CUSP.md:34 + review/EMERGING_QUESTIONS.md (EQ-026 stage-fixing)")
    {
        if (N < 1) throw new ArgumentOutOfRangeException(nameof(N), $"N must be ≥ 1; got {N}");
        this.N = N;
    }

    public override string DisplayName =>
        $"w_XY = {WXY} ({(IsHalfIntegerRegime ? "half-integer" : "integer")}-mirror regime)";

    public override string Summary => IsHalfIntegerRegime
        ? $"odd N = {N}: w_XY = {WXY} half-integer; no Pauli strings on mirror axis; balanced 50/50 below/above"
        : $"even N = {N}: w_XY = {WXY:F0} integer; Pauli strings ON mirror axis exist";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("regime",
                summary: IsHalfIntegerRegime ? "half-integer mirror" : "integer mirror");
            yield return new InspectableNode("w_XY",
                summary: WXY.ToString("F1"));
            yield return new InspectableNode("framework anchor",
                summary: "EQ-026 stage-fixing: N=3 introduces half-integer mirror, N=4 fixes 15/46/59 with integer mirror, N=5 inherits half-integer regime");
        }
    }
}
