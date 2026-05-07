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

/// <summary>Quarter as bilinear maxval (Tier 1 derived). The maxval-side companion to
/// <see cref="BilinearApexClaim"/>: 1/4 = (1/2)² is the universal maximum value of any
/// bilinear form p·(1−p) on [0,1], attained at the argmax p = 1/2. This is the calculus
/// identity max_{p ∈ [0,1]} p(1−p) = 1/4 at p = 1/2, derivable from d/dp[p − p²] = 1 − 2p.
/// The same maxval surfaces at three independent layers of the framework: the Mandelbrot
/// cardioid cusp where R = CΨ²'s discriminant 1 − 4CΨ vanishes (PROOF_ROADMAP_QUARTER_BOUNDARY
/// Layer 1), the Theorem 2 ceiling on the (popcount-n, popcount-n+1) coherence-block CΨ
/// over any density matrix on 2^N (PROOF_BLOCK_CPSI_QUARTER), and the squared-dimension
/// Pauli-normalization 1/d² = 1/4 at d = 2.</summary>
public sealed class QuarterAsBilinearMaxvalClaim : Claim
{
    public QuarterAsBilinearMaxvalClaim()
        : base("Quarter as bilinear maxval: max p·(1−p) = 1/4 at p = 1/2",
               Tier.Tier1Derived,
               "docs/proofs/PROOF_BLOCK_CPSI_QUARTER.md (Theorem 2 + Reading 1/4 is half of half) + docs/proofs/PROOF_ROADMAP_QUARTER_BOUNDARY.md (Layer 1 Mandelbrot cardioid)")
    { }

    public override string DisplayName => "Maxval 1/4 of any bilinear form p·(1−p)";

    public override string Summary =>
        "1/4 = (1/2)² is the universal maximum value of any bilinear form p·(1−p) on [0,1], attained at p = 1/2; surfaces at the Mandelbrot cardioid cusp, the c-block CΨ ceiling, and the d=2 Pauli-normalization";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("calculus identity",
                summary: "d/dp[p − p²] = 1 − 2p = 0 → p = 1/2; value (1/2)·(1/2) = 1/4; off-peak p(1−p) < 1/4 strictly");
            yield return new InspectableNode("Mandelbrot cardioid cusp",
                summary: "discriminant 1 − 4CΨ = 0 at CΨ = 1/4 = (1/2)²; factor 4 = 2² from completing the square (PROOF_ROADMAP_QUARTER_BOUNDARY Layer 1)");
            yield return new InspectableNode("Theorem 2 c-block ceiling",
                summary: "C_block ≤ p_n · p_{n+1} ≤ 1/4 over any density matrix on 2^N, AM-GM-saturated at p_n = p_{n+1} = 1/2; chromaticity-universal at c ≥ 2 (PROOF_BLOCK_CPSI_QUARTER Theorem 2)");
            yield return new InspectableNode("d = 2 Pauli normalization",
                summary: "1/d² = 1/4 at d = 2; the squared-dimension anchor for Pauli-basis decomposition (vec[α] = Tr(σ_α·ρ)/2^N)");
            yield return new InspectableNode("companion to BilinearApexClaim",
                summary: "argmax/maxval pair: BilinearApexClaim names p = 1/2 as the argmax; this claim names 1/4 = (1/2)² as the maxval (ON_THE_HALF Coda 2026-05-07)");
        }
    }
}

/// <summary>The argmax/maxval pair of the bilinear apex closes (Tier 1 derived). The
/// framework's recurring 1/2 and 1/4 are not two independent constants but the argmax and
/// maxval of one parabola p·(1−p): 1/2 is where the bilinear is maximised
/// (<see cref="BilinearApexClaim"/>), 1/4 = (1/2)² is the value at that argmax
/// (<see cref="QuarterAsBilinearMaxvalClaim"/>). The pair is invariant: every layer that
/// instances the bilinear form inherits both numbers together, never one without the other.
///
/// <para>Tom Wicht's coda (2026-05-07, after PROOF_BLOCK_CPSI_QUARTER landed): <em>"1/4 ist
/// die Hälfte von 0.5"</em>. The quarter is the half's quadratic shadow, not a separate
/// constant. Wherever the framework finds 1/4, it is finding 1/2 squared once through the
/// natural quadratic; wherever it finds 1/2, the corresponding 1/4 ceiling comes with it.
/// </para>
///
/// <para>Inheritance evidence: the 1/2 lineage (<see cref="QubitDimensionalAnchorClaim"/>
/// six layers) and the 1/4 lineage (Mandelbrot cardioid + Theorem 2 c-block ceiling) are
/// not parallel but paired. The same parabola, two readings: the half is the axis, the
/// quarter is the height. See reflections/ON_THE_HALF.md "Coda: her quadratic shadow" for
/// the synthesis statement and docs/proofs/PROOF_BLOCK_CPSI_QUARTER.md "Reading: 1/4 is
/// half of half" for the proof-internal articulation.</para>
/// </summary>
public sealed class ArgmaxMaxvalPairClaim : Claim
{
    public ArgmaxMaxvalPairClaim()
        : base("Argmax/maxval pair of p·(1−p) closes: 1/2 and 1/4 are two readings of one parabola",
               Tier.Tier1Derived,
               "reflections/ON_THE_HALF.md (Coda her quadratic shadow, 2026-05-07) + docs/proofs/PROOF_BLOCK_CPSI_QUARTER.md (Reading 1/4 is half of half)")
    { }

    public override string DisplayName => "Argmax/maxval pair of bilinear apex (1/2, 1/4) closes";

    public override string Summary =>
        "1/2 (BilinearApexClaim, argmax) and 1/4 (QuarterAsBilinearMaxvalClaim, maxval) are the argmax/maxval pair of one parabola p·(1−p); inseparable invariant of every layer instancing the bilinear apex";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("argmax side: BilinearApexClaim",
                summary: "p = 1/2 is the universal argmax of p·(1−p); F81 50/50, F83 anti-fraction r=0, balanced Π² partition, slow-mode apex (Pi2KnowledgeBaseClaims.cs:130 BilinearApexClaim)");
            yield return new InspectableNode("maxval side: QuarterAsBilinearMaxvalClaim",
                summary: "1/4 = (1/2)² is the universal maxval of p·(1−p); Mandelbrot cardioid cusp, Theorem 2 c-block ceiling, 1/d² at d=2 (Pi2KnowledgeBaseClaims.cs QuarterAsBilinearMaxvalClaim)");
            yield return new InspectableNode("pair invariance",
                summary: "every framework layer instancing the bilinear apex inherits both 1/2 and 1/4 together; the half is the axis, the quarter is the height; the parabola does not separate them");
            yield return new InspectableNode("Tom's coda (2026-05-07)",
                summary: "'1/4 ist die Hälfte von 0.5' — the quarter is the half's quadratic shadow; the recurring 1/2 (six lineage layers) and the recurring 1/4 (Mandelbrot + Theorem 2) are one structure read two ways");
            yield return new InspectableNode("inheritance reading",
                summary: "1/4 inherits across layers because every layer instances the same quadratic form, and that form's argmax/maxval pair is invariant; the synthesis claim closes the pair, complementing HalfAsStructuralFixedPointClaim's three-faces closure of the 1/2 side");
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

/// <summary>The three faces of 1/2 close on a structural fixed point (Tier 1 derived).
/// The framework's recurring 0.5 has three readings (where we sit, where we asymptote,
/// and what we are) that are not three accidents but one structure read three ways. This
/// claim makes the closure explicit so the next reader does not have to rediscover it.
///
/// <para>Face 1 (where we are): V-Effect bridge at C = 1/2. The d = 2 off-diagonal content
/// between static sectors P_n is where dynamics happens; at C = 1 no V-Effect bridge can
/// form because there are no boundary modes to orphan; at C = 1/2 the V-Effect always
/// works (HEISENBERG_RELOADED, V_EFFECT_AS_OBSERVATION_OF_INCOMPLETENESS).</para>
///
/// <para>Face 2 (where we go): observation horizon ends at ρ_mm = I/d = I/2. The slowest
/// mode the initial state has overlap with sets the felt-time horizon; beyond that the
/// envelope flattens to the stationary equilibrium, the maximally mixed state, whose
/// diagonal entries are exactly 1/2 (ON_TWO_TIMES, PROOF_ASYMPTOTIC_SECTOR_PROJECTION).</para>
///
/// <para>Face 3 (what we are): qubit dimensional anchor 1/d at d = 2. Setting Ψ = 0 and
/// C = 1/2 in R = C(Ψ + R)² collapses to R(R − 2) = 0; the qubit IS the framework
/// polynomial at C = 1/2 (EXCLUSIONS:251, PRIMORDIAL_QUBIT §9).</para>
///
/// <para>Closure: bridge = horizon = substrate = 1/2. The horizon is not somewhere else;
/// it is the same value as the substrate which is the same value as the bridge which is
/// where we are. The inside observer, at her own asymptote, at her own substrate, looks
/// at herself and sees the same number. She is not over the horizon; she IS the horizon,
/// and we live in her (ON_THE_HALF).</para>
/// </summary>
public sealed class HalfAsStructuralFixedPointClaim : Claim
{
    public HalfAsStructuralFixedPointClaim()
        : base("Three faces of 1/2 close: bridge = horizon = substrate (structural fixed point)",
               Tier.Tier1Derived,
               "reflections/ON_THE_HALF.md + reflections/ON_TWO_TIMES.md + hypotheses/HEISENBERG_RELOADED.md + docs/EXCLUSIONS.md:251 + docs/proofs/PROOF_ASYMPTOTIC_SECTOR_PROJECTION.md")
    { }

    public override string DisplayName => "1/2 as structural fixed point (three faces close)";

    public override string Summary =>
        "C = 1/2 V-Effect bridge (where we are) = ρ_mm = I/2 horizon (where we go) = 1/d qubit anchor (what we are); the framework is self-referential at 1/2";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("face 1 (where we are): V-Effect bridge at C = 1/2",
                summary: "d = 2 off-diagonal content between static sectors P_n is where dynamics happens; at C = 1 no V-Effect bridge forms (no orphaned boundary modes); at C = 1/2 the V-Effect always works (HEISENBERG_RELOADED, V_EFFECT_AS_OBSERVATION_OF_INCOMPLETENESS)");
            yield return new InspectableNode("face 2 (where we go): observation horizon = ρ_mm",
                summary: "felt-time horizon ends at the maximally mixed state ρ_mm = I/d = I/2; diagonal entries of the most-uncertain qubit state are exactly 1/2 (ON_TWO_TIMES, PROOF_ASYMPTOTIC_SECTOR_PROJECTION)");
            yield return new InspectableNode("face 3 (what we are): qubit dimensional anchor",
                summary: "1/d = 1/2 at d = 2 (only non-zero solution of d² − 2d = 0); R = CΨ² collapses to R(R − 2) = 0 at Ψ = 0, C = 1/2; the qubit IS the framework polynomial at C = 1/2 (EXCLUSIONS:251, PRIMORDIAL_QUBIT §9)");
            yield return new InspectableNode("closure (self-referential at 1/2)",
                summary: "bridge = horizon = substrate = 1/2; the inside observer at her own asymptote at her own substrate looks at herself and sees the same number; she IS the horizon, and we live in her (ON_THE_HALF)");
        }
    }
}

/// <summary>The 90° rotation back to the mirror is the framework's first primitive memory
/// (Tier 1 derived). The factor i in F80's Spec(M) = ±2i · Spec_nontrivial(H_non-truly) is
/// a 90° rotation in the complex plane; the mathematical fact that the mirror projects
/// H's spectrum onto M by this rotation is what makes memory possible across the two
/// sides of the framework. H lives on the real axis (energies, what our side measures);
/// M lives on the imaginary axis (frequencies and decay rates, what governs time
/// evolution); the 90° rotation is the channel between them, verified bit-exact at
/// N = 3..7 plus k-body extensions.
///
/// <para>The memory mechanism: 90° from one side (+i) and 90° from the other (−i) sum to
/// an imaginary spectrum on which both readings agree. The operators differ by the F81
/// shift Π·M·Π⁻¹ = M − 2·L_{H_odd}; the spectra coincide. What the framework remembers,
/// it remembers by sharing eigenvalues across the 90° turn (ON_BOTH_SIDES_OF_THE_MIRROR,
/// 2026-04-30).</para>
///
/// <para>Lineage: Π² = I (root, Π is the involution whose square is identity, the
/// structural √I) → F1 palindrome Π·L·Π⁻¹ = −L − 2σ·I (two Π applications give a sign
/// flip = i² = two 90° rotations) → F80 Spec(M) = ±2i·Spec(H_non-truly) (the explicit
/// 90° rotation, bit-exact N=3..7 + k-body) → F81 Π·M·Π⁻¹ = M − 2·L_{H_odd} (the
/// 90°-channel made explicit at operator level) → M_anti = L_{H_odd} (the
/// Π-antisymmetric half of M IS the Π²-odd dynamics generator).</para>
///
/// <para>Companion to <see cref="HalfAsStructuralFixedPointClaim"/>: 1/2 is the framework's
/// number-anchor (where we sit, where we asymptote, what we are); 90° is the framework's
/// angle-anchor (what we can remember across the mirror). Both are d = 2 read from two
/// sides, the dimension as a number and the dimension as a rotation.</para>
/// </summary>
public sealed class NinetyDegreeMirrorMemoryClaim : Claim
{
    public NinetyDegreeMirrorMemoryClaim()
        : base("90° rotation back to the mirror (i in F80's 2i; first primitive memory)",
               Tier.Tier1Derived,
               "docs/ANALYTICAL_FORMULAS.md F80 (Spec(M) = ±2i·Spec(H_non-truly)) + docs/proofs/PROOF_F80_BLOCH_SIGNWALK.md + docs/proofs/PROOF_F81_PI_CONJUGATION_OF_M.md + reflections/ON_BOTH_SIDES_OF_THE_MIRROR.md")
    { }

    public override string DisplayName => "90° rotation back to the mirror (memory channel)";

    public override string Summary =>
        "the i in F80's Spec(M) = ±2i·Spec(H_non-truly) is a 90° rotation that maps H (real axis, energies) to M (imaginary axis, time/decay); the mirror projects everything 90° onto itself so that it does not forget";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("layer 0 (root): Π² = I",
                summary: "Π is the involution whose square is identity; structurally √I; two Π applications close back to identity (90° + 90° = 180° = i²)");
            yield return new InspectableNode("layer 1: F1 palindrome",
                summary: "Π·L·Π⁻¹ = −L − 2σ·I; at γ = 0 truly the residual closes (Π·L·Π⁻¹ = −L); the −1 is two 90° rotations summed");
            yield return new InspectableNode("layer 2 (F80): the i is explicit 90°",
                summary: "Spec(M) = ±2i · Spec_nontrivial(H_non-truly); H's real spectrum (energies) is rotated 90° onto M's imaginary spectrum (rates); bit-exact verified N = 3..7 + k-body extensions (PROOF_F80_BLOCH_SIGNWALK)");
            yield return new InspectableNode("layer 3 (F81): operator shift between sides",
                summary: "Π·M·Π⁻¹ = M − 2·L_{H_odd}; same Frobenius norm, same spectrum, different matrix elements; the shift is exactly the unitary commutator of the Π²-odd part of H (PROOF_F81_PI_CONJUGATION_OF_M)");
            yield return new InspectableNode("layer 3: M_anti = L_{H_odd}",
                summary: "Π-antisymmetric half of M IS the Π²-odd dynamics generator; the 90°-channel materialises as M = M_sym + M_anti with 50/50 norm split at pure Π²-odd N = 3");
            yield return new InspectableNode("memory channel",
                summary: "90° on one side (+i) and 90° on the other (−i) sum to an imaginary spectrum where both readings agree; what the framework remembers, it remembers by sharing eigenvalues across this 90° turn");
            yield return new InspectableNode("companion to 1/2 anchor",
                summary: "1/2 = number-anchor (where we sit, where we asymptote, what we are); 90° = angle-anchor (what we can remember across the mirror); both are readings of d = 2");
            yield return new InspectableNode("ontological anchor",
                summary: "ON_BOTH_SIDES_OF_THE_MIRROR (2026-04-30): 'the mirror projects everything 90° onto itself so that it does not forget'; Tom's recognition months before the typed claim");
        }
    }
}

/// <summary>The framework's foundational polynomial (Tier 1 derived; deepest root).
/// The dimension equation d² − 2d = 0 is the same polynomial as R = CΨ² (set Ψ = 0
/// and C = 1/2 in R = C(Ψ + R)² to collapse to R(R − 2) = 0); its two solutions are
/// d = 0 (the substrate axis) and d = 2 (the count of ±0.5 signed positions around
/// the axis). d = 1 is algebraically excluded; the polynomial does not allow it. The
/// polynomial IS the minimum-memory equation: memory requires distinction, distinction
/// requires at least two values, and the polynomial selects exactly that minimum.
///
/// <para>The pair at d = 2 is not {|0⟩, |1⟩} (conventional basis labels) but
/// {−0.5, +0.5}: signed positions around d = 0. 0.5 is the unsigned basis magnitude;
/// the ±-sign around 0 is the pair-maker. Bloch makes it concrete: ρ = (I + r·σ)/2,
/// diagonal (1 ± r_z)/2; at the maximally mixed state both diagonals are 1/2 (the
/// axis); at the pure states |0⟩ and |1⟩ they are 1 and 0, i.e. 1/2 + 0.5 and
/// 1/2 − 0.5.</para>
///
/// <para>This is the trunk that generates both framework anchors:
/// <see cref="HalfAsStructuralFixedPointClaim"/> (1/2 number-anchor, the unsigned
/// magnitude basis) and <see cref="NinetyDegreeMirrorMemoryClaim"/> (90° angle-anchor,
/// the rotation that flips +0.5 ↔ −0.5 across the d = 0 axis). The formula R = CΨ²
/// encoded its own qubit-dimensionality as the Ψ = 0 fixed point
/// (THE_BRIDGE_WAS_ALWAYS_OPEN, ON_WHAT_THE_FORMULA_KNEW): the formula knows the
/// dimension, the dimension knows the formula.</para>
/// </summary>
public sealed class PolynomialFoundationClaim : Claim
{
    public PolynomialFoundationClaim()
        : base("d²−2d=0 ↔ R=CΨ² (foundational polynomial; minimum-memory equation)",
               Tier.Tier1Derived,
               "docs/EXCLUSIONS.md:251 + docs/THE_BRIDGE_WAS_ALWAYS_OPEN.md + reflections/ON_WHAT_THE_FORMULA_KNEW.md + reflections/ON_THE_HALF.md + hypotheses/ZERO_IS_THE_MIRROR.md")
    { }

    public override string DisplayName => "d²−2d=0 ↔ R=CΨ² (foundational polynomial)";

    public override string Summary =>
        "memory requires at least two values; the polynomial d²−2d=0 selects exactly the minimum dimension that supports a pair (d=0 axis + d=2 count of ±0.5 positions); d=1 is algebraically excluded; R=CΨ² is the same polynomial with C=1/2";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("layer −2 (semantic root): minimum-memory requirement",
                summary: "memory needs distinction; distinction needs at least two values; one value alone cannot remember another; this is the bottom of the framework in words, not in algebra");
            yield return new InspectableNode("layer −1: polynomial enforces the minimum",
                summary: "d²−2d=0 has solutions {0, 2}; d=1 is not a solution; the polynomial IS the minimum-dimension equation, the algebraic version of 'at least two'");
            yield return new InspectableNode("solution d = 0: substrate axis",
                summary: "the half we conjugate around; not a dimension we inhabit; ZERO_IS_THE_MIRROR (March 2026); operator-level: the kernel of L (static sectors P_n)");
            yield return new InspectableNode("solution d = 2: count of ±0.5 positions",
                summary: "two signed positions around the d=0 axis; the qubit dimension; where we live; the structural pair is {−0.5, +0.5}, not {0, 1}");
            yield return new InspectableNode("the pair-maker: ± around 0",
                summary: "0.5 is the unsigned basis magnitude; ±-sign around 0 is the operation that turns the basis into a pair; Bloch (1 ± r_z)/2 at pure states gives 1/2 ± 0.5");
            yield return new InspectableNode("Π shifted by 0.5 = ±0.5 (spectrum-to-pair map)",
                summary: "Π has natural spectrum {+1, −1} (since Π² = I); shifting each eigenvalue by 0.5 toward 0 gives {+0.5, −0.5}, the structural pair at d=2; this is why the 1/2 factor appears everywhere in Bloch form ρ = (I + r·σ)/2 and Pauli-basis normalization a_α = Tr(σ_α·ρ)/d — the 1/2 transports Π's natural eigenvalues onto the memory axis");
            yield return new InspectableNode("R = CΨ² is the same polynomial",
                summary: "R = C(Ψ+R)² with Ψ=0 and C=1/2 collapses to R(R−2) = 0; the framework's defining identity IS the dimension equation (THE_BRIDGE_WAS_ALWAYS_OPEN, EXCLUSIONS.md:251)");
            yield return new InspectableNode("generates 1/2 number-anchor",
                summary: "1/d = 1/2 at d=2 is the unsigned magnitude of ±0.5; HalfAsStructuralFixedPoint is the synthesis (three faces close at 1/2)");
            yield return new InspectableNode("generates 90° angle-anchor",
                summary: "the i in F80's 2i is the rotation that flips +0.5 ↔ −0.5 across the d=0 axis; NinetyDegreeMirrorMemory is the synthesis (mirror does not forget)");
            yield return new InspectableNode("self-reference",
                summary: "R = CΨ² encoded its own qubit-dimensionality as the Ψ=0 fixed point of R(R−2)=0; the formula knows the dimension, the dimension knows the formula (ON_WHAT_THE_FORMULA_KNEW)");
        }
    }
}

/// <summary>The +0 / 0 / −0 polarity layer at d = 0 is the origin of the {−0.5, +0.5}
/// pair at d = 2 (Tier 1 derived). The 0.5-shift is the explicit operation that maps
/// Π's natural eigenvalue spectrum {+1, −1} through the Bloch form ρ = (I + r·σ)/2 onto
/// the qubit's diagonal coordinates 1/2 ± r/2; at the X-eigenstates |+⟩, |−⟩ this gives
/// 1/2 ± 0.5 = {1, 0}, and at the d=0 axis (r = 0) exactly 1/2.
///
/// <para>The polarity layer has no dimension of its own; it IS the polarity differentiation
/// that turns the d = 0 substrate into a layer that supports the d = 2 minimum-memory pair.
/// Without differentiation there is the pre-polarized vacuum substrate (the Stromanschluss);
/// with it there is a layer.</para>
///
/// <para>Multi-axis structure (verified at full enumeration k=2 and sample k=3):
/// at k = 2 the polarity is Z₂² = Klein-Vierergruppe with axes (bit_a, bit_b); at k ≥ 3 a
/// third axis (Y-parity) becomes independent and the polarity is Z₂³ with 8 sectors.
/// The +0/−0 reading lives on the bit_a axis; the Klein view of <see cref="KleinFourCellClaim"/>
/// is the k = 2 collapse of this Z₂³ structure.</para>
///
/// <para>Two readings, both mathematically consistent with the F-chain:
/// <list type="bullet">
///   <item><b>Palindromic-over-the-layer:</b> Π couples +0 (one site) to −0 (its reflected
///         site); the layer pre-exists as a scaffold for the polarities.</item>
///   <item><b>ARE-the-layer:</b> the +0/−0 differentiation IS what constitutes a layer;
///         layers are emergent from polarity, not pre-existing scaffolds.</item>
/// </list>
/// Both consistent; (b) is the more radical reading.</para>
///
/// <para>Companion to <see cref="PolynomialFoundationClaim"/> (the trunk d²−2d=0 selecting
/// the minimum-memory dimension), <see cref="HalfAsStructuralFixedPointClaim"/> (the
/// 1/2 number-anchor where bridge = horizon = substrate close), and
/// <see cref="NinetyDegreeMirrorMemoryClaim"/> (the 90° angle-anchor that flips
/// +0.5 ↔ −0.5). Where the trio lives at the algebraic surface, this claim names the
/// polarity-layer beneath: where the ±0.5 comes from before the 0.5-shift makes it
/// a pair around the d=0 axis.</para>
/// </summary>
public sealed class PolarityLayerOriginClaim : Claim
{
    public PolarityLayerOriginClaim()
        : base("+0/0/−0 polarity layer at d=0 (origin of the ±0.5 pair at d=2 via the 0.5-shift)",
               Tier.Tier1Derived,
               "hypotheses/THE_POLARITY_LAYER.md + hypotheses/ZERO_IS_THE_MIRROR.md + reflections/ON_THE_HALF.md + docs/EXCLUSIONS.md:251 + Core.States.PolarityState")
    { }

    public override string DisplayName => "+0/0/−0 polarity layer (origin of the 0.5-shift)";

    public override string Summary =>
        "the +0/0/−0 polarity differentiation at d=0 generates the {−0.5, +0.5} pair at d=2 via the 0.5-shift ρ = (I + r·σ)/2; multi-axis Z₂² (k=2) / Z₂³ (k≥3); two readings (palindromic-over / ARE the layer)";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("layer −1 (substrate)",
                summary: "d=0 is the active Stromanschluss, not a passive midpoint; THE_POLARITY_LAYER.md: 'the active substrate from which the system draws coherence', analogous to QFT vacuum's zero-point activity");
            yield return new InspectableNode("layer 0 (Π involution)",
                summary: "Π² = I; Π has natural eigenvalues {+1, −1}; the involution is the structural √I that the polarity layer hangs on");
            yield return new InspectableNode("layer 1 (X-basis diagonalization)",
                summary: "|+⟩ = (|0⟩+|1⟩)/√2 carries +0 (X-eigenstate +1); |−⟩ = (|0⟩−|1⟩)/√2 carries −0 (X-eigenstate −1); operationally Core.States.PolarityState.Build(N, signs) builds tensor products on this axis (e.g. |+,−,+⟩ = X-Néel, the canonical hardware initial state)");
            yield return new InspectableNode("layer 2 (the 0.5-shift)",
                summary: "ρ = (I + r·σ)/2 maps Π-spectrum {+1, −1} onto Bloch-diagonal (1±r)/2 = 1/2 ± r/2; at pure X-eigenstates (r = ±1) this is 1/2 ± 0.5 = {1, 0}; at d=0 axis (r = 0) exactly 1/2; the 1/2 in the denominator is d, the 0.5-shift around 1/2 is r/2");
            yield return new InspectableNode("layer 3 (multi-axis refinement)",
                summary: "at k = 2 the polarity is Z₂² = Klein-Vierergruppe (bit_a, bit_b); at k ≥ 3 Y-parity becomes a third independent axis; +0/−0 lives on bit_a; the 4 Klein cells (Pp/Pm/Mp/Mm) collapse from Z₂³ at k = 2 via Y-par = bit_a XOR bit_b");
            yield return new InspectableNode("layer 4 (bridges to the trio)",
                summary: "PolynomialFoundationClaim: R = C(Ψ+R)² with C=1/2 collapses to R(R−2)=0 (the qubit IS the polynomial at C=1/2); HalfAsStructuralFixedPointClaim Face 1: V-Effect bridge at C=1/2 is the d=2 off-diagonal content; NinetyDegreeMirrorMemoryClaim: the i in F80's 2i flips +0.5 ↔ −0.5 across the d=0 axis");
            yield return new InspectableNode("reading (a): palindromic-over-the-layer",
                summary: "Π couples +0 (site i) to −0 (site N−1−i, F71-mirror); the layer pre-exists as scaffold; polarities live on it");
            yield return new InspectableNode("reading (b): ARE-the-layer",
                summary: "the +0/−0 differentiation IS what constitutes the layer; without differentiation there is only the pre-polarized substrate; layers are emergent from polarity, not pre-existing scaffolds");
            yield return new InspectableNode("operational anchor (X-Néel)",
                summary: "Core.States.PolarityState.Build(N=3, [+1, -1, +1]) builds |+,−,+⟩, the canonical Marrakesh hardware initial state; projects onto the bit_a polarity axis palindromically across the chain");
            yield return new InspectableNode("Tom's reading (2026-04-30)",
                summary: "'Wir SIND das Selbst-koppelnde System. 0 ist der Stromanschluss. +0/−0 ist der Layer.' The polarity layer is not somewhere we observe from outside; we ARE the system whose bra and ket indices Π conjugates, and the +0/−0 axis is what gives our memory pair its sign");
        }
    }
}
