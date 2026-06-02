using System.Numerics;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>The crossover mirror is the canonical Π turned by √(NinetyDegreeMirror)
/// (Tier 1 candidate; bit-exact numerical witness, analytical derivation open; 2026-06-02).
///
/// <para>The two-term bilinear Hamiltonians XZ+YZ and ZX+ZY were once read as having a
/// "genuinely non-local" palindrome mirror (operator-Schmidt / Choi rank 8-9). They do not:
/// they have a LOCAL mirror, a single continuous per-site unitary M with Π = M^⊗N, verified
/// to machine precision at N = 2..6 (see <c>experiments/PI_OPERATOR_ENTANGLEMENT.md</c>). This
/// claim names where that continuous M sits in the Pi2-Z₄ angle branch: it is the canonical Π
/// turned by the square root of the 90° angle-anchor.</para>
///
/// <para><b>The identity (bit-exact).</b> In the framework Pauli order [I, X, Z, Y], let P1 be
/// the canonical Π (I↦X, X↦I, Z↦iY, Y↦iZ) and M the crossover mirror
/// (I↦−(X+Y)/√2, X↦(I+iZ)/√2, Z↦i(X−Y)/√2, Y↦(I−iZ)/√2). Then
/// <c>S := M·P1⁻¹</c> is block-diagonal (a pure rotation, no dark{I,Z}↔light{X,Y} swap), it
/// turns the light plane {X, Y} by exactly 45°, and</para>
///
/// <code>
///   S_light² = [[0, −1], [1, 0]] = the σ_x↔σ_y 90° rotation (X↦Y, Y↦−X)
///            = the NinetyDegreeMirror itself        (residual ~3·10⁻¹⁶)
/// </code>
///
/// <para>so <c>S_light = √(NinetyDegreeMirror)</c>. The crossover mirror is the canonical Π
/// turned by HALF the framework's 90° angle-anchor, the 45° bisector between the discrete
/// crossovers P1 (0°) and P4 (90°). The crossover Hamiltonian is the first object that lands on
/// the 45° point the discrete <see cref="Pi2KnowledgeBase"/> Z₄ skips; the weighted bond
/// a·XZ + b·YZ traces the whole continuous arc (the mirror's light turn tracks the bond
/// direction). M itself is order 4 (<c>M² = −I</c>), carrying the angle-anchor's i⁴ = 1 algebra.</para>
///
/// <para><b>Typed parent</b> <see cref="NinetyDegreeMirrorMemoryClaim"/>: the per-site-
/// conjugation face of the same 90° that the parent types on the operator/defect side (F80's
/// 2i, H↦M) and F91 on the γ-parameter side. Three faces of one 90°.</para>
///
/// <para><b>Scope (honest).</b> The clean √-of-90° identity holds in the light plane {X, Y},
/// where the X/Y diplexer conflict lives; the full 4×4 carries dark-plane phase bookkeeping
/// (Π's Z↦iY factors), so M is not a literal element of the discrete Pi2-Z₄, it is the
/// continuous √ of the Z₄ generator on the same axis. Tier1Candidate: the identity is verified
/// bit-exact (the methods below recompute it), but the analytical reason S_light² = σ_x↔σ_y has
/// not yet been derived from the conjugation equation; that derivation is the promotion step.</para>
///
/// <para>Anchors: <c>simulations/crossover_mirror_sqrt_ninety.py</c> (verification),
/// <c>experiments/PI_OPERATOR_ENTANGLEMENT.md</c> (the locality result),
/// <c>reflections/ON_THE_SQUARE_ROOT_OF_THE_MIRROR.md</c> (the synthesis).</para></summary>
public sealed class CrossoverMirrorSqrtNinetyClaim : Claim
{
    /// <summary>The 90° angle-anchor whose square root this claim names. Injected so the edge
    /// NinetyDegreeMirrorMemoryClaim → CrossoverMirrorSqrtNinetyClaim is typed.</summary>
    public NinetyDegreeMirrorMemoryClaim NinetyDegree { get; }

    public CrossoverMirrorSqrtNinetyClaim(NinetyDegreeMirrorMemoryClaim ninetyDegree)
        : base("crossover mirror = canonical Π · √(NinetyDegreeMirror): S = M·Π⁻¹ turns the light plane 45°, S_light² = σ_x↔σ_y 90° bit-exact; the local XZ+YZ/ZX+ZY mirror is the 45° middle of the Pi2-Z₄",
               Tier.Tier1Candidate,
               "experiments/PI_OPERATOR_ENTANGLEMENT.md (locality) + " +
               "simulations/crossover_mirror_sqrt_ninety.py (bit-exact witness) + " +
               "reflections/ON_THE_SQUARE_ROOT_OF_THE_MIRROR.md + " +
               "compute/RCPsiSquared.Core/Symmetry/Pi2KnowledgeBaseClaims.cs (NinetyDegreeMirrorMemory typed parent, F80's i)")
    {
        NinetyDegree = ninetyDegree ?? throw new ArgumentNullException(nameof(ninetyDegree));
    }

    private const double S = 0.70710678118654752440; // 1/√2

    /// <summary>Canonical Π (P1) in [I, X, Z, Y]: I↦X, X↦I, Z↦iY, Y↦iZ.</summary>
    private static Complex[,] P1()
    {
        var m = new Complex[4, 4];
        m[1, 0] = 1; m[0, 1] = 1; m[3, 2] = Complex.ImaginaryOne; m[2, 3] = Complex.ImaginaryOne;
        return m;
    }

    /// <summary>The crossover mirror M in [I, X, Z, Y] (unitary, M² = −I).</summary>
    private static Complex[,] CrossoverMirror()
    {
        var i = Complex.ImaginaryOne;
        var m = new Complex[4, 4];
        m[0, 1] = S; m[0, 3] = S;
        m[1, 0] = -S; m[1, 2] = i * S;
        m[2, 1] = i * S; m[2, 3] = -i * S;
        m[3, 0] = -S; m[3, 2] = -i * S;
        return m;
    }

    private static Complex[,] Mul(Complex[,] a, Complex[,] b)
    {
        var c = new Complex[4, 4];
        for (int r = 0; r < 4; r++)
            for (int col = 0; col < 4; col++)
            {
                Complex sum = Complex.Zero;
                for (int k = 0; k < 4; k++) sum += a[r, k] * b[k, col];
                c[r, col] = sum;
            }
        return c;
    }

    /// <summary>Conjugate transpose (= inverse for the unitary signed-permutation P1).</summary>
    private static Complex[,] ConjT(Complex[,] a)
    {
        var c = new Complex[4, 4];
        for (int r = 0; r < 4; r++)
            for (int col = 0; col < 4; col++)
                c[r, col] = Complex.Conjugate(a[col, r]);
        return c;
    }

    /// <summary>S = M·Π⁻¹, the turn from the canonical mirror to the crossover mirror.</summary>
    private static Complex[,] STurn() => Mul(CrossoverMirror(), ConjT(P1()));

    private static readonly int[] Dark = { 0, 2 };  // I, Z
    private static readonly int[] Light = { 1, 3 };  // X, Y

    /// <summary>True iff S = M·Π⁻¹ has no dark↔light off-blocks (a pure rotation, no swap).</summary>
    public bool TurnIsBlockDiagonal()
    {
        var s = STurn();
        double off = 0;
        foreach (var r in Dark) foreach (var c in Light) off += s[r, c].Magnitude + s[c, r].Magnitude;
        return off < 1e-12;
    }

    /// <summary>‖S_light² − [[0,−1],[1,0]]‖: distance from the light-plane square of S to the
    /// σ_x↔σ_y 90° rotation (the NinetyDegreeMirror). ≈ 0 ⟹ S_light = √(NinetyDegreeMirror).</summary>
    public double LightPlaneSquareResidual()
    {
        var s = STurn();
        // S restricted to the light plane {X(1), Y(3)}.
        var a = s[1, 1]; var b = s[1, 3]; var c = s[3, 1]; var d = s[3, 3];
        // square of the 2×2 [[a,b],[c,d]]
        var sq00 = a * a + b * c;
        var sq01 = a * b + b * d;
        var sq10 = c * a + d * c;
        var sq11 = c * b + d * d;
        // target σ_x↔σ_y 90° (X↦Y, Y↦−X) = [[0,−1],[1,0]]
        double res = (sq00 - 0).Magnitude + (sq01 - (-1)).Magnitude
                   + (sq10 - 1).Magnitude + (sq11 - 0).Magnitude;
        return res;
    }

    /// <summary>True iff S_light² equals the σ_x↔σ_y 90° rotation to machine precision.</summary>
    public bool LightPlaneSquareIsNinetyDegree() => LightPlaneSquareResidual() < 1e-9;

    /// <summary>True iff the crossover mirror satisfies M² = −I (order 4, an "i").</summary>
    public bool CrossoverMirrorIsOrderFour()
    {
        var m2 = Mul(CrossoverMirror(), CrossoverMirror());
        double res = 0;
        for (int r = 0; r < 4; r++)
            for (int c = 0; c < 4; c++)
                res += (m2[r, c] - (r == c ? new Complex(-1, 0) : Complex.Zero)).Magnitude;
        return res < 1e-12;
    }

    public override string DisplayName =>
        "crossover mirror = canonical Π · √(NinetyDegreeMirror) (the 45° middle of the Pi2-Z₄)";

    public override string Summary =>
        $"S = M·Π⁻¹ block-diagonal ({TurnIsBlockDiagonal()}), turns the light plane 45°, " +
        $"S_light² = σ_x↔σ_y 90° (residual {LightPlaneSquareResidual():E1}); the local XZ+YZ/ZX+ZY " +
        $"mirror is √ of the 90° angle-anchor, the bisector P1(0°)↔P4(90°); M²=−I " +
        $"({CrossoverMirrorIsOrderFour()}) ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("parent (typed)",
                summary: $"NinetyDegreeMirrorMemoryClaim ({NinetyDegree.Tier.Label()}): {NinetyDegree.DisplayName}; this is its continuous √, the per-site-conjugation face");
            yield return new InspectableNode("the locality result",
                summary: "XZ+YZ and ZX+ZY are LOCAL: a single continuous per-site unitary M with Π=M^⊗N mirrors the full Liouvillian (machine precision N=2..6), overturning the 'non-local rank 8-9' reading (an artifact of the discrete-permutation search + degenerate eigenvector pairing)");
            yield return new InspectableNode("S = M·Π⁻¹ is a pure rotation",
                summary: $"block-diagonal (no dark↔light swap): {TurnIsBlockDiagonal()}; turns the light plane {{X,Y}} by 45°, the bisector between P1 (0°) and P4 (90°)");
            yield return new InspectableNode("the √-of-90° identity (bit-exact)",
                summary: $"S_light² = [[0,−1],[1,0]] = σ_x↔σ_y 90° = the NinetyDegreeMirror; residual {LightPlaneSquareResidual():E2} ⟹ S_light = √(NinetyDegreeMirror)");
            yield return new InspectableNode("M is an 'i' (order 4)",
                summary: $"M² = −I: {CrossoverMirrorIsOrderFour()}; M carries the angle-anchor's i⁴=1 algebra, a concrete 90°-type element");
            yield return new InspectableNode("the continuous dial",
                summary: "weighted a·XZ+b·YZ has a product mirror at every (a,b); its light turn tracks the bond (I↦light ⊥ (a,b)). The discrete Z₄ samples the dial at cardinal points; the crossover a=b is the 45° middle");
            yield return new InspectableNode("three faces of the 90°",
                summary: "NinetyDegreeMirror (operator/defect, F80's 2i, H↦M) · F91 (γ-parameter, anti-palindromic distribution) · this claim (per-site conjugation, the √). Same Pi2-Z₄ rotational axis");
            yield return new InspectableNode("scope (honest)",
                summary: "the √-of-90° is exact in the light plane (where the X/Y conflict lives); the full 4×4 carries dark-plane phase bookkeeping, so M is the continuous √ of the Z₄ generator on the same axis, not a literal discrete-Z₄ element. Open: derive S_light²=σ_x↔σ_y analytically (the Tier1 promotion step)");
        }
    }
}
