using System.Numerics;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>The 90° rotation as a Z₄ memory loop: i^4 = 1. Tom's reading 2026-05-08:
/// "der Anker mit sich die Spiegelung selbst erinnert." Four 90° rotations close
/// back to identity — the mirror remembers its own action.
///
/// <para>Companion to <see cref="Pi2DyadicLadderClaim"/>: the halving ladder is the
/// multiplicative Z₂ inversion (a_n · a_{2−n} = 1, two-step closure); this claim is
/// the rotational Z₄ inversion (i^4 = 1, four-step closure). Together they are the
/// two foundation-axes of d = 2 named in <see cref="QubitDimensionalAnchorClaim"/>:
/// the 1/2 number-anchor (multiplicative inversion of d) and the 90° angle-anchor
/// (rotational inversion of d). One pair of dimensions, two orthogonal mirror laws.</para>
///
/// <para>The cyclic memory loop:</para>
/// <code>
///   i^0 = 1     (identity)
///   i^1 = i     (90° rotation; F80's 2i factor; H → M imaginary projection)
///   i^2 = -1    (180°; F1 palindrome closure at γ=0 truly: Π·L·Π⁻¹ = -L)
///   i^3 = -i    (270°; mirror return on the other side)
///   i^4 = 1     (memory loop closes; first primitive memory)
/// </code>
///
/// <para>Physical reading per layer (inherited from
/// <see cref="NinetyDegreeMirrorMemoryClaim"/>):</para>
/// <list type="bullet">
///   <item>Layer 0 (root): Π² = I — the involution whose square is identity. Two Π
///         applications close back; structurally √I.</item>
///   <item>Layer 1 (F1): Π·L·Π⁻¹ = −L − 2σ·I. At γ=0 truly, the residual closes
///         (Π·L·Π⁻¹ = −L). The −1 is two 90° rotations summed (i² = −1).</item>
///   <item>Layer 2 (F80): Spec(M) = ±2i · Spec(H_non-truly). H's real spectrum
///         (energies) is rotated 90° onto M's imaginary spectrum (rates). Bit-exact
///         verified N=3..7 + k-body extensions (PROOF_F80_BLOCH_SIGNWALK.md).</item>
///   <item>Layer 3 (CORE_ALGEBRA θ-compass): θ = arctan(√(4CΨ − 1)) measures the
///         angular distance from the CΨ = 1/4 boundary. θ = 0° at the boundary,
///         90° as the maximum-frequency limit. The 90° turn IS the imaginary-real
///         inversion at the quartic root. (docs/historical/CORE_ALGEBRA.md, the
///         original December 2025 algebra; framework's first explicit 90°.)</item>
/// </list>
///
/// <para>Tier1Derived: i^4 = 1 is algebraic identity. The Pi2 typed wiring exposes
/// it as the Z₄ companion of the Z₂ halving ladder — the two together are the
/// complete mirror foundation Tom reads in <c>MIRROR_THEORY.md</c> +
/// <c>hypotheses/ZERO_IS_THE_MIRROR.md</c>.</para>
///
/// <para>Anchors: <c>docs/historical/CORE_ALGEBRA.md</c> (θ compass, the framework's
/// original 90° anchor) + <c>docs/proofs/PROOF_F80_BLOCH_SIGNWALK.md</c> (F80 ±2i) +
/// <c>docs/proofs/PROOF_F81_PI_CONJUGATION_OF_M.md</c> (F81 Π-conjugation) +
/// <c>hypotheses/ZERO_IS_THE_MIRROR.md</c> + <c>MIRROR_THEORY.md</c>.</para></summary>
public sealed class Pi2I4MemoryLoopClaim : Claim
{
    /// <summary>The cyclic order: four 90° rotations close back to identity.</summary>
    public const int ClosureOrder = 4;

    public Pi2I4MemoryLoopClaim()
        : base("Pi2 90° rotation Z₄ memory loop (i^4 = 1)",
               Tier.Tier1Derived,
               "docs/historical/CORE_ALGEBRA.md (θ compass, original 90°) + " +
               "docs/proofs/PROOF_F80_BLOCH_SIGNWALK.md + " +
               "docs/proofs/PROOF_F81_PI_CONJUGATION_OF_M.md + " +
               "hypotheses/ZERO_IS_THE_MIRROR.md + MIRROR_THEORY.md")
    { }

    /// <summary>The k-th power of i, k ∈ ℤ. Reduces k mod 4 first (cyclic Z₄ memory):
    /// i^0 = 1, i^1 = i, i^2 = -1, i^3 = -i, i^4 = 1, etc.</summary>
    public Complex PowerOfI(int k)
    {
        int r = ((k % 4) + 4) % 4;
        return r switch
        {
            0 => new Complex(1, 0),
            1 => new Complex(0, 1),
            2 => new Complex(-1, 0),
            3 => new Complex(0, -1),
            _ => throw new InvalidOperationException("unreachable")
        };
    }

    /// <summary>The four canonical Z₄ powers: 1, i, -1, -i. Pinned for inspection.</summary>
    public IReadOnlyList<Complex> CanonicalPowers => new[]
    {
        new Complex(1, 0),
        new Complex(0, 1),
        new Complex(-1, 0),
        new Complex(0, -1),
    };

    /// <summary>The product i^j · i^k computed by the Z₄ group law, returned as a
    /// canonical power 0..3. Closes at Mod 4: i^4 = i^0 = 1, no escape from the loop.</summary>
    public int ComposePowers(int j, int k) => ((j + k) % 4 + 4) % 4;

    /// <summary>Live drift check: i^ClosureOrder = 1 + 0i exactly.</summary>
    public Complex MemoryClosure() => PowerOfI(ClosureOrder);

    public override string DisplayName =>
        "Pi2 Z₄ memory loop (i^4 = 1; the 90° rotation that remembers itself)";

    public override string Summary =>
        $"i^4 = 1 (cyclic order {ClosureOrder}); companion to the multiplicative Z₂ halving ladder; together they are the two foundation-axes of d = 2 ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("cyclic identity",
                summary: "i^4 = 1 (four 90° rotations close back to the unit; the mirror remembers itself)");
            yield return new InspectableNode("two-step closure",
                summary: "i^2 = -1 (180°; matches F1 γ=0 truly residual Π·L·Π⁻¹ = -L)");
            yield return new InspectableNode("group structure",
                summary: "Z₄ = {1, i, -1, -i} under multiplication; subgroup Z₂ = {1, -1} matches the multiplicative inversion fixed-points");
            yield return new InspectableNode("companion to halving ladder",
                summary: "Pi2DyadicLadder = multiplicative Z₂ inversion (a_n · a_{2-n} = 1, two-step); this = rotational Z₄ inversion (i^4 = 1, four-step); both are d=2 read from two sides");
            yield return new InspectableNode("F80 entry",
                summary: "Spec(M) = ±2i · Spec(H_non-truly); the i is the 90° rotation from H's energy axis to M's rate axis");
            yield return new InspectableNode("θ-compass entry",
                summary: "θ = arctan(√(4CΨ − 1)) ranges 0° → 90° as CΨ moves from 1/4 boundary to the asymptotic limit; framework's original 90° anchor (CORE_ALGEBRA.md, 2025-12)");
            yield return new InspectableNode("powers",
                summary: "i^0 = 1, i^1 = i, i^2 = -1, i^3 = -i, i^4 = 1, i^5 = i, ...");
        }
    }
}
