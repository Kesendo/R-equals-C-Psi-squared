using System.Numerics;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>F80's <c>Spec(M) = ±2i · Spec(H_non-truly)</c> closed form has its
/// "±2i" factor sitting on the Pi2-Foundation: the 2 is <c>a_0</c> on the dyadic
/// halving ladder (= <c>d</c>, the qubit dimension), the i is <c>i^1</c> on the Z₄
/// memory loop (= 90° rotation), and the ± is the mirror-partner inversion between
/// <c>i^1</c> and <c>i^3</c> (= <c>±i</c>). One algebra, both Pi2-axes participating.
///
/// <list type="bullet">
///   <item><b>2</b> = <see cref="Pi2DyadicLadderClaim.Term"/>(0) = qubit dimension d.
///         The number-anchor side of d=2.</item>
///   <item><b>i</b> = <see cref="Pi2I4MemoryLoopClaim.PowerOfI"/>(1). The angle-anchor
///         side; the 90° rotation that maps real H-spectrum onto imaginary M-spectrum.</item>
///   <item><b>±</b> = <see cref="Pi2I4MemoryLoopClaim.MirrorPartnerIndex"/>(1) gives 1
///         (self-mirror would require the multiplicative inversion, but here we mean
///         the Z₄ partner of i: that is i^3 = −i). Both signs i and −i sit on the loop.</item>
/// </list>
///
/// <para>This is not a new identity — F80 is already Tier1Derived in
/// <c>docs/proofs/PROOF_F80_BLOCH_SIGNWALK.md</c> and verified bit-exact at N=3..7.
/// What this claim makes typed is the <b>inheritance</b>: the "2i" lives in the
/// Pi2-Foundation, both pieces (2 and i) are typed Schicht-1-facts on their own
/// axes. F80 inherits them rather than introducing a new constant.</para>
///
/// <para>Tier1Derived: pure composition. The TwoIFactor = 2 · i = (0, 2) (Complex
/// imaginary axis at magnitude 2) is computed live from the parent claims; drift
/// in either parent surfaces here.</para>
///
/// <para>Anchors: <c>docs/proofs/PROOF_F80_BLOCH_SIGNWALK.md</c> +
/// <c>docs/ANALYTICAL_FORMULAS.md</c> F80 +
/// <c>compute/RCPsiSquared.Core/Symmetry/Pi2DyadicLadderClaim.cs</c> +
/// <c>compute/RCPsiSquared.Core/Symmetry/Pi2I4MemoryLoopClaim.cs</c>.</para></summary>
public sealed class F80FactorPi2Inheritance : Claim
{
    private readonly Pi2DyadicLadderClaim _ladder;
    private readonly Pi2I4MemoryLoopClaim _loop;

    /// <summary>The "2" in F80's "±2i" — exactly <c>a_0</c> on the Pi2 dyadic ladder.</summary>
    public double TwoFactor => _ladder.Term(0);

    /// <summary>The "i" in F80's "±2i" — exactly <c>i^1</c> on the Z₄ memory loop.</summary>
    public Complex IFactor => _loop.PowerOfI(1);

    /// <summary>The "−i" partner — exactly <c>i^3</c> on the Z₄ memory loop. Together
    /// with <see cref="IFactor"/> they form the "±" of F80's Spec(M).</summary>
    public Complex MinusIFactor => _loop.PowerOfI(3);

    /// <summary>The "+2i" half of F80's spectrum factor: (0, 2). Live composition of
    /// <see cref="TwoFactor"/> · <see cref="IFactor"/>.</summary>
    public Complex PlusTwoIFactor => new Complex(TwoFactor, 0) * IFactor;

    /// <summary>The "−2i" half: (0, −2). Live composition.</summary>
    public Complex MinusTwoIFactor => new Complex(TwoFactor, 0) * MinusIFactor;

    public F80FactorPi2Inheritance(Pi2DyadicLadderClaim ladder, Pi2I4MemoryLoopClaim loop)
        : base("F80's ±2i factor inherits from Pi2-Foundation (2 = a_0, i = i^1)",
               Tier.Tier1Derived,
               "docs/proofs/PROOF_F80_BLOCH_SIGNWALK.md + docs/ANALYTICAL_FORMULAS.md F80 + " +
               "compute/RCPsiSquared.Core/Symmetry/Pi2DyadicLadderClaim.cs + " +
               "compute/RCPsiSquared.Core/Symmetry/Pi2I4MemoryLoopClaim.cs")
    {
        _ladder = ladder ?? throw new ArgumentNullException(nameof(ladder));
        _loop = loop ?? throw new ArgumentNullException(nameof(loop));
    }

    public override string DisplayName =>
        "F80 ±2i factor as Pi2-Foundation composition";

    public override string Summary =>
        $"Spec(M) = ±2i · Spec(H_non-truly): 2 = a_0 (dyadic ladder), i = i^1 (Z₄ loop), ±2i = ({PlusTwoIFactor}, {MinusTwoIFactor}) ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("F80 closed form",
                summary: "Spec(M) = ±2i · Spec(H_non-truly); Tier1Derived in PROOF_F80_BLOCH_SIGNWALK; bit-exact verified N=3..7");
            yield return InspectableNode.RealScalar("TwoFactor (= a_0 = d)", TwoFactor);
            yield return new InspectableNode("IFactor (= i^1)",
                summary: $"({IFactor.Real}, {IFactor.Imaginary}) — the 90° rotation");
            yield return new InspectableNode("MinusIFactor (= i^3)",
                summary: $"({MinusIFactor.Real}, {MinusIFactor.Imaginary}) — the −90° partner");
            yield return new InspectableNode("PlusTwoIFactor",
                summary: $"({PlusTwoIFactor.Real}, {PlusTwoIFactor.Imaginary}) = +2i");
            yield return new InspectableNode("MinusTwoIFactor",
                summary: $"({MinusTwoIFactor.Real}, {MinusTwoIFactor.Imaginary}) = −2i");
            yield return new InspectableNode("inheritance",
                summary: "the constant 2i is not introduced by F80; both 2 and i live on Pi2-Foundation typed axes; F80 inherits from Pi2DyadicLadder (number-anchor) and Pi2I4MemoryLoop (angle-anchor)");
        }
    }
}
