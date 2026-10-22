using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Numerics;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>The one closed 2 × 2 behind many of the repository's exceptional points (Tier 1 derived,
/// <c>docs/THE_DOUBLE_ROOT.md</c>, section "One 2 × 2 in many settings"). Any 2 × 2 is −c·1 + M₀ with M₀ traceless,
/// −c its half trace, and M₀² = d·1 (Cayley-Hamilton), so the eigenvalues are −c ± √d and meet where d = 0; there M₀,
/// when non-zero, is nilpotent: one 2 × 2 Jordan block (a 2 × 2 with a double root and a non-zero traceless part is
/// defective, the rule in <c>docs/CAUGHT_ERRORS.md</c>'s coherence-horizon entry). Near any isolated meeting of two
/// eigenvalues a matrix reduces to such a 2 × 2 (Kato's local picture); what the settings this claim gathers share is
/// that the 2 × 2 is exact and closed, written out by its page.
///
/// <para>Two writings recur: real entries of opposite sign [[a, g], [−g, b]] and the complex-symmetric
/// [[a, ig], [ig, b]], carried into each other by the similarity diag(1, i). In a complex-symmetric writing the
/// coalesced vector is self-orthogonal, vᵀv = 0; the similarity does not keep vᵀv, so self-orthogonality belongs to the
/// bilinear form of the writing, not to the block.</para>
///
/// <para>Its children in the typed graph are the claims that own a setting of it, most of them for one part of their
/// content: <see cref="F86.QEpLaw"/> (the F86 toy's Q_EP = 2/g_eff, a definition in the two-level reduction, not the full
/// generator's exceptional point), <see cref="F86.LEffMirrorAxisClaim"/> (the toy's mirror axis),
/// <see cref="CoherenceHorizonClaim"/> (Q* = 2/√κ_N at N = 2, 3 only, where the horizon is the clean 2 × 2),
/// <see cref="TransitionBridgeF95SiblingClaim"/> (the toy's angle), <see cref="Spin32QuadrupolarMirrorsClaim"/> (the sodium
/// bright pairs), <see cref="SeedHolonomyClaim"/> (the seed toy) and the F50 exceptional couplings
/// (<c>ExceptionalCouplingSetClaim</c>, the (1,1) pair at N = 2). The toy's Q_EP and the horizon's Q* are critical damping
/// of a 2 × 2 whose diagonal splitting is 4γ, the same shape; they are different modes (docs/Q_REGIME_ANCHORS.md).
/// Related instruments outside the graph: <see cref="Resonance.EpAlgebra"/> (the toy's closed forms), the Diagnostics'
/// <c>ExceptionalPointClock</c> and <c>EpCharacterWitness</c> (the known-defective toy), MirrorWorld's <c>FlavorSplit</c>.
/// Where an owner hands a row's point out (TwoLevelEpModel, the sodium witness's EP points), the witness reads it there; the
/// remaining rows are the pages' blocks written into the witness.</para>
///
/// <para>Live witness: <c>inspect --root twobytwo</c>
/// (<c>compute/RCPsiSquared.Diagnostics/Foundation/TwoByTwoExceptionalPointWitness.cs</c>).</para></summary>
public sealed class TwoByTwoExceptionalPointClaim : Claim
{
    /// <summary>The half trace (M₀₀ + M₁₁)/2, which is −c in M = −c·1 + M₀.</summary>
    public static GaussianRational HalfTrace(GaussianRational[,] m) => (m[0, 0] + m[1, 1]) / 2;

    /// <summary>M₀ = M − (tr M / 2)·1.</summary>
    public static GaussianRational[,] TracelessPart(GaussianRational[,] m)
    {
        var h = HalfTrace(m);
        return new[,] { { m[0, 0] - h, m[0, 1] }, { m[1, 0], m[1, 1] - h } };
    }

    /// <summary>d with M₀² = d·1: ((M₀₀ − M₁₁)/2)² + M₀₁·M₁₀.</summary>
    public static GaussianRational D(GaussianRational[,] m)
    {
        var half = (m[0, 0] - m[1, 1]) / 2;
        return half * half + m[0, 1] * m[1, 0];
    }

    /// <summary>The two eigenvalues meet and form one 2 × 2 Jordan block: d = 0 and M₀ ≠ 0.</summary>
    public static bool IsJordanBlock(GaussianRational[,] m)
    {
        var m0 = TracelessPart(m);
        bool nonZero = !(m0[0, 0].IsZero && m0[0, 1].IsZero && m0[1, 0].IsZero && m0[1, 1].IsZero);
        return D(m).IsZero && nonZero;
    }

    /// <summary>A non-zero vector in the kernel of M₀ at a Jordan point (the coalesced eigenvector).</summary>
    public static (GaussianRational X, GaussianRational Y) CoalescedVector(GaussianRational[,] m)
    {
        var m0 = TracelessPart(m);
        return !(m0[0, 1].IsZero && m0[0, 0].IsZero) ? (m0[0, 1], -m0[0, 0]) : (-m0[1, 1], m0[1, 0]);
    }

    /// <summary>vᵀv, the bilinear (not Hermitian) square, zero exactly for a self-orthogonal vector.</summary>
    public static GaussianRational BilinearSquare((GaussianRational X, GaussianRational Y) v) => v.X * v.X + v.Y * v.Y;

    public TwoByTwoExceptionalPointClaim()
        : base(
            "Any 2 × 2 is −c·1 + M₀ with M₀² = d·1: eigenvalues −c ± √d, meeting at d = 0, and there M₀ ≠ 0 is one " +
            "Jordan block; the repository's settings where this 2 × 2 is exact and closed are one object written in " +
            "several books; self-orthogonality of the coalesced vector belongs to the complex-symmetric writing",
            Tier.Tier1Derived,
            "docs/THE_DOUBLE_ROOT.md (One 2 × 2 in many settings) + simulations/two_by_two_ep_family.py")
    {
    }

    public override string DisplayName => "One 2 × 2 in many settings: the exact closed block behind many exceptional points";

    public override string Summary =>
        "M = −c·1 + M₀, M₀² = d·1, a Jordan block where d = 0 and M₀ ≠ 0; gathered over the horizon, F50, the ring, the " +
        "proton, the edge block, the sodium pairs, the neural pair, the F86 toy, the seed toy and the Route-B control; beside " +
        $"them F155's dimer, an exceptional point of a Hamiltonian ({Tier.Label()}; live witness: inspect --root twobytwo)";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("the law",
                summary: "M = −c·1 + M₀, −c the half trace; M₀² = d·1 for every traceless 2 × 2, d = ((M₀₀ − M₁₁)/2)² + M₀₁M₁₀; " +
                         "eigenvalues −c ± √d; at d = 0 with M₀ ≠ 0 one 2 × 2 Jordan block. A companion matrix with a double " +
                         "root is always one, so for a page that writes a polynomial the verdict is the physical block's");
            yield return new InspectableNode("two writings",
                summary: "[[a, g], [−g, b]] and [[a, ig], [ig, b]], carried into each other by diag(1, i); vᵀv = 0 for the " +
                         "coalesced vector only in the complex-symmetric writing (the ring's (1, −1), vᵀv = 2, becomes (1, i))");
            yield return new InspectableNode("one shape, two modes",
                summary: "Q_EP = 2/g_eff (QEpLaw, the F86 toy, a definition in the two-level reduction) and Q* = 2/√κ_N " +
                         "(CoherenceHorizonClaim at N = 2, 3) are both d = 0 of a 2 × 2 whose diagonal splitting is 4γ, " +
                         "d = 4γ² − J²g² and d = 4γ² − κJ²: critical damping of the same shape, read on different modes " +
                         "(docs/Q_REGIME_ANCHORS.md)");
        }
    }
}
