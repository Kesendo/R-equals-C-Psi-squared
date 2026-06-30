using System;
using System.Collections.Generic;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Diagnostics.Ptf;

/// <summary>The empirical Q=20 instance refinement of <see cref="DefectReadingEquivarianceClaim"/> (the
/// structural parent): grounded by feeding the decoder across N (play session 2026-06-30,
/// DictionaryParityInvestigationTests). The decoder projects an observed reading onto the PAINTED per-site
/// α-dictionary {F_b} (painted = the per-site profile read off a PROPAGATED run at Q=20, as opposed to the
/// BARE algebraic dictionary; least-squares plain Euclidean cosine, distinct from the parent's mode-basis
/// Gram metric), and that dictionary's worst ANTI-collinear pair (min SIGNED cos, the sign-location
/// confuser) follows the PARITY of N: STRONG at odd N, WEAK at even N, measured at the canonical γ=0.05
/// (Q=20): N=3 −0.976, N=4 −0.541, N=5 −0.965, N=6 −0.378.
///
/// <para>Two honest qualifiers that keep this Tier2Empirical, not a law:
/// (1) it is a PAINTING effect, not structural. The BARE single-excitation dictionary M[b,k] is flat for
/// N≥4 (worst anti ≈ −0.5 to −0.74, no odd/even alternation; N=3 is the −1.0 single-pair degenerate case,
/// 2 bonds → one pair; verified N=3..9). Propagation at Q=20 AMPLIFIES the anti-collinearity (the parent
/// attributes the amplification, not the parity, to weight concentrating on the R-odd seesaw mode); but
/// WHY the amplification is parity-selective (odd strong, even weak) and lands on the distance-2 pair is
/// OPEN (see below).
/// (2) the realizing pair is the DISTANCE-2 bond pair (i, i+2) (e.g. (1,3) at N=5, exactly the decode
/// demo's "bond 3 weakened reads like bond 1 strengthened"), NOT the mirror pair (b, N−2−b) the parent's
/// closed form cos = Σ(−1)^{k−1}w_k addresses. So the decoder trips on the distance-2 confuser, while the
/// parent explains the (weaker) mirror-pair anti-collinearity; sibling faces of one dictionary.</para>
///
/// <para>Mechanism for the PARITY (why painting amplifies odd N specifically, why distance-2) is OPEN, the
/// concrete next move. A second, distinct channel, collinear adjacent EDGE bonds (max positive cos, pure
/// location confusion), GROWS with N (+0.846 at N=5, +0.982 at N=6) rather than alternating; it is the
/// location-only counterpart, not this sign-confusion channel. The de-loss
/// (DefectDecoder.DecodeDeviation) resolves the odd-N pair by preserving the defect SIGN, not by escaping
/// the angle (the deviation dictionary is just as anti-collinear). Anchor:
/// hypotheses/HANDSHAKE_GEOMETRY.md +
/// compute/RCPsiSquared.Diagnostics/Foundation/DefectDecoder.cs (DictionaryWorstAntiCos, the live quantity) +
/// compute/RCPsiSquared.Diagnostics/Foundation/ReadingPowerWitness.cs (the live node, inspect --root decoder) +
/// compute/RCPsiSquared.Diagnostics.Tests/Foundation/DictionaryParityInvestigationTests.cs (the verifier:
/// painted N=3..6 + bare N=3..9).</para></summary>
public sealed class DecoderAntiCollinearityParityClaim : Claim
{
    private readonly DefectReadingEquivarianceClaim _equivariance;

    public DecoderAntiCollinearityParityClaim(DefectReadingEquivarianceClaim equivariance)
        : base("Decoder α-dictionary anti-collinearity follows odd-N parity (a Q=20 painting effect): the " +
               "painted per-site location dictionary {F_b} the decoder least-squares-projects onto has worst " +
               "SIGNED cos −0.976/−0.541/−0.965/−0.378 at N=3/4/5/6 (STRONG at odd N, WEAK at even N); the worst " +
               "confuser is the distance-2 bond pair (i,i+2), e.g. (1,3) at N=5, NOT the mirror pair. The bare " +
               "single-excitation dictionary is flat for N≥4 (≈ −0.5..−0.74, no parity; N=3 is the −1.0 single-pair " +
               "degenerate case; verified N=3..9), so the parity is a propagation amplification of the anti-collinearity " +
               "(the parent attributes the amplification to the R-odd seesaw; the parity itself is OPEN): the empirical " +
               "Q=20 instance of the structural DefectReadingEquivarianceClaim (mirror-pair R-equivariance, Gram metric)",
               Tier.Tier2Empirical,
               "hypotheses/HANDSHAKE_GEOMETRY.md + " +
               "compute/RCPsiSquared.Diagnostics/Foundation/DefectDecoder.cs (DictionaryWorstAntiCos, the live quantity) + " +
               "compute/RCPsiSquared.Diagnostics/Foundation/ReadingPowerWitness.cs (the live node, inspect --root decoder) + " +
               "compute/RCPsiSquared.Diagnostics.Tests/Foundation/DictionaryParityInvestigationTests.cs (the verifier)")
    {
        _equivariance = equivariance ?? throw new ArgumentNullException(nameof(equivariance));
    }

    public string Parity =>
        "Painted α-dictionary worst anti-collinearity (min signed cos), γ=0.05 / Q=20: " +
        "N=3 −0.976, N=4 −0.541, N=5 −0.965, N=6 −0.378. STRONG at odd N (the sign-location confusion the " +
        "decoder reports ambiguous), WEAK at even N (comparatively clean).";

    public string PaintingNotStructural =>
        "The bare single-excitation dictionary M[b,k] is flat for N≥4 (worst anti ≈ −0.5..−0.74, no odd/even " +
        "alternation; N=3 is the −1.0 single-pair degenerate case; verified N=3..9, no propagation). The " +
        "parity appears only after PAINTING (propagation at Q=20), which AMPLIFIES the anti-collinearity (the " +
        "parent attributes the amplification, not the parity, to weight concentrating on the R-odd seesaw). WHY the amplification is " +
        "parity-selective (odd vs even) is OPEN. Hence Tier2Empirical, not a law.";

    public string DistanceTwoConfuser =>
        "The worst confuser is the DISTANCE-2 bond pair (i, i+2): (0,2)/(1,3) at N=5 (the decode demo's N=5, " +
        "bond 3 weakened, reads like bond 1 strengthened), not the mirror pair (b, N−2−b). The parent " +
        "DefectReadingEquivarianceClaim's closed-form cos = Σ(−1)^{k−1}w_k addresses the (weaker) mirror pair; " +
        "this claim records the actual decoder confuser.";

    public override string DisplayName =>
        "Decoder α-dictionary anti-collinearity: odd-N parity, distance-2 confuser (Q=20 painting, Tier2Empirical)";

    public override string Summary =>
        "The decoder's painted α-dictionary worst anti-collinearity follows N-parity (strong odd, weak even: " +
        "−0.976/−0.541/−0.965/−0.378 at N=3..6, Q=20); the confuser is the distance-2 bond pair, not the mirror " +
        "pair; a propagation/painting effect (the bare dictionary is flat for N≥4), the empirical Q=20 instance " +
        $"of the structural parent, with the parity mechanism OPEN ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("The parity (painted α-dictionary, Q=20)", summary: Parity);
            yield return new InspectableNode("A painting effect, not structural (bare flat for N≥4; N=3 is −1.0)", summary: PaintingNotStructural);
            yield return new InspectableNode("The confuser is the distance-2 pair, not the mirror pair", summary: DistanceTwoConfuser);
            yield return _equivariance; // the typed parent edge (the structural cousin)
        }
    }
}
