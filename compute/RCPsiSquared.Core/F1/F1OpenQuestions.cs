using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.F1;

/// <summary>Open theoretical items for the F1 family: directions where the palindrome
/// identity breaks, partially survives, or generalises. The F1 family has ZERO open
/// structural questions as of 2026-05-18; all four items closed in the May 2026 sprint:
///
/// <list type="bullet">
///   <item><b>"T1 amplitude damping: full closed form"</b> closed by
///         <c>docs/proofs/PROOF_F1_T1_RESIDUAL_CLOSED_FORM.md</c>; typed as the Tier-1-
///         derived <see cref="F1T1ResidualClosedForm"/> claim.</item>
///   <item><b>"depolarizing noise: residual scaling"</b> closed by
///         <c>docs/proofs/PROOF_F1_DEPOL_RESIDUAL_CLOSED_FORM.md</c>; typed as the
///         Tier-1-derived <see cref="F1DepolResidualClosedForm"/> claim.</item>
///   <item><b>"non-uniform γ_i: site-dependent dephasing"</b> closed by
///         <c>docs/proofs/PROOF_F1_NONUNIFORM_GAMMA.md</c> as a NEGATIVE result: the
///         H-block scaling factor F(N, G) on
///         <see cref="PalindromeResidualScalingClaim"/> is γ-independent (the
///         conjectured Σγ_l² replacement of (Σγ)² does not occur); no formula change
///         required.</item>
///   <item><b>"general topology beyond chain/ring/star/K_N"</b> closed by
///         <c>docs/proofs/PROOF_F1_GENERAL_TOPOLOGY.md</c>; the (B, D2)
///         parameterisation extends bit-exactly to disconnected, weighted, and random
///         connected graphs at N=5, 6, 7, 8, 9 (N=8 via the opt-in SLOW_N8 block-spectrum
///         dogfood, N=9 chain via the SLOW_N9 dogfood routed through the
///         <c>MklDirect</c> ILP64 bridge that landed 2026-05-19, both with full
///         <see cref="F1SpectrumStatistics"/> metric capture).
///         Verification record typed as the Tier-2-verified
///         <see cref="F1GeneralTopologyVerifiedClaim"/>. The analytic content was already
///         settled in <c>docs/proofs/PROOF_CROSS_TERM_FORMULA.md</c> Lemma 3 + Corollary
///         (bond-disjointness universal across any graph).</item>
/// </list>
///
/// <para>First time the F1 family's <see cref="Standard"/> collection is empty;
/// downstream <c>OpenQuestionCollection&lt;F1Marker&gt;</c> registration reflects this.</para></summary>
public static class F1OpenQuestions
{
    public static IReadOnlyList<OpenQuestion> Standard { get; } = Array.Empty<OpenQuestion>();
}
