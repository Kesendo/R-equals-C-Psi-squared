using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Diagnostics.F87;

/// <summary>Registers the parameterless F87 trichotomy claims into the typed-knowledge
/// runtime, with the structural chain F87Trichotomy → DissipatorResonance →
/// DissipatorAxisSelectsPolarity made explicit through typed parent edges.
///
/// <list type="bullet">
///   <item><see cref="F87TrichotomyClassification"/>: the F87 truly/soft/hard split via
///         F1 residual (Tier1Derived). No typed parents — the F87 trichotomy is the
///         foundation of this family.</item>
///   <item><see cref="DissipatorResonanceLaw"/>: SU(2)-symmetric (bit_a, bit_b)-cell
///         alignment law (Tier1Derived). Typed parent: F87TrichotomyClassification —
///         the resonance law is the empirical statement about WHERE F87-hardness lives
///         (in the cell matching the dephase letter), so F87 is its source.</item>
///   <item><see cref="DissipatorAxisSelectsPolarityClaim"/>: the typed bridge claim
///         (Tier1Derived) that declares the dissipator letter as polarity-axis
///         selector. Two typed parents:
///         <list type="bullet">
///           <item><see cref="PolarityLayerOriginClaim"/> in <c>Pi2KnowledgeBase</c>
///                 — what is being differentiated.</item>
///           <item><see cref="DissipatorResonanceLaw"/> in this family — the witness
///                 table (the docstring explicitly says "witnesses live upstream").</item>
///         </list></item>
/// </list>
///
/// <para>The DissipatorAxisSelectsPolarityClaim cross-KB edge to PolarityLayerOriginClaim
/// requires the polarity family registered first. F87CanonicalWitness is parameterised
/// by chain + terms + expected class and is left for a per-witness audit iteration.</para>
///
/// <para>Layer note: this extension lives in <c>RCPsiSquared.Diagnostics</c> because the
/// F87 typed claims live there. Diagnostics gains a project reference to Runtime to
/// reach <see cref="ClaimRegistryBuilder"/>; the production-side asymmetry remains
/// Core ← Runtime + Core ← Diagnostics → Runtime, no circular reference.</para></summary>
public static class F87FamilyRegistration
{
    public static ClaimRegistryBuilder RegisterF87Family(this ClaimRegistryBuilder builder) =>
        builder
            .Register<F87TrichotomyClassification>(_ => new F87TrichotomyClassification())
            .Register<DissipatorResonanceLaw>(b =>
            {
                _ = b.Get<F87TrichotomyClassification>();
                return new DissipatorResonanceLaw();
            })
            .Register<DissipatorAxisSelectsPolarityClaim>(b =>
            {
                _ = b.Get<PolarityLayerOriginClaim>();
                _ = b.Get<DissipatorResonanceLaw>();
                return new DissipatorAxisSelectsPolarityClaim();
            });
}
