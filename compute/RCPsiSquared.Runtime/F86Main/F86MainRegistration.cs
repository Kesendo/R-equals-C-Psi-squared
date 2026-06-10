using RCPsiSquared.Core.F86;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.F86Main;

/// <summary>Registers a subset of the F86 main family: parameter-free Tier-1-derived facts
/// (<see cref="ChiralAiiiClassification"/>, <see cref="DressedModeWeightClaim"/>,
/// <see cref="F71MirrorInvariance"/>) plus the two parameterised laws
/// (<see cref="TPeakLaw"/>, <see cref="QEpLaw"/>) that take γ₀ and g_eff respectively.
///
/// <para>This is an audit-style registration: it exercises the architecture's anchor-file
/// existence check and Tier inheritance against existing Core code. Bugs surfaced and fixed
/// during the first audit pass: F71MirrorInvariance previously anchored at a non-existent
/// <c>docs/proofs/PROOF_F71.md</c>; corrected to <c>PROOF_C1_MIRROR_SYMMETRY.md</c>.</para>
///
/// <para>F86KnowledgeBase aggregator's instances are unrelated; the Runtime registry holds
/// parallel stateless instances. Out of scope here: PerBlockQPeakClaim,
/// PerBondQPeakWitnessTable, PerF71OrbitObservation, UniversalShapePrediction,
/// ShapeFunctionWitnesses (all parameterised at Tier 1 candidate or Tier 2),
/// SigmaZeroChromaticityScaling, C2UniversalShapeDerivation, RetractedClaim
/// (constructor-protected). These can be added in later iterations once their
/// parameterisation strategies are settled.</para>
///
/// <para>Requires (since 2026-06-10):
/// <see cref="PolarityArchitecture.AbsorptionTheoremClaimRegistration.RegisterAbsorptionTheoremClaim"/>
/// (and through it Pi2Family + Pi2DyadicLadder), because the
/// <see cref="TPeakLaw"/> factory declares the rung-2 parent edge.</para></summary>
public static class F86MainRegistration
{
    public static ClaimRegistryBuilder RegisterF86Main(
        this ClaimRegistryBuilder builder,
        double gammaZero,
        double gEff) =>
        builder
            .Register<ChiralAiiiClassification>(_ => new ChiralAiiiClassification())
            .Register<DressedModeWeightClaim>(_ => new DressedModeWeightClaim())
            .Register<F71MirrorInvariance>(_ => new F71MirrorInvariance())
            .Register<TPeakLaw>(b =>
            {
                // t_peak = 1/(4γ₀): the rung-2 four of the Absorption Theorem
                // (two absorption quanta, 2γ·2; per-site −2(γ₁+γ₂), exact), NOT the
                // discriminant four a₋₁ = d². See PROOF_ABSORPTION_THEOREM.md,
                // Remark (two different fours). Edge added 2026-06-10.
                _ = b.Get<AbsorptionTheoremClaim>();
                return new TPeakLaw(gammaZero);
            })
            .Register<QEpLaw>(_ => new QEpLaw(gEff));
}
