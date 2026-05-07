using RCPsiSquared.Core.F86;
using RCPsiSquared.Core.Resonance;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.F86Main;

/// <summary>Registers an additional batch of F86 claims that exercise the architecture's
/// multi-Tier validation: <see cref="PerBlockQPeakClaim"/> (Tier1Candidate when given a
/// caveat), <see cref="UniversalShapePrediction"/> (Tier1Candidate),
/// <see cref="SigmaZeroChromaticityScaling"/> (Tier2Empirical), and a representative
/// <see cref="RetractedClaim"/>. Together with <c>RegisterF86Main</c> these cover Tier1Derived,
/// Tier1Candidate, Tier2Empirical, and Retracted in one registry pass.
///
/// <para>The Builder's DuplicateRegistration invariant is per-Type, so the four-instance
/// PerBlockQPeakClaim.Standard list and two-instance RetractedClaim.Standard list cannot
/// be registered as multiple separate entries today. We register the first of each as a
/// representative; multi-instance "list claim" registration is a separate architectural
/// question for a future iteration.</para></summary>
public static class F86ExtendedRegistration
{
    public static ClaimRegistryBuilder RegisterF86Extended(
        this ClaimRegistryBuilder builder,
        double gammaZero) =>
        builder
            .Register<PerBlockQPeakClaim>(_ => PerBlockQPeakClaim.Standard[0])
            .Register<UniversalShapePrediction>(_ =>
                new UniversalShapePrediction(
                    bondClass: BondClass.Interior,
                    expectedRatio: 0.756,
                    tolerance: 0.005,
                    witnesses: Array.Empty<UniversalShapeWitness>()))
            .Register<SigmaZeroChromaticityScaling>(_ =>
                new SigmaZeroChromaticityScaling(gammaZero))
            .Register<RetractedClaim>(_ => RetractedClaim.Standard[0]);
}
