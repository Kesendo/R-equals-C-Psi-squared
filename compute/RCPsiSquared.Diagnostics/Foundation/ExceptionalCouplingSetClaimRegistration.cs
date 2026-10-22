using RCPsiSquared.Core.BlockSpectrum;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Diagnostics.Foundation;

/// <summary>Wiring of <see cref="ExceptionalCouplingSetClaim"/> (F50 at exceptional couplings: the count 2N holds
/// off a finite set E(N, G) of γ/J). Four typed parents: <see cref="AbsorptionTheoremClaim"/> (the floor that
/// puts a real −2γ mode at ⟨n_XY⟩ = 1), <see cref="F50WeightOneDegeneracyPi2Inheritance"/> (the count whose
/// converse is fenced), <see cref="JointPopcountSectors"/> (the grading the localisation is stated in) and, for one part of its content,
/// <see cref="TwoByTwoExceptionalPointClaim"/> (the (1,1) pair at N = 2 is a setting of the one closed 2 × 2). All four
/// must be registered.</summary>
public static class ExceptionalCouplingSetClaimRegistration
{
    public static ClaimRegistryBuilder RegisterExceptionalCouplingSetClaim(
        this ClaimRegistryBuilder builder) =>
        builder.Register<ExceptionalCouplingSetClaim>(b =>
        {
            _ = b.Get<TwoByTwoExceptionalPointClaim>();   // partial scope: the (1,1) pair at N = 2 is a setting of the one closed block
            return new ExceptionalCouplingSetClaim(
                b.Get<AbsorptionTheoremClaim>(),
                b.Get<F50WeightOneDegeneracyPi2Inheritance>(),
                b.Get<JointPopcountSectors>());
        });
}
