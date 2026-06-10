using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.PolarityArchitecture;

/// <summary>Schicht-1 wiring of <see cref="F25CPsiBellPlusPi2Inheritance"/>:
/// F25's <c>CΨ(t) = f(1+f²)/6</c> Bell+ Z-dephasing closed form. F25 is the
/// mother claim of F57's prefactor 1.080088 = 2 / 1.851701. Three typed parent
/// edges:
///
/// <list type="bullet">
///   <item><see cref="Pi2DyadicLadderClaim"/>: provides <c>a_{−1} = 4</c>
///         (decay rate), <c>a_0 = 2</c> (dCΨ/dt coefficient), <c>a_3 = 1/4</c>
///         (CΨ crossing threshold).</item>
///   <item><see cref="QuarterAsBilinearMaxvalClaim"/>: typed parent grounding
///         CrossingThreshold = 1/4 on the bilinear-apex maxval axis. Same
///         anchor as F57, Dicke, F60 fold, F62 fold. Added 2026-05-16.</item>
///   <item><see cref="AbsorptionTheoremClaim"/>: the Bell+ rate 4γ is the
///         rung-2 four (two absorption quanta), split per-site by the exact
///         |00⟩⟨11| eigenmode at −2(γ₁+γ₂). Added 2026-06-10.</item>
/// </list>
///
/// <para>Tier consistency: F25 is Tier 1 proven (PROOF_MONOTONICITY_CPSI);
/// O(1) evaluation. All four claims Tier1Derived (5 ≥ 5).</para>
///
/// <para>Requires: <see cref="Pi2DyadicLadderRegistration.RegisterPi2DyadicLadder"/>
/// + <see cref="Pi2FamilyRegistration.RegisterPi2Family"/> (which registers
/// QuarterAsBilinearMaxvalClaim)
/// + <see cref="AbsorptionTheoremClaimRegistration.RegisterAbsorptionTheoremClaim"/>
/// in the builder pipeline.</para></summary>
public static class F25CPsiBellPlusPi2InheritanceRegistration
{
    public static ClaimRegistryBuilder RegisterF25CPsiBellPlusPi2Inheritance(
        this ClaimRegistryBuilder builder) =>
        builder.Register<F25CPsiBellPlusPi2Inheritance>(b =>
        {
            var ladder = b.Get<Pi2DyadicLadderClaim>();
            var quarter = b.Get<QuarterAsBilinearMaxvalClaim>();
            // The carrier 4γ = rung 2 of the absorption ladder, split per-site into
            // two quanta 2γ_l (|00⟩⟨11| exact at −2(γ₁+γ₂), residual 0.0; see
            // simulations/at_rung2_per_site_split.py). Edge added 2026-06-10.
            _ = b.Get<AbsorptionTheoremClaim>();
            return new F25CPsiBellPlusPi2Inheritance(ladder, quarter);
        });
}
