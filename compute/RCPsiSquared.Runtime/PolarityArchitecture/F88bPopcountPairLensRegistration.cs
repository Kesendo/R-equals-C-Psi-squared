using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.PolarityArchitecture;

/// <summary>Schicht-1 wiring of <see cref="F88bPopcountPairLens"/>: a parameterised
/// (N, n_p, n_q) Tier1Derived primitive exposing F88b's α + Π²-odd memory closed form
/// for one specific popcount-pair coherence configuration.
///
/// <para>Edge declared: <see cref="F88bPopcountPairLens"/> ←
/// <see cref="PopcountCoherenceClaim"/>. The general F88b closed form (parameter-free)
/// is the parent; the lens is the parameterised application of the closed form to a
/// concrete (N, n_p, n_q) point.</para>
///
/// <para>Tier consistency: both are Tier1Derived (5 ≥ 5 trivially). The lens does not
/// downgrade the parent's tier — it is a typed projection of the same closed form.</para>
///
/// <para>Requires: <see cref="F88bPopcountCoherenceRegistration.RegisterF88bPopcountCoherence"/>
/// must have been called first so <see cref="PopcountCoherenceClaim"/> is available as a
/// parent. The (N, n_p, n_q) tuple is independent of the parent registration's parameters
/// (the closed form is universal).</para>
///
/// <para>Multiple lens instances for different (N, n_p, n_q) tuples cannot coexist in the
/// same registry (single-type-key registration). Consumers that need a sweep of points
/// either build a separate registry per tuple, or use the static helpers in
/// <see cref="PopcountCoherencePi2Odd"/> directly.</para></summary>
public static class F88bPopcountPairLensRegistration
{
    public static ClaimRegistryBuilder RegisterF88bPopcountPairLens(
        this ClaimRegistryBuilder builder,
        int N,
        int np,
        int nq) =>
        builder.Register<F88bPopcountPairLens>(b =>
        {
            _ = b.Get<PopcountCoherenceClaim>();
            return new F88bPopcountPairLens(N, np, nq);
        });
}
