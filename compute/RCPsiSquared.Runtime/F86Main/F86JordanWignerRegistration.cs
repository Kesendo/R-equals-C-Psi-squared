using RCPsiSquared.Core.CoherenceBlocks;
using RCPsiSquared.Core.F86.JordanWigner;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.F86Main;

/// <summary>Schicht-1 wiring for the JW free-fermion track's lightweight Tier1Derived
/// primitives (Item 1' Direction (b'')). Mirrors the audit-style of
/// <see cref="F86Item1Registration"/>: only fast, parameter-free / fast-block-derived
/// primitives are eager-registered; heavier composers live in their <c>Build()</c> factories.
///
/// <para><b>Eager-registered (Tier1Derived foundation):</b></para>
/// <list type="bullet">
///   <item><see cref="XyJordanWignerModes"/> — sine-mode basis ψ_k(j) + dispersion ε_k.</item>
///   <item><see cref="BondHdChannelWeights"/> — per-bond HD=1/HD=3 column-Frobenius² split.</item>
///   <item><see cref="JwBlockBasis"/> — JW basis transformation U on the c=2 block.</item>
///   <item><see cref="JwDispersionStructure"/> — dispersion-degenerate (k, k₁, k₂) triples.</item>
///   <item><see cref="JwClusterDEigenstructure"/> — per-cluster D-eigenvalue spectra.</item>
/// </list>
///
/// <para><b>Edges declared (parent ← child):</b> JwBlockBasis ← XyJordanWignerModes;
/// JwDispersionStructure ← XyJordanWignerModes; JwClusterDEigenstructure ← JwBlockBasis +
/// JwDispersionStructure. BondHdChannelWeights is independent (uses
/// <c>BlockLDecomposition.MhPerBond</c>, not the JW modes).</para>
///
/// <para><b>Out of scope here</b> (deferred for the same reason as
/// <see cref="F86Item1Registration"/>'s heavy primitives): <see cref="JwDispersionDProjection"/>
/// (Tier2Verified), <see cref="JwBondClusterPairAffinity"/> (Tier2Verified),
/// <see cref="C2BlockJwDecomposition"/> (Tier2Verified, full per-bond bilinear),
/// <see cref="JwBondQPeakPrediction"/> + <see cref="JwBondQPeakUnified"/> (Tier2Verified
/// composers; the latter pulls <c>C2HwhmRatio</c> with its 153 Q-grid eigendecompositions).
/// A future iteration can introduce a Lazy registration path, or these can wire as their
/// own follow-up Schicht-1 plan.</para>
/// </summary>
public static class F86JordanWignerRegistration
{
    public static ClaimRegistryBuilder RegisterF86JordanWignerLight(
        this ClaimRegistryBuilder builder,
        int N,
        int n,
        double gammaZero,
        double J = 1.0)
    {
        // Block ctor itself enforces c=2 via downstream BondHdChannelWeights / JwBlockBasis;
        // we let those throw with their canonical messages rather than re-checking here.
        var block = new CoherenceBlock(N, n, gammaZero);
        return builder
            .Register<XyJordanWignerModes>(_ => XyJordanWignerModes.Build(N, J))
            .Register<BondHdChannelWeights>(_ => BondHdChannelWeights.Build(block))
            .Register<JwBlockBasis>(b =>
            {
                _ = b.Get<XyJordanWignerModes>();
                return JwBlockBasis.Build(block);
            })
            .Register<JwDispersionStructure>(b =>
            {
                _ = b.Get<XyJordanWignerModes>();
                return JwDispersionStructure.Build(N);
            })
            .Register<JwClusterDEigenstructure>(b =>
            {
                _ = b.Get<JwBlockBasis>();
                _ = b.Get<JwDispersionStructure>();
                return JwClusterDEigenstructure.Build(block);
            });
    }
}
