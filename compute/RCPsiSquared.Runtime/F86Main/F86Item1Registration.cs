using RCPsiSquared.Core.CoherenceBlocks;
using RCPsiSquared.Core.F86.Item1Derivation;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.F86Main;

/// <summary>Registers the lightweight c=2 Item 1 derivation primitives:
/// <see cref="C2BlockShape"/> (Stage A: closed-form HD=1 and HD=3 pair counts) and
/// <see cref="C2ChannelUniformAnalytical"/> (Stage A2: closed-form |c_1⟩, |c_3⟩ vectors).
/// Both are constructed from a <see cref="CoherenceBlock"/> at c=2 and complete in
/// microseconds (no SVD, no eigendecomposition, no Q-scan).
///
/// <para>The heavier c=2 primitives are intentionally skipped here:</para>
/// <list type="bullet">
///   <item><see cref="C2InterChannelAnalytical"/>: Stage A3 SVD-top vectors,
///         Tier2Verified due to even-N σ_0 degeneracy.</item>
///   <item><see cref="C2BondCoupling"/>: Stages B1-B3 per-bond V_b matrices.</item>
///   <item><see cref="C2EffectiveSpectrum"/>: Stage C 4-mode L_eff eigendecomposition.</item>
///   <item><see cref="C2KShape"/>: Stage D Duhamel evaluation.</item>
///   <item><see cref="C2HwhmRatio"/>: Stage D2 full Q-scan
///         (153 Q-grid points × eigendecompositions per N).</item>
/// </list>
///
/// <para>These all cache via Lazy{T} inside F86KnowledgeBase. Registering them in the
/// runtime would force eager construction at Build() time and add seconds to the registry
/// build cost. A future iteration can introduce a Lazy registration path that defers their
/// construction until first Get{T}() access.</para></summary>
public static class F86Item1Registration
{
    public static ClaimRegistryBuilder RegisterF86Item1Light(
        this ClaimRegistryBuilder builder,
        int N,
        int n,
        double gammaZero)
    {
        var block = new CoherenceBlock(N, n, gammaZero);
        if (block.C != 2)
            throw new ArgumentException(
                $"F86 Item 1 derivation applies only to the c=2 stratum; got c={block.C} from (N={N}, n={n}).",
                nameof(n));

        return builder
            .Register<C2BlockShape>(_ => new C2BlockShape(block))
            .Register<C2ChannelUniformAnalytical>(b =>
            {
                _ = b.Get<C2BlockShape>();
                return new C2ChannelUniformAnalytical(block);
            });
    }
}
