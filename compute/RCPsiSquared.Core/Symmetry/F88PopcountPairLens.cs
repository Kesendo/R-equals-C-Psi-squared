using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>F88 Π²-odd memory closed-form lens for a specific popcount-pair coherence
/// configuration (N, n_p, n_q). Wraps the static helpers in
/// <see cref="PopcountCoherencePi2Odd"/> as a typed parameterised Claim so the runtime
/// can register one F88 fact per coherence block (e.g. c=2 stratum at N=6 with
/// n_p=1, n_q=2) and expose it queryable through the registry.
///
/// <para>Tier1Derived: the α anchors are proven analytically via Krawtchouk
/// reflection-orthogonality (see <see cref="PopcountCoherencePi2Odd"/> docs) and
/// bit-exact verified against the universal Krawtchouk verifier for 213 configurations
/// at N=2..7 (max deviation 8.88e−16).</para>
///
/// <para>Classification of the (n_p, n_q) configuration:</para>
/// <list type="bullet">
///   <item><b>Popcount-mirror</b> (n_p + n_q = N): α = 0.</item>
///   <item><b>K-intermediate</b> (even N, exactly one of n_p, n_q is N/2):
///         α = C(N, N/2) / (2·(C(N, n_other) + C(N, N/2))).</item>
///   <item><b>Generic</b> (everything else): α = 1/2.</item>
/// </list>
///
/// <para>Π²-odd-fraction-within-memory at Hamming distance hd then reads
/// (1/2 − α·s) / (1 − s) with s = static-fraction; HD = N (full bit-flip / Π²-classical
/// anchor: GHZ_N, Bell, intra-complements) collapses to 0.</para>
///
/// <para>Anchors: <c>docs/ANALYTICAL_FORMULAS.md</c> F88 +
/// <c>docs/proofs/PROOF_F86_QPEAK.md</c> (Structural inheritance from F88) +
/// <c>compute/RCPsiSquared.Core/Symmetry/PopcountCoherencePi2Odd.cs</c> (the static helpers).</para></summary>
public sealed class F88PopcountPairLens : Claim
{
    public int N { get; }
    public int Np { get; }
    public int Nq { get; }

    /// <summary>α anchor for the (N, n_p, n_q) configuration: 0 (mirror), closed-form
    /// (K-intermediate), or 1/2 (generic). Always equals
    /// <see cref="PopcountCoherencePi2Odd.AlphaKrawtchouk"/> at machine precision.</summary>
    public double Alpha { get; }

    /// <summary>Bit-exact universal Krawtchouk re-computation of α; should match
    /// <see cref="Alpha"/>. Drift between the two indicates a regression in the closed-form
    /// classification logic.</summary>
    public double AlphaVerifier { get; }

    /// <summary>Static Frobenius² fraction s of the kernel projection of |ψ⟩⟨ψ| with
    /// |ψ⟩ = (|p⟩ + |q⟩)/√2. HD/bit-position invariant.</summary>
    public double StaticFraction { get; }

    public bool IsPopcountMirror { get; }
    public bool IsKIntermediate { get; }
    public ConfigurationKind Kind { get; }

    /// <summary>The 1/2 generic α fallback.</summary>
    public const double GenericAlpha = 0.5;

    public F88PopcountPairLens(int N, int np, int nq)
        : base($"F88 popcount-pair lens (N={N}, n_p={np}, n_q={nq})",
               Tier.Tier1Derived,
               "docs/ANALYTICAL_FORMULAS.md F88 + " +
               "docs/proofs/PROOF_F86_QPEAK.md (Structural inheritance from F88) + " +
               "compute/RCPsiSquared.Core/Symmetry/PopcountCoherencePi2Odd.cs")
    {
        if (N < 2)
            throw new ArgumentOutOfRangeException(nameof(N), N, "F88PopcountPairLens requires N ≥ 2.");
        if (np < 0 || np > N)
            throw new ArgumentOutOfRangeException(nameof(np), np, $"n_p must be in [0, {N}]; got {np}.");
        if (nq < 0 || nq > N)
            throw new ArgumentOutOfRangeException(nameof(nq), nq, $"n_q must be in [0, {N}]; got {nq}.");

        this.N = N;
        Np = np;
        Nq = nq;

        IsPopcountMirror = PopcountCoherencePi2Odd.IsPopcountMirror(N, np, nq);
        IsKIntermediate = PopcountCoherencePi2Odd.IsKIntermediate(N, np, nq);
        Kind = IsPopcountMirror ? ConfigurationKind.PopcountMirror
            : IsKIntermediate ? ConfigurationKind.KIntermediate
            : ConfigurationKind.Generic;

        Alpha = PopcountCoherencePi2Odd.AlphaAnchor(N, np, nq);
        AlphaVerifier = PopcountCoherencePi2Odd.AlphaKrawtchouk(N, np, nq);
        StaticFraction = PopcountCoherencePi2Odd.StaticFraction(N, np, nq);
    }

    /// <summary>Π²-odd-fraction-within-memory at Hamming distance <paramref name="hd"/>:
    /// (1/2 − α·s) / (1 − s); collapses to 0 at hd = N (Π²-classical anchor).</summary>
    public double Pi2OddInMemory(int hd) =>
        PopcountCoherencePi2Odd.Pi2OddInMemory(N, Np, Nq, hd);

    public override string DisplayName =>
        $"F88 lens (N={N}, n_p={Np}, n_q={Nq}, {Kind})";

    public override string Summary =>
        $"α = {Alpha:G6} ({Kind}); s = {StaticFraction:G6}; verifier-drift = {Math.Abs(Alpha - AlphaVerifier):E2} ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return InspectableNode.RealScalar("N", N);
            yield return InspectableNode.RealScalar("n_p", Np);
            yield return InspectableNode.RealScalar("n_q", Nq);
            yield return new InspectableNode("kind", summary: Kind.ToString());
            yield return InspectableNode.RealScalar("Alpha (anchor)", Alpha);
            yield return InspectableNode.RealScalar("Alpha (Krawtchouk verifier)", AlphaVerifier);
            yield return InspectableNode.RealScalar("StaticFraction", StaticFraction);
            yield return new InspectableNode("Pi2-odd anchors",
                summary: "popcount-mirror α=0 → Pi2-odd memory = 1/2; K-intermediate closed form; generic α=1/2 → Pi2-odd memory = (1/2 − s/2)/(1 − s)");
            yield return new InspectableNode("HD = N anchor",
                summary: "Pi2-odd memory = 0 at full bit-flip (GHZ_N, Bell, intra-complements)");
        }
    }
}

/// <summary>Three-way classification of an (N, n_p, n_q) configuration under F88's
/// Krawtchouk reflection-orthogonality lemma. Determines which closed-form branch
/// produces α.</summary>
public enum ConfigurationKind
{
    /// <summary>n_p + n_q = N. α = 0 (X-flip conjugation cancels all odd-|S| Paulis).</summary>
    PopcountMirror,
    /// <summary>Even N, exactly one of n_p, n_q equals N/2. α = C(N, N/2) / (2·(C(N, n_other) + C(N, N/2))).</summary>
    KIntermediate,
    /// <summary>Everything else. α = 1/2.</summary>
    Generic,
}
