using System.Numerics;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>F89 path-3 (SE, DE) S_2-symmetric sub-block characteristic
/// polynomial factorisation (Tier 1 derived; sympy verified):
///
/// <code>
///   char_(S_2-sym)(λ) = F_a(λ) · F_b(λ) · F_8(λ)
///
///   F_a(λ) = λ² + (2iq+4)·λ + (4q²+4iq+4)
///         roots: λ = −2γ + iJ·(−1±√5)
///
///   F_b(λ) = λ² + (2iq+12)·λ + (4q²+12iq+36)
///         roots: λ = −6γ + iJ·(−1±√5)
///
///   F_8(λ) = degree-8 irreducible over Q[i, √5]
/// </code>
///
/// <para>The 12-dim S_2-sym sub-block factors as deg-2 · deg-2 · deg-8.
/// Both quadratics give 4 closed-form eigenvalues at AT rates 2γ (overlap)
/// and 6γ (no-overlap) with frequencies J·(−1±√5). The octic is irreducible
/// over Q[i, √5] and conjecturally Galois-non-solvable (separate Claim).</para>
///
/// <para>Anchors: <c>simulations/_f89_path3_se_de_symbolic.py</c> (sympy
/// charpoly + factor), <c>simulations/_f89_path3_octic_factor_test.py</c>
/// (irreducibility), <c>experiments/F89_TOPOLOGY_ORBIT_CLOSURE.md</c>
/// § "Path-3 (SE, DE) S_2-symmetric sub-block: deg-2·deg-2·deg-8 factorisation".</para></summary>
public sealed class F89Path3SeDeFactorisationClaim : Claim
{
    private readonly F89TopologyOrbitClosure _f89;
    private readonly F89PathKAtLockMechanismClaim _atLock;

    /// <summary>Dimension of the S_2-symmetric sub-block of the (SE, DE) sector
    /// for path-3 (N_block=4): 12.</summary>
    public const int S2SymSubBlockDimension = 12;

    private static readonly IReadOnlyList<int> _factorDegrees = new[] { 2, 2, 8 };

    /// <summary>Factor degrees of char_(S_2-sym)(λ) = F_a · F_b · F_8: {2, 2, 8}.
    /// Sum equals <see cref="S2SymSubBlockDimension"/> = 12.</summary>
    public static IReadOnlyList<int> FactorDegrees => _factorDegrees;

    /// <summary>F_a quadratic roots at given (γ, J): λ = −2γ + iJ·(−1±√5).
    /// At J=0 both roots collapse to the degenerate double root λ = −2γ
    /// (F_a becomes (λ+2γ)² when H_B vanishes — no XY coupling regime).</summary>
    public static Complex[] FaRoots(double gamma, double j)
    {
        if (gamma < 0) throw new ArgumentOutOfRangeException(nameof(gamma), gamma, "γ must be ≥ 0.");
        if (j < 0) throw new ArgumentOutOfRangeException(nameof(j), j, "J must be ≥ 0.");
        double sqrt5 = Math.Sqrt(5);
        return new[]
        {
            new Complex(-2 * gamma, j * (-1 + sqrt5)),
            new Complex(-2 * gamma, j * (-1 - sqrt5)),
        };
    }

    /// <summary>F_b quadratic roots at given (γ, J): λ = −6γ + iJ·(−1±√5).
    /// At J=0 both roots collapse to the degenerate double root λ = −6γ.</summary>
    public static Complex[] FbRoots(double gamma, double j)
    {
        if (gamma < 0) throw new ArgumentOutOfRangeException(nameof(gamma), gamma, "γ must be ≥ 0.");
        if (j < 0) throw new ArgumentOutOfRangeException(nameof(j), j, "J must be ≥ 0.");
        double sqrt5 = Math.Sqrt(5);
        return new[]
        {
            new Complex(-6 * gamma, j * (-1 + sqrt5)),
            new Complex(-6 * gamma, j * (-1 - sqrt5)),
        };
    }

    /// <summary>F_8 octic factor irreducibility over Q[i, √5]: true (sympy
    /// verified per <c>_f89_path3_octic_factor_test.py</c>).</summary>
    public const bool OcticIsIrreducibleOverQiSqrt5 = true;

    public F89Path3SeDeFactorisationClaim(F89TopologyOrbitClosure f89, F89PathKAtLockMechanismClaim atLock)
        : base("F89 path-3 (SE,DE) S_2-sym sub-block char(λ) = F_a · F_b · F_8 (deg 2·2·8); F_a, F_b roots in closed form λ = −{2,6}γ + iJ(−1±√5); F_8 octic irreducible over Q[i, √5]",
               Tier.Tier1Derived,
               "experiments/F89_TOPOLOGY_ORBIT_CLOSURE.md + " +
               "simulations/_f89_path3_se_de_symbolic.py + " +
               "simulations/_f89_path3_octic_factor_test.py + " +
               "compute/RCPsiSquared.Core/Symmetry/F89PathKAtLockMechanismClaim.cs")
    {
        _f89 = f89 ?? throw new ArgumentNullException(nameof(f89));
        _atLock = atLock ?? throw new ArgumentNullException(nameof(atLock));
    }

    public override string DisplayName =>
        "F89 path-3 (SE,DE) S_2-sym factorisation: deg-2 · deg-2 · deg-8";

    public override string Summary =>
        $"12-dim sub-block factors as F_a (deg 2) · F_b (deg 2) · F_8 (deg 8 irreducible over Q[i,√5]); F_a/F_b roots at λ = -{{2,6}}γ + iJ(-1±√5) ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            var faRoots = FaRoots(0.05, 0.075);
            yield return new InspectableNode("F_a roots at γ=0.05, J=0.075 (q=1.5)",
                summary: $"λ ≈ {faRoots[0]} and λ ≈ {faRoots[1]}");
            var fbRoots = FbRoots(0.05, 0.075);
            yield return new InspectableNode("F_b roots at γ=0.05, J=0.075",
                summary: $"λ ≈ {fbRoots[0]} and λ ≈ {fbRoots[1]}");
            yield return new InspectableNode("Octic F_8 over Q[i, √5]",
                summary: $"Irreducible: {OcticIsIrreducibleOverQiSqrt5}");
        }
    }
}
