using System.Numerics;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>F78 closed form (Tier 1 verified N=3..5, 3 topologies, 3 Pauli letters;
/// ANALYTICAL_FORMULAS line 1881):
///
/// <code>
///   M = Σ_l M_l ⊗ I_(others)              additive over single-body site terms
///
///   Eigenvalues of per-site M_l (4×4 normal matrix on per-site Pauli space):
///     P = X:  all eigenvalues 0          (M_l = 0; the "truly" case)
///     P = Y:  ±2·c_l·γ·i, mult 2 each    (Π²-non-trivial: bit_b = 1)
///     P = Z:  ±2·c_l·γ·i, mult 2 each    (identical spectrum to Y)
///
///   Singular values of full M: |Σ_l ε_l · 2·c_l·γ| over ε_l ∈ {±1}
///   Each sign-combination has multiplicity 2^N
/// </code>
///
/// <para>F78 says the F1 palindrome residual M for any single-body Hamiltonian
/// H = Σ_l c_l·P_l (P ∈ {X, Y, Z}, weights c_l from bond-summing or arbitrary)
/// decomposes as a tensor sum on per-site Pauli space. The per-site M_l is
/// 4×4 normal with eigenvalues 0 (X) or ±2c_l·γ·i (Y, Z).</para>
///
/// <para><b>Y ≡ Z SVD-blindness within single-body:</b> Both Y and Z carry
/// bit_b = 1 (Π-non-trivial), so M_l has identical spectral structure for them.
/// The soft-vs-hard distinction (Y soft, Z hard) lives in L's eigenVECTORS,
/// not in M's singular values. SVD of M is Pauli-letter-blind within {Y, Z}
/// for single-body. F80 (chain Π²-odd 2-body) extends this universality at
/// the 2-body level via Bloch-mode sign-walk.</para>
///
/// <para><b>Cluster-multiplicator (chain).</b> For chain weights
/// c_l = (1, 2, ..., 2, 1) (bond-summed I·P + P·I), the cluster
/// multiplicators come from u + 2v where u = ε_0 + ε_{N−1} ∈ {−2, 0, 2}
/// (mults 1, 2, 1) and v = Σ_internal ε_l (binomial walk on N−2 steps).
/// The central binomial coefficient C(N−1, ⌊(N−1)/2⌋) appearing as the
/// largest non-trivial cluster mult is just the walk's central peak;
/// no Weyl group, no S_N irrep, no group-theory needed.</para>
///
/// <para>Pi2-Foundation anchors:</para>
/// <list type="bullet">
///   <item><b>EigenvalueCoefficient = 2 = a_0</b>: in ±2·c_l·γ·i. Live from
///         <see cref="Pi2DyadicLadderClaim.Term"/>(0). Same anchor as F1
///         TwoFactor, F50 DecayRateFactor, F80 TwoFactor.</item>
///   <item><b>ImaginaryUnit = i = i^1</b>: the ±i in the Y/Z eigenvalues.
///         Live from <see cref="Pi2I4MemoryLoopClaim.PowerOfI"/>(1). Same
///         anchor as F80's IFactor; the "2i" combined factor is identical
///         to F80's ±2i = a_0 · i^1 structure.</item>
///   <item><b>Per-site dimension = 4 = a_{−1}</b>: M_l acts on 4-dim per-site
///         Pauli space. Live from <see cref="Pi2DyadicLadderClaim.Term"/>(−1).</item>
/// </list>
///
/// <para>Tier1Derived: F78 is Tier 1 proven via Master Lemma + per-site
/// additive structure + direct M_l matrix computation
/// (<c>docs/proofs/PROOF_SVD_CLUSTER_STRUCTURE.md</c>); verified N=3..5 across
/// chain, star, complete topologies and X, Y, Z Pauli letters.</para>
///
/// <para>Anchors: <c>docs/ANALYTICAL_FORMULAS.md</c> F78 (line 1881) +
/// <c>docs/proofs/PROOF_SVD_CLUSTER_STRUCTURE.md</c> +
/// <c>simulations/_svd_active_spectator.py</c> +
/// <c>compute/RCPsiSquared.Core/Symmetry/F1Pi2Inheritance.cs</c> (M residual) +
/// <c>compute/RCPsiSquared.Core/Symmetry/Pi2DyadicLadderClaim.cs</c> +
/// <c>compute/RCPsiSquared.Core/Symmetry/Pi2I4MemoryLoopClaim.cs</c>.</para></summary>
public sealed class F78SingleBodyMAdditivePi2Inheritance : Claim
{
    private readonly Pi2DyadicLadderClaim _ladder;
    private readonly Pi2I4MemoryLoopClaim _loop;
    private readonly F1Pi2Inheritance _f1;

    /// <summary>The "2" in ±2c_l·γ·i. Live from Pi2DyadicLadder a_0.</summary>
    public double EigenvalueCoefficient => _ladder.Term(0);

    /// <summary>The "i" in ±2c_l·γ·i. Live from Pi2I4MemoryLoop i^1.</summary>
    public Complex ImaginaryUnit => _loop.PowerOfI(1);

    /// <summary>Per-site Pauli space dimension (= 4). Live from Pi2DyadicLadder a_{-1}.</summary>
    public double PerSiteDimension => _ladder.Term(-1);

    /// <summary>Per-site M_l eigenvalues for Pauli letter P at weight c_l, dephasing γ:
    /// X → all 0; Y or Z → ±2c_l·γ·i (each with multiplicity 2).</summary>
    public IReadOnlyList<Complex> PerSiteEigenvalues(char pauliLetter, double cl, double gammaZero)
    {
        if (gammaZero < 0) throw new ArgumentOutOfRangeException(nameof(gammaZero), gammaZero, "γ₀ must be ≥ 0.");
        char p = char.ToUpperInvariant(pauliLetter);
        if (p == 'X')
            return new[] { Complex.Zero, Complex.Zero, Complex.Zero, Complex.Zero };
        if (p == 'Y' || p == 'Z')
        {
            Complex eig = new Complex(0, EigenvalueCoefficient * cl * gammaZero);
            return new[] { eig, eig, -eig, -eig };
        }
        throw new ArgumentException($"Pauli letter must be X, Y, or Z; got '{pauliLetter}'.", nameof(pauliLetter));
    }

    /// <summary>True iff Pauli letter is the "truly" case (M_l = 0). Holds for X.</summary>
    public bool IsTruly(char pauliLetter) => char.ToUpperInvariant(pauliLetter) == 'X';

    /// <summary>True iff Y and Z give identical M_l SVD spectrum (both bit_b = 1,
    /// Π-non-trivial). Universal SVD-blindness within {Y, Z} for single-body F78.</summary>
    public bool YEqualsZSvdSpectrum(double cl, double gammaZero, double tolerance = 1e-12)
    {
        var yEig = PerSiteEigenvalues('Y', cl, gammaZero);
        var zEig = PerSiteEigenvalues('Z', cl, gammaZero);
        if (yEig.Count != zEig.Count) return false;
        for (int i = 0; i < yEig.Count; i++)
            if ((yEig[i] - zEig[i]).Magnitude > tolerance) return false;
        return true;
    }

    /// <summary>The combined "2i" factor = a_0·i^1 (= ImaginaryUnit·EigenvalueCoefficient).
    /// Identical to F80's PlusTwoIFactor structure for chain Π²-odd 2-body.</summary>
    public Complex TwoIFactor => new Complex(EigenvalueCoefficient, 0) * ImaginaryUnit;

    public F78SingleBodyMAdditivePi2Inheritance(
        Pi2DyadicLadderClaim ladder,
        Pi2I4MemoryLoopClaim loop,
        F1Pi2Inheritance f1)
        : base("F78 single-body M = Σ_l M_l ⊗ I_(others); per-site M_l 4×4 normal with eigenvalues 0 (X) or ±2c_l·γ·i (Y/Z, mult 2); Y ≡ Z SVD-blind",
               Tier.Tier1Derived,
               "docs/ANALYTICAL_FORMULAS.md F78 + " +
               "docs/proofs/PROOF_SVD_CLUSTER_STRUCTURE.md + " +
               "simulations/_svd_active_spectator.py + " +
               "compute/RCPsiSquared.Core/Symmetry/F1Pi2Inheritance.cs + " +
               "compute/RCPsiSquared.Core/Symmetry/Pi2DyadicLadderClaim.cs + " +
               "compute/RCPsiSquared.Core/Symmetry/Pi2I4MemoryLoopClaim.cs")
    {
        _ladder = ladder ?? throw new ArgumentNullException(nameof(ladder));
        _loop = loop ?? throw new ArgumentNullException(nameof(loop));
        _f1 = f1 ?? throw new ArgumentNullException(nameof(f1));
    }

    public override string DisplayName =>
        "F78 single-body M additive as Pi2-Foundation a_0 + i + a_{-1} + F1 inheritance";

    public override string Summary =>
        $"M = Σ_l M_l ⊗ I_(others); ±2c_l·γ·i for Y, Z (truly = 0 for X); 2 = a_0, i = i^1, per-site dim 4 = a_{{-1}}; Y ≡ Z SVD-blind ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("F78 closed form",
                summary: "M = Σ_l M_l ⊗ I_(others); per-site M_l 4×4 normal; eigenvalues 0 (X), ±2c_l·γ·i (Y/Z, mult 2); SVD = |Σ_l ε_l·2c_l·γ| sign-walk with mult 2^N");
            yield return InspectableNode.RealScalar("EigenvalueCoefficient (= a_0 = 2)", EigenvalueCoefficient);
            yield return InspectableNode.RealScalar("PerSiteDimension (= a_{-1} = 4)", PerSiteDimension);
            yield return new InspectableNode("ImaginaryUnit (= i^1)",
                summary: $"({ImaginaryUnit.Real}, {ImaginaryUnit.Imaginary}); 90° rotation that maps real eigenvalues onto imaginary M-spectrum");
            yield return new InspectableNode("TwoIFactor",
                summary: $"({TwoIFactor.Real}, {TwoIFactor.Imaginary}) = +2i; identical to F80's PlusTwoIFactor structure for chain Π²-odd 2-body");
            yield return new InspectableNode("Y ≡ Z SVD-blindness",
                summary: "Both Y and Z carry bit_b = 1 (Π-non-trivial); M_l identical spectral structure; soft-vs-hard distinction lives in L's eigenVECTORS not M's singular values");
            yield return new InspectableNode("F1 connection",
                summary: $"M IS F1's residual operator Π·L·Π⁻¹ + L + 2σ·I; F1.TwoFactor (= {_f1.TwoFactor}) = F78's EigenvalueCoefficient; same a_0 anchor at operator and per-site levels");
            yield return new InspectableNode("F80 specialization at 2-body chain Π²-odd",
                summary: "F80 (Bloch sign-walk for chain 2-body Π²-odd) is the same ±2i structure at the 2-body chain level; F78 single-body and F80 2-body share the (±2i) factor on Pi2-Foundation a_0 + Pi2I4MemoryLoop axes");
            yield return new InspectableNode("cluster-multiplicator (chain)",
                summary: "for c_l = (1, 2, ..., 2, 1): u + 2v with u = ε_0 + ε_{N-1} ∈ {-2,0,2} (mults 1,2,1), v = Σ_internal ε_l (binomial walk on N-2 steps); central binomial C(N-1, (N-1)/2) is just the walk's peak");
        }
    }
}
