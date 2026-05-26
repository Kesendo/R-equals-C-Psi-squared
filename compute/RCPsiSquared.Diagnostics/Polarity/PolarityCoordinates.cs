using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.ChainSystems;
using RCPsiSquared.Core.Lindblad;
using RCPsiSquared.Core.Pauli;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Diagnostics.F81;
using ComplexMatrix = MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex>;

namespace RCPsiSquared.Diagnostics.Polarity;

/// <summary>Three-way polarity decomposition of the F1 residual
/// M = Π·L·Π⁻¹ + L + 2σ·I into the polarity-triple {−1/2, 0, +1/2} at d=2:
///
/// <code>
///   M_zero       = (M + Π·M·Π⁻¹) / 2                      (0-axis, Π²-symmetric = F81 M_sym)
///   M_plus_half  = (M_anti − i · Π·M_anti·Π⁻¹) / 2        (+1/2 polarity, Π eigenvalue +i)
///   M_minus_half = (M_anti + i · Π·M_anti·Π⁻¹) / 2        (−1/2 polarity, Π eigenvalue −i)
/// </code>
///
/// <para>where M_anti = (M − Π·M·Π⁻¹) / 2 is the F81 antisymmetric part. Refinement of F81:
/// F81 splits M = M_sym + M_anti by Π-conjugation parity (eigenvalues ±1 of the linear map
/// X ↦ Π·X·Π⁻¹). Π is order-4 on Liouville space (Π⁴ = I), so the full Π-eigenvalue
/// spectrum is {+1, −1, +i, −i}. The +1 / −1 eigenspaces together form M_sym (Π²-even),
/// and the +i / −i eigenspaces together form M_anti (Π²-odd). This primitive refines the
/// ±1 sub-split into the explicit +i / −i projections, giving the typed polarity triple.</para>
///
/// <para><b>Frobenius-orthogonal invariant</b>:
/// ‖M‖² = ‖M_zero‖² + ‖M_plus_half‖² + ‖M_minus_half‖². The
/// <see cref="PolarityCoordinatesResult.OrthogonalityResidual"/> field records
/// |‖M‖² − (‖M_zero‖² + ‖M_plus_half‖² + ‖M_minus_half‖²)| as a numerical sanity check
/// (machine precision when the invariant holds).</para>
///
/// <para><b>Connection to F81</b>:
/// F81 M_sym = M_zero;
/// F81 M_anti = M_plus_half + M_minus_half (further split by Π ±i eigenvalue).</para>
///
/// <para><b>F112 typed scope (Tier1Derived)</b>: for any Lindblad-form Liouvillian
/// L = -i[H, ·] + Σ_k γ_k · D[c_k] with Hermitian H and each c_k bit_b-homogeneous,
/// <see cref="PolarityCoordinatesResult.Asymmetry"/> = ‖M_plus_half‖² − ‖M_minus_half‖²
/// is exactly 0 bit-exact. See <c>docs/proofs/PROOF_F112_LINDBLAD_BIT_B_PI_BALANCE.md</c>
/// and <see cref="LindbladBitBPiBalance"/>. The diagnostic asymmetry is the precise
/// witness for c with cross-bit_b Pauli support (outside F108's closure regime) or for
/// non-Hermitian H combined with bit_b-mixed c.</para>
///
/// <para>Python sibling: <c>simulations/framework/diagnostics/polarity_coordinates.py</c>
/// (<c>polarity_coordinates_from_L</c>, <c>polarity_coordinates_from_hc</c>); this C#
/// primitive matches the Python output bit-exactly via the same Π construction
/// (<see cref="PiOperator.BuildFull"/>) and the same vec→Pauli transform
/// (<see cref="PauliBasis.VecToPauliBasisTransform"/>).</para>
/// </summary>
public sealed record PolarityCoordinatesResult(
    ComplexMatrix M,
    ComplexMatrix MZero,
    ComplexMatrix MPlusHalf,
    ComplexMatrix MMinusHalf,
    double MNormSquared,
    double MZeroNormSquared,
    double MPlusHalfNormSquared,
    double MMinusHalfNormSquared,
    double Asymmetry,
    double OrthogonalityResidual,
    double F81Violation);

/// <summary>Matrix-level compute primitive for the F112 polarity decomposition.
/// Mirrors the Python <c>polarity_coordinates_from_L</c> family by delegating M /
/// M_sym / M_anti / L_HOdd construction to <see cref="PiDecomposition.Decompose"/>
/// and adding the +i / −i Π-eigenvalue refinement of M_anti on top.</summary>
public static class PolarityCoordinates
{
    /// <summary>Decompose M for a chain-built Hamiltonian from bilinear bond terms.
    /// Delegates to <see cref="PiDecomposition.Decompose"/> for M / M_sym / M_anti /
    /// L_HOdd / F81Violation, then computes the +i / −i Π-eigenvalue refinement
    /// of M_anti.</summary>
    /// <param name="chain">N-qubit chain providing N, bonds, J, γ₀.</param>
    /// <param name="terms">Bilinear bond Pauli-pair terms; see
    /// <see cref="PiDecomposition.Decompose"/>.</param>
    /// <param name="gammaT1">Optional uniform per-site σ⁻ amplitude-damping rate.
    /// When non-null and non-zero, the Lindbladian uses the
    /// <see cref="T1Dissipator"/> path; the F112 typed Tier1Derived scope no longer
    /// covers this case (σ⁻ is bit_b-mixed), and the asymmetry need not be zero.</param>
    /// <param name="dephaseLetter">Dephasing letter for the Π operator; default Z.</param>
    public static PolarityCoordinatesResult Decompose(
        ChainSystem chain,
        IReadOnlyList<PauliPairBondTerm> terms,
        double? gammaT1 = null,
        PauliLetter dephaseLetter = PauliLetter.Z)
    {
        if (chain is null) throw new ArgumentNullException(nameof(chain));
        if (terms is null) throw new ArgumentNullException(nameof(terms));

        IReadOnlyList<double>? gammaT1PerSite = null;
        if (gammaT1.HasValue && gammaT1.Value != 0.0)
            gammaT1PerSite = Enumerable.Repeat(gammaT1.Value, chain.N).ToArray();

        var pi = PiDecomposition.Decompose(chain, terms, gammaT1PerSite);
        var piOp = PiOperator.BuildFull(chain.N, dephaseLetter);
        return RefineMAnti(pi.M, pi.MSym, pi.MAnti, piOp, pi.F81Violation);
    }

    /// <summary>Decompose M for a general k-body Pauli Hamiltonian. Builds H from
    /// raw <see cref="PauliTerm"/> entries (supports single-site h_x·X_l, 3-body
    /// X⊗X⊗X, etc., that the bilinear bond overload cannot express), constructs L
    /// via <see cref="PauliDephasingDissipator.BuildZ"/> or
    /// <see cref="T1Dissipator.Build"/>, then runs the same Π-decomposition +
    /// ±i refinement pipeline.</summary>
    /// <param name="chain">N-qubit chain providing N, γ₀ (uniform Z-dephasing rate).</param>
    /// <param name="kBodyTerms">Pauli-string terms with arbitrary k-body weight.
    /// Each term must have <see cref="PauliTerm.N"/> = <see cref="ChainSystem.N"/>.
    /// Coefficients are <see cref="Complex"/>; for Hermitian H pass real coefficients
    /// on Hermitian Pauli strings (the canonical case for F112 Tier1Derived scope).</param>
    /// <param name="gammaT1">Optional uniform per-site σ⁻ amplitude-damping rate.
    /// When non-null and non-zero, routes through <see cref="T1Dissipator"/>; otherwise
    /// pure Z-dephasing via <see cref="PauliDephasingDissipator"/>. F112 typed
    /// Tier1Derived scope covers the gammaT1=null branch (single-Pauli dephasing is
    /// bit_b-homogeneous); the T1 branch is bit_b-mixed and out of F112 typed scope.</param>
    /// <param name="dephaseLetter">Dephasing letter for the Π operator; default Z.</param>
    /// <remarks>No F81 inner check on this path: the F81 identity is a 2-body
    /// statement built on bond-term H_odd; for k-body or single-site H the analogous
    /// identity is open, and <see cref="PolarityCoordinatesResult.F81Violation"/> is
    /// reported as 0 (sentinel, no closed-form prediction to compare against).</remarks>
    public static PolarityCoordinatesResult Decompose(
        ChainSystem chain,
        IReadOnlyList<PauliTerm> kBodyTerms,
        double? gammaT1 = null,
        PauliLetter dephaseLetter = PauliLetter.Z)
    {
        if (chain is null) throw new ArgumentNullException(nameof(chain));
        if (kBodyTerms is null) throw new ArgumentNullException(nameof(kBodyTerms));
        foreach (var term in kBodyTerms)
        {
            if (term.N != chain.N)
                throw new ArgumentException(
                    $"term has N={term.N} letters; expected N={chain.N}", nameof(kBodyTerms));
        }

        var H = new PauliHamiltonian(chain.N, kBodyTerms).ToMatrix();
        var gammaZ = Enumerable.Repeat(chain.GammaZero, chain.N).ToArray();

        ComplexMatrix L;
        if (gammaT1.HasValue && gammaT1.Value != 0.0)
        {
            var gammaT1PerSite = Enumerable.Repeat(gammaT1.Value, chain.N).ToArray();
            L = T1Dissipator.Build(H, gammaZ, gammaT1PerSite);
        }
        else
        {
            L = PauliDephasingDissipator.BuildZ(H, gammaZ);
        }

        var M = PalindromeResidual.Build(L, chain.N, chain.SigmaGamma, dephaseLetter);
        var piOp = PiOperator.BuildFull(chain.N, dephaseLetter);
        var piInv = piOp.ConjugateTranspose();
        var piMpi = piOp * M * piInv;
        var mSym = (M + piMpi) / 2.0;
        var mAnti = (M - piMpi) / 2.0;

        // No F81 inner check on this path (k-body / single-site H_odd identity is open).
        return RefineMAnti(M, mSym, mAnti, piOp, f81Violation: 0.0);
    }

    /// <summary>Shared core: given M, M_sym, M_anti and Π, compute the ±i
    /// Π-eigenvalue projectors on M_anti and assemble the result record.</summary>
    private static PolarityCoordinatesResult RefineMAnti(
        ComplexMatrix M, ComplexMatrix mSym, ComplexMatrix mAnti, ComplexMatrix piOp, double f81Violation)
    {
        var piInv = piOp.ConjugateTranspose(); // Π is a unitary signed permutation; Π⁻¹ = Π†.
        var piMAntiPiInv = piOp * mAnti * piInv;
        var mPlusHalf = (mAnti - Complex.ImaginaryOne * piMAntiPiInv) / 2.0;
        var mMinusHalf = (mAnti + Complex.ImaginaryOne * piMAntiPiInv) / 2.0;
        var mZero = mSym; // F81 M_sym is the 0-axis component by definition.

        double normM = Math.Pow(M.FrobeniusNorm(), 2);
        double normZero = Math.Pow(mZero.FrobeniusNorm(), 2);
        double normPlus = Math.Pow(mPlusHalf.FrobeniusNorm(), 2);
        double normMinus = Math.Pow(mMinusHalf.FrobeniusNorm(), 2);

        double orthogonalityResidual = Math.Abs(normM - (normZero + normPlus + normMinus));
        double asymmetry = normPlus - normMinus;

        return new PolarityCoordinatesResult(
            M: M,
            MZero: mZero,
            MPlusHalf: mPlusHalf,
            MMinusHalf: mMinusHalf,
            MNormSquared: normM,
            MZeroNormSquared: normZero,
            MPlusHalfNormSquared: normPlus,
            MMinusHalfNormSquared: normMinus,
            Asymmetry: asymmetry,
            OrthogonalityResidual: orthogonalityResidual,
            F81Violation: f81Violation);
    }
}
