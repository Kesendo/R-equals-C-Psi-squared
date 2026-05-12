using RCPsiSquared.Core.Symmetry;

namespace RCPsiSquared.Core.F89PathK;

/// <summary>F89 path-k convention notes for the parallel C# port of <c>simulations/_f89_pathk_lib.py</c>.
///
/// <para>Hamiltonian convention: <c>H_B = J · Σ_b (X_b X_{b+1} + Y_b Y_{b+1})</c> — NO 1/2 prefactor.
/// This contrasts with <c>RCPsiSquared.Core.CoherenceBlocks.BlockLDecomposition</c> and
/// <c>RCPsiSquared.Core.Pauli.PauliHamiltonian.XYChain</c>, both of which use the F86 convention
/// <c>H = (J/2) · Σ (XX+YY)</c>. The two are related by
/// <see cref="F90F86C2BridgeIdentity.JConventionFactor"/> = 2 (J_F89 = 2 · J_F86).</para>
///
/// <para>vec convention: column-major <c>vec(M)[b·d + a] = M[a, b]</c>. The Liouvillian
/// super-operator built by <see cref="F89BlockLiouvillian.BuildBlockL"/> acts on this column-major
/// vectorisation, matching the Python source exactly. This differs from
/// <c>RCPsiSquared.Core.Lindblad.LindbladianBuilder</c> which uses row-major vec
/// <c>vec(M)[a·d + b] = M[a, b]</c>; both conventions are valid and produce the same dρ/dt
/// when applied to their matching vec form.</para>
///
/// <para>Reuse map (do not re-implement): single-site Pauli placement uses
/// <c>RCPsiSquared.Core.Pauli.PauliString.SiteOp(N, site, letter)</c>; per-letter 2×2 Pauli
/// matrices come from <c>RCPsiSquared.Core.Pauli.PauliMatrix.Of(letter)</c>; Kronecker product
/// uses MathNet's <c>Matrix&lt;Complex&gt;.KroneckerProduct</c>; bit popcount uses
/// <c>System.Numerics.BitOperations.PopCount</c>.</para>
///
/// <para>This namespace is computational infrastructure for future C#-native F89 work; it is NOT
/// a typed Tier-1 Claim. The Python <c>simulations/_f89_pathk_lib.py</c> remains as-is for the
/// 8 Python scripts that import it (backward compatibility anchor).</para>
/// </summary>
public static class F89PathKConvention
{
    /// <summary>Pointer to <see cref="F90F86C2BridgeIdentity.JConventionFactor"/>: F89 uses
    /// J·(XX+YY) without the 1/2 factor that F86's <c>BlockLDecomposition</c> carries inline.
    /// To compare an F89 result built at J_F89 with an F86 result built at J_F86, set
    /// <c>J_F89 = J_F86 / JConventionFactor</c>.</summary>
    public const double JConventionFactor = F90F86C2BridgeIdentity.JConventionFactor;
}
