using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using ComplexMatrix = MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex>;

namespace RCPsiSquared.Core.Pauli;

/// <summary>Sliding-window k-body chain Hamiltonian builder (parallel to Python
/// framework's <c>_build_kbody_chain</c>).
///
/// <para>For each template (Letters = a₁..aₖ, Coefficient = c), and each chain
/// position l ∈ [0, N-k], sums c · σ_{a₁}^l ⊗ σ_{a₂}^{l+1} ⊗ ... ⊗ σ_{aₖ}^{l+k-1}
/// (with identity at all other positions). Generalises <see cref="PauliHamiltonian.Bilinear"/>
/// from k=2 (bond list) to arbitrary k via chain sliding window.</para>
///
/// <para>Mixed-k templates in the same list are allowed; each template uses its
/// own k independently. Throws <see cref="ArgumentException"/> if any template
/// has more letters than N.</para>
///
/// <para>See <see cref="PauliString.SiteOp"/> for the per-site identity-padded
/// Pauli operator and <see cref="PauliString.Build"/> for the full N-letter tensor.
/// F104 uses this to feed k≥3 Hamiltonians into <c>PauliPairTrichotomy.Classify</c>.</para></summary>
public static class PauliKBodyChainExtensions
{
    public static ComplexMatrix ChainKBody(this IReadOnlyList<PauliTerm> templates, int N)
    {
        if (templates is null) throw new ArgumentNullException(nameof(templates));
        if (N <= 0) throw new ArgumentOutOfRangeException(nameof(N), $"N must be positive; got {N}");

        int d = 1 << N;
        var H = Matrix<Complex>.Build.Dense(d, d);

        foreach (var template in templates)
        {
            int k = template.Letters.Count;
            if (k > N)
                throw new ArgumentException(
                    $"k-body template of length {k} cannot fit on N={N} chain");

            for (int l = 0; l <= N - k; l++)
            {
                var op = PauliString.SiteOp(N, l, template.Letters[0]);
                for (int i = 1; i < k; i++)
                    op = op * PauliString.SiteOp(N, l + i, template.Letters[i]);
                H = H + template.Coefficient * op;
            }
        }

        return H;
    }
}
