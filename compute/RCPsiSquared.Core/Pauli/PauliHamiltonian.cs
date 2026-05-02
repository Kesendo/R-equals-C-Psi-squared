using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.ChainSystems;
using ComplexMatrix = MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex>;

namespace RCPsiSquared.Core.Pauli;

/// <summary>An N-qubit Hamiltonian as a sum of weighted Pauli-string terms.
///
/// The dense 2^N × 2^N matrix is built lazily via <see cref="ToMatrix"/>. For typical
/// project N (≤ 8), this is up to 256 × 256 — tractable.
/// </summary>
public sealed record PauliHamiltonian(int N, IReadOnlyList<PauliTerm> Terms)
{
    public ComplexMatrix ToMatrix()
    {
        int d = 1 << N;
        var H = Matrix<Complex>.Build.Dense(d, d);
        foreach (var term in Terms)
        {
            if (term.N != N)
                throw new ArgumentException($"term has N={term.N} letters; expected N={N}");
            H = H + term.Coefficient * PauliString.Build(term.Letters);
        }
        return H;
    }

    /// <summary>Uniform XY chain: H = (J/2) Σ_b (X_b X_{b+1} + Y_b Y_{b+1}).</summary>
    public static PauliHamiltonian XYChain(int N, double J)
    {
        var terms = new List<PauliTerm>(2 * (N - 1));
        Complex c = J / 2.0;
        for (int b = 0; b < N - 1; b++)
        {
            terms.Add(PauliTerm.TwoSite(N, b, PauliLetter.X, b + 1, PauliLetter.X, c));
            terms.Add(PauliTerm.TwoSite(N, b, PauliLetter.Y, b + 1, PauliLetter.Y, c));
        }
        return new PauliHamiltonian(N, terms);
    }

    /// <summary>Uniform Heisenberg chain: H = (J/4) Σ_b (X_b X_{b+1} + Y_b Y_{b+1} + Z_b Z_{b+1}).</summary>
    public static PauliHamiltonian HeisenbergChain(int N, double J)
    {
        var terms = new List<PauliTerm>(3 * (N - 1));
        Complex c = J / 4.0;
        for (int b = 0; b < N - 1; b++)
        {
            terms.Add(PauliTerm.TwoSite(N, b, PauliLetter.X, b + 1, PauliLetter.X, c));
            terms.Add(PauliTerm.TwoSite(N, b, PauliLetter.Y, b + 1, PauliLetter.Y, c));
            terms.Add(PauliTerm.TwoSite(N, b, PauliLetter.Z, b + 1, PauliLetter.Z, c));
        }
        return new PauliHamiltonian(N, terms);
    }

    /// <summary>Bilinear bond Hamiltonian H = Σ_bond Σ_term coeff · σ_la^i σ_lb^j on the given bonds.</summary>
    public static PauliHamiltonian Bilinear(int N, IReadOnlyList<Bond> bonds,
        IReadOnlyList<(PauliLetter La, PauliLetter Lb, Complex Coeff)> terms)
    {
        var allTerms = new List<PauliTerm>(bonds.Count * terms.Count);
        foreach (var bond in bonds)
            foreach (var (la, lb, coeff) in terms)
                allTerms.Add(PauliTerm.TwoSite(N, bond.Site1, la, bond.Site2, lb, coeff));
        return new PauliHamiltonian(N, allTerms);
    }
}
