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
    /// <summary>Set of distinct Klein indices (bit_a, bit_b) across all non-identity terms.
    /// Identity terms are excluded (their Klein index (0,0) is trivial).</summary>
    public IReadOnlySet<(int BitA, int BitB)> KleinSet =>
        Terms.Where(t => t.KBody > 0).Select(t => t.KleinIndex).ToHashSet();

    /// <summary>True if all non-identity terms share the same Klein index. Empirical
    /// structural fact (k=2 full enumeration; k=3 sample): Klein-homogeneous Hamiltonians
    /// under matched-Klein dephasing are always F87 soft or truly, never F87 hard.</summary>
    public bool IsKleinHomogeneous => KleinSet.Count <= 1;

    /// <summary>Set of distinct full Z₂³ signatures (bit_a, bit_b, Y-par) across non-identity
    /// terms. At k=2 this has the same cardinality as <see cref="KleinSet"/>; at k≥3 it
    /// can be strictly finer (Y-parity becomes independent).</summary>
    public IReadOnlySet<(int BitA, int BitB, int YParity)> FullZ2SignatureSet =>
        Terms.Where(t => t.KBody > 0).Select(t => t.FullZ2Signature).ToHashSet();

    /// <summary>True if all non-identity terms share the same full Z₂³ signature. Strictly
    /// finer than <see cref="IsKleinHomogeneous"/> at k≥3; equivalent at k=2.</summary>
    public bool IsZ2Homogeneous => FullZ2SignatureSet.Count <= 1;

    /// <summary>Per-term Klein index list (in term order, non-identity terms only).</summary>
    public IReadOnlyList<(int BitA, int BitB)> PerTermKleinIndices =>
        Terms.Where(t => t.KBody > 0).Select(t => t.KleinIndex).ToList();

    /// <summary>Per-term full Z₂³ signature list (in term order, non-identity terms only).</summary>
    public IReadOnlyList<(int BitA, int BitB, int YParity)> PerTermFullZ2Signatures =>
        Terms.Where(t => t.KBody > 0).Select(t => t.FullZ2Signature).ToList();

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
