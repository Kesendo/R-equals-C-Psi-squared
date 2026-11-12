using System.Numerics;
using RCPsiSquared.Core.Numerics;
using ComplexMatrix = MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex>;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>The lemma of F173's Theorem H (<c>docs/proofs/PROOF_SIGNED_MASK_BALANCE.md</c>): L = −i[J·H, ·] +
/// Σ_l γ_l (Z_l·Z_l − ·) at second order in J on its two edges. On the populations |a⟩⟨a| the 2^N eigenvalues near 0
/// are −J²·spec(𝓛) + O(J³), and on the far cells |a⟩⟨ā| the 2^N near −2σ are −2σ + J²·spec(𝓜) + O(J³), with γ(m) the
/// sum of the rates over the bits of m and
/// <code>
///   𝓛_ba = −|H_ba|²/γ(a ⊕ b),        𝓛_aa = Σ_b |H_ba|²/γ(a ⊕ b),
///   𝓜_ba = −H_ba·H_āb̄/γ(a ⊕ b),      𝓜_aa = Σ_b (|H_ba|² + |H_b̄ā|²)/(2γ(a ⊕ b)).
/// </code>
/// 𝓛 is the Zeno generator of PROOF_WEIGHT1_DEGENERACY's Theorem E (a) on every coset; 𝓜 is, edge by edge, |h||h′|/γ times a
/// magnetic Laplacian's block carrying Theorem 2's ratio as a phase, plus (|h| − |h′|)²/(2γ) on both cells, h = H_ba and
/// h′ = H_āb̄, the potential present exactly where condition (i) fails. Both are Hermitian and positive semidefinite, and their
/// kernels have Theorem 2's two dimensions, the components of Γ_H and the good components of Γ. The diagonal of H must
/// satisfy Theorem 2's condition (ii), H_aa = H_āā, at every bitstring; otherwise the far edge moves at first order.</summary>
public static class StrongDephasingEdges
{
    /// <summary>𝓛 and 𝓜 from H's columns (<see cref="ComplementConnectionGraph.Columns"/>, bit l = site l) and the rates.</summary>
    public static (ComplexMatrix Near, ComplexMatrix Far) Build(int n, Dictionary<int, GaussianInteger>[] h, IReadOnlyList<double> gamma)
    {
        if (gamma.Count != n) throw new ArgumentException($"{gamma.Count} rates for {n} sites", nameof(gamma));
        if (gamma.Any(g => !(g > 0))) throw new ArgumentException("every rate strictly positive", nameof(gamma));
        int dim = 1 << n, full = dim - 1;
        Complex At(int r, int c) { var e = ComplementConnectionGraph.Entry(h, r, c); return new Complex((double)e.Re, (double)e.Im); }
        double Rate(int m) { double s = 0; for (int l = 0; l < n; l++) if (((m >> l) & 1) == 1) s += gamma[l]; return s; }

        var near = ComplexMatrix.Build.Dense(dim, dim);
        var far = ComplexMatrix.Build.Dense(dim, dim);
        for (int a = 0; a < dim; a++)
        {
            if (ComplementConnectionGraph.Entry(h, a, a) != ComplementConnectionGraph.Entry(h, a ^ full, a ^ full))
                throw new ArgumentException($"H_aa ≠ H_āā at a = {a}: the far edge moves at first order", nameof(h));
            // the neighbours of a in Γ: the entries of column a, and the complements of column ā's
            var neighbours = new HashSet<int>(h[a].Keys);
            foreach (int k in h[a ^ full].Keys) neighbours.Add(k ^ full);
            neighbours.Remove(a);
            foreach (int b in neighbours)
            {
                double s = Rate(a ^ b);
                Complex hba = At(b, a), hfar = At(a ^ full, b ^ full);
                double m2 = hba.Magnitude * hba.Magnitude, f2 = hfar.Magnitude * hfar.Magnitude;
                near[b, a] -= m2 / s;
                near[a, a] += m2 / s;
                far[b, a] -= hba * hfar / s;
                far[a, a] += (m2 + f2) / (2 * s);
            }
        }
        return (near, far);
    }
}
