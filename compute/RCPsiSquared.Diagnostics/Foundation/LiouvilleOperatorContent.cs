using System.Collections.Generic;
using System.Linq;
using System.Numerics;
using ComplexVector = MathNet.Numerics.LinearAlgebra.Vector<System.Numerics.Complex>;

namespace RCPsiSquared.Diagnostics.Foundation;

/// <summary>Operator content of a Liouville-space vector. Reshape vec[a·d+b] → ρ[a,b] (the
/// <c>PauliDephasingDissipator.BuildZ</c> convention <see cref="Symphony"/> and
/// <c>ChainSystem.BuildLiouvillian</c> use), then bucket its weight by n_diff = popcount(a⊕b).
/// Histogram {0:½, 2:½} marks a population/antisymmetric coherence (the freezing block of the
/// coherence horizon); {1:1} marks a |vac⟩⟨ψ_k| band coherence (the γ-protected survivor).</summary>
public static class LiouvilleOperatorContent
{
    /// <summary>The n_diff histogram (normalized weight per popcount(a⊕b)) and its mean, for an
    /// eigenvector of the d²×d² Liouvillian at N qubits (d = 2^N).</summary>
    public static (double MeanNDiff, IReadOnlyDictionary<int, double> Histogram) NDiffHistogram(
        ComplexVector vec, int n)
    {
        int d = 1 << n;
        var weight = new Dictionary<int, double>();
        double total = 0.0;
        for (int a = 0; a < d; a++)
        {
            for (int b = 0; b < d; b++)
            {
                Complex c = vec[a * d + b];
                double w = c.Real * c.Real + c.Imaginary * c.Imaginary;
                if (w < 1e-12) continue;
                int nd = System.Numerics.BitOperations.PopCount((uint)(a ^ b));
                weight[nd] = weight.GetValueOrDefault(nd) + w;
                total += w;
            }
        }
        if (total < 1e-30) return (0.0, weight);
        foreach (int k in weight.Keys.ToList()) weight[k] /= total;
        double mean = weight.Sum(kv => kv.Key * kv.Value);
        return (mean, weight);
    }
}
