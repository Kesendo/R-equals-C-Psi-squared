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

    /// <summary>The TOTAL Pauli-weight histogram (normalized mass per total weight = number of
    /// non-identity Pauli letters) and its mean, for the same d²-length Liouville vector
    /// (vec[a·d+b] = ρ[a,b]). Companion to <see cref="NDiffHistogram"/>, and the load-bearing
    /// distinction behind the survival_incompleteness_mirror "V-Effect identity" question:
    ///
    /// <para>Each density component |i⟩⟨j| factorizes per site, agreeing sites (i_s=j_s) carry
    /// {I,Z}, disagreeing sites carry {X,Y}. So the XY-weight (count of X,Y letters) EQUALS
    /// n_diff = popcount(a⊕b) exactly, while the TOTAL weight = n_diff + (#Z letters on the
    /// agreeing sites). After the tensor-factorized single-site change of basis |i⟩⟨j| → Pauli
    /// coefficients (O(N·4^N), no 4^N×4^N projector), the output index y carries letter ≠ I at a
    /// site iff that site's (ket-bit, bra-bit) ≠ (0,0), so the total weight collapses to
    /// popcount(a|b) and the XY/Z split to popcount(a⊕b)/popcount(a&b). The two histograms
    /// therefore differ EXACTLY by the Z-shadow popcount(a&b): a "dark" {0,2}-coherence survivor
    /// (n_diff ∈ {0,2}) is NOT a half-Pauli-weight (w=N/2) object, because its Z-shadow spreads the
    /// total weight high. (Single-site B: cI=(m00+m11)/2, cX=(m01+m10)/2, cY=i(m01−m10)/2,
    /// cZ=(m00−m11)/2.)</para></summary>
    public static (double MeanWeight, IReadOnlyDictionary<int, double> Histogram) PauliWeightHistogram(
        ComplexVector vec, int n)
    {
        int d = 1 << n;
        long d2 = (long)d * d;
        var c = new Complex[d2];
        for (long x = 0; x < d2; x++) c[x] = vec[(int)x];

        // Apply the single-site Pauli change of basis on each site in place (separable tensor map).
        // The site-p bits live at x-bit (n+p) for the ket-half a and x-bit p for the bra-half b.
        for (int p = 0; p < n; p++)
        {
            long hi = 1L << (n + p);   // ket bit p  (high half of x = a·d + b)
            long lo = 1L << p;         // bra bit p  (low half)
            for (long x = 0; x < d2; x++)
            {
                if ((x & hi) != 0 || (x & lo) != 0) continue;   // each 4-sibling quad once, from its base
                Complex m00 = c[x], m01 = c[x | lo], m10 = c[x | hi], m11 = c[x | hi | lo];
                c[x]           = 0.5 * (m00 + m11);                  // I -> (ket,bra)=(0,0)
                c[x | lo]      = 0.5 * (m01 + m10);                  // X -> (0,1)
                c[x | hi]      = new Complex(0, 0.5) * (m01 - m10);  // Y -> (1,0)
                c[x | hi | lo] = 0.5 * (m00 - m11);                  // Z -> (1,1)
            }
        }

        var weight = new Dictionary<int, double>();
        double total = 0.0;
        for (long y = 0; y < d2; y++)
        {
            double w = c[y].Real * c[y].Real + c[y].Imaginary * c[y].Imaginary;
            if (w < 1e-12) continue;
            int a = (int)(y / d), b = (int)(y % d);
            int wt = System.Numerics.BitOperations.PopCount((uint)(a | b));   // total Pauli weight
            weight[wt] = weight.GetValueOrDefault(wt) + w;
            total += w;
        }
        if (total < 1e-30) return (0.0, weight);
        foreach (int k in weight.Keys.ToList()) weight[k] /= total;
        double mean = weight.Sum(kv => kv.Key * kv.Value);
        return (mean, weight);
    }
}
