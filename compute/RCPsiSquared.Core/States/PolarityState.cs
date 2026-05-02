using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using ComplexVector = MathNet.Numerics.LinearAlgebra.Vector<System.Numerics.Complex>;

namespace RCPsiSquared.Core.States;

/// <summary>X-basis polarity states on the +0/−0 axis (THE_POLARITY_LAYER.md).
///
/// |+⟩ = (|0⟩+|1⟩)/√2 carries the +0 polarity; |−⟩ = (|0⟩−|1⟩)/√2 the −0 polarity.
/// Tensor products over N sites — |s_0, s_1, …, s_{N−1}⟩ with s_j ∈ {+, −} — span the
/// polarity sublattice. The d=0 axis lives at the centre, between +0 and −0.
///
/// Big-endian convention (site 0 = leftmost Kronecker factor) matches the rest of Core.
/// </summary>
public static class PolarityState
{
    private static readonly Complex[] XPlus = { 1.0 / Math.Sqrt(2.0), 1.0 / Math.Sqrt(2.0) };
    private static readonly Complex[] XMinus = { 1.0 / Math.Sqrt(2.0), -1.0 / Math.Sqrt(2.0) };

    /// <summary>Uniform polarity (+1 or −1) across all N sites.</summary>
    public static ComplexVector Uniform(int N, int sign) => Build(N, Enumerable.Repeat(sign, N).ToArray());

    /// <summary>Per-site polarity signs (+1 or −1), in site order.</summary>
    public static ComplexVector Build(int N, IReadOnlyList<int> signs)
    {
        if (signs.Count != N)
            throw new ArgumentException($"signs has length {signs.Count}, expected N={N}");
        foreach (var s in signs)
            if (s != 1 && s != -1) throw new ArgumentException($"signs must each be +1 or -1; got {s}");

        Complex[] state = signs[0] == 1 ? XPlus : XMinus;
        for (int i = 1; i < N; i++)
        {
            var next = signs[i] == 1 ? XPlus : XMinus;
            var combined = new Complex[state.Length * 2];
            for (int j = 0; j < state.Length; j++)
                for (int kk = 0; kk < 2; kk++)
                    combined[j * 2 + kk] = state[j] * next[kk];
            state = combined;
        }
        return MathNet.Numerics.LinearAlgebra.Vector<Complex>.Build.DenseOfArray(state);
    }
}
