namespace RCPsiSquared.Core.Symmetry;

/// <summary>Π²-odd anchor classification for the Dicke superposition
/// (|D_n⟩+|D_{n+1}⟩)/√2 on N qubits. Three structurally distinct cases give
/// three closed-form total Π²-odd Frobenius² values (Tier 1 derived via
/// X⊗N-eigenbasis decomposition; see PROOF_F86B_UNIVERSAL_SHAPE.md §Statement 2):
///
/// <list type="bullet">
///   <item><b>Mirror</b> (N odd, 2n+1=N): ψ is an X⊗N-eigenstate, γ = 1, α_total = 0.</item>
///   <item><b>KIntermediate</b> (N even, n ∈ {N/2−1, N/2}): γ = 1/2, α_total = 3/8.</item>
///   <item><b>Generic</b> (everywhere else): γ = 0, α_total = 1/2.</item>
/// </list>
///
/// Companion to <see cref="F86.Item1Derivation.BondSubClass"/> and
/// <see cref="Resonance.BondClass"/> in the F86 orbit-classification family.
/// </summary>
public enum DickeAnchor
{
    Mirror,
    KIntermediate,
    Generic,
}

public static class DickeAnchorExtensions
{
    /// <summary>Classify the Dicke superposition (|D_n⟩+|D_{n+1}⟩)/√2 on N qubits.</summary>
    public static DickeAnchor Classify(int N, int n)
    {
        if (2 * n + 1 == N) return DickeAnchor.Mirror;
        if (N % 2 == 0 && (n == N / 2 - 1 || n == N / 2)) return DickeAnchor.KIntermediate;
        return DickeAnchor.Generic;
    }

    /// <summary>Total Π²-odd Frobenius² value α_total at the anchor.</summary>
    public static double AlphaTotal(this DickeAnchor anchor) => anchor switch
    {
        DickeAnchor.Mirror => 0.0,
        DickeAnchor.KIntermediate => 3.0 / 8.0,
        DickeAnchor.Generic => 0.5,
        _ => throw new ArgumentOutOfRangeException(nameof(anchor)),
    };

    /// <summary>X⊗N overlap γ = ⟨ψ|X⊗N|ψ⟩ at the anchor. Drives the closed
    /// form α_total = (1 − γ²)/2.</summary>
    public static double Gamma(this DickeAnchor anchor) => anchor switch
    {
        DickeAnchor.Mirror => 1.0,
        DickeAnchor.KIntermediate => 0.5,
        DickeAnchor.Generic => 0.0,
        _ => throw new ArgumentOutOfRangeException(nameof(anchor)),
    };
}
