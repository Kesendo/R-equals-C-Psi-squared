namespace RCPsiSquared.Core.Symmetry;

/// <summary>Classification of a Pi²-Inheritance Claim by which Z₂ axis (or axes) of
/// the cubic Z₂³ polarity structure it lives on.
///
/// Per F34/QUBIT_NECESSITY the operator algebra of a single qubit admits at most
/// two independent Π² operators (Π²_Z = X⊗N and Π²_X = Z⊗N); the third Z₂
/// classifier (Y-parity = #Y mod 2) lives at the Pauli-term level and becomes
/// independent only at k≥3-body terms. So the cubic Z₂³ is term-level, not
/// operator-level.</summary>
public enum Z2Axis
{
    /// <summary>Π²_Z = X⊗N axis. bit_b parity (n_Y + n_Z mod 2). The F1² family.</summary>
    BitB,

    /// <summary>Π²_X = Z⊗N axis. bit_a parity (n_X + n_Y mod 2). The F61 family.</summary>
    BitA,

    /// <summary>Uses both Π²_Z and Π²_X axes actively (Klein-Vierergruppe Z₂ × Z₂).</summary>
    Klein2,

    /// <summary>Independent term-level Y-parity refinement (relevant only at k≥3-body).</summary>
    YParity,

    /// <summary>Uses all three Z₂ classifiers (full Z₂³ = 8 sectors at k≥3).</summary>
    Cubic3,

    /// <summary>Utility / non-Π²-axis Claim (e.g., chiral K, F71 spatial mirror, dyadic ladder anchors).</summary>
    NotApplicable
}
