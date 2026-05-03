using RCPsiSquared.Core.Pauli;

namespace RCPsiSquared.Diagnostics.F87;

/// <summary>The spectral 3-way classification produced by <see cref="PauliPairTrichotomy.Classify"/>:
/// truly (‖M‖ = 0), soft (M ≠ 0 but spectrum still pairs), hard (no pairing).
///
/// <para>The 3-way spectral cut and the algebraic 4-way <see cref="Pi2Class"/> describe the
/// same Hamiltonian from two angles. Use <see cref="Pi2ClassExtensions.PredictedTrichotomy"/>
/// to map the algebraic class to its expected spectral class without an L-build.</para>
/// </summary>
public enum TrichotomyClass
{
    Truly,
    Soft,
    Hard,
}

public static class Pi2ClassExtensions
{
    /// <summary>Mapping from the algebraic 4-way <see cref="Pi2Class"/> to the spectral 3-way
    /// <see cref="TrichotomyClass"/>. Verified empirically across all canonical witnesses
    /// (XX+YY → Truly, XY+YX → Soft, YZ+ZY → Soft, XX+XY → Hard) but not formally proven;
    /// callers needing certainty must invoke <see cref="PauliPairTrichotomy.Classify"/>
    /// against the actual L spectrum.</summary>
    public static TrichotomyClass PredictedTrichotomy(this Pi2Class pi2) => pi2 switch
    {
        Pi2Class.Truly => TrichotomyClass.Truly,
        Pi2Class.Pi2OddPure => TrichotomyClass.Soft,
        Pi2Class.Pi2EvenNonTruly => TrichotomyClass.Soft,
        Pi2Class.Mixed => TrichotomyClass.Hard,
        _ => throw new ArgumentOutOfRangeException(nameof(pi2)),
    };
}
