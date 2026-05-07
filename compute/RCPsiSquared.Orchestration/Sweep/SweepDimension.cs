using RCPsiSquared.Core.Lindblad;

namespace RCPsiSquared.Orchestration.Sweep;

/// <summary>One sweep dimension. Each subtype defines a parameter grid and a target
/// claim type. The F73Frobenius dimension sweeps (N, HamiltonianClass, ChainOnly) and
/// targets PalindromeResidualScalingClaim.</summary>
public abstract record SweepDimension
{
    /// <summary>F73 closed-form sweep. The closed form depends only on N and HamiltonianClass
    /// (chain default); ChainOnly is a flag list that reflects the spec's γ_T1 grid surrogate,
    /// but the closed form is γ-independent so each ChainOnly entry produces the same
    /// numeric output: this proves γ-independence is a structural property of the closed
    /// form rather than a mere numeric coincidence.</summary>
    public sealed record F73Frobenius(
        IReadOnlyList<int> NValues,
        IReadOnlyList<HamiltonianClass> HClasses,
        IReadOnlyList<bool> ChainOnly
    ) : SweepDimension;
}
