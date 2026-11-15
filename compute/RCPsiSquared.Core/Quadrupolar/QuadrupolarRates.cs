using System;
using RCPsiSquared.Core.Numerics;

namespace RCPsiSquared.Core.Quadrupolar;

/// <summary>The spectral densities J₀, J₁, J₂ in the normalization of the claims: the dissipator is
/// Σ_m J_|m|·D[T₂ₘ], so at spin 3/2 R₁ = 6J₁ and 6J₂ and R₂f = 3(J₀ + J₁). The papers read on the data
/// pages write R₂f = J₀ + J₁; <see cref="FromLiterature"/> divides by 3. Signs are not checked: the
/// witnesses' controls use negative densities.</summary>
public sealed record QuadrupolarRates(BigRational J0, BigRational J1, BigRational J2)
{
    public static QuadrupolarRates FromLiterature(BigRational j0, BigRational j1, BigRational j2) =>
        new(j0 / 3, j1 / 3, j2 / 3);

    public BigRational ByOrder(int q) => Math.Abs(q) switch
    {
        0 => J0,
        1 => J1,
        2 => J2,
        _ => throw new ArgumentOutOfRangeException(nameof(q)),
    };
}
