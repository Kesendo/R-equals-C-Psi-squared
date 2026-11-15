using RCPsiSquared.Core.Numerics;

namespace RCPsiSquared.Core.Quadrupolar;

/// <summary>H = (q/3)·(3I_z² − I(I + 1)) + ω·I_z + K₁(T₂,₁†T₂,₁ − T₂,₋₁†T₂,₋₁)
/// + K₂(T₂,₂†T₂,₂ − T₂,₋₂†T₂,₋₂) + ν·I_x + an optional exact Hermitian matrix in the physical basis.
/// At spin 3/2 the first term is q·diag(1, −1, −1, 1), the witnesses' q = ω_Q·√6/2. The extra matrix
/// carries a field that is a selected part of I_x (the sodium witness's selective lock) and is
/// honoured by the physical views only.</summary>
public sealed record SecularHamiltonian(BigRational Q, BigRational Omega, BigRational K1, BigRational K2,
                                        BigRational Nu, GaussianRational[,]? ExtraPhysical)
{
    public static SecularHamiltonian Of(BigRational q, BigRational omega) =>
        new(q, omega, BigRational.Zero, BigRational.Zero, BigRational.Zero, null);

    public static SecularHamiltonian None => Of(BigRational.Zero, BigRational.Zero);

    public SecularHamiltonian WithLambShift(BigRational k1, BigRational k2) => this with { K1 = k1, K2 = k2 };

    public SecularHamiltonian WithLock(BigRational nu) => this with { Nu = nu };

    public SecularHamiltonian WithExtraPhysical(GaussianRational[,] h) => this with { ExtraPhysical = h };
}
