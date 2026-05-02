namespace RCPsiSquared.Core.ChainSystems;

public enum HamiltonianType
{
    /// <summary>H = (J/2) Σ_b (X_b X_{b+1} + Y_b Y_{b+1}). Hopping bilinear, popcount-conserving.</summary>
    XY,

    /// <summary>H = (J/4) Σ_b (X_b X_{b+1} + Y_b Y_{b+1} + Z_b Z_{b+1}). Heisenberg adds Z·Z.</summary>
    Heisenberg,
}
