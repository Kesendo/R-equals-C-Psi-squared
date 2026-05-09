namespace RCPsiSquared.Core.Numerics;

/// <summary>Shared Shannon entropy helpers in bits (log₂). Used by F75 mirror-pair MI
/// and F76 t-decay envelope; promoted out of those classes when the duplicate was
/// flagged in the 2026-05-09 simplify pass.</summary>
public static class Entropy
{
    private const double Log2 = 0.6931471805599453;

    /// <summary>Binary entropy in bits: <c>h(x) = -x·log₂(x) - (1-x)·log₂(1-x)</c>.
    /// Returns 0 at x = 0 or x = 1 (no log singularity).</summary>
    public static double Binary(double x)
    {
        if (x <= 0.0 || x >= 1.0) return 0.0;
        return -(x * System.Math.Log(x) + (1.0 - x) * System.Math.Log(1.0 - x)) / Log2;
    }

    /// <summary>F76 joint pair entropy in bits for eigenvalues
    /// <c>{1 − 2p, p(1+λ), p(1−λ), 0}</c>:
    /// <c>S_ab(p, λ) = −(1−2p) log₂(1−2p) − p(1+λ) log₂(p(1+λ)) − p(1−λ) log₂(p(1−λ))</c>.
    /// Each term skipped if its eigenvalue is below 10⁻¹⁵.</summary>
    public static double JointPair(double p, double lambda)
    {
        double s = 0.0;
        double v;
        v = 1.0 - 2.0 * p;     if (v > 1e-15) s -= v * System.Math.Log(v);
        v = p * (1.0 + lambda); if (v > 1e-15) s -= v * System.Math.Log(v);
        v = p * (1.0 - lambda); if (v > 1e-15) s -= v * System.Math.Log(v);
        return s / Log2;
    }
}
