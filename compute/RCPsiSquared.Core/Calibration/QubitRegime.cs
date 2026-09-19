namespace RCPsiSquared.Core.Calibration;

/// <summary>Position of the calibration ratio r = T2/(2*T1) relative to the
/// named R* threshold of one normalized-purity proxy.</summary>
public enum Regime
{
    /// <summary>r &lt; R*.</summary>
    BelowRStar,

    /// <summary>|r - R*| &lt;= epsilon; returned only for explicit epsilon &gt; 0.</summary>
    NearRStar,

    /// <summary>r &gt;= R* when no explicit near band captures the value.</summary>
    AtOrAboveRStar,
}

/// <summary>Calibration-level classifier for r = T2/(2*T1).
///
/// <para>R* belongs to the free-single-transmon |+> proxy
/// P_norm(y;r) = 1 - 2y + y^2 + y^(1/r), y = exp(-t/T1). It is the value of r
/// for which that proxy's stationary minimum touches the scalar 1/4. This is
/// not a universal phase label and it is not the radial |c| = 1/4 boundary of
/// the quadratic map.</para>
///
/// <para>IBM calibration rows provide T1 and T2, so this type reports a
/// reproducible normalized-purity proxy reading without changing the measured
/// columns. Operationality remains a separate property of the calibration
/// row.</para></summary>
public static class QubitRegime
{
    /// <summary>Binary64 value of the proxy-touch ratio. Bits:
    /// 0x3fcb3b8c72b39bd9.</summary>
    public const double R_STAR = 0.21275477982200533;

    /// <summary>r = T2/(2*T1). A non-positive T1 maps to positive infinity;
    /// callers must use the separate operationality flag when validity matters.</summary>
    public static double RParam(double t1Us, double t2Us) =>
        t1Us > 0 ? t2Us / (2.0 * t1Us) : double.PositiveInfinity;

    /// <summary>Free-single-transmon |+> normalized-purity proxy
    /// P_norm(y;r) = 1 - 2y + y^2 + y^(1/r).</summary>
    public static double NormalizedPurityProxy(double y, double r) =>
        1.0 - 2.0 * y + y * y + Math.Pow(y, 1.0 / r);

    /// <summary>Derivative of <see cref="NormalizedPurityProxy"/> with respect
    /// to y: -2 + 2y + (1/r)y^(1/r-1).</summary>
    public static double NormalizedPurityProxyStationarity(double y, double r) =>
        -2.0 + 2.0 * y + (1.0 / r) * Math.Pow(y, 1.0 / r - 1.0);

    /// <summary>Classify a T1/T2 pair. With epsilon = 0 the result is binary,
    /// and exact equality belongs to <see cref="Regime.AtOrAboveRStar"/>.
    /// With epsilon &gt; 0, the closed interval |r-R*| &lt;= epsilon is
    /// <see cref="Regime.NearRStar"/>.</summary>
    public static Regime Classify(double t1Us, double t2Us, double epsilon = 0.0)
    {
        if (epsilon < 0)
            throw new ArgumentOutOfRangeException(nameof(epsilon), "epsilon must be non-negative");

        double r = RParam(t1Us, t2Us);
        if (epsilon > 0 && Math.Abs(r - R_STAR) <= epsilon)
            return Regime.NearRStar;
        return r < R_STAR ? Regime.BelowRStar : Regime.AtOrAboveRStar;
    }

    /// <summary>True exactly when r &lt; R*.</summary>
    public static bool IsBelowRStar(double t1Us, double t2Us) =>
        RParam(t1Us, t2Us) < R_STAR;

    /// <summary>Counts the three R* bands on a path.</summary>
    public static (int BelowRStar, int NearRStar, int AtOrAboveRStar) PathComposition(
        IReadOnlyList<QubitData> qubits,
        IReadOnlyList<int> path,
        double epsilon = 0.0)
    {
        var byId = qubits.ToDictionary(q => q.Qubit);
        int below = 0, near = 0, atOrAbove = 0;

        foreach (int qid in path)
        {
            if (!byId.TryGetValue(qid, out var data))
                throw new ArgumentException($"qubit {qid} not in calibration", nameof(path));

            switch (Classify(data.T1Us, data.T2Us, epsilon))
            {
                case Regime.BelowRStar: below++; break;
                case Regime.NearRStar: near++; break;
                case Regime.AtOrAboveRStar: atOrAbove++; break;
            }
        }

        return (below, near, atOrAbove);
    }
}
