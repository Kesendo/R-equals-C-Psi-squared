namespace MirrorWorld;

// Closed forms adopted verbatim from the F-registry (docs/ANALYTICAL_FORMULAS.md). Each replaces a
// matrix computation. No interpretation, just the formula and its tier.
public static class Formulas
{
    // F2 (T1, D10): w=1 Liouvillian dispersion, Heisenberg chain. omega_k = 4J(1 - cos(pi k/N)), k=1..N-1.
    public static double[] F2_Dispersion(int n, double j)
    {
        var w = new double[n - 1];
        for (int k = 1; k <= n - 1; k++) w[k - 1] = 4.0 * j * (1.0 - Math.Cos(Math.PI * k / n));
        return w;
    }

    // F2b (T1): XY chain single-excitation spectrum. E_k = 2J cos(pi k/(N+1)), k=1..N.
    public static double[] F2b_SingleExcitation(int n, double j)
    {
        var e = new double[n];
        for (int k = 1; k <= n; k++) e[k - 1] = 2.0 * j * Math.Cos(Math.PI * k / (n + 1));
        return e;
    }

    // F1 residual norms (T1, H-independent, gamma_Z-independent closed forms). sg = Sigma gamma, sg2 = Sigma gamma^2.
    public static double F1_T1Residual(int n, double sg, double sg2) => Math.Pow(4, n - 1) * (3.0 * sg2 + 4.0 * sg * sg);
    public static double F1_DepolResidual(int n, double sg, double sg2) => Math.Pow(4, n - 1) * (16.0 / 9.0 * sg2 + 16.0 * sg * sg);
}
