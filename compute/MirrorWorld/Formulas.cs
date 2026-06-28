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

    // F2b corollary (T1): the coherence hand. omega_mem = 2J cos(pi/(N+1)) for N>=3 (sqrt2, phi, sqrt3
    // at N=3,4,5; gamma-independent); 2 sqrt(J^2 - gamma^2) at N=2 (-> 0 at the EP Q=1).
    public static double OmegaMem(int n, double j, double gamma) =>
        n >= 3 ? 2.0 * j * Math.Cos(Math.PI / (n + 1)) : 2.0 * Math.Sqrt(Math.Max(0.0, j * j - gamma * gamma));

    // Coherence horizon Q*(N) (T1): exact values. N=2,3 clean (1, sqrt2); N>=4 transcendental SE-EP
    // (1.8787, 2.3737); asymptotic slope exactly 2/pi (Q*(N) -> 2N/pi).
    public static double Qstar(int n) => n switch
    {
        2 => 1.0,
        3 => Math.Sqrt(2.0),
        4 => 1.8787,
        5 => 2.3737,
        _ => 2.0 * n / Math.PI,
    };

    // F3 (T1, AT corollary): decay rate bounds. min=2γ (w=1), max=2(N-1)γ (w=N-1), bw=2(N-2)γ.
    public static (double Min, double Max, double Bw) F3_RateBounds(int n, double gamma) =>
        (2.0 * gamma, 2.0 * (n - 1) * gamma, 2.0 * (n - 2) * gamma);

    // F4 (T1): kernel dim = N+1 for one connected component (identity + N magnetization projectors).
    public static int F4_KernelDim(int n) => n + 1;

    // F5 (T1): depolarizing palindrome error = (2/3)Σγ = γ·2N/3.
    public static double F5_DepolError(int n, double gamma) => gamma * 2.0 * n / 3.0;

    // F23 (T1): XOR-drain fraction = (N+1)/4^N (GHZ fragility vanishes at large N).
    public static double F23_XorFraction(int n) => (n + 1.0) / Math.Pow(4, n);

    // F33 (T1): N=3 exact decay rates {2γ, 8γ/3, 10γ/3} (<n_XY> = 1, 4/3, 5/3).
    public static double[] F33_N3Rates(double gamma) => new[] { 2.0 * gamma, 8.0 * gamma / 3.0, 10.0 * gamma / 3.0 };

    // F50 (T1 lower bound): weight-1 degeneracy d_real(-2γ) = 2N (chain); 8 for the K_3 triangle (N=3).
    public static int F50_Weight1Degeneracy(int n, bool triangleK3 = false) => triangleK3 && n == 3 ? 8 : 2 * n;

    // F7 (T1): w=1 Q-factor spectrum. Q_max=2J/g(1+cos pi/N), Q_min=2J/g(1-cos pi/N), Q_mean=2J/g,
    // Q_spread = Q_max/Q_min = cot^2(pi/(2N)).
    public static (double Max, double Min, double Mean, double Spread) F7_QSpectrum(int n, double j, double gamma)
    {
        double s = 2.0 * j / gamma, c = Math.Cos(Math.PI / n), cot = 1.0 / Math.Tan(Math.PI / (2.0 * n));
        return (s * (1 + c), s * (1 - c), s, cot * cot);
    }

    // F8 (T1): 2x universal decay law. unpaired = 2Ng (<n_XY>=N), paired mean = Ng, ratio = 2 exactly.
    public static (double Unpaired, double PairedMean) F8_DecayLaw(int n, double gamma) => (2.0 * n * gamma, n * gamma);

    // F12 (T2): single-qubit universal crossing fraction t*/T2 = 0.858367, the root of x^3 + x = 1/2.
    public const double F12_CrossingFraction = 0.858367;
}
