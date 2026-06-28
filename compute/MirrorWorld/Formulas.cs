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

    // F16 (T1): the fold normal form R = C(Psi+R)^2 = Mandelbrot u->u^2+c, c=C*Psi; boundary at C*Psi=1/4.
    public const double F16_FoldBoundary = 0.25;

    // F25 (T1): CPsi(t) for Bell+ under Z-dephasing. CPsi = f(1+f^2)/6, f = e^{-4 gamma t}; crossing
    // f* = 0.8612 (f(1+f^2)=3/2), K = gamma*t_cross = 0.0374.
    public static double F25_CPsi(double f) => f * (1.0 + f * f) / 6.0;
    public const double F25_CrossingF = 0.8612;
    public const double F25_K = 0.0374;

    // F26 (T1): CPsi for Bell+ under general Pauli noise. CPsi = u(1+u^2+v^2+w^2)/12.
    public static double F26_CPsi(double u, double v, double w) => u * (1.0 + u * u + v * v + w * w) / 12.0;

    // F27 (T1): K per noise channel. K_X=K_Y=ln(2)/8=0.0867, K_Z=0.0374, K_depol=0.0440.
    public static readonly double F27_KX = Math.Log(2) / 8.0;
    public const double F27_KZ = 0.0374;
    public const double F27_KDepol = 0.0440;

    // F15 (T2): theta compass, angular distance from CPsi=1/4. theta = arctan(sqrt(4 C Psi - 1)); 0 at crossing.
    public static double F15_ThetaDeg(double cpsi) => Math.Atan(Math.Sqrt(4.0 * cpsi - 1.0)) * 180.0 / Math.PI;

    // F34 (T1, proven): qubit necessity. d^2 - 2d = 0 -> d = 0 (nothing) or d = 2 (the qubit). Palindromic
    // dephasing needs exactly 2 immune (I,Z) and 2 decaying (X,Y) per site, fixing d=2. The polarity root.
    public static int[] F34_QubitNecessity() => new[] { 0, 2 };

    // D1 (T1, from F2): w=1 bandwidth = omega_{N-1} - omega_1 = 8J cos(pi/N) -> 8J at large N.
    public static double D1_Bandwidth(int n, double j) => 8.0 * j * Math.Cos(Math.PI / n);

    // D4 (T1): the crossing condition scales with Hilbert dimension as (d-1)/2. d=2: f*(1+f*^2)=1/2; d=4: =3/2.
    public static double D4_CrossingRhs(int d) => (d - 1) / 2.0;

    // D6 (T1, AT): spectral gap = 2γ (min nonzero rate); mixing time <= N ln(4)/(2γ).
    public static double D6_Gap(double gamma) => 2.0 * gamma;
    public static double D6_MixingTime(int n, double gamma) => n * Math.Log(4.0) / (2.0 * gamma);

    // F38 (T1): Pi^2 = (-1)^{w_YZ} = (-1)^{n_Y+n_Z} on a Pauli string (order 4, Pi^4=I); = conjugation by X^N.
    public static int F38_PiSquared(int nY, int nZ) => (nY + nZ) % 2 == 0 ? +1 : -1;
}
