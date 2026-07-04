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

    // F18 (T2, N-independent): fold threshold Sg_crit/J. Below: CPsi oscillates forever; above: crosses
    // 1/4 irreversibly. Bell 0.00249, product 0.00497 (max/min across N=2-5: 1.5%).
    public const double F18_FoldThresholdBell = 0.00249;
    public const double F18_FoldThresholdProduct = 0.00497;

    // F36/F37 (T1, neural): the Wilson-Cowan palindrome Q*J*Q + J + 2S = 0 (structural analog of
    // Pi*L*Pi^-1 = -L - 2 Sg), eigenvalue pairing mu_k + mu_k' = -(1/tau_E + 1/tau_I). C.elegans
    // connectome residual 0.013 vs random 0.108 (8x more palindromic than chance).
    public static double F37_NeuralPairSum(double tauE, double tauI) => -(1.0 / tauE + 1.0 / tauI);

    // F61 (T1): the bit_a parity Pi^2_X = Z^{tensor N} (the global Z-string), companion to F38's X^N.
    // On a Pauli string it is (-1)^{n_X+n_Y} = (-1)^k, the disagreement-count parity the Hamiltonian
    // conserves in its even-step mixing. F63: [L, Pi^2] = 0 (Pi^2 a conserved quantum number, all N).
    public static int F61_PiSquaredX(int nX, int nY) => (nX + nY) % 2 == 0 ? +1 : -1;

    // F39 (T1): det(Pi) = (-1)^{N 4^{N-1}}. N=1: -1; N>=2: +1 (4^{N-1} is even).
    public static int F39_DetPi(int n) => n == 1 ? -1 : +1;

    // F41 (T1, D10): palindromic time t_Pi = 2pi/omega_min = pi/(4J sin^2(pi/(2N))) ~ N^2/(pi J).
    public static double F41_PalindromicTime(int n, double j) => Math.PI / (4.0 * j * Math.Pow(Math.Sin(Math.PI / (2.0 * n)), 2));

    // F44 (T1, D08): Crooks-like rate identity ln(d_fast/d_slow) = 2 artanh(Delta_d/(2 Sg)) for a
    // palindromic pair d_fast + d_slow = 2 Sg (algebraic, NOT a thermodynamic Crooks theorem).
    public static double F44_LogRatio(double dFast, double dSlow, double sg) => 2.0 * Math.Atanh((dFast - dSlow) / (2.0 * sg));

    // F49 (T1, proven): cross-term ratio R(N) = sqrt((N-2)/(N 4^{N-1})). N=2: 0 (exact Pythagorean);
    // N=3: 1/sqrt48; N=4: 1/sqrt128. gamma/J/topology-independent, depends only on N.
    public static double F49_CrossTerm(int n) => Math.Sqrt((n - 2.0) / (n * Math.Pow(4, n - 1)));

    // F49b (T1, proven): centered dissipator norm ||L_Dc||^2 = gamma^2 4^N N (uniform Z-dephasing).
    public static double F49b_CenteredDissipatorNormSq(int n, double gamma) => gamma * gamma * Math.Pow(4, n) * n;

    // F49c (T1, proven): cross-term for shadow-crossing couplings (one bond Pauli in {X,Y}, the other in
    // {I,Z}): R(N) = sqrt((N-1)/(N 4^(N-1))). Companion to F49 (bond-site variance 1 not 0, so N-2 -> N-1).
    public static double F49c_CrossTermCrossing(int n) => Math.Sqrt((n - 1.0) / (n * Math.Pow(4, n - 1)));

    // F55 (T1, from D6): universal absorption dose K_death = ln(10) = 2.303 (99% absorption of the slowest
    // mortal mode, rate 2gamma). Immortal modes = N+1 (pure {I,Z}, zero absorption, invisible to the light).
    public const double F55_KDeath = 2.302585092994046;   // ln(10)
    public static int F55_ImmortalModes(int n) => n + 1;

    // F56 (T1, closed form): critical-slowing iteration count near the cardioid cusp c = 1/4 - eps.
    // K(eps,tol) = (1/2) ln(4 eps/tol) + alpha(tol) sqrt(eps), alpha(tol) = -4 + (1/2) ln(16 tol).
    public static double F56_CriticalSlowing(double eps, double tol) =>
        0.5 * Math.Log(4.0 * eps / tol) + (-4.0 + 0.5 * Math.Log(16.0 * tol)) * Math.Sqrt(eps);

    // F57 (T1): trajectory dwell at CPsi=1/4, t_dwell = 2 delta/|dCPsi/dt|. Bell+ K_dwell/delta =
    // 1.080088 = 2/1.851701 (gamma-independent to machine precision).
    public const double F57_DwellPrefactorBell = 1.080088;

    // F59 (T1): generalized dwell prefactor for a two-sector state (stationary W_0, coherent W_k at
    // XY-weight k): (4/k)(W_0 + W_k)/(W_0 + 3 W_k). Independent of N, d. Reduces to F57 at W_0=1/2,k=2.
    public static double F59_DwellPrefactor(int k, double w0, double wk) => (4.0 / k) * (w0 + wk) / (w0 + 3.0 * wk);

    // F60 (T1): GHZ_N born below the fold. CPsi(0) = 1/(2^N - 1); < 1/4 for all N >= 3 (gamma-independent).
    public static double F60_GhzCPsi0(int n) => 1.0 / ((1 << n) - 1);

    // F62 (T1): W_N initial CPsi. CPsi(0) = 2(N^2 - 4N + 8)/(3 N^3).
    public static double F62_WstateCPsi0(int n) => 2.0 * (n * n - 4.0 * n + 8.0) / (3.0 * n * n * n);

    // F63 (T1, proven): [L, Pi^2] = 0. With F61 (n_XY parity), L has two independent Z2 symmetries; the
    // d=2 Pauli algebra splits the operator space into 4 blocks of dim 4^(N-1). Per Pi^2-sector conserved
    // mode count (boundary Z-dephasing): even = floor(N/2)+1, odd = ceil(N/2) (the e_d(Z) by parity).
    public static long F63_BlockDim(int n) => 1L << (2 * (n - 1));   // 4^(N-1)
    public static (int Even, int Odd) F63_ConservedPerSector(int n) => (n / 2 + 1, (n + 1) / 2);

    // F65 (T1, proven): single-excitation dissipation spectrum, uniform open XY chain, endpoint Z-dephasing.
    // alpha_k/gamma0 = (4/(N+1)) sin^2(k pi/(N+1)), k=1..N. All in [0, 2]; alpha_k = alpha_{N+1-k}.
    public static double[] F65_SingleExcitationRates(int n)
    {
        var a = new double[n];
        for (int k = 1; k <= n; k++) a[k - 1] = 4.0 / (n + 1) * Math.Pow(Math.Sin(k * Math.PI / (n + 1)), 2);
        return a;
    }

    // F65 Niven rationality: all alpha_k/gamma0 rational iff N+1 in {1,2,3,4,6}, i.e. N in {0,1,2,3,5}
    // (N=4 is the first golden-irrational, sin^2(pi/5)). Niven's theorem on cos(2k pi/(N+1)).
    public static bool F65_RatesRational(int n) => n is 0 or 1 or 2 or 3 or 5;

    // F66 (T1): the dissipation interval [0, 2 gamma0] has poles at both endpoints, multiplicity N+1 each
    // (the N+1 elementary symmetric polynomials e_d(Z) at alpha=0, their Pi-partners at alpha=2 gamma0).
    public static int F66_PoleMultiplicity(int n) => n + 1;

    // F68 (T1): palindromic partner rate of the bonding mode, alpha_p = 2 gamma0 - alpha_b (from F1).
    public static double F68_PartnerRate(double alphaB, double gamma0) => 2.0 * gamma0 - alphaB;

    // F69 (T1): GHZ_3 + W_3 sector mix lifts pair-CPsi(0) above the fold; the optimum is a degree-6
    // algebraic number (irreducible sextic), pair-CPsi = 0.320412, ratio 1.281646 to 1/4. N>=4: stays below.
    public const double F69_N3Optimum = 0.320411541127025;
    public const double F69_RatioToQuarter = 1.281646;

    // F70 (T1, kinematic): a k-local partial trace annihilates coherence blocks with |Delta N| >= k+1, so a
    // k-local observable sees only |Delta N| <= k content (single-site <= 1, pair <= 2).
    public static int F70_MaxVisibleDeltaN(int kLocal) => kLocal;

    // F71 (T1, kinematic): the closure-breaking c1 bond profile is mirror-symmetric c1(b)=c1(N-2-b), so it
    // has ceil((N-1)/2) = floor(N/2) independent components; the SE reflection R|psi_k> = (-1)^(k+1)|psi_k>.
    public static int F71_C1IndependentComponents(int n) => n / 2;
    public static int F71_ReflectionParity(int k) => k % 2 == 1 ? +1 : -1;

    // F98 (T1): K-intermediate Dicke long-time Pi^2-odd asymptote (even N): (N+2)/(4(N+1)) -> 1/4.
    public static double F98_DickeAsymptote(int n) => (n + 2.0) / (4.0 * (n + 1.0));

    // F121 (T1, combinatorial): the qudit partial palindrome. Under full-Cartan dephasing |i><j| decays at
    // -2gamma*Hamming(i,j) (same ladder as the qubit); multiplicity c_k = d^N C(N,k) (d-1)^k, Sum=d^(2N).
    // The dissipator pairs rung k <-> N-k; paired(d,N) = Sum_k d^N C(N,k) (d-1)^(min(k,N-k)), = d^(2N)
    // iff d=2 (the d^2-2d=0 necessity re-seen). d=3,N=2: c=[9,36,36], paired=54/81.
    public static long F121_CoherenceCount(int d, int n, int k) => IntPow(d, n) * Block.Binomial(n, k) * IntPow(d - 1, k);
    public static long F121_PairedCeiling(int d, int n)
    {
        long s = 0;
        for (int k = 0; k <= n; k++) s += IntPow(d, n) * Block.Binomial(n, k) * IntPow(d - 1, Math.Min(k, n - k));
        return s;
    }

    // F122 (T1): the structural ceiling g2 = <n_XY> of the slowest mode. Chain g2=1 (band edge); more
    // connected graphs grow a darker [H,A]=0 coherence: g2(K_N)=4/N (N>=5), g2(star_N)=4/(N-1) (N>=6),
    // g2(K_4)=2-2/sqrt3. The ring has no ceiling (g2=1); its (1,1) commutant = 2(N-2)/N (even), 2(N-1)/N (odd).
    public static double F122_CompleteCeiling(int n) => 4.0 / n;
    public static double F122_StarCeiling(int n) => 4.0 / (n - 1);
    public static double F122_K4Ceiling() => 2.0 - 2.0 / Math.Sqrt(3.0);
    public static double F122_RingCommutant(int n) => n % 2 == 0 ? 2.0 * (n - 2) / n : 2.0 * (n - 1) / n;

    // F85 (T1): the k-body generalization of the F49 Frobenius scaling. Per k-body Pauli term:
    // truly iff #Y even AND #Z even (contributes M = 0); else the Pi^2-class decides the factor,
    // c(Pi^2-odd) = 1, c(Pi^2-even non-truly) = 2, and ||M||^2_F per term = 4 c ||H_k||^2_F 2^N.
    // For k >= 3, n_YZ is NOT the determining quantity (YYY has n_YZ=3 but c=1); only the class is.
    // Pi^2-odd count over {X,Y,Z}^k: (3^k - (-1)^k)/2.
    public static int F85_FrobeniusFactor(string letters)
    {
        int nY = letters.Count(c => c == 'Y'), nZ = letters.Count(c => c == 'Z');
        if (nY % 2 == 0 && nZ % 2 == 0) return 0;                   // truly
        return (nY + nZ) % 2 == 1 ? 1 : 2;                          // Pi^2-odd : Pi^2-even non-truly
    }
    public static long F85_Pi2OddCount(int k) => (IntPow(3, k) - (k % 2 == 0 ? 1 : -1)) / 2;
    public static double F85_ResidualNormSqPerTerm(int n, double hNormSq, int c) => 4.0 * c * hNormSq * (1L << n);

    // F97 (T1): the Mandelbrot main cardioid at the framework anchor b = 1/2. The period-1 fixed
    // point of z^2 + c sits at magnitude exactly b on the marginally-stable boundary:
    // z*(phi) = b e^{i phi}, c(phi) = b e^{i phi} - b^2 e^{2i phi} = z*(1 - z*). Cusp c(0) = 1/4
    // (the F16 fold boundary re-seen), tail c(pi) = -3/4, top c(pi/2) = 1/4 + i/2.
    public static System.Numerics.Complex F97_FixedPoint(double phi)
        => 0.5 * System.Numerics.Complex.Exp(System.Numerics.Complex.ImaginaryOne * phi);
    public static System.Numerics.Complex F97_Cardioid(double phi)
        => 0.5 * System.Numerics.Complex.Exp(System.Numerics.Complex.ImaginaryOne * phi)
         - 0.25 * System.Numerics.Complex.Exp(System.Numerics.Complex.ImaginaryOne * (2.0 * phi));

    // F124 (T1): the band-edge transition invariant. For the open chain's full single-excitation
    // bond-transition matrix M (band-edge carrier), ||M||_F^2 + lambda_min(M M^T) = 2 (the
    // coordination number), split as (2 - E) + E where E = (4/(N+1)) sin^2(pi/(N+1)) is the
    // carrier's weight on the two free ends -- exactly the k=1 rung of the F65 ladder. The floor
    // vanishes as E (N+1)^3 -> 4 pi^2 (the resolution limit reading).
    public static double F124_EndWeight(int n) => 4.0 / (n + 1) * Math.Pow(Math.Sin(Math.PI / (n + 1)), 2);
    public static double F124_FrobeniusNormSq(int n) => 2.0 - F124_EndWeight(n);
    public static double F124_SpectralFloor(int n) => F124_EndWeight(n);

    private static long IntPow(int b, int e) { long r = 1; for (int i = 0; i < e; i++) r *= b; return r; }
}
