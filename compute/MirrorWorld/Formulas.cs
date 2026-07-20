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

    // Coherence horizon Q*(N) (T1): the single-excitation EP (coherence_horizon_se_block.py, qstar_se;
    // PROOF_COHERENCE_HORIZON_SLOPE). N=2,3 clean closed forms (1, sqrt2); N=4..8 the transcendental
    // SE-EP values, each strictly BELOW the 2N/pi asymptote it approaches from below. Only the slope
    // 2/pi is the clean N->inf limit, so N>=9 falls back to the asymptote -- a documented approximation
    // that OVERSHOOTS the true finite-N EP (still ~+29% at N=8), so extend the exact table if a larger
    // N matters. (The N=6,7,8 exact values landed 2026-07-12 after an empty review caught N>=6 silently
    // using the asymptote, ~30% high, as if exact.)
    public static double Qstar(int n) => n switch
    {
        2 => 1.0,
        3 => Math.Sqrt(2.0),
        4 => 1.8787,
        5 => 2.3737,
        6 => 2.889253,
        7 => 3.419782,
        8 => 3.961618,
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

    // F33 (T1): N=3 rate ladder. Rungs {0, 2γ, 4γ, 6γ} are exact at every J; the two values
    // below are the J/γ → ∞ limit of bands that split at finite coupling, NOT exact rationals (<n_XY> = 1, 4/3, 5/3).
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

    // D6 (T1 above the coupling threshold): spectral gap = 2γ (min nonzero rate);
    // mixing time <= N ln(4)/(2γ). PRECONDITION for both: Q = J/γ above Q*_gap(N), which is
    // 0.50, 0.80, 1.34, 1.82 at N = 2..5 on the Heisenberg chain (Pauli-J units). Below
    // that the gap is Zeno-suppressed, these return values that are too large (gap) and
    // too small (mixing time). The signatures take no J, so the caller carries the
    // precondition.
    public static double D6_Gap(double gamma) => 2.0 * gamma;
    public static double D6_MixingTime(int n, double gamma) => n * Math.Log(4.0) / (2.0 * gamma);

    // NO predicate is offered here on purpose. The thresholds are known in PAULI-J units
    // (Heisenberg chain: 0.500000, 0.800243, 1.342243, 1.819350 at N = 2..5; XY chain:
    // 0.500000, 0.707107, 0.939271, 1.186087), while this module's own coupling is the
    // hopping amplitude, J_MirrorWorld = 2 * J_Pauli (see Restless: "zz=1 with the hopping
    // J=2 is the isotropic Heisenberg bond"), and Restless defaults to zz=0, i.e. XY. A
    // predicate taking this module's J and comparing it against the Pauli-J Heisenberg
    // table answers wrongly, and wrongly in the unsafe direction: at the canonical regime
    // (gamma=0.05, J_Pauli=0.075, so J_MirrorWorld=0.15) it would report the gap law
    // holding at N=5, where the true gap is 1.2754*gamma rather than 2*gamma.
    // Convert to Pauli-J and pick the right family before comparing.

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

    // F55 (T1, from D6): absorption dose K_death = ln(10) = 2.303 (99% absorption of the slowest
    // mortal mode, rate 2gamma). The rate 2gamma is D6's, and holds only above the coupling
    // threshold Q*_gap(N) in Q = J/gamma; below it the slowest mortal mode is Zeno-suppressed and the dose
    // is larger. Immortal modes = N+1 (zero absorption, invisible to the light).
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

    // F74 (T1, combinatorial; adopted 2026-07-12): chromaticity of the (n, n+1) popcount coherence
    // block at J = 0 -- exactly min(n, N-1-n) + 1 distinct pure dephasing rates, the odd ladder
    // 2*gamma0*{1, 3, ..., 2c-1}. A pair (x, y) with popcounts (n, n+1) differs at HD = 2n+1-2*match
    // sites (match = ones shared); match ranges over [max(0, 2n+1-N), n], so HD runs the odd values
    // 1..min(2n+1, 2N-2n-1) and the Pair rate -2*gamma0*HD takes exactly c(n, N) values. The block
    // generalization of Pair's single-coherence rate: mono-chromatic at the ends, maximal at centre.
    public static int F74_Chromaticity(int n, int nSites) => Math.Min(n, nSites - 1 - n) + 1;
    public static double[] F74_RateLadder(int n, int nSites, double gamma0)
        => Enumerable.Range(0, F74_Chromaticity(n, nSites)).Select(i => 2.0 * gamma0 * (2 * i + 1)).ToArray();

    // F94 + F96 (T1 derived, bit-exact Dyson; adopted 2026-07-12 together): the Born-deviation
    // table of ONE setup -- |0+0+> on the N = 4 Heisenberg ring (H = (J/4) sum XX+YY+ZZ) under
    // uniform Z-dephasing, pair (0, 2). SETUP-SPECIFIC scalars, not an N-family law. The single
    // anchor is F94's 4/3 = (sym3 element 8) / 3!, a pure COUNT (32 surviving Dyson diagrams,
    // each 1/4, no cancellations): the dominant outcome |00> pays Delta = +(4/3) Q^2 K^3
    // (Q = J/gamma, K = gamma t); F96 makes the three subdominant outcomes algebra in the same
    // anchor via the universal slope M_{2k+1} / ((2k+1) U_{2k}): |01> = |10> = -(4/3)^2 K
    // (= (-4)/(3 * 3/4)) and |11> = -2 (4/3) K (= (-20)/(5 * 3/2)), linear in K, Q-independent.
    // The qudit lift c(d) = 4(d+2)(d-1)/(3d^2) (the dynamics is the d-independent (J/2) SWAP)
    // refutes the "4 = d^2" reading: c(2) = 4/3 is both the qubit value and the d -> inf limit,
    // with the finite-d bump peaking at c(4) = 3/2.
    public static double F94_DominantDeviation(double q, double k) => 4.0 / 3.0 * q * q * k * k * k;
    public static double F94_QuditCoefficient(int d) => 4.0 * (d + 2) * (d - 1) / (3.0 * d * d);
    public static double F96_SlopeFromDyson(double m, int k, double u) => m / ((2 * k + 1) * u);
    public static double F96_SubdominantSlopeSingle() => -(4.0 / 3.0) * (4.0 / 3.0);   // |01>, |10>
    public static double F96_SubdominantSlopeDouble() => -2.0 * (4.0 / 3.0);           // |11>

    // F72 (T1, corollary of F70; adopted 2026-07-12): the block-diagonal DD + CC split of per-site
    // purity. Valid for any Hamiltonian conserving excitation number, any sector-preserving
    // dissipator, any rho0; purely kinematic. Tr(rho_i^2) = 1/2
    // + P_DD[rho_diag] + P_CC[rho_coh] with NO cross term: the site marginal's z Bloch component
    // reads only the DeltaN = 0 (diagonal) blocks of rho, its x and y read only the |DeltaN| = 1
    // blocks (F70: a single site sees at most one excitation step), and squaring keeps each
    // contribution in its own sector. Blocks with |DeltaN| >= 2 are invisible to any single site.
    // Generalizes at k = 2 sites to three sub-blocks DD + DC + CC.
    public static double F72_PurityDD(double z) => 0.5 * z * z;
    public static double F72_PurityCC(double x, double y) => 0.5 * (x * x + y * y);
    public static double F72_SitePurity(double x, double y, double z)
        => 0.5 + F72_PurityDD(z) + F72_PurityCC(x, y);

    // F73 (T1, proven; adopted 2026-07-12): the spatial-sum coherence closure. Any Hermitian H with
    // [H, N_total] = 0 and uniform Z-dephasing gamma0, on the vac-SE coherent probe
    // rho0 = (|vac><alpha| + |alpha><vac|)/2 (|alpha> any normalized single-excitation state), obeys
    // sum_i 2*|(rho_i)_{01}(t)|^2 = (1/2) e^{-4 gamma0 t} EXACTLY -- blind to H. The (vac, SE) block
    // is F74's n = 0 mono-chromatic block (c(0, N) = 1, the single Pair rate -2*gamma0): the block
    // vector evolves as x(t) = e^{-2 gamma0 t} U_SE(t) x(0) with U_SE unitary, so the summed squared
    // magnitude pays twice the rate and the H-rotation drops out of the norm. Breaks for non-uniform
    // gamma_l, non-U(1) H, dissipators that shift the d_H = 1 rate (mixed X/Z, amplitude damping),
    // and probes with two-or-more-excitation admixture (the rate leaves 2*gamma0).
    public static double F73_SpatialSumClosure(double gamma0, double t) => 0.5 * Math.Exp(-4.0 * gamma0 * t);

    // F75 (T1): mirror-pair mutual information for a single-excitation mirror-symmetric state
    // |psi> = sum c_j |1_j> with c_{N-1-j} = +-c_j: MI(l, N-1-l) = 2 h(p) - h(2p), p = |c_l|^2,
    // h the binary entropy. Sign-independent; saturates at 2 bits (a Bell pair) at p = 1/2.
    // Bonding:k populations are the F65 amplitudes squared, p_l = (2/(N+1)) sin^2(pi k (l+1)/(N+1));
    // the mirror-pair sum MM(0) over pairs l = 0..floor(N/2)-1 is O(N), no propagation.
    public static double F75_MirrorPairMI(double p) => 2.0 * H2(p) - H2(2.0 * p);
    public static double F75_BondingSitePopulation(int n, int k, int l)
        => 2.0 / (n + 1) * Math.Pow(Math.Sin(Math.PI * k * (l + 1) / (n + 1)), 2);
    public static double F75_MirrorPairSum(int n, int k)
    {
        double mm = 0;
        for (int l = 0; l < n / 2; l++) mm += F75_MirrorPairMI(F75_BondingSitePopulation(n, k, l));
        return mm;
    }

    // F77 (T1, asymptotic proven; adopted 2026-07-12): the multi-drop MM(0) saturates at 1 bit --
    // MM(0)(N, k*) = 1 + 3/(4(N+1) ln 2) + O(N^-2). Per-pair information ~ 4/(N+1) shrinks exactly as
    // the ~ N/2 pair count grows (the probability normalisation matches the two scalings); the
    // correction is the entropy non-linearity via sum sin^4 = 3(N+1)/8 at generic k. Sits ON F75:
    // the from-below pin is convergence of the exact F75 mirror-pair sum, (MM-1)(N+1) -> 3/(4 ln 2)
    // = 1.0820. The resonant k = (N+1)/2 carries the enhanced 1/((N+1) ln 2) deviation, rescaled
    // 1/ln 2 = 1.4427 (the registry's 1.445; its symbolic line carried a stray *2, fixed at this
    // adoption); isolated, density zero in the limit.
    public static double F77_MMSaturation(int n) => 1.0 + 3.0 / (4.0 * (n + 1) * Math.Log(2.0));
    public static double F77_RescaledDeviationLimit() => 3.0 / (4.0 * Math.Log(2.0));

    // F76 (T1): pure-dephasing decay of the mirror-pair MI. The pair coherence decays at 4 gamma0
    // (lambda = e^{-4 gamma0 t}) while the populations stay, so the pair eigenvalues become
    // {1-2p, p(1+lambda), p(1-lambda), 0} and MI(p, t) = 2 h(p) - S_ab(p, lambda). lambda = 1
    // recovers F75 (S_ab = h(2p)); lambda = 0 gives S_ab = h(1-2p) + 2p. The 0.93 envelope at
    // gamma0 = 0.05, t = 0.1 is the gamma0 signature, not a hidden constant (0.964 at gamma0 =
    // 0.025, 0.888 at 0.10, both at (5,2)); Heisenberg mixing is second-order small (< 0.5%).
    public static double F76_PairEntropy(double p, double lambda)
        => -XLog2(1.0 - 2.0 * p) - XLog2(p * (1.0 + lambda)) - XLog2(p * (1.0 - lambda));
    public static double F76_MirrorPairMI(double p, double lambda) => 2.0 * H2(p) - F76_PairEntropy(p, lambda);
    public static double F76_Envelope(int n, int k, double gamma0, double t)
    {
        double lambda = Math.Exp(-4.0 * gamma0 * t);
        double now = 0, start = 0;
        for (int l = 0; l < n / 2; l++)
        {
            double p = F75_BondingSitePopulation(n, k, l);
            now += F76_MirrorPairMI(p, lambda);
            start += F75_MirrorPairMI(p);
        }
        return now / start;
    }

    // F95 (T1): the theta-compass at the quadratic discriminant zero. For z^2 - 2bz + c = 0 the
    // complex-root angle above the threshold c = b^2 is theta = arctan(sqrt(c/b^2 - 1)); zero at
    // the degenerate double root, undefined (NaN) below it. At b = 1/2 the threshold is 1/4 and
    // the Februar compass arctan(sqrt(4c - 1)) (= F15) is recovered; the Lindblad specialization
    // (lambda^2 + 2 gamma lambda + gamma^2 + J^2, b = -gamma, c = gamma^2 + J^2) gives
    // theta = arctan(J/gamma) = arctan(Q), the Clock's angle: the compass and the clock are one.
    public static double F95_Theta(double c, double b)
        => c < b * b ? double.NaN : Math.Atan(Math.Sqrt(c / (b * b) - 1.0));
    public static double F95_ThetaHalf(double c) => F95_Theta(c, 0.5);

    // F99 (T1): the five canonical trigonometric anchors. The F86b alpha-formula
    // alpha(theta) = sin^2(theta)/2 at {0, 30, 45, 60, 90} degrees produces the five Pi2 dyadic
    // anchors {0, 1/8, 1/4, 3/8, 1/2}; the non-uniform Dicke weight realizing gamma = cos(theta)
    // is c^2 = cos(theta)/(2 sin^2(theta/2)) (2 sqrt3 + 3 at 30, the silver ratio 1 + sqrt2 at 45,
    // the uniform Dicke c = 1 at 60, 0 at 90). The standard trig triangles ARE the anchor triangles.
    public static double F99_Alpha(double theta) => Math.Pow(Math.Sin(theta), 2) / 2.0;
    public static double F99_DickeWeightSq(double theta)
        => Math.Cos(theta) / (2.0 * Math.Pow(Math.Sin(theta / 2.0), 2));

    // F88b (T1): the popcount-coherence Pi^2-odd / memory closed form for (|p> + |q>)/sqrt2 with
    // popcounts n_p, n_q and Hamming distance HD: 0 at HD = N (Pi^2-classical, GHZ), otherwise
    // (1/2 - alpha s)/(1 - s). Three alpha anchors from one Krawtchouk identity: 0 at the
    // popcount-mirror n_p + n_q = N; C(N,N/2)/(2(C(N,n_other) + C(N,N/2))) at K-intermediate
    // (even N, exactly one popcount at N/2; the adjacent case IS F98's (N+2)/(4(N+1)));
    // 1/2 generic. Static fraction s: inter-sector 1/(4 C(N,n_p)) + 1/(4 C(N,n_q)), intra 1/C(N,n).
    // Multi-state Dicke extension: alpha_total = (1 - gamma^2)/2, anchors {1/2, 3/8, 0}.
    public static double F88b_Alpha(int n, int np, int nq)
    {
        if (np + nq == n) return 0.0;                                       // popcount-mirror
        if (n % 2 == 0 && (np == n / 2) != (nq == n / 2))                   // K-intermediate
        {
            double half = Block.Binomial(n, n / 2);
            double other = Block.Binomial(n, np == n / 2 ? nq : np);
            return half / (2.0 * (other + half));
        }
        return 0.5;                                                         // generic
    }
    public static double F88b_StaticFraction(int n, int np, int nq)
        => np == nq
            ? 1.0 / Block.Binomial(n, np)
            : 1.0 / (4.0 * Block.Binomial(n, np)) + 1.0 / (4.0 * Block.Binomial(n, nq));
    public static double F88b_Pi2OddInMemory(int n, int np, int nq, int hd)
    {
        if (hd == n) return 0.0;                                            // Pi^2-classical
        double s = F88b_StaticFraction(n, np, nq);
        return (0.5 - F88b_Alpha(n, np, nq) * s) / (1.0 - s);
    }
    public static double F88b_DickeAlphaTotal(double gamma) => (1.0 - gamma * gamma) / 2.0;

    // F116 (T1): the metallic mean of the router family, r(c) = (c + sqrt(c^2+4))/2 -- the
    // positive root of r^2 = c r + 1 and the frame ratio of the period-4 router that
    // palindromizes the weighted Z-middle soft line t2 = t3 (c = t1/t2). Golden phi at c = 1,
    // the silver ratio 1 + sqrt2 at c = 2, bronze at c = 3, the 45-degree frame r = 1 at c = 0;
    // r(-c) = 1/r(c). The frame directions are the roots of the locus alpha^2 - c alpha beta
    // - beta^2 = 0 -- the identity-column determinant factors as c times exactly this locus.
    public static double F116_MetallicMean(double c) => (c + Math.Sqrt(c * c + 4.0)) / 2.0;

    private static double H2(double x) => -XLog2(x) - XLog2(1.0 - x);
    private static double XLog2(double x) => x <= 0.0 ? 0.0 : x * Math.Log2(x);

    private static long IntPow(int b, int e) { long r = 1; for (int i = 0; i < e; i++) r *= b; return r; }
}
