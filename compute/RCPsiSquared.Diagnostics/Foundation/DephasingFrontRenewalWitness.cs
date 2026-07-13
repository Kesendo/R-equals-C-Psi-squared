using System;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;
using System.Numerics;
using MathNet.Numerics;                    // SpecialFunctions: BesselJ, BesselIScaled, AiryAi
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.Inspection;

namespace RCPsiSquared.Diagnostics.Foundation;

/// <summary>The live witness of the dephasing-front renewal representation
/// (<c>docs/proofs/PROOF_DEPHASING_FRONT_RENEWAL.md</c>, <c>inspect --root renewal</c>): the single
/// excitation released on a chain and WATCHED evolves by the exact renewal
///
/// <code>
///   P_n(t) = e^{−Γt}·S_n(t),   S_n(t) = |G_{n0}(t)|² + Γ∫₀ᵗ ds Σ_m |G_{nm}(t−s)|²·S_m(s),   Γ = 4γ
/// </code>
///
/// with G the clean single-particle propagator (the F2b band, on the infinite chain
/// G_{nm}(τ) = (−i)^{|n−m|}·J_{|n−m|}(2Jτ)); the j=0 term is the coherent front, the j ≥ 1 refill
/// ladder the incoherent halo, and in momentum-Laplace space the ladder closes to
/// Ŝ(p, z) = 1/(√(z² + a²) − Γ), a = 4J·sin(p/2). The witness recomputes six from-below checks at
/// inspect time (all sub-second, small N×N sectors): renewal-vs-RK4 agreement, probability
/// conservation, the coherent-front Bessel identity, the Γ=0 clean-wave limit, the Haken-Strobl
/// diffusive plateau, and the closed refill constant I₁ (the Airy corollary).</summary>
public sealed class DephasingFrontRenewalWitness : IInspectable
{
    private static readonly CultureInfo Inv = CultureInfo.InvariantCulture;

    // ---- renewal-vs-RK4 sector (cases a, b): a small interior-seeded chain, t small enough that the
    // ballistic front (speed 2J) stays interior; agreement is exact up to the discretization only.
    public int N { get; }
    public int Seed { get; }
    public double J { get; }
    public double Gamma { get; }
    public double Tmax { get; }
    public double Dt { get; }
    public double GammaPhi => 4.0 * Gamma;   // Γ = 4γ (the Absorption-Theorem sector rate)

    // ---- the analytic corollary anchors (cases e, f).
    private const int PlateauN = 200;                       // the diffusive plateau read at n = 200
    private const double PlateauTarget = 0.24197072451914337;   // e^{−1/2}/√(2π)
    private const double AiryZeroAlpha = 1.0187929716;      // |first zero of Ai′|; 2c = 2^{2/3}·α
    private const double I1Target = 0.27694424;             // 1/12 + ¼∫₀^{2c} Ai(−w) dw

    public DephasingFrontRenewalWitness(
        int n = 27, int? seed = null, double j = 1.0, double gamma = 0.15,
        double tMax = 2.5, double dt = 0.005)
    {
        if (n < 9) throw new ArgumentOutOfRangeException(nameof(n), n, "N >= 9 (an interior seed with room for the front)");
        if (j <= 0) throw new ArgumentOutOfRangeException(nameof(j), j, "J must be > 0");
        if (gamma < 0) throw new ArgumentOutOfRangeException(nameof(gamma), gamma, "gamma must be >= 0");
        if (dt <= 0 || tMax <= 0) throw new ArgumentOutOfRangeException(nameof(dt), "dt, tMax must be > 0");
        N = n; Seed = seed ?? n / 2; J = j; Gamma = gamma; Tmax = tMax; Dt = dt;
    }

    private IReadOnlyList<BatteryCase>? _cases;
    public IReadOnlyList<BatteryCase> Cases => _cases ??= BuildBattery();
    public int PassCount => Cases.Count(c => c.Passes);

    /// <summary>One from-below check of the renewal representation: a named claim, its computed detail,
    /// and the expected-vs-actual verdict tokens (equal ⟺ PASS).</summary>
    public readonly record struct BatteryCase(string Name, string Detail, string Expected, string Actual)
    {
        public bool Passes => string.Equals(Expected, Actual, StringComparison.Ordinal);
    }

    // ==================================================================
    // battery
    // ==================================================================

    private IReadOnlyList<BatteryCase> BuildBattery()
    {
        var (renewalDev, traceDrift) = RenewalVersusRk4();
        double frontDev = CoherentFrontIdentity();
        double cleanDev = CleanWaveLimit();
        double plateau = HakenStroblPlateau();
        double i1 = RefillConstantI1();

        return new List<BatteryCase>
        {
            new("renewal representation reproduces the RK4 Lindblad sector",
                $"N={N} seed={Seed}, t<= {Tmax.ToString("0.##", Inv)}, dt={Dt.ToString("0.###", Inv)}: " +
                $"max |P_renewal − P_RK4| = {renewalDev.ToString("E2", Inv)}",
                "match(<=1e-4)", renewalDev <= 1e-4 ? "match(<=1e-4)" : $"dev={renewalDev.ToString("E2", Inv)}"),

            new("the RK4 Lindblad run conserves probability",
                $"max_t |Σ_n P_n − 1| = {traceDrift.ToString("E2", Inv)}",
                "conserved(<=1e-12)", traceDrift <= 1e-12 ? "conserved(<=1e-12)" : $"drift={traceDrift.ToString("E2", Inv)}"),

            new("the j=0 term is the coherent front e^{−Γt}·J_n(2Jt)² = |<a_n>|²",
                $"amplitude ODE (decay Γ/2) vs Bessel: max dev = {frontDev.ToString("E2", Inv)}",
                "coherent-front(<=1e-8)", frontDev <= 1e-8 ? "coherent-front(<=1e-8)" : $"dev={frontDev.ToString("E2", Inv)}"),

            new("the Γ=0 limit returns the clean Bessel wave J_n(2Jt)²",
                $"P_n(Γ=0) vs J_n(2Jt)²: max dev = {cleanDev.ToString("E2", Inv)}",
                "clean-wave(<=1e-10)", cleanDev <= 1e-10 ? "clean-wave(<=1e-10)" : $"dev={cleanDev.ToString("E2", Inv)}"),

            new("the Haken-Strobl diffusive plateau max_t e^{−2Dt}I_n(2Dt)·n -> e^{−1/2}/√(2π)",
                $"n={PlateauN}: plateau·n = {plateau.ToString("0.00000", Inv)} (target {PlateauTarget.ToString("0.00000", Inv)})",
                "plateau~0.24197",
                Math.Abs(plateau - PlateauTarget) / PlateauTarget < 0.01 ? "plateau~0.24197" : plateau.ToString("0.00000", Inv)),

            new("the closed refill constant I₁ = 1/12 + ¼∫₀^{2c} Ai(−w) dw",
                $"2c = 2^(2/3)·α, α={AiryZeroAlpha.ToString("0.#######", Inv)}: I₁ = {i1.ToString("0.000000", Inv)} (target {I1Target.ToString("0.000000", Inv)})",
                "I1~0.276944",
                Math.Abs(i1 - I1Target) < 1e-5 ? "I1~0.276944" : i1.ToString("0.000000", Inv)),
        };
    }

    // ==================================================================
    // (a) + (b): the watched-walk sector, two ways
    // ==================================================================

    /// <summary>Solve the sector both ways on the SAME finite chain h and compare the site populations:
    /// (i) direct RK4 of ρ̇ = −i[h, ρ] − Γ(ρ − diag ρ); (ii) the Volterra renewal (★) marched forward
    /// with the trapezoidal rule (the b=a self-term is the diagonal K(0)=I, solved implicitly). Returns
    /// (max population deviation, max probability-conservation drift of the RK4 run).</summary>
    private (double renewalDev, double traceDrift) RenewalVersusRk4()
    {
        int steps = (int)Math.Round(Tmax / Dt);
        double gPhi = GammaPhi;

        // clean propagator on the grid: G(τ)_{nm} = Σ_k U_nk e^{−iE_k τ} U_mk, from the real symmetric
        // hopping eigendecomposition. K(Δ) = |G(Δ·dt)|² is the one-segment population kernel.
        var (evals, evecs) = HoppingEigen();
        var K = new double[steps + 1][];          // K[d] flattened N*N: |G(d·dt)_{nm}|²
        for (int d = 0; d <= steps; d++)
            K[d] = PopulationKernel(evals, evecs, d * Dt);

        // (i) RK4 of the N×N density matrix (flat row-major, index n*N+m).
        int nn = N * N;
        var rho = new Complex[nn];
        rho[Seed * N + Seed] = Complex.One;
        var rk4P = new double[steps + 1][];
        rk4P[0] = Diagonal(rho);
        double traceDrift = 0.0;
        var k1 = new Complex[nn]; var k2 = new Complex[nn]; var k3 = new Complex[nn]; var k4 = new Complex[nn];
        var tmp = new Complex[nn];
        for (int a = 1; a <= steps; a++)
        {
            Deriv(rho, k1, gPhi);
            Axpy(rho, k1, Dt / 2, tmp); Deriv(tmp, k2, gPhi);
            Axpy(rho, k2, Dt / 2, tmp); Deriv(tmp, k3, gPhi);
            Axpy(rho, k3, Dt, tmp); Deriv(tmp, k4, gPhi);
            for (int i = 0; i < nn; i++)
                rho[i] += Dt / 6.0 * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]);
            rk4P[a] = Diagonal(rho);
            double tr = 0.0;
            for (int n = 0; n < N; n++) tr += rk4P[a][n];
            traceDrift = Math.Max(traceDrift, Math.Abs(tr - 1.0));
        }

        // (ii) the renewal march. S[a][n]; source = K[a][n, seed]; the b=a term (K(0)=I) is implicit.
        var S = new double[steps + 1][];
        S[0] = new double[N];
        for (int n = 0; n < N; n++) S[0][n] = K[0][n * N + Seed];
        double denom = 1.0 - 0.5 * gPhi * Dt;
        var acc = new double[N];
        for (int a = 1; a <= steps; a++)
        {
            Array.Clear(acc);
            // b = 0 carries the trapezoidal endpoint weight 1/2; 0 < b < a carry weight 1.
            AddKtimesS(K[a], S[0], acc, 0.5);
            for (int b = 1; b < a; b++)
                AddKtimesS(K[a - b], S[b], acc, 1.0);
            var Sa = new double[N];
            for (int n = 0; n < N; n++)
                Sa[n] = (K[a][n * N + Seed] + gPhi * Dt * acc[n]) / denom;
            S[a] = Sa;
        }

        double renewalDev = 0.0;
        for (int a = 0; a <= steps; a++)
        {
            double e = Math.Exp(-gPhi * a * Dt);
            for (int n = 0; n < N; n++)
                renewalDev = Math.Max(renewalDev, Math.Abs(e * S[a][n] - rk4P[a][n]));
        }
        return (renewalDev, traceDrift);
    }

    /// <summary>K(Δ)·S accumulated into acc with weight w (K flat N*N row-major).</summary>
    private void AddKtimesS(double[] K, double[] Sb, double[] acc, double w)
    {
        for (int n = 0; n < N; n++)
        {
            double s = 0.0;
            int row = n * N;
            for (int m = 0; m < N; m++) s += K[row + m] * Sb[m];
            acc[n] += w * s;
        }
    }

    private double[] Diagonal(Complex[] rho)
    {
        var p = new double[N];
        for (int n = 0; n < N; n++) p[n] = rho[n * N + n].Real;
        return p;
    }

    private static void Axpy(Complex[] baseM, Complex[] k, double scale, Complex[] outM)
    {
        for (int i = 0; i < baseM.Length; i++) outM[i] = baseM[i] + scale * k[i];
    }

    /// <summary>ρ̇ = −i[h, ρ] − Γ(ρ − diag ρ), h the tridiagonal J-hopping (flat row-major).</summary>
    private void Deriv(Complex[] r, Complex[] d, double gPhi)
    {
        for (int n = 0; n < N; n++)
            for (int m = 0; m < N; m++)
            {
                // (h r)_{nm} = J(r_{n-1,m} + r_{n+1,m}); (r h)_{nm} = J(r_{n,m-1} + r_{n,m+1}).
                Complex hr = Complex.Zero, rh = Complex.Zero;
                if (n > 0) hr += r[(n - 1) * N + m];
                if (n < N - 1) hr += r[(n + 1) * N + m];
                if (m > 0) rh += r[n * N + (m - 1)];
                if (m < N - 1) rh += r[n * N + (m + 1)];
                Complex comm = J * (hr - rh);
                Complex val = -Complex.ImaginaryOne * comm;
                if (n != m) val -= gPhi * r[n * N + m];
                d[n * N + m] = val;
            }
    }

    // ==================================================================
    // (c) + (d): the coherent front and the clean-wave limit vs Bessel
    // ==================================================================

    /// <summary>The j=0 identity e^{−Γt}·J_{|n−seed|}(2Jt)² = |⟨a_n⟩|²: the noise-averaged amplitude
    /// obeys ȧ = (−ih − Γ/2)a, so |a_n(t)|² = e^{−Γt}|G_{n,seed}(t)|². Compared against the infinite-chain
    /// Bessel form on a large chain with an interior seed and t small enough that the front stays interior
    /// (the finite-vs-infinite propagator difference is then boundary-reflection-small).</summary>
    private double CoherentFrontIdentity()
    {
        // A genuinely independent route to the j=0 term: integrate the amplitude ODE
        // a-dot = -i h a - (Gamma/2) a by RK4 (the noise-averaged amplitude damps at HALF the
        // coherence rate) and compare |a_n(t)|^2 against e^{-Gamma t} J_{|n-seed|}(2Jt)^2.
        // This is not the closed form restated: the Gamma/2 decay content is exercised, not cancelled.
        const int nBig = 41, seedBig = 20;
        const double dtAmp = 0.0025;
        double gPhi = GammaPhi;
        var a = new Complex[nBig];
        a[seedBig] = Complex.One;

        Complex[] Rhs(Complex[] x)
        {
            var r = new Complex[nBig];
            for (int n = 0; n < nBig; n++)
            {
                Complex hop = Complex.Zero;
                if (n > 0) hop += J * x[n - 1];
                if (n < nBig - 1) hop += J * x[n + 1];
                r[n] = -Complex.ImaginaryOne * hop - 0.5 * gPhi * x[n];
            }
            return r;
        }

        Complex[] Axpy(Complex[] u, Complex[] v, double s)
        {
            var r = new Complex[nBig];
            for (int n = 0; n < nBig; n++) r[n] = u[n] + s * v[n];
            return r;
        }

        double maxDev = 0.0;
        double t = 0.0;
        var checkTimes = new[] { 1.0, 2.0, 3.0 };
        int next = 0;
        while (next < checkTimes.Length)
        {
            var k1 = Rhs(a);
            var k2 = Rhs(Axpy(a, k1, dtAmp / 2));
            var k3 = Rhs(Axpy(a, k2, dtAmp / 2));
            var k4 = Rhs(Axpy(a, k3, dtAmp));
            for (int n = 0; n < nBig; n++)
                a[n] += (dtAmp / 6) * (k1[n] + 2 * k2[n] + 2 * k3[n] + k4[n]);
            t += dtAmp;
            if (t >= checkTimes[next] - 1e-9)
            {
                for (int n = 0; n < nBig; n++)
                {
                    double amp2 = a[n].Real * a[n].Real + a[n].Imaginary * a[n].Imaginary;
                    double bj = SpecialFunctions.BesselJ(Math.Abs(n - seedBig), 2.0 * J * t);
                    double bessel = Math.Exp(-gPhi * t) * bj * bj;
                    maxDev = Math.Max(maxDev, Math.Abs(amp2 - bessel));
                }
                next++;
            }
        }
        return maxDev;
    }

    /// <summary>The Γ=0 limit of the renewal (no refill): P_n(t) = |G_{n,seed}(t)|² = J_{|n−seed|}(2Jt)².</summary>
    private double CleanWaveLimit()
    {
        const int nBig = 41, seedBig = 20;
        var (evals, evecs) = HoppingEigen(nBig);
        double maxDev = 0.0;
        foreach (double t in new[] { 1.0, 2.0, 3.0 })
        {
            var col = PropagatorColumn(evals, evecs, nBig, seedBig, t);
            for (int n = 0; n < nBig; n++)
            {
                double p = col[n].Real * col[n].Real + col[n].Imaginary * col[n].Imaginary;
                double bj = SpecialFunctions.BesselJ(Math.Abs(n - seedBig), 2.0 * J * t);
                maxDev = Math.Max(maxDev, Math.Abs(p - bj * bj));
            }
        }
        return maxDev;
    }

    // ==================================================================
    // (e) + (f): the two analytic corollary anchors
    // ==================================================================

    /// <summary>The diffusive plateau: the (☆)-pole long-wavelength Green's function is the discrete
    /// heat kernel P_n(t) = e^{−2Dt}I_n(2Dt), whose global-in-time maximum times n approaches
    /// e^{−1/2}/√(2π) (peak at 2Dt = n²). D drops out under the max over t; scan x = 2Dt across n².</summary>
    private static double HakenStroblPlateau()
    {
        // the peak sits at x = n²; scan a broad window with the scaled modified Bessel (overflow-safe).
        double n = PlateauN;
        double best = 0.0;
        double lo = 0.5 * n * n, hi = 1.5 * n * n, step = (hi - lo) / 800.0;
        for (double x = lo; x <= hi; x += step)
        {
            double v = SpecialFunctions.BesselIScaled(n, x);   // e^{−x}·I_n(x)
            if (v > best) best = v;
        }
        return best * n;
    }

    /// <summary>The single-refill front integral saturates to I₁ = 1/12 + ¼∫₀^{2c} Ai(−w) dw with
    /// 2c = 2^{2/3}·α, α the first zero of Ai′ (Simpson quadrature of the Airy function).</summary>
    private static double RefillConstantI1()
    {
        double twoC = Math.Pow(2.0, 2.0 / 3.0) * AiryZeroAlpha;
        const int m = 2000;                       // Simpson intervals (even)
        double h = twoC / m;
        double sum = SpecialFunctions.AiryAi(0.0) + SpecialFunctions.AiryAi(-twoC);
        for (int i = 1; i < m; i++)
        {
            double w = i * h;
            sum += (i % 2 == 1 ? 4.0 : 2.0) * SpecialFunctions.AiryAi(-w);
        }
        double integral = sum * h / 3.0;
        return 1.0 / 12.0 + 0.25 * integral;
    }

    // ==================================================================
    // hopping propagator helpers
    // ==================================================================

    private (double[] E, double[,] U) HoppingEigen() => HoppingEigen(N);

    /// <summary>Real symmetric tridiagonal J-hopping eigendecomposition h = U diag(E) Uᵀ (open chain).</summary>
    private (double[] E, double[,] U) HoppingEigen(int n)
    {
        var h = Matrix<double>.Build.Dense(n, n);
        for (int i = 0; i < n - 1; i++) { h[i, i + 1] = J; h[i + 1, i] = J; }
        var evd = h.Evd(Symmetricity.Symmetric);
        var E = evd.EigenValues.Select(z => z.Real).ToArray();
        var vec = evd.EigenVectors;
        var U = new double[n, n];
        for (int i = 0; i < n; i++)
            for (int k = 0; k < n; k++) U[i, k] = vec[i, k];
        return (E, U);
    }

    /// <summary>|G(τ)_{nm}|² flattened row-major, G(τ) = U diag(e^{−iEτ}) Uᵀ.</summary>
    private double[] PopulationKernel(double[] E, double[,] U, double tau)
    {
        int n = E.Length;
        var cos = new double[n]; var sin = new double[n];
        for (int k = 0; k < n; k++) { cos[k] = Math.Cos(E[k] * tau); sin[k] = Math.Sin(E[k] * tau); }
        var K = new double[n * n];
        for (int a = 0; a < n; a++)
            for (int b = 0; b < n; b++)
            {
                double re = 0.0, im = 0.0;
                for (int k = 0; k < n; k++)
                {
                    double w = U[a, k] * U[b, k];
                    re += w * cos[k];        // Re e^{−iEτ} = cos
                    im -= w * sin[k];        // Im e^{−iEτ} = −sin
                }
                K[a * n + b] = re * re + im * im;
            }
        return K;
    }

    /// <summary>The seed column of G(t) = U diag(e^{−iEt}) Uᵀ (the clean amplitude from the seed site).</summary>
    private static Complex[] PropagatorColumn(double[] E, double[,] U, int n, int seed, double t)
    {
        var col = new Complex[n];
        for (int a = 0; a < n; a++)
        {
            Complex s = Complex.Zero;
            for (int k = 0; k < n; k++)
                s += U[a, k] * U[seed, k] * Complex.Exp(-Complex.ImaginaryOne * E[k] * t);
            col[a] = s;
        }
        return col;
    }

    // ==================================================================
    // IInspectable
    // ==================================================================

    public string DisplayName => $"DephasingFrontRenewalWitness (N={N}, seed={Seed}, γ={Gamma.ToString("0.###", Inv)})";

    public string Summary =>
        "the dephasing-front renewal representation, recomputed live: the watched single excitation is the " +
        "unwatched wave repeatedly caught and released, P_n(t) = e^{−Γt}·S_n(t) with S the Volterra refill " +
        $"ladder (★), Γ = 4γ (the Absorption-Theorem sector rate). Battery {PassCount}/{Cases.Count} pass: " +
        "the renewal reproduces the RK4 Lindblad sector, probability is conserved, the j=0 term is the " +
        "coherent front e^{−Γt}J_n(2Jt)², the Γ=0 limit is the clean Bessel wave, the diffusive plateau " +
        "hits e^{−1/2}/√(2π), and the refill constant I₁ closes on the Airy integral. Momentum-Laplace form " +
        "Ŝ(p,z) = 1/(√(z²+a²)−Γ), a = 4J·sin(p/2). Proof: docs/proofs/PROOF_DEPHASING_FRONT_RENEWAL.md.";

    public IEnumerable<IInspectable> Children
    {
        get
        {
            yield return new InspectableNode("the representation",
                summary: "P_n(t) = e^{−Γt}·S_n(t), S_n = |G_{n0}|² + Γ∫₀ᵗ Σ_m |G_{nm}(t−s)|²·S_m(s), Γ=4γ. " +
                         "Every refill order carries the same e^{−Γt}; the j=0 term is the coherent front, the " +
                         "j≥1 halo the incoherent refill. Closed form Ŝ(p,z) = 1/(√(z²+a²)−Γ), a = 4J·sin(p/2): " +
                         "conserves probability at p=0, returns the clean Bessel wave at Γ=0, and the small-p pole " +
                         "√(Γ²−a²) is the diffusive branch with D = 2J²/Γ (the F123 SurvivorDiffusionGradient sibling).");
            foreach (var c in Cases)
                yield return new InspectableNode(c.Name,
                    summary: $"{c.Detail}; expected {c.Expected}, got {c.Actual}, " + (c.Passes ? "PASS" : "FAIL"));
        }
    }

    public InspectablePayload Payload => InspectablePayload.Empty;
}
