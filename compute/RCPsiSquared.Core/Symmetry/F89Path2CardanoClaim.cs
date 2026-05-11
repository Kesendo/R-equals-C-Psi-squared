using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>F89 path-2 (SE, DE) S_2-symmetric closed-form characteristic
/// polynomial (Tier 1 derived; sympy verified):
///
/// <code>
///   char(λ) = −(λ + 2γ)·(λ + 6γ)·[λ³ + 10γ·λ² + (28γ² + 32J²)·λ + (24γ³ + 96γJ²)]
///   dimensionless: μ³ + 10μ² + (28 + 32q²)μ + 24(1 + 4q²) = 0   (μ=λ/γ, q=J/γ)
/// </code>
///
/// <para>Cubic factor solvable in radicals via Cardano. λ=−6γ is unpopulated
/// (orthogonal to ρ_block(0) S_3-sym content); λ=−2γ + 3 cubic roots are
/// populated. Path-2 thus fully analytically tractable.</para>
///
/// <para>Derivation, Cardano closure detail, generalisation note, and
/// numerical verification at q=1.5: see
/// <c>experiments/F89_TOPOLOGY_ORBIT_CLOSURE.md</c> § "Path-2 (SE, DE) S_2-sym
/// sub-block" and <c>simulations/_f89_path2_se_de_symbolic.py</c>.</para></summary>
public sealed class F89Path2CardanoClaim : Claim
{
    private readonly F89TopologyOrbitClosure _f89;

    /// <summary>The two linear-factor eigenvalues at path-2 (SE, DE) S_2-sym
    /// sub-block: <c>(−2γ, −6γ)</c>. Pure-AT rates per F89c; no J-dependence.
    /// At our ρ_block(0) projection, only −2γ is populated (the −6γ eigenvector
    /// is orthogonal to ρ_block(0)'s S_3-symmetric content).</summary>
    public static (double LowRate, double HighRate) LinearFactorRates(double gammaZero)
    {
        if (gammaZero < 0) throw new ArgumentOutOfRangeException(nameof(gammaZero), gammaZero, "γ₀ must be ≥ 0.");
        return (-2.0 * gammaZero, -6.0 * gammaZero);
    }

    /// <summary>The cubic-factor coefficients in dimensionless form <c>μ = λ/γ</c>,
    /// <c>q = J/γ</c>: returns the cubic <c>μ³ + a·μ² + b·μ + c = 0</c> with
    /// a = 10, b = 28 + 32q², c = 24·(1 + 4q²).</summary>
    public static (double A, double B, double C) CubicCoefficients(double q)
    {
        if (q < 0) throw new ArgumentOutOfRangeException(nameof(q), q, "q = J/γ must be ≥ 0.");
        double a = 10.0;
        double b = 28.0 + 32.0 * q * q;
        double c = 24.0 * (1.0 + 4.0 * q * q);
        return (a, b, c);
    }

    /// <summary>The 3 cubic eigenvalues μ_k = λ_k/γ at given q = J/γ, computed
    /// via numerical polynomial roots (numpy-equivalent). Returns 1 real root
    /// and 1 complex-conjugate pair (for q in the path-2 typical range; the
    /// cubic discriminant determines real-vs-complex behaviour).</summary>
    public static System.Numerics.Complex[] CubicEigenvaluesNumerical(double q)
    {
        var (a, b, c) = CubicCoefficients(q);
        // Cubic μ³ + a·μ² + b·μ + c = 0; depressed via μ = ν − a/3:
        // ν³ + p·ν + r = 0 with p = b − a²/3, r = 2a³/27 − a·b/3 + c
        double p = b - a * a / 3.0;
        double r = 2.0 * a * a * a / 27.0 - a * b / 3.0 + c;
        // Cardano discriminant: Δ = −4p³ − 27r²
        // Δ < 0: 1 real + 2 complex; Δ > 0: 3 real
        double disc = -4.0 * p * p * p - 27.0 * r * r;
        var roots = new System.Numerics.Complex[3];
        if (disc < 0)
        {
            // Standard Cardano (one real root + complex pair)
            double sqrt_term = Math.Sqrt(-disc / 108.0);
            double u_cube = -r / 2.0 + sqrt_term;
            double v_cube = -r / 2.0 - sqrt_term;
            double u = Math.Sign(u_cube) * Math.Pow(Math.Abs(u_cube), 1.0 / 3.0);
            double v = Math.Sign(v_cube) * Math.Pow(Math.Abs(v_cube), 1.0 / 3.0);
            double real_root = u + v - a / 3.0;
            // Complex pair: −(u + v)/2 ± i·(u − v)·√3/2 (in ν), shift by −a/3
            double real_pair = -(u + v) / 2.0 - a / 3.0;
            double imag_pair = (u - v) * Math.Sqrt(3.0) / 2.0;
            roots[0] = new System.Numerics.Complex(real_root, 0);
            roots[1] = new System.Numerics.Complex(real_pair, imag_pair);
            roots[2] = new System.Numerics.Complex(real_pair, -imag_pair);
        }
        else
        {
            // Trigonometric solution for 3 real roots (q very small or specific values)
            double rho = Math.Sqrt(-p * p * p / 27.0);
            double theta = Math.Acos(-r / (2.0 * rho)) / 3.0;
            double scale = 2.0 * Math.Pow(rho, 1.0 / 3.0);
            for (int k = 0; k < 3; k++)
                roots[k] = new System.Numerics.Complex(
                    scale * Math.Cos(theta - 2.0 * Math.PI * k / 3.0) - a / 3.0, 0);
        }
        return roots;
    }

    /// <summary>Drift check: at q = 1.5, the cubic real root should be μ ≈
    /// −3.0448 and the complex pair μ ≈ −3.4776 ± 8.169i. Returns true if
    /// numerical roots match these reference values within 1·10⁻³.</summary>
    public static bool VerifyAtStandardQ()
    {
        var roots = CubicEigenvaluesNumerical(1.5);
        // Expected: 1 real near -3.0448; complex pair near -3.4776 ± 8.169i
        bool foundReal = false, foundComplex = false;
        foreach (var r in roots)
        {
            if (Math.Abs(r.Imaginary) < 1e-9 && Math.Abs(r.Real - (-3.0448)) < 1e-3)
                foundReal = true;
            if (Math.Abs(r.Real - (-3.4776)) < 1e-3 && Math.Abs(Math.Abs(r.Imaginary) - 8.169) < 1e-2)
                foundComplex = true;
        }
        return foundReal && foundComplex;
    }

    public F89Path2CardanoClaim(F89TopologyOrbitClosure f89)
        : base("F89 path-2 (SE,DE) S_2-sym closed-form factorisation: char(λ) = −(λ+2γ)(λ+6γ)·[λ³ + 10γλ² + (28γ²+32J²)λ + (24γ³+96γJ²)]; cubic solvable in radicals via Cardano; path-2 fully analytically tractable",
               Tier.Tier1Derived,
               "experiments/F89_TOPOLOGY_ORBIT_CLOSURE.md + " +
               "simulations/_f89_path2_se_de_symbolic.py + " +
               "compute/RCPsiSquared.Core/Symmetry/F89TopologyOrbitClosure.cs")
    {
        _f89 = f89 ?? throw new ArgumentNullException(nameof(f89));
    }

    public override string DisplayName =>
        "F89 path-2 (SE,DE) S_2-sym characteristic polynomial: linear×linear×cubic via Cardano";

    public override string Summary =>
        $"5×5 char(λ) factors as −(λ+2γ)(λ+6γ)·[cubic]; cubic μ³ + 10μ² + (28+32q²)μ + 24(1+4q²) = 0 in dimensionless μ=λ/γ, q=J/γ; Cardano-solvable; verified at q=1.5 ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("Cubic coefficients at q=1.5",
                summary: $"a={CubicCoefficients(1.5).A}, b={CubicCoefficients(1.5).B}, c={CubicCoefficients(1.5).C}");
            yield return new InspectableNode("Sample roots at q=1.5",
                summary: $"μ ≈ −3.0448 (real), μ ≈ −3.4776 ± 8.169i (complex pair); VerifyAtStandardQ = {VerifyAtStandardQ()}");
        }
    }
}
