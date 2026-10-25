using System.Numerics;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Numerics;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>F168, the cat under a turned field (Tier 1 derived,
/// <c>experiments/THE_CAT_UNDER_A_TURNED_FIELD.md</c>): N bits under Z dephasing at rates γ_l
/// (σ = Σγ_l &gt; 0), a field h_l·(cos φ·X + sin φ·Y) at one common angle with every h_l ≠ 0, and a
/// Hermitian coupling H_c commuting with X^⊗N and with Σ_l Z_l. For the cat with phase θ the
/// centred odd moments m_k = Tr(ρ·(L + σ)^k ρ) vanish for odd k &lt; 2N + 1, every odd moment is
/// sin²(Nφ − θ) times a factor free of the angle, and
///
/// <code>
///     m_{2N+1} = (σ/3)·(2N + 1)!·∏_l h_l²·sin²(Nφ − θ),
/// </code>
///
/// free of H_c; so the cat meets itself the same way forwards and, at the price e^{−2σt},
/// backwards exactly when sin(Nφ − θ) = 0.
///
/// <para>The proof is read on the letter cube: the turn about Z carries the field to X and the cat
/// to the cat with phase θ − Nφ; the parity of the count of Y and Z letters, F63's turn by X⊗N,
/// never changes along a step; the first move returns the cat at sin = 0 (F99's α = 0 endpoint),
/// which with the parity makes every odd moment sin² times an angle-free factor; the charge
/// d = p − q puts the two far corners 2N field steps apart, and without the watching every odd
/// moment vanishes, so the first is the (2N + 1)-th and contains no step of H_c; one bit's round
/// trip gives the third.</para>
///
/// <para>Live witness: <c>inspect --root catfield</c>
/// (<c>compute/RCPsiSquared.Diagnostics/Foundation/CatTurnedFieldWitness.cs</c>), which reads the
/// moments exactly in the Pauli-string basis from the mask algebra, a different road from the
/// page's sympy gate in <c>simulations/the_wave_between.py</c> section 10.</para>
///
/// <para>Typed parents: <see cref="F63LCommutesPi2Pi2Inheritance"/> (step 3, the parity that keeps
/// the two kinds of words apart) and <see cref="XGlobalEigenstateMirrorPi2Inheritance"/> (the
/// X⊗N-eigenstate endpoint, the two points where the property holds).</para></summary>
public sealed class CatTurnedFieldSelfEncounterClaim : Claim
{
    public F63LCommutesPi2Pi2Inheritance F63 { get; }
    public XGlobalEigenstateMirrorPi2Inheritance XGlobalEigenstate { get; }

    /// <summary>The order of the first odd moment that can be nonzero, 2N + 1.</summary>
    public static int FirstOddOrder(int n)
    {
        if (n < 1) throw new ArgumentOutOfRangeException(nameof(n), n, "N must be ≥ 1");
        return 2 * n + 1;
    }

    /// <summary>m_{2N+1} = (σ/3)·(2N + 1)!·∏h_l²·sin²(Nφ − θ), exact, from σ, the product of the
    /// squared field strengths and the squared sine.</summary>
    public static BigRational LeadingOddMoment(int n, BigRational sigma, BigRational productOfSquaredStrengths,
                                               BigRational sinSquared)
    {
        int order = FirstOddOrder(n);
        BigInteger factorial = BigInteger.One;
        for (int j = 2; j <= order; j++) factorial *= j;
        return sigma / new BigRational(3) * new BigRational(factorial) * productOfSquaredStrengths * sinSquared;
    }

    /// <summary>sin(Nφ) = Im((c + i·s)^N) for a rational point (c, s) of the unit circle, exact.</summary>
    public static BigRational SinOfMultiple(int n, BigRational c, BigRational s)
    {
        if (n < 0) throw new ArgumentOutOfRangeException(nameof(n), n, "N must be ≥ 0");
        if (c * c + s * s != BigRational.One)
            throw new ArgumentException("(c, s) must lie on the unit circle");
        BigRational re = BigRational.One, im = BigRational.Zero;
        for (int j = 0; j < n; j++)
            (re, im) = (re * c - im * s, re * s + im * c);
        return im;
    }

    public CatTurnedFieldSelfEncounterClaim(F63LCommutesPi2Pi2Inheritance f63,
                                            XGlobalEigenstateMirrorPi2Inheritance xGlobalEigenstate)
        : base(
            "The cat under a field in the equator plane meets itself the same way forwards and, at the price, backwards " +
            "exactly when sin(Nφ − θ) = 0: its centred odd moments vanish below order 2N + 1 and the first is " +
            "(σ/3)·(2N + 1)!·∏h_l²·sin²(Nφ − θ), free of any coupling commuting with X^⊗N and Σ Z_l",
            Tier.Tier1Derived,
            "experiments/THE_CAT_UNDER_A_TURNED_FIELD.md + docs/ANALYTICAL_FORMULAS.md (F168, the index entry)")
    {
        F63 = f63 ?? throw new ArgumentNullException(nameof(f63));
        XGlobalEigenstate = xGlobalEigenstate ?? throw new ArgumentNullException(nameof(xGlobalEigenstate));
    }

    public override string DisplayName =>
        "F168: the cat under a turned field, the first odd moment at order 2N + 1";

    public override string Summary =>
        "For the cat (|0…0⟩ + e^{iθ}|1…1⟩)/√2 under Z dephasing (σ = Σγ_l > 0), a field h_l·(cos φ·X + sin φ·Y) at one " +
        "angle and any coupling commuting with X^⊗N and Σ Z_l: odd m_k = 0 for k < 2N + 1, every odd m_k = " +
        "sin²(Nφ − θ)·G_k with G_k free of the angle, m_{2N+1} = (σ/3)·(2N + 1)!·∏h_l²·sin²(Nφ − θ); the property " +
        "S(t) = e^{−2σt}S(−t) holds exactly when sin(Nφ − θ) = 0 (Tier 1 derived; live witness: inspect --root catfield).";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("the turn",
                summary: "U = ⊗diag(1, e^{iφ}) carries X onto cos φ·X + sin φ·Y and leaves Z and H_c alone, so the field at φ is the field along X seen through a copy, and the cat with phase θ becomes the cat with phase θ − Nφ");
            yield return new InspectableNode("the parity (F63)",
                summary: "Y and Z anticommute with X; a step carries P onto T·P and the count of Y and Z letters adds modulo 2, so every term with an even count keeps a word's parity: the cat's quiet and cosine words never meet its sine words");
            yield return new InspectableNode("the two points (F99)",
                summary: "the first move ρ ↦ ρ·X^⊗N returns the cat with phase 0 and with phase π up to the sign, so every odd moment vanishes there; with the parity every odd moment is sin²(Nφ − θ) times a factor free of the angle; sin² is 2α of the turned cat");
            yield return new InspectableNode("the order",
                summary: "the two far corners |0…0⟩⟨1…1| and |1…1⟩⟨0…0| sit at d = −N and +N; the field changes d by one per step, H_c and the watching keep it, and without the watching every odd moment vanishes: the first odd moment is the (2N + 1)-th, 2N field steps and one rate, no step of H_c");
            yield return new InspectableNode("the value",
                summary: "with H_c gone every bit runs alone; one bit's doubly flipped pair carries, times e^{γt}, h²t²e^{−2iφ}(1 − γt/3) + O(t⁴), the third read on the cube as one net charge −γ over the three slots of a two-step round trip, times 2!/3!; so m_{2N+1} = (σ/3)·(2N + 1)!·∏h_l²·sin²(Nφ − θ)");
            yield return new InspectableNode("the coupling at identical bits",
                summary: "with one strength, one rate and one angle on every bit the cat stays invariant under every permutation of the bits, a swap commutes with it, and an isotropic coupling J_b(2·SWAP_b − 1) never enters the self-encounter on any graph; a swap needs a Z_aZ_b step on each side, so m_{2N+3} sees each bond only through Δ_b − J_b");
        }
    }
}
