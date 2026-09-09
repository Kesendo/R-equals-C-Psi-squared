using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>F164: the N=5 Route-B A₂ locus that sits at a REAL coupling and a REAL eigenvalue carries a
/// SEMISIMPLE double eigenvalue, and the reason is a field argument rather than a contour reading
/// (Tier 1 derived; proof <c>docs/proofs/PROOF_N5_REAL_Q_DIABOLIC.md</c>).
///
/// <para><b>Why this locus and no other.</b> One A₂ root in the layer parameter w has both its coupling
/// and its eigenvalue real; it lifts to the two q loci the inventory counts, q = ±1.129 250 970 874 767 1
/// at the same λ = −4.791 960 365 179 641 0, in the (1,2) coherence block of the open N=5 XY chain under
/// uniform Z dephasing, R-odd sector, Δ = 0, γ = 1. One root, two loci. They are the only points in the
/// Route-B corpus at N=5 and N=6 an experiment could be set to: at N=6 none of the 266 loci is real,
/// counting the stored atlas by a float read of a field that artifact calls a display midpoint rather than
/// an algebraic root, with 118 at purely imaginary q where the block is Hermitian and 148 generic complex.</para>
///
/// <para><b>The argument.</b> Let w₀ be the positive real root of the R-odd A₂ layer and t = √(−w).
/// A₂_O is irreducible over ℚ of degree 13, so ℚ(w₀) is a degree-13 subfield of ℝ, while t₀ = √(−w₀) is
/// nonzero purely imaginary and therefore lies outside it. Hence [ℚ(t₀) : ℚ] = 26, the minimal polynomial
/// of t₀ is A₂_O(−t²), and its 26 roots form a SINGLE Galois orbit: 1 positive w gives 2 imaginary t,
/// 6 negative w give 12 real t, 6 nonreal w give 12 nonreal t. At each of the 12 real t the block
/// L_O = D_O + t·K_O is real symmetric, so geometric multiplicity equals algebraic multiplicity there.
/// Rank is preserved by every embedding, because a minor is zero or nonzero before and after it, and both
/// multiplicities are Galois invariant, so the Hermitian-axis verdict travels the orbit to the physical
/// point. The chain never uses the value 2, which is what keeps it clear of χ = AT·F_res, whose AT factor
/// has real roots at negative w.</para>
///
/// <para><b>What this upgrades, and what it leaves alone.</b> The stores classified this locus without
/// certifying it: the artifact's <c>exactRankCertificates</c> array is empty and the reading came from a
/// floating-point contour classifier. This supplies the exact half for the 2 real-q loci of those 34; the
/// other 32 keep the numerical reading. It is not an all-N theorem and it says nothing metrological.</para>
///
/// <para><b>The sibling at N=6.</b> F163 makes the same move one chain longer and for the whole layer: the
/// R-even degree-133 A₂ factor is irreducible, one real Hermitian embedding carries semisimplicity to each
/// of its conjugates, and a diagonal sign operator exchanging the two reflection sectors then transports
/// the verdict to the 133 odd loci, covering all 266. The two are one method at two sizes, which is why the
/// N=6 inventory's <c>ExactRankExecuted = 0</c> is a statement about that classifier and not about what is
/// known.</para>
///
/// <para><b>The gate beside it, and its error model.</b> <c>simulations/n5_a2_rodd_character_gate.py</c> is
/// numerical corroboration and says so. Its coupling response is σ₂ = c·(δJ/J) with the ratio flat to
/// 1.000 18 over eight decades. The constant is direction-dependent and no parity fixes it: three
/// PALINDROMIC directions give c = 5.8391e-2, 8.3097e-2 and 1.1004e-1, a factor 1.885 within one parity.
/// What parity decides is whether a response exists at all, and that is a theorem: on a
/// REFLECTION-SYMMETRIC base an ANTI-palindromic detune obeys R V R = −V, so U_Oᵀ V U_O vanishes
/// identically at every N, exactly 0.0 in the integer orbit basis; without that symmetry R does not
/// commute with L and there is no R-odd sector to miss. It is a BLIND SPOT of the sector reading rather
/// than a null response: PVP = 0 keeps σ₂ on the float floor over the ten decades δJ/J = 1e-12 to 1e-2,
/// the restricted block is unchanged bit for bit at every detune size, while the physical pair
/// splits through the EVEN sector at second order, 9.69·(δJ/J)², already 9.7e-4 apart at δJ/J = 1e-2.</para>
///
/// <para><b>Typed parent.</b> <see cref="F89Path3OcticEpClaim"/>: the same (1,2) coherence block one chain
/// shorter, the same character question, and the same verdict reached the other way, by a scalar
/// compression and live diagnostics rather than by exact arithmetic. This claim is that reading's
/// exact-arithmetic form at N=5.</para>
///
/// <para><b>Live witness.</b> <c>inspect --root n5diabolic</c>
/// (<c>compute/RCPsiSquared.Diagnostics/Foundation/RouteBN5RealQSemisimpleWitness.cs</c>), which pins the
/// coupling book by spectrum rather than assuming it, reads the real-symmetry residual on the Hermitian
/// axis against 0.0 exactly, and recomputes the direction-parity law as a projection norm with no
/// eigensolver in it.</para></summary>
public sealed class RouteBN5RealQSemisimpleClaim : Claim
{
    /// <summary>Parent: the same (1,2) coherence block at N=4, whose on-axis diabolic character was
    /// established by compression and live diagnostics. This claim is the exact-arithmetic form of that
    /// same reading one chain longer.</summary>
    public F89Path3OcticEpClaim OcticSibling { get; }

    public RouteBN5RealQSemisimpleClaim(F89Path3OcticEpClaim octicSibling)
        : base("The N=5 Route-B A2 locus at q = +-1.1292509708747671, lambda = -4.7919603651796410 (the (1,2) " +
               "coherence block of the open XY chain under uniform Z dephasing, R-odd sector, Delta = 0, gamma = 1) " +
               "carries a SEMISIMPLE double eigenvalue, proved by Galois conjugacy and not read off a contour. " +
               "A2_O is irreducible over Q of degree 13, so Q(w0) is a degree-13 subfield of R while t0 = sqrt(-w0) " +
               "is nonzero purely imaginary and lies outside it; hence [Q(t0):Q] = 26 and the 26 lifts form a SINGLE " +
               "Galois orbit, 1 positive w giving 2 imaginary t, 6 negative w giving 12 real t, 6 nonreal w giving " +
               "12 nonreal t. At each of the 12 real t the block L_O = D_O + t*K_O is real symmetric, so geometric " +
               "multiplicity equals algebraic multiplicity there, and rank is preserved by every embedding because a " +
               "minor is zero or nonzero before and after it. The orbit therefore carries the Hermitian-axis verdict " +
               "to the physical point. This is the only Route-B point at N=5 or N=6 with both q and lambda real, so " +
               "the only ones an experiment could be set to; at N=6 none of the 266 is real by the atlas's own float " +
               "midpoint read. It supplies the exact half " +
               "for the 2 real-q loci of the 34 that carried only a stable numerical EpCharacter reading, and the " +
               "other 32 keep the reading. The chain never uses the value 2, which keeps it clear of chi = AT*F_res. " +
               "Not an all-N theorem and no metrological claim.",
               Tier.Tier1Derived,
               "docs/ANALYTICAL_FORMULAS.md F164 + docs/proofs/PROOF_N5_REAL_Q_DIABOLIC.md + " +
               "simulations/n5_a2_rodd_character_gate.py (numerical corroboration, five steps) + " +
               "simulations/o2b_gcd_certificate.py (the exact irreducibility and layer inputs) + " +
               "compute/RCPsiSquared.Diagnostics/Foundation/RouteBN5RealQSemisimpleWitness.cs " +
               "(RouteBN5RealQSemisimpleWitness, inspect --root n5diabolic)")
    {
        OcticSibling = octicSibling ?? throw new ArgumentNullException(nameof(octicSibling));
    }

    public override string DisplayName => "F164: the N=5 real-q A₂ locus is semisimple, by Galois conjugacy";

    public override string Summary =>
        $"the only Route-B operating point with real q and real λ (q = ±1.1292509708747671, λ = −4.7919603651796410) " +
        $"carries a semisimple double eigenvalue: A₂_O irreducible of degree 13 puts it in one Galois orbit with the " +
        $"12 real t where the block is real symmetric, and rank is embedding-invariant. Exact half for 2 of the 34 " +
        $"classified-not-certified loci; the other 32 keep the numerical reading ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("the locus, and why it is the only one",
                summary: "q = ±1.129 250 970 874 767 1, λ = −4.791 960 365 179 641 0, both real, so an operating " +
                         "point one can set rather than an analytic continuation. It is the only such point in the " +
                         "Route-B corpus at N=5 and N=6: at N=5, 2 of the 58 q-loci are strictly real and they are " +
                         "±the same point; at N=6, 0 of 266 are, with 118 at purely imaginary q where the block is " +
                         "Hermitian and 148 generic complex.");

            yield return new InspectableNode("step 1: one Galois orbit",
                summary: "A₂_O is irreducible over ℚ of degree 13, so [ℚ(w₀):ℚ] = 13 and ℚ(w₀) ⊂ ℝ. Then " +
                         "t₀ = √(−w₀) is nonzero purely imaginary (w₀ > 0 and A₂_O does not vanish at w = 0), and no " +
                         "nonzero purely imaginary number lies in a subfield of ℝ, so [ℚ(t₀):ℚ(w₀)] = 2 and " +
                         "[ℚ(t₀):ℚ] = 26. The minimal polynomial of t₀ has degree 26 and divides A₂_O(−t²), which " +
                         "also has degree 26, so they agree up to a constant and the 26 roots form a single orbit. " +
                         "The orbit accounts for every root: 2 imaginary t + 12 real t + 12 nonreal t = 26.");

            yield return new InspectableNode("step 2: the orbit meets the Hermitian axis",
                summary: "A₂_O has 6 negative real roots, and for each of them −w > 0 makes t real. There " +
                         "L_O = D_O + t·K_O is real symmetric, hence Hermitian, hence geometric multiplicity equals " +
                         "algebraic multiplicity for every one of its eigenvalues. Twelve of the orbit's 26 members " +
                         "are such t. Live and exact: inspect --root n5diabolic reads the residual " +
                         "max(|Im L|, |L − Lᵀ|) against 0.0 rather than against a tolerance, at three real t, and the " +
                         "block being affine in t with real coefficients carries that to every real t.");

            yield return new InspectableNode("step 3: the transfer",
                summary: "rank over a number field is preserved by every embedding, since a minor is zero or nonzero " +
                         "before and after it, and algebraic multiplicity is Galois invariant because the " +
                         "characteristic polynomial has rational coefficients. Both properties are therefore constant " +
                         "on the orbit, and the orbit contains both the Hermitian-axis points and the physical one. " +
                         "The transfer is per (t, λ) pair, which is the granularity the conclusion needs.");

            yield return new InspectableNode("what it upgrades, and what stays a reading",
                summary: "the stores classified this locus and did not certify it: the artifact's " +
                         "exactRankCertificates array is empty and the character came from a floating-point contour " +
                         "classifier whose own gates are 1e-4 and 1e-6. This supplies the exact half for the 2 " +
                         "real-q loci of those 34; the other 32 keep only the numerical reading. Not an all-N " +
                         "theorem, and no metrological claim follows from it.");

            yield return new InspectableNode("the sibling at N=6 (F163)",
                summary: "the same method one chain longer and for the whole layer: the R-even degree-133 A₂ " +
                         "factor is irreducible, one real Hermitian embedding carries semisimplicity to each of " +
                         "its conjugates, and a diagonal sign operator exchanging the two reflection sectors " +
                         "transports the verdict to the 133 odd loci, covering all 266. So the N=6 inventory's " +
                         "ExactRankExecuted = 0 is a statement about that classifier's route and not about what " +
                         "is known.");

            yield return new InspectableNode("the gate's error model, and where it stops applying",
                summary: "σ₂ = c·(δJ/J) at w₀, the ratio flat to 1.000 18 over eight decades. The constant is " +
                         "direction-dependent and no parity fixes it: three palindromic directions give " +
                         "c = 5.8391e-2, 8.3097e-2 and 1.1004e-1 (spreads 1.00018, 1.00043, 1.00067), a factor " +
                         "1.885 within one parity. Parity decides only whether a response EXISTS, and that is a " +
                         "theorem: on a reflection-symmetric base an anti-palindromic detune obeys R V R = −V, " +
                         "so U_Oᵀ V U_O vanishes identically at every N, exactly 0.0 in the integer orbit basis; " +
                         "without that symmetry there is no R-odd sector to miss. A blind spot of the sector " +
                         "reading, not a null response: the restricted block is unchanged bit for bit at every " +
                         "detune size, while the physical pair splits through the EVEN sector at second order, " +
                         "9.69·(δJ/J)². The three constants are the gate's readings, not this claim's; what " +
                         "inspect --root n5diabolic recomputes is the zero-or-not projection, exactly, with no " +
                         "eigensolver in it.");

            yield return new InspectableNode("what the gate cannot decide",
                summary: "the nullity reading corroborates and does not carry: a defective point whose Jordan " +
                         "coupling happened to be weak would read nullity 2 as well, so the one-sidedness of that " +
                         "failure is measured at one defective point and is not a bound. The gate is also blind to " +
                         "the commutator sign: flipping it conjugates the R-odd spectrum, and λ₀ is real, so the " +
                         "corroboration fixes the operator only up to conjugation. Neither touches the field " +
                         "argument, whose two inputs come from the exact certificate rather than from that builder.");

            yield return OcticSibling;   // typed parent edge
        }
    }
}
