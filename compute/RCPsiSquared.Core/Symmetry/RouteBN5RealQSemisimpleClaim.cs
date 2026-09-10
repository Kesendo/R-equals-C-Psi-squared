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
/// uniform Z dephasing, R-odd sector, Δ = 0, γ = 1. One root, two loci. They are the only real-q points
/// in the N=5/N=6 A2 inventories; q is the settable parameter and λ is the resulting spectral eigenvalue. At N=6 none of the 266 loci is real,
/// certified by the exact rational <c>tBox.real</c> intervals, all of which exclude zero (q is real iff
/// Re(t)=0), with 118 at purely imaginary q where the block is Hermitian and 148 generic complex.</para>
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
/// <para><b>What is exact, and what remains numerical.</b> Of the N=5 inventory's 58 q-loci, all 24
/// imaginary-q loci are directly Hermitian-certified. F164 certifies the complete 26-member R-odd orbit,
/// overlapping 12 Hermitian loci and adding the 2 real-q plus 12 nonreal-q loci. Thus 38 distinct loci are
/// exact-certified semisimple and only the 20 nonreal-q R-even loci remain numerical-only. The classifier
/// still reads all 32 nonreal-q loci numerically. The artifact's empty <c>exactRankCertificates</c> array
/// records only that producer route; it is not the current knowledge verdict. It is not an all-N theorem
/// and says nothing metrological.</para>
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
/// What parity decides is whether a response exists at all, and that is a theorem: an
/// ANTI-palindromic detune obeys R V R = −V, so P_- V P_- vanishes identically at every N, independently
/// of the base point and exactly 0.0 in the integer orbit basis. Reflection symmetry of the base is needed
/// only so the R-even and R-odd subspaces are invariant under L. It is a BLIND SPOT of the sector reading rather
/// than a null response: PVP = 0 keeps σ₂ on the float floor over the ten decades δJ/J = 1e-12 to 1e-2,
/// the restricted block is unchanged bit for bit at every detune size, while the physical pair
/// splits through the EVEN sector at second order, 9.69·(δJ/J)², already 9.7e-4 apart at δJ/J = 1e-2.</para>
///
/// <para><b>Live witness.</b> <c>inspect --root n5diabolic</c>
/// (<c>compute/RCPsiSquared.Diagnostics/Foundation/RouteBN5RealQSemisimpleWitness.cs</c>), which pins the
/// coupling book by spectrum rather than assuming it, reads the real-symmetry residual on the Hermitian
/// axis against 0.0 exactly, and recomputes the direction-parity law as a projection norm with no
/// eigensolver in it.</para></summary>
public sealed class RouteBN5RealQSemisimpleClaim : Claim
{
    public RouteBN5RealQSemisimpleClaim()
        : base("The N=5 Route-B A2 locus at q = +-1.1292509708747671, lambda = -4.7919603651796410 (the (1,2) " +
               "coherence block of the open XY chain under uniform Z dephasing, R-odd sector, Delta = 0, gamma = 1) " +
               "carries a SEMISIMPLE double eigenvalue, proved by Galois conjugacy and not read off a contour. " +
               "A2_O is irreducible over Q of degree 13, so Q(w0) is a degree-13 subfield of R while t0 = sqrt(-w0) " +
               "is nonzero purely imaginary and lies outside it; hence [Q(t0):Q] = 26 and the 26 lifts form a SINGLE " +
               "Galois orbit, 1 positive w giving 2 imaginary t, 6 negative w giving 12 real t, 6 nonreal w giving " +
               "12 nonreal t. At each of the 12 real t the block L_O = D_O + t*K_O is real symmetric, so geometric " +
               "multiplicity equals algebraic multiplicity there, and rank is preserved by every embedding because a " +
               "minor is zero or nonzero before and after it. The orbit therefore carries the Hermitian-axis verdict " +
               "to the physical point. These are the only real-q points in the N=5/N=6 A2 inventories; q is the " +
               "settable parameter and lambda is the resulting spectral eigenvalue. At N=6 every exact rational tBox.real interval " +
               "excludes zero. Of the N=5 inventory's 58 q-loci, all 24 imaginary-q loci are directly Hermitian-certified. " +
               "F164 certifies the complete 26-member R-odd orbit, overlapping 12 Hermitian loci and adding the 2 real-q " +
               "plus 12 nonreal-q loci. Thus 38 distinct loci are exact-certified semisimple and only the 20 nonreal-q " +
               "R-even loci remain numerical-only; the classifier still reads all 32 nonreal-q loci. An empty " +
               "exactRankCertificates array is classifier-local provenance, not the current verdict. The chain never " +
               "uses the value 2, which keeps it clear of chi = AT*F_res. " +
               "Not an all-N theorem and no metrological claim.",
               Tier.Tier1Derived,
               "docs/ANALYTICAL_FORMULAS.md F164 + docs/proofs/PROOF_N5_REAL_Q_DIABOLIC.md + " +
               "simulations/n5_a2_rodd_character_gate.py (numerical corroboration, five steps) + " +
               "simulations/o2b_gcd_certificate.py (the exact irreducibility and layer inputs) + " +
               "compute/RCPsiSquared.Diagnostics/Foundation/RouteBN5RealQSemisimpleWitness.cs " +
               "(RouteBN5RealQSemisimpleWitness, inspect --root n5diabolic)")
    { }

    public override string DisplayName => "F164: the N=5 real-q A₂ locus is semisimple, by Galois conjugacy";

    public override string Summary =>
        $"the only real-q points in the N=5/N=6 A2 inventories have q = ±1.1292509708747671 and resulting spectral " +
        $"λ = −4.7919603651796410; q is the settable parameter and lambda is the resulting spectral eigenvalue. They carry " +
        $"a semisimple double eigenvalue: A₂_O irreducible of degree 13 puts them in one Galois orbit with the " +
        $"12 real t where the block is real symmetric, and rank is embedding-invariant. N=5 current truth: 24 " +
        $"imaginary-q directly Hermitian-certified; F164 certifies all 26 R-odd loci and yields 38 distinct " +
        $"exact-certified loci in total, leaving 20 nonreal-q R-even loci numerical-only " +
        $"({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("the locus, and why it is the only one",
                summary: "q = ±1.129 250 970 874 767 1 is the settable parameter; λ = −4.791 960 365 179 641 0 " +
                         "is the resulting spectral eigenvalue, not a second control. These are the only real-q points in the " +
                         "N=5/N=6 A2 inventories: at N=5, 2 of the 58 q-loci are strictly real and they are " +
                         "±the same point; at N=6, 0 of 266 are, because every exact rational tBox.real interval " +
                         "excludes zero (q real iff Re(t)=0), with 118 at purely imaginary q where the block is " +
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
                summary: "current truth for the 58 N=5 q-loci: 24 imaginary-q loci are directly Hermitian-certified. " +
                         "F164 certifies the complete 26-member R-odd orbit, overlapping 12 Hermitian loci and adding " +
                         "2 real-q plus 12 nonreal-q loci. Hence 38 distinct loci are exact-certified and only 20 " +
                         "nonreal-q R-even loci remain numerical-only; all 32 nonreal-q loci retain the classifier's " +
                         "numerical reading. The artifact's empty exactRankCertificates array is classifier-local " +
                         "provenance, not the current knowledge verdict. Not an all-N " +
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
                         "theorem: an anti-palindromic detune obeys R V R = −V, so P_- V P_- vanishes identically " +
                         "at every N, independently of the base point and exactly 0.0 in the integer orbit basis. " +
                         "Reflection symmetry of the base is needed only to make the R-even and R-odd subspaces " +
                         "invariant under L. A blind spot of the sector " +
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
        }
    }
}
