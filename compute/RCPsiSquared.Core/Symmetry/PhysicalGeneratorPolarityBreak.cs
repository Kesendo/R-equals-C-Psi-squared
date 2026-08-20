using System.Numerics;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Pauli;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>F155 (Tier1Derived, universal N): the polarity break of the PHYSICAL (no-jump)
/// generator is a DIAGONAL bilinear form in the Pauli basis.
///
/// <para>For H = A + iB with A and B both Hermitian, the un-normalised conditional generator
/// G_H : ρ ↦ −i(Hρ − ρH†) has</para>
///
/// <code>
///   asymmetry(G_H) = 4^(N+1) · Σ_{σ : bit_b(σ) odd} (−1)^(#Z(σ)) · a_σ · b_σ
/// </code>
///
/// <para>where a_σ and b_σ are the Pauli coefficients of A and B (A = Σ_σ a_σ σ over the 4^N
/// strings, likewise B), bit_b(σ) = (#Y(σ) + #Z(σ)) mod 2, and #Z counts the Z letters. The form
/// is DIAGONAL: distinct strings never pair, so only strings carried by BOTH A and B contribute,
/// and bit_b-even content is silent against every B. The engine is that Ad_Π exchanges LEFT with
/// RIGHT multiplication, paying a per-letter sign; see
/// <c>docs/proofs/PROOF_F155_PHYSICAL_GENERATOR_POLARITY_BREAK.md</c> §(d).</para>
///
/// <para><b>Which generator, and it is the distinction the repo's costliest scope error turned
/// on.</b>
/// G_H is NOT the commutator −i[H, ·] that <see cref="LindbladBitBPiBalance"/> (F112) is about,
/// and it is NOT the full Lindbladian: a physical gain medium pumps incoherently, and the
/// recycling term it adds changes the asymmetry outside the σ⁻/σ⁺ family. G_H is the operator
/// PT and gain-loss models are usually written with, and it is not trace-preserving. Conflating
/// the two is recorded twice in <c>docs/CAUGHT_ERRORS.md</c>: the 2026-06-20 entry and its
/// 2026-08-19 sequel, which are one entry with a sequel rather than two separate ones.</para>
///
/// <para><b>The parent edge.</b> <see cref="LindbladBitBPiBreakMagnitude"/> (F113) is the
/// Z-diagonal special case: substituting a_{Z_l} = ω_l/2 and b_{Z_l} = (γ_T1,l − γ_pump,l)/4
/// turns the sum above into F113's (4^N/2)·Σ_l ω_l·(γ_pump,l − γ_T1,l). The two prefactors
/// 4^(N+1) and 4^N/2 differ by exactly the factor 8 that change of variables carries; they do
/// not contradict. Proof corollary 10, gate G7.</para>
///
/// <para><b>The label does not decide; the value does.</b> Put a detuning on ONE site of a
/// gain/loss profile and the asymmetry is nonzero; put the SAME drive on every site and it is
/// exactly 0, since the gain end then carries +ωγ and the loss end −ωγ and the two cancel. Both
/// are live configurations of the witness. (The registry entry's worked instances, 0.64 at N=2
/// and 2.56 at N=3, belong to F113's own fixture at ω₀ = 0.8 and γ = 0.1 on an open XY chain;
/// they are quoted there with those parameters and are not this file's numbers.) The textbook PT
/// dimer is balanced at every coupling INCLUDING its exceptional point (corollary 7): the
/// quantity is a detuning × gain-imbalance meter, not a PT meter.</para>
///
/// <para><b>The value is not basis-invariant</b>, and that is the scope fact a reader most needs:
/// it moves under a per-site Hadamard and under a local rotation (gate G13). The law reads the
/// Pauli content of A and B in the frame the Z-palindrome convention fixes, exactly as F113's
/// number is a lab-frame reading. The SIGN is additionally fixed by the stacking pairing
/// documented on <see cref="PauliBasis"/>; the magnitude and the zero-versus-nonzero verdict are
/// convention-free.</para>
///
/// <para>Live witness: <c>dotnet run --project compute/RCPsiSquared.Cli -- inspect --root f155</c>.</para>
/// </summary>
public sealed class PhysicalGeneratorPolarityBreak : Claim, IZ2AxisClaim
{
    /// <summary>Typed parent: F113, whose closed form is this law's Z-diagonal special case.</summary>
    public LindbladBitBPiBreakMagnitude BreakMagnitude { get; }

    public PhysicalGeneratorPolarityBreak(LindbladBitBPiBreakMagnitude breakMagnitude)
        : base("F155: for H = A + iB with A, B Hermitian, the physical (no-jump) generator "
             + "ρ ↦ −i(Hρ − ρH†) has polarity asymmetry = 4^(N+1) · Σ over bit_b-ODD Pauli "
             + "strings of (−1)^#Z · a_σ · b_σ, bilinear in (A, B) and DIAGONAL in the Pauli "
             + "basis; bit_b-odd content in A is necessary and NOT sufficient. Universal N.",
               Tier.Tier1Derived,
               "docs/ANALYTICAL_FORMULAS.md F155 + "
             + "docs/proofs/PROOF_F155_PHYSICAL_GENERATOR_POLARITY_BREAK.md + "
             + "simulations/f155_polarity_break_bilinear.py + "
             + "simulations/f155_dephase_siblings.py")
    {
        BreakMagnitude = breakMagnitude ?? throw new ArgumentNullException(nameof(breakMagnitude));
    }

    public override string DisplayName =>
        "F155 (physical generator polarity break: one diagonal bilinear form)";

    /// <summary>The theorem in one line.</summary>
    public string Theorem =>
        "asymmetry(−i(Hρ − ρH†)) = 4^(N+1) · Σ_{bit_b(σ) odd} (−1)^#Z(σ) · a_σ · b_σ "
      + "for H = A + iB, A and B Hermitian";

    public override string Summary => $"{Theorem} ({Tier.Label()})";

    /// <summary>The Π²-axis this claim's palindrome convention sits on: Π_Z grades by bit_b.</summary>
    public Z2Axis Z2Axis => Z2Axis.BitB;

    /// <summary>No bespoke typed bit_a twin is owed, so the pointer stays null; see
    /// <see cref="BitATwinStatus"/> for why.</summary>
    public Claim? BitATwin => null;

    /// <summary>The bit_a twin is a corollary of the global Hadamard X ↔ Z duality
    /// (<c>docs/proofs/PROOF_BIT_A_TWIN_VIA_HADAMARD.md</c>, case 1: the asymmetry is a
    /// Frobenius-norm statement on the 4^N Pauli space, and Ad_{Q_zx} is Frobenius-unitary with
    /// Q_zx·Π_Z·Q_zx⁻¹ = Π_X). Proof §(g) states the transport as an identity,
    /// asymmetry_X(A, B) = asymmetry_Z(hAh, hBh) for h = H^⊗N with EVERY site rotated, and gate
    /// S10 measures it.
    ///
    /// <para><b>Why this differs from F113's BitBSpecific, which is not a contradiction.</b> That
    /// classification is about a physical CHANNEL: σ⁺/σ⁻ jumps single out the Z-basis energy
    /// ladder, and the Hadamard relabelling maps them to (Z ± iY)/2, a different channel outside
    /// the standard damping family. F155 makes no channel assumption at all: A and B are arbitrary
    /// Hermitian operators, so the statement is operator-space and transports. The two verdicts
    /// are about two different objects.</para>
    ///
    /// <para>The Y leg is stronger than "no twin" and is not a slot at all: Π_Y = Π_Z⁻¹ exactly,
    /// so asymmetry_Y(M) = −asymmetry_Z(M) for EVERY superoperator M, with no hypothesis (proof
    /// §(g), gate S4, bit-exact on fifteen arbitrary superoperators).</para></summary>
    public BitATwinClassification BitATwinStatus =>
        BitATwinClassification.CoveredByHadamardDuality;

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("theorem", summary: Theorem);

            yield return new InspectableNode(
                "which generator",
                summary: "−i(Hρ − ρH†), the un-normalised conditional (no-jump) generator. NOT "
                       + "the commutator −i[H, ·] of F112, and NOT the full Lindbladian: the "
                       + "recycling term c ⊗ c* is asymmetry-inert only for the σ⁻/σ⁺ family "
                       + "(measured N = 2, 3, 4) and the X/Y-dephase break of the full "
                       + "Lindbladian remains open.");

            yield return new InspectableNode(
                "diagonality",
                summary: "distinct Pauli strings never pair, so only strings carried by BOTH A "
                       + "and B contribute; measured on all 210 ordered pairs of DISTINCT "
                       + "non-identity strings at N=2 (15·14) and all 3906 at N=3 (63·62), "
                       + "gate G4.");

            yield return new InspectableNode(
                "typed parent",
                summary: $"LindbladBitBPiBreakMagnitude ({BreakMagnitude.Tier.Label()}): F113 is "
                       + "the Z-diagonal special case, reached by a_{Z_l} = ω_l/2 and "
                       + "b_{Z_l} = (γ_T1,l − γ_pump,l)/4; the prefactors 4^(N+1) and 4^N/2 "
                       + "differ by exactly the factor 8 that substitution carries (corollary 10, "
                       + "gate G7).");

            yield return new InspectableNode(
                "the sorting is one-directional",
                summary: "bit_b-odd content in A is NECESSARY and NOT sufficient: B = 0 (Hermitian "
                       + "H), B ∝ I, and disjoint odd supports each give exactly 0 with A as odd "
                       + "as one likes (corollaries 1, 3, 4).");

            yield return new InspectableNode(
                "the PT dimer is balanced at its own exceptional point",
                summary: "J(XX+YY)/2 + i·g(Z₀ − Z₁)/2 reads exactly 0 at every g including g = 1, "
                       + "and the same dimer with a +0.4·Z₀ detuning reads exactly −12.8·g. The "
                       + "quantity is a detuning × gain-imbalance meter, not a PT meter, and it "
                       + "is blind to the PT-breaking threshold (corollary 7, gate G11).");

            yield return new InspectableNode(
                "the bit_a reading, and the third answer",
                summary: "X is this theorem in the Hadamard-rotated frame "
                       + "(asymmetry_X(A, B) = asymmetry_Z(hAh, hBh), every site rotated, gate "
                       + "S10); Y is not a sibling at all, since Π_Y = Π_Z⁻¹ gives "
                       + "asymmetry_Y = −asymmetry_Z for every superoperator (gate S4). The F113 "
                       + "configuration read under Π_X is SILENT rather than balanced: both ±i "
                       + "norms are individually and exactly zero, as is any bit_a-even content "
                       + "at any N (gate S8).");

            yield return new InspectableNode(
                "frame",
                summary: "the value is NOT basis-invariant (gate G13) and must be quoted in the "
                       + "frame the Z-palindrome convention fixes; the sign is additionally fixed "
                       + "by the stacking pairing PauliBasis documents. Magnitude and the "
                       + "zero-versus-nonzero verdict are convention-free.");
        }
    }

    /// <summary>The closed form, evaluated from the Pauli coefficients of A and B directly:
    /// 4^(N+1) · Σ over strings σ whose Π_ℓ²-parity is ODD of (−1)^(#ℓ(σ)) · a_σ · b_σ, where ℓ is
    /// <paramref name="dephaseLetter"/>. This is the theorem's own arithmetic and touches no
    /// matrix: the sum runs over the strings <paramref name="bTerms"/> carries, because the form
    /// is diagonal and a string absent from either operator contributes nothing.
    ///
    /// <para>The support is Π_ℓ²-odd (bit_b for ℓ = Y and Z, bit_a for ℓ = X) and the weight
    /// counts the dephasing letter itself, both read off proof §(g)'s table through
    /// <see cref="PiOperator.SquaredEigenvalue"/> rather than hard-coded here.</para>
    ///
    /// <para><b>Coefficients must be exactly real.</b> A Hermitian operator has real Pauli
    /// coefficients, so a nonzero imaginary part means the caller has not split H into its
    /// Hermitian halves and the theorem's hypothesis does not hold. That is refused rather than
    /// silently truncated.</para>
    ///
    /// <para>Repeated strings are summed, so a caller may pass an unreduced term list. The sum is
    /// accumulated in the order <paramref name="bTerms"/> is given, so the float result is a
    /// deterministic function of the input order.</para></summary>
    /// <param name="aTerms">Pauli terms of A, the Hermitian part of H.</param>
    /// <param name="bTerms">Pauli terms of B, the Hermitian part of −i·(H − H†)/2.</param>
    /// <param name="n">Qubit count; every term must carry exactly N letters.</param>
    /// <param name="dephaseLetter">The palindrome convention; Z is the theorem of §(a).</param>
    public static double PredictAsymmetry(
        IReadOnlyList<PauliTerm> aTerms,
        IReadOnlyList<PauliTerm> bTerms,
        int n,
        PauliLetter dephaseLetter = PauliLetter.Z)
    {
        if (aTerms is null) throw new ArgumentNullException(nameof(aTerms));
        if (bTerms is null) throw new ArgumentNullException(nameof(bTerms));
        if (n < 1) throw new ArgumentOutOfRangeException(nameof(n), $"N must be ≥ 1; got {n}");
        if (dephaseLetter == PauliLetter.I)
            throw new ArgumentException(
                $"dephaseLetter must be X, Y, or Z; got {dephaseLetter}", nameof(dephaseLetter));

        var aByString = new Dictionary<long, double>();
        foreach (var term in aTerms)
        {
            double coefficient = RealCoefficient(term, n, nameof(aTerms));
            long key = PauliIndex.ToFlat(term.Letters);
            aByString[key] = aByString.TryGetValue(key, out double running)
                ? running + coefficient
                : coefficient;
        }

        double sum = 0.0;
        foreach (var term in bTerms)
        {
            double bCoefficient = RealCoefficient(term, n, nameof(bTerms));
            if (!aByString.TryGetValue(PauliIndex.ToFlat(term.Letters), out double aCoefficient))
                continue;
            if (PiOperator.SquaredEigenvalue(term.Letters, dephaseLetter) != -1) continue;

            int letterCount = CountLetter(term.Letters, dephaseLetter);
            double weight = (letterCount & 1) == 0 ? 1.0 : -1.0;
            sum += weight * aCoefficient * bCoefficient;
        }

        return Math.Pow(4.0, n + 1) * sum;
    }

    private static double RealCoefficient(PauliTerm term, int n, string parameterName)
    {
        if (term is null)
            throw new ArgumentException("term list contains a null term", parameterName);
        if (term.N != n)
            throw new ArgumentException(
                $"term '{term.Label}' has N={term.N} letters; expected N={n}", parameterName);
        if (term.Coefficient.Imaginary != 0.0)
            throw new ArgumentException(
                $"term '{term.Label}' has coefficient {term.Coefficient}, whose imaginary part is "
              + "nonzero; A and B are Hermitian and therefore carry exactly real Pauli "
              + "coefficients. Split H into A = (H + H†)/2 and B = −i·(H − H†)/2 first.",
                parameterName);
        return term.Coefficient.Real;
    }

    private static int CountLetter(IReadOnlyList<PauliLetter> letters, PauliLetter letter)
    {
        int count = 0;
        foreach (var l in letters) if (l == letter) count++;
        return count;
    }
}
