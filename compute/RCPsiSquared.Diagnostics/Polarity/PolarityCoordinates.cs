using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.ChainSystems;
using RCPsiSquared.Core.Lindblad;
using RCPsiSquared.Core.Pauli;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Diagnostics.F81;
using ComplexMatrix = MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex>;

namespace RCPsiSquared.Diagnostics.Polarity;

/// <summary>Three-way polarity decomposition of the F1 residual
/// M = Π·L·Π⁻¹ + L + 2σ·I into the polarity-triple {−1/2, 0, +1/2} at d=2:
///
/// <code>
///   M_zero       = (M + Π·M·Π⁻¹) / 2                      (0-axis, = F81 M_sym)
///   M_plus_half  = (M_anti − i · Π·M_anti·Π⁻¹) / 2        (+1/2 polarity)
///   M_minus_half = (M_anti + i · Π·M_anti·Π⁻¹) / 2        (−1/2 polarity)
/// </code>
///
/// <para>where M_anti = (M − Π·M·Π⁻¹) / 2 is the F81 antisymmetric part.</para>
///
/// <para><b>These three are NOT the eigenspaces of the conjugation map.</b> A := Ad_Π has
/// order 4 on Liouville space, so its spectrum is {+1, −1, +i, −i} and (1 ± A)/2 are not
/// eigenprojections. Writing M = u + v + p + m for the true A-eigencomponents (+1, −1, +i,
/// −i), what is built above is M_zero = u + (1+i)p/2 + (1−i)m/2, M_plus_half = (1+i)v/2 +
/// (1−i)p/2, M_minus_half = (1−i)v/2 + (1+i)m/2. So M_zero is not the Π²-even part (for
/// XY+YX at N=3 it is entirely Π²-ODD), and the three pieces are not mutually
/// Frobenius-orthogonal. The Π²-even projector is (1 + A²)/2.</para>
///
/// <para>What survives untouched is the quantity everything downstream reads:
/// <see cref="PolarityCoordinatesResult.Asymmetry"/> = ½(‖P_{+i}M‖² − ‖P_{−i}M‖²), since
/// the v contamination enters both halves equally and cancels. F112's zero and F113's
/// (4^N/2) are unaffected.</para>
///
/// <para><b>The norm-sum identity is not an orthogonality check</b>:
/// ‖M‖² = ‖M_zero‖² + ‖M_plus_half‖² + ‖M_minus_half‖² holds algebraically for ANY unitary
/// Π and ANY M, the cross terms cancelling identically, so
/// <see cref="PolarityCoordinatesResult.OrthogonalityResidual"/> is machine zero by
/// construction and cannot fail. It checks the arithmetic wiring, not Π.</para>
///
/// <para><b>Connection to F81</b>:
/// F81 M_sym = M_zero;
/// F81 M_anti = M_plus_half + M_minus_half (further split by Π ±i eigenvalue).</para>
///
/// <para><b>F112 typed scope (Tier1Derived)</b>: for any Lindblad-form Liouvillian
/// L = -i[H, ·] + Σ_k γ_k · D[c_k] with Hermitian H and each c_k bit_b-homogeneous,
/// <see cref="PolarityCoordinatesResult.Asymmetry"/> = ‖M_plus_half‖² − ‖M_minus_half‖²
/// is exactly 0 bit-exact. See <c>docs/proofs/PROOF_F112_LINDBLAD_BIT_B_PI_BALANCE.md</c>
/// and <see cref="LindbladBitBPiBalance"/>. The diagnostic asymmetry is the precise
/// witness for c with cross-bit_b Pauli support (outside F108's closure regime) or for
/// non-Hermitian H combined with bit_b-mixed c.</para>
///
/// <para>Python sibling: <c>simulations/framework/diagnostics/polarity_coordinates.py</c>
/// (<c>polarity_coordinates_from_L</c>, <c>polarity_coordinates_from_hc</c>); this C#
/// primitive matches the Python output bit-exactly via the same Π construction
/// (<see cref="PiOperator.BuildFull"/>) and the same vec→Pauli transform
/// (<see cref="PauliBasis.VecToPauliBasisTransform"/>).</para>
/// </summary>
/// <param name="F81Violation">F81 closed-form ‖M_anti − L_{H_odd}‖_F on Z-dephase; numerical diagnostic only on X / Y (see <see cref="PiDecomposition.Decompose"/>).</param>
/// <param name="IsStructurallySilent">The EXACT verdict of
/// <see cref="PolarityCoordinates.IsStructurallyDegenerate"/> on the inputs this result was
/// built from, carried in because the record holds only matrices and cannot re-derive it from
/// Pauli letters it never saw. <c>null</c> on any path that had no term list, and the float
/// <see cref="PolarityCoordinatesResult.MAntiNormSquared"/> == 0.0 test governs there.</param>
public sealed record PolarityCoordinatesResult(
    ComplexMatrix M,
    ComplexMatrix MZero,
    ComplexMatrix MPlusHalf,
    ComplexMatrix MMinusHalf,
    double MNormSquared,
    double MZeroNormSquared,
    double MPlusHalfNormSquared,
    double MMinusHalfNormSquared,
    double Asymmetry,
    double OrthogonalityResidual,
    double F81Violation,
    double? LHOddNormSquared = null,
    bool? IsStructurallySilent = null)
{
    /// <summary>‖M_anti‖² = ‖M_plus_half‖² + ‖M_minus_half‖², the total polarity content:
    /// the whole of what <see cref="Asymmetry"/> is a difference of, and equal to M_anti's own
    /// norm². With A = M_anti and B = Π·A·Π⁻¹ the halves are (A ∓ iB)/2, so the parallelogram
    /// law gives ‖(A−iB)/2‖² + ‖(A+iB)/2‖² = (‖A‖² + ‖B‖²)/2, and Π unitary gives ‖B‖² = ‖A‖²;
    /// the sum is therefore ‖A‖² in exact arithmetic. In FLOAT it is not bit-identical to a
    /// direct ‖M_anti‖² and must not be claimed as such: each half costs a separate
    /// <c>Math.Pow(FrobeniusNorm(), 2)</c> (two roundings), so the two routes agree to a few ulp
    /// (measured 0 to 12 ulp over the Heisenberg / YZ+ZY / XY / Heisenberg+T1 / Z-drive+T1 cases
    /// at N=2 and N=3), and agree exactly only where both sides are 0.0 or land on a power of
    /// two. An earlier version of this comment said "bit-identical"; that was read off fixtures
    /// whose halves are exact powers of two, which cannot break it.</summary>
    public double MAntiNormSquared => MPlusHalfNormSquared + MMinusHalfNormSquared;

    /// <summary>True when L carries no Π²-odd content at all, so there is no polarity content
    /// and the balance <see cref="RelativeAsymmetry"/> would report is vacuous: it compares 0
    /// against 0.
    ///
    /// <para><b>STRUCTURAL first, float second, and the order is the whole point.</b>
    /// <see cref="IsStructurallySilent"/> is decided from the Pauli letters
    /// (<see cref="PolarityCoordinates.IsStructurallyDegenerate"/>), is exact and
    /// N-independent, and a <c>true</c> from it settles the question on its own. The float
    /// ‖M_anti‖² == 0.0 stays in the disjunction beneath it (and is the only test left where
    /// the flag is <c>null</c>, i.e. on a path with no term list), because the structural test
    /// has known false negatives, cancelling +c/−c pairs among them. That float test is
    /// sufficient but NOT necessary: ‖M_anti‖² is a difference of two rounded matrices and
    /// carries rounding where the physics is zero, so on Heisenberg + Z-dephasing it reaches
    /// exact 0.0 at N = 2 alone (the per-N figures are quoted once, on
    /// <see cref="PolarityCoordinates.IsStructurallyDegenerate"/>, rather than twice here).
    /// Read alone it is an N=2 guard. Mirrors the Python <c>polarity_degenerate</c> in
    /// <c>simulations/framework/workflows/polarity_fingerprint.py</c>, which composes the
    /// same two tests in the same order.</para>
    ///
    /// <para>This is the generic case for Π²-even H, not an edge case: for <i>truly</i> H
    /// (Heisenberg, XX, ZZ), where M itself vanishes by the F83 closed form
    /// (<c>docs/ANALYTICAL_FORMULAS.md</c>, F83 special-cases table: "truly | undefined
    /// (M=0)"), and equally for <i>non-truly Π²-even</i> H such as YZ+ZY, where ‖M‖² is large
    /// (256 at N=2, J=1, γ-independent) while M_anti still vanishes. That second family is why
    /// ‖M‖² is the wrong scale for <see cref="RelativeAsymmetry"/>: it can be far from zero
    /// while the object the asymmetry measures is empty.</para>
    ///
    /// <para><b>Not the same question as "H has no odd part"</b>, which is
    /// <see cref="HasNoPi2OddPart"/>: M_anti = L_{H_odd} is a Z-dephase, no-T1 identity.
    /// Heisenberg WITH γ_T1 = 0.1 measures ‖M_anti‖² = 8.0e-2 and is NOT degenerate although
    /// H_odd = 0, and the cross-dephase (Π_X / Π_Y) paths have no closed-form F81
    /// prediction.</para></summary>
    public bool IsPolarityDegenerate =>
        (IsStructurallySilent ?? false) || MAntiNormSquared == 0.0;

    /// <summary>The EXACT structural companion to <see cref="IsPolarityDegenerate"/>, and a
    /// statement about H alone: <c>true</c> when H carries no Π²-odd part, i.e.
    /// <see cref="LHOddNormSquared"/> is exactly 0.0. <c>null</c> where no H_odd projection was
    /// made, which is every path except the bilinear-bond one.
    ///
    /// <para>Why both exist. This one is exact but answers only "does H have an odd part";
    /// <see cref="IsPolarityDegenerate"/> is a float reading of the M_anti that was actually
    /// built, and the two come apart in BOTH directions. At N=3, γ_z=0.1, YZ+ZY: this is
    /// <c>true</c> while <see cref="IsPolarityDegenerate"/> is <c>false</c>, because M_anti
    /// carries 5.8e-33 of rounding where the physics is zero. With γ_T1 = 0.1 on Heisenberg:
    /// this is <c>true</c> while ‖M_anti‖² = 8.0e-2 is genuinely non-zero, because
    /// M_anti = L_{H_odd} is an F81 identity that holds on the Z-dephase, no-T1 path and
    /// nowhere else. So a <c>true</c> here does NOT license skipping the division; read the two
    /// together, and read this one as "H is the reason" rather than "the number is zero".</para></summary>
    public bool? HasNoPi2OddPart => LHOddNormSquared is double lhOdd ? lhOdd == 0.0 : null;

    /// <summary>Contrast ratio |Asymmetry| / (‖M_plus_half‖² + ‖M_minus_half‖²), in [0, 1]:
    /// the imbalance of the two polarity halves relative to their own total. No tolerance
    /// constant enters the denominator.
    ///
    /// <para><b>It is NOT scale-free in H.</b> The numerator is linear in the drive by F113
    /// while the denominator is quadratic, so the ratio goes like 1/ω: measured on the
    /// Z-drive + T1 fixture, scaling ω by 0.5 / 1 / 2 / 4 gives 1.538279e-2, 7.692080e-3,
    /// 3.846125e-3, 1.923073e-3. A WEAKER physical break therefore reads as a LARGER ratio.
    /// Use it to compare against a threshold at fixed parameters, never to rank two drives.
    /// (The retired ‖M‖² denominator was not scale-free either; this is a property to know,
    /// not a regression.)</para>
    ///
    /// <para><b>The degenerate branch returns NaN, not 0.0</b>, since 2026-08-07. Where
    /// <see cref="IsPolarityDegenerate"/> holds there is nothing for the two halves to be a
    /// contrast BETWEEN, and the quotient is 0/0: NaN is the entry that says so. A 0.0 there
    /// reads as "perfectly balanced", which is precisely the vacuous confirmation the
    /// DEGENERATE verdict exists to refuse, and it is the same move as the retired
    /// max(‖M‖², 1e-15) floor, a comfortable number in place of an uninterpretable one. NaN
    /// also cannot be aggregated away in silence: it loses every comparison, so a mean or a
    /// max over mixed rows surfaces it instead of hiding it behind a finite neighbour.
    /// <b>Branch on <see cref="IsPolarityDegenerate"/>, not on <c>double.IsNaN</c> of this
    /// value</b>: the implication runs one way only. A NaN rate or a NaN entry in H also lands
    /// a NaN here, on a row that is NOT silent and whose verdict is a confident "BROKEN"
    /// (<c>NaN &lt;= tolerance</c> is false), and IsNaN cannot tell the two apart while the bool
    /// can. The Python sibling carries the same warning and refuses a NaN gamma_T1 or
    /// gamma_pump at its entry point, though not a NaN gamma_z; on this side there is no guard
    /// at all. (Underflow caveat, unchanged:
    /// <c>Math.Pow(x, 2)</c> is 0.0 for x below about 2.2e-162, so a preposterously small J
    /// would set the float half of the flag on a non-zero M_anti.)</para>
    ///
    /// <para>Replaces the earlier |Asymmetry| / max(‖M‖², 1e-15), which was saturated: on the
    /// flagship Heisenberg + Z-dephasing case ‖M‖² is exactly 0.0 at N=2 and ~1e-31 at N=3..5,
    /// so the 1e-15 floor was the divisor and a threshold of 1e-10 admitted any |Asymmetry|
    /// below 1e-25. The two denominators are related in closed form by F83's anti-fraction:
    /// ‖M‖² / ‖M_anti‖² = 2 + 4r with r = ‖H_even_nontruly‖² / ‖H_odd‖², so
    /// <c>new = old · (2 + 4r)</c> converts a legacy number without a re-run: measured exactly
    /// 2.0 at r=0 and 6.0 at r=1 for N=2,3,4. <b>Scope</b>: that is the pure-dephasing
    /// bilinear family F83 models. With T1 the identity does not hold, because the damping
    /// puts content into M that the anti-fraction does not describe; on the Z-drive + T1
    /// witness the measured factor is 2.000266, not 2.</para>
    ///
    /// <para>Prior art, and it predates this property: <c>simulations/polarity_step5_stress.py</c>
    /// line 86 has always divided by <c>norm_plus_i + norm_minus_i</c>, the same polarity
    /// content. That script is the Step-5 verification cited by
    /// <c>docs/proofs/PROOF_F112_LINDBLAD_BIT_B_PI_BALANCE.md</c>.</para></summary>
    public double RelativeAsymmetry =>
        IsPolarityDegenerate ? double.NaN : Math.Abs(Asymmetry) / MAntiNormSquared;
}

/// <summary>Matrix-level compute primitive for the F112 polarity decomposition.
/// Mirrors the Python <c>polarity_coordinates_from_L</c> family by delegating M /
/// M_sym / M_anti / L_HOdd construction to <see cref="PiDecomposition.Decompose"/>
/// and adding the +i / −i Π-eigenvalue refinement of M_anti on top.</summary>
public static class PolarityCoordinates
{
    /// <summary>EXACT structural test for "this configuration has no polarity content to
    /// read", decided from Pauli letters and integer parities alone. No float, no tolerance.
    ///
    /// <para>M_anti is the Π²-odd part of L, M_anti = (L − Π²LΠ⁻²)/2, so it vanishes exactly
    /// when L carries no Π²-odd content, and that needs BOTH halves: H must be Π²-even on the
    /// axis AND every collapse operator must be Π²-homogeneous on it. Where this returns true
    /// the asymmetry is 0 − 0 and the relative asymmetry is a ratio of two zeros: the
    /// configuration is not balanced, it is SILENT, and a "BALANCED" reading there is vacuous.</para>
    ///
    /// <para>Why this and not <see cref="PolarityCoordinatesResult.IsPolarityDegenerate"/>: that
    /// one reads the float M_anti actually built, and the structurally-zero case only reaches
    /// exact 0.0 at N = 2. Measured on Heisenberg + Z-dephasing through the PYTHON sibling at
    /// its own defaults (γ₀ = 0.05, Pauli-J convention): 0.0 at N=2, but 1.4e-33 at N=3,
    /// 8.5e-33 at N=4 and 2.2e-31 at N=5. Those are the sibling's numbers, quoted for the
    /// SHAPE (exact only at N=2, rounding above it) and not as this path's measurements, which
    /// build H at J/4 with a different γ default; what is pinned on this side is the shape
    /// itself, by Decompose_SilentAtNAboveTwo_... at N = 2, 3 and 4. Either way the float guard
    /// is silent from N=3 up and a noise/noise ratio can exceed any threshold. This predicate
    /// is N-independent.</para>
    ///
    /// <para>Why not <see cref="PolarityCoordinatesResult.HasNoPi2OddPart"/>: that is a statement
    /// about H ALONE (‖L_HOdd‖² == 0) and cannot see the c side. Heisenberg has H_odd = 0 either
    /// way, yet ‖M_anti‖² is 0 without T1 and 8.0e-2 with γ_T1 = 0.1, because σ⁻ is bit_b-mixed.</para>
    ///
    /// <para>The axis matters and the asymmetry between the two axes is real, not a convention:
    /// σ⁻ = (X ∓ iY)/2 has Pauli support {X, Y}, which is bit_b-MIXED (X and Y disagree) but
    /// bit_a-HOMOGENEOUS (both are 1). So amplitude damping rescues a Π²-even H from silence on
    /// the Π_Z / Π_Y axis and does NOT on the Π_X axis. Both facts fall out of
    /// <see cref="PiOperator.SquaredEigenvalue"/> rather than being hard-coded here.</para>
    ///
    /// <para><b>Detailed balance is a third silent case this cannot reach, and deliberately
    /// carries no parameter for it.</b> σ⁻ and σ⁺ are each Π²-MIXED on the Z axis, yet their
    /// Π²-odd contributions are equal and opposite, so at γ_T1 = γ_pump SITE BY SITE they
    /// cancel and M_anti vanishes after all. The Python sibling
    /// (<c>framework.diagnostics.polarity_coordinates.is_structurally_degenerate</c>) takes
    /// both rate profiles and tests that elementwise. No path INTO this predicate can build
    /// such an L, and that narrower statement is the whole justification. This method is
    /// public and the witnesses call it directly, so the claim is not "every caller comes
    /// through <see cref="PolarityCoordinates.Decompose"/>"; it is that every place a
    /// LIOUVILLIAN is built on this axis does, and Decompose reaches only
    /// <see cref="PauliDephasingDissipator.BuildZ"/> and <see cref="T1Dissipator.Build"/>, the
    /// latter taking γ_Z and γ_T1 only. A caller passing <c>hasAmplitudeDamping: true</c> is
    /// describing a σ⁻, because that is the only damping the pipeline can have produced. A γ_pump parameter here would therefore be
    /// unreachable code that reads as coverage. The REPO is not without a σ⁺: Core carries a
    /// typed <c>HardwareDissipators.T1Pump</c> (σ⁺ = (X − iY)/2) and
    /// <see cref="LindbladianBuilder.Build"/> takes arbitrary collapse operators, so the pump
    /// arrives here the moment a Decompose overload accepts one, and this predicate gains the
    /// parameter in that same change.</para>
    /// </summary>
    /// <param name="hTerms">The Hamiltonian's Pauli terms.</param>
    /// <param name="hasAmplitudeDamping">Whether a σ⁻ (or σ⁺) channel is present. The
    /// Z-dephasing collapse operators are single-letter and therefore homogeneous on every
    /// axis; only amplitude damping can break c-homogeneity.</param>
    /// <param name="dephaseLetter">The Π axis: X grades by bit_a, Y and Z by bit_b.</param>
    public static bool IsStructurallyDegenerate(
        IReadOnlyList<PauliTerm> hTerms,
        bool hasAmplitudeDamping,
        PauliLetter dephaseLetter = PauliLetter.Z)
    {
        if (hTerms is null) throw new ArgumentNullException(nameof(hTerms));

        // H side: every term that actually CONTRIBUTES must be Π²-even on this axis. A term
        // with a zero coefficient is not in H at all, and skipping it is not a nicety: a first
        // version of this test read the letters only, so appending a zero-weighted Π²-odd term
        // flipped the answer to "not silent" and handed the verdict straight back to the
        // noise/noise ratio it exists to bypass. Cancelling pairs (+c and −c on the same
        // string) are NOT handled here and remain a known false negative; they fail safe, in
        // that the ratio is consulted rather than a silence being asserted wrongly.
        foreach (var term in hTerms)
        {
            if (term is null || term.Letters is null) return false;
            if (term.Coefficient == Complex.Zero) continue;
            if (PiOperator.SquaredEigenvalue(term.Letters, dephaseLetter) != +1) return false;
        }

        // c side: Z-dephasing is single-letter, hence homogeneous on every axis. Amplitude
        // damping carries support {X, Y} and is homogeneous only where those two agree.
        if (hasAmplitudeDamping)
        {
            int sx = PiOperator.SquaredEigenvalue(new[] { PauliLetter.X }, dephaseLetter);
            int sy = PiOperator.SquaredEigenvalue(new[] { PauliLetter.Y }, dephaseLetter);
            if (sx != sy) return false;
        }

        return true;
    }

    /// <summary>The same exact predicate for the bilinear bond-term representation; see the
    /// <see cref="IsStructurallyDegenerate(IReadOnlyList{PauliTerm}, bool, PauliLetter)"/>
    /// overload for what it decides and why.
    ///
    /// <para>A <see cref="PauliPairBondTerm"/> carries no coefficient of its own (every bond
    /// term enters at the chain's J), so <paramref name="coupling"/> supplies the one thing the
    /// letters cannot: at J = 0 the Hamiltonian is the zero operator and therefore silent, and
    /// without the parameter this overload would answer "not silent" on a Π²-odd letter pair
    /// where the <see cref="PauliTerm"/> overload (which skips zero-coefficient terms) answers
    /// "silent". Two public entry points must not give opposite exact verdicts on the same
    /// physics, so the coupling is required rather than defaulted.</para></summary>
    /// <param name="bondTerms">The Hamiltonian's bilinear bond terms.</param>
    /// <param name="coupling">The coupling every bond term enters at, i.e. the chain's J.
    /// A zero here means H = 0, which is silent whatever the letters say.</param>
    /// <param name="hasAmplitudeDamping">As in the <see cref="PauliTerm"/> overload.</param>
    /// <param name="dephaseLetter">As in the <see cref="PauliTerm"/> overload.</param>
    public static bool IsStructurallyDegenerate(
        IReadOnlyList<PauliPairBondTerm> bondTerms,
        double coupling,
        bool hasAmplitudeDamping,
        PauliLetter dephaseLetter = PauliLetter.Z)
    {
        if (bondTerms is null) throw new ArgumentNullException(nameof(bondTerms));

        if (coupling != 0.0)
        {
            foreach (var term in bondTerms)
            {
                if (term is null) return false;
                if (PiOperator.SquaredEigenvalue(
                        new[] { term.LetterA, term.LetterB }, dephaseLetter) != +1) return false;
            }
        }

        if (hasAmplitudeDamping)
        {
            int sx = PiOperator.SquaredEigenvalue(new[] { PauliLetter.X }, dephaseLetter);
            int sy = PiOperator.SquaredEigenvalue(new[] { PauliLetter.Y }, dephaseLetter);
            if (sx != sy) return false;
        }

        return true;
    }

    /// <summary>Decompose M for a chain-built Hamiltonian from bilinear bond terms.
    /// Delegates to <see cref="PiDecomposition.Decompose"/> for M / M_sym / M_anti /
    /// L_HOdd / F81Violation, then computes the +i / −i Π-eigenvalue refinement
    /// of M_anti.</summary>
    /// <param name="chain">N-qubit chain providing N, bonds, J, γ₀.</param>
    /// <param name="terms">Bilinear bond Pauli-pair terms; see
    /// <see cref="PiDecomposition.Decompose"/>.</param>
    /// <param name="gammaT1">Optional uniform per-site σ⁻ amplitude-damping rate.
    /// When non-null and non-zero, the Lindbladian uses the
    /// <see cref="T1Dissipator"/> path; the F112 typed Tier1Derived scope no longer
    /// covers this case (σ⁻ is bit_b-mixed), and the asymmetry need not be zero.</param>
    /// <param name="dephaseLetter">Dephasing letter for the Π operator; default Z.</param>
    public static PolarityCoordinatesResult Decompose(
        ChainSystem chain,
        IReadOnlyList<PauliPairBondTerm> terms,
        double? gammaT1 = null,
        PauliLetter dephaseLetter = PauliLetter.Z)
    {
        if (chain is null) throw new ArgumentNullException(nameof(chain));
        if (terms is null) throw new ArgumentNullException(nameof(terms));

        IReadOnlyList<double>? gammaT1PerSite = null;
        if (gammaT1.HasValue && gammaT1.Value != 0.0)
            gammaT1PerSite = Enumerable.Repeat(gammaT1.Value, chain.N).ToArray();

        var pi = PiDecomposition.Decompose(chain, terms, gammaT1PerSite, dephaseLetter: dephaseLetter);
        var piOp = PiOperator.BuildFull(chain.N, dephaseLetter);
        return RefineMAnti(pi.M, pi.MSym, pi.MAnti, piOp, pi.F81Violation,
            lHOddNormSquared: pi.LHOddNormSquared,
            isStructurallySilent: IsStructurallyDegenerate(
                terms, chain.J, gammaT1PerSite is not null, dephaseLetter));
    }

    /// <summary>Decompose M for a general k-body Pauli Hamiltonian. Builds H from
    /// raw <see cref="PauliTerm"/> entries (supports single-site h_x·X_l, 3-body
    /// X⊗X⊗X, etc., that the bilinear bond overload cannot express), constructs L
    /// via <see cref="PauliDephasingDissipator.BuildZ"/> or
    /// <see cref="T1Dissipator.Build"/>, then runs the same Π-decomposition +
    /// ±i refinement pipeline.</summary>
    /// <param name="chain">N-qubit chain providing N, γ₀ (uniform Z-dephasing rate).</param>
    /// <param name="kBodyTerms">Pauli-string terms with arbitrary k-body weight.
    /// Each term must have <see cref="PauliTerm.N"/> = <see cref="ChainSystem.N"/>.
    /// Coefficients are <see cref="Complex"/>; for Hermitian H pass real coefficients
    /// on Hermitian Pauli strings (the canonical case for F112 Tier1Derived scope).</param>
    /// <param name="gammaT1">Optional uniform per-site σ⁻ amplitude-damping rate.
    /// When non-null and non-zero, routes through <see cref="T1Dissipator"/>; otherwise
    /// pure Z-dephasing via <see cref="PauliDephasingDissipator"/>. F112 typed
    /// Tier1Derived scope covers the gammaT1=null branch (single-Pauli dephasing is
    /// bit_b-homogeneous); the T1 branch is bit_b-mixed and out of F112 typed scope.</param>
    /// <param name="dephaseLetter">Dephasing letter for the Π operator; default Z.</param>
    /// <remarks>No F81 inner check on this path: it takes H and c directly, so there is no
    /// term list to project H_odd from, and <see cref="PolarityCoordinatesResult.F81Violation"/>
    /// is reported as 0 (sentinel, no closed-form prediction to compare against). The identity
    /// itself is not 2-body-only; it holds at k = 3 and k = 4 (Python
    /// test_F85_kbody_F81_identity_at_k3).</remarks>
    public static PolarityCoordinatesResult Decompose(
        ChainSystem chain,
        IReadOnlyList<PauliTerm> kBodyTerms,
        double? gammaT1 = null,
        PauliLetter dephaseLetter = PauliLetter.Z)
    {
        if (chain is null) throw new ArgumentNullException(nameof(chain));
        if (kBodyTerms is null) throw new ArgumentNullException(nameof(kBodyTerms));
        foreach (var term in kBodyTerms)
        {
            if (term.N != chain.N)
                throw new ArgumentException(
                    $"term has N={term.N} letters; expected N={chain.N}", nameof(kBodyTerms));
        }

        var H = new PauliHamiltonian(chain.N, kBodyTerms).ToMatrix();
        var gammaZ = Enumerable.Repeat(chain.GammaZero, chain.N).ToArray();

        ComplexMatrix L;
        if (gammaT1.HasValue && gammaT1.Value != 0.0)
        {
            var gammaT1PerSite = Enumerable.Repeat(gammaT1.Value, chain.N).ToArray();
            L = T1Dissipator.Build(H, gammaZ, gammaT1PerSite);
        }
        else
        {
            L = PauliDephasingDissipator.BuildZ(H, gammaZ);
        }

        var M = PalindromeResidual.Build(L, chain.N, chain.SigmaGamma, dephaseLetter);
        var piOp = PiOperator.BuildFull(chain.N, dephaseLetter);
        var piInv = piOp.ConjugateTranspose();
        var piMpi = piOp * M * piInv;
        var mSym = (M + piMpi) / 2.0;
        var mAnti = (M - piMpi) / 2.0;

        // No F81 inner check on this path (k-body / single-site H_odd identity is open), so
        // LHOddNormSquared stays null rather than a zero that would read as "no polarity
        // content". The STRUCTURAL predicate is available here regardless: it reads the Pauli
        // letters of kBodyTerms and never needs the H_odd projection.
        return RefineMAnti(M, mSym, mAnti, piOp, f81Violation: 0.0, lHOddNormSquared: null,
            isStructurallySilent: IsStructurallyDegenerate(
                kBodyTerms, gammaT1.HasValue && gammaT1.Value != 0.0, dephaseLetter));
    }

    /// <summary>Shared core: given M, M_sym, M_anti and Π, compute the ±i
    /// Π-eigenvalue projectors on M_anti and assemble the result record.</summary>
    private static PolarityCoordinatesResult RefineMAnti(
        ComplexMatrix M, ComplexMatrix mSym, ComplexMatrix mAnti, ComplexMatrix piOp, double f81Violation,
        double? lHOddNormSquared, bool? isStructurallySilent = null)
    {
        var piInv = piOp.ConjugateTranspose(); // Π is a unitary signed permutation; Π⁻¹ = Π†.
        var piMAntiPiInv = piOp * mAnti * piInv;
        var mPlusHalf = (mAnti - Complex.ImaginaryOne * piMAntiPiInv) / 2.0;
        var mMinusHalf = (mAnti + Complex.ImaginaryOne * piMAntiPiInv) / 2.0;
        var mZero = mSym; // F81 M_sym is the 0-axis component by definition.

        double normM = Math.Pow(M.FrobeniusNorm(), 2);
        double normZero = Math.Pow(mZero.FrobeniusNorm(), 2);
        double normPlus = Math.Pow(mPlusHalf.FrobeniusNorm(), 2);
        double normMinus = Math.Pow(mMinusHalf.FrobeniusNorm(), 2);

        double orthogonalityResidual = Math.Abs(normM - (normZero + normPlus + normMinus));
        double asymmetry = normPlus - normMinus;

        return new PolarityCoordinatesResult(
            M: M,
            MZero: mZero,
            MPlusHalf: mPlusHalf,
            MMinusHalf: mMinusHalf,
            MNormSquared: normM,
            MZeroNormSquared: normZero,
            MPlusHalfNormSquared: normPlus,
            MMinusHalfNormSquared: normMinus,
            Asymmetry: asymmetry,
            OrthogonalityResidual: orthogonalityResidual,
            F81Violation: f81Violation,
            LHOddNormSquared: lHOddNormSquared,
            IsStructurallySilent: isStructurallySilent);
    }
}
