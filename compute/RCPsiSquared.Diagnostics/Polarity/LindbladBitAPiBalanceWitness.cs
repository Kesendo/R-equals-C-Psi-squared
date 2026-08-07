using System.Numerics;
using RCPsiSquared.Core.ChainSystems;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Pauli;

namespace RCPsiSquared.Diagnostics.Polarity;

/// <summary>An F112-X polarity-balance witness: a named (chain + H + c_ops + γ_T1) triple
/// paired with the expected verdict (<c>"BALANCED"</c> or <c>"BROKEN"</c>) that the
/// computed F112-X polarity asymmetry should produce under <b>X-dephase Π_X polarity</b>
/// decomposition. Self-computes the actual asymmetry lazily via
/// <see cref="PolarityCoordinates.Decompose(ChainSystem, IReadOnlyList{PauliTerm}, double?, PauliLetter)"/>
/// with <c>dephaseLetter = PauliLetter.X</c>, so listing many witnesses does not trigger
/// N L-builds upfront.
///
/// <para>Mirrors <see cref="LindbladBitBPiBalanceWitness"/> exactly except for the
/// X-dephase routing. The 5-witness <see cref="StandardSet"/> is bit_a-axis adapted from
/// the BitB witness set: every term-list is independently checked for bit_a-homogeneity
/// to land inside the F112-X Tier1Derived scope (per
/// <c>docs/proofs/PROOF_F112_CROSS_DEPHASE_VIA_KLEIN_V4.md</c> Theorem statement).</para>
///
/// <para><b>F112-X typed scope (Tier1Derived)</b>: for any Lindblad-form Liouvillian
/// with any H (Hermitian or non-Hermitian) and each c_k bit_a-homogeneous, the asymmetry
/// is exactly 0 bit-exact. See <see cref="LindbladBitAPiBalance"/>.</para>
///
/// <para><b>Why all 5 witnesses are BALANCED (no BROKEN counterexample in StandardSet)</b>:
/// the structural F112-X breakage requires bit_a-MIXED c. Under X-dephase, the dephasing
/// dissipator is c = X per site (bit_a = 1), and the optional T1 σ⁻ has bit_a = 1 as
/// well (σ⁻ = (X − i Y) / 2; both X and Y have bit_a = 1). Combining X-dephase + T1
/// keeps c bit_a-homogeneous, so the F112-Z's f95-style "Z-drive + T1" BROKEN
/// counterexample has NO analog reachable through the current PolarityCoordinates
/// plumbing on the X-dephase axis (which only exposes the (dephase-letter, optional
/// uniform T1) c-construction knobs). The 5-witness StandardSet substantively spans
/// F112-X's in-scope half across Heisenberg / Π²_X-even / Π²_X-odd bilinears + single-
/// site X-drive ± T1 envelopes; each is expected BALANCED bit-exact.</para>
///
/// <para>Anchor: <c>docs/ANALYTICAL_FORMULAS.md</c> F112 +
/// <c>docs/proofs/PROOF_F112_CROSS_DEPHASE_VIA_KLEIN_V4.md</c> +
/// <see cref="LindbladBitAPiBalance"/> + <see cref="PolarityCoordinates"/>.</para>
/// </summary>
public sealed class LindbladBitAPiBalanceWitness : Claim
{
    private readonly Lazy<PolarityCoordinatesResult> _polarity;

    /// <summary>Human-readable witness name.</summary>
    public string WitnessName { get; }

    /// <summary>Chain providing N, γ₀, bond layout.</summary>
    public ChainSystem Chain { get; }

    /// <summary>k-body Pauli term list representing the Hamiltonian.</summary>
    public IReadOnlyList<PauliTerm> Terms { get; }

    /// <summary>Optional uniform per-site σ⁻ T1 rate. Null or zero means pure X-dephasing.</summary>
    public double? GammaT1 { get; }

    /// <summary>Expected verdict: <c>"BALANCED"</c> or <c>"BROKEN"</c>.</summary>
    public string ExpectedVerdict { get; }

    /// <summary>True iff <see cref="ExpectedVerdict"/> = <c>"BALANCED"</c>.</summary>
    public bool ExpectedBalanced => ExpectedVerdict == "BALANCED";

    /// <summary>Relative-asymmetry threshold below which the witness reads as
    /// <c>"BALANCED"</c>, applied to
    /// <see cref="PolarityCoordinatesResult.RelativeAsymmetry"/>, the contrast ratio
    /// |asymmetry| / ‖M_anti‖². Default 1e-10 is not a bit-exact scale and must not be called one: it is a
    /// SEPARATION, and on THIS axis it currently separates nothing: <see cref="StandardSet"/>
    /// has no BROKEN witness at all (three DEGENERATE, two BALANCED), the two BALANCED ones
    /// measure rel asym exactly 0.0 and the three DEGENERATE ones measure no ratio and read
    /// NaN. An earlier version of this paragraph quoted "the counterexample measures
    /// 7.692080e-3"; that number is the Pi_Z axis's, transplanted from
    /// <see cref="LindbladBitBPiBalanceWitness"/>, and F112-X has no analogue because the
    /// F113-style break does not trigger under X-dephase (witness 5 below is the sighting).
    /// The tolerance is therefore a guard against a future break, not a measured separation. The theorem it serves (asymmetry exactly 0 in F112-X's typed scope) IS exact,
    /// but a threshold on a float difference of two Frobenius norms cannot witness that;
    /// docs/GLOSSARY.md:304 is the governing rule. Where the witness is polarity-degenerate
    /// (<see cref="IsDegenerate"/>) the threshold decides nothing and is not consulted: the
    /// ratio is 0/0 and <see cref="ActualRelativeAsymmetry"/> is NaN.</summary>
    public double Tolerance { get; }

    /// <summary>Lazily computed polarity decomposition result. First access triggers the
    /// L-build and Π_X-decomposition via <see cref="PolarityCoordinates.Decompose"/>.</summary>
    public Lazy<PolarityCoordinatesResult> Polarity => _polarity;

    /// <summary>Actual relative asymmetry, the contrast ratio
    /// <see cref="PolarityCoordinatesResult.RelativeAsymmetry"/>. Forces lazy evaluation of
    /// <see cref="Polarity"/>.</summary>
    public double ActualRelativeAsymmetry => _polarity.Value.RelativeAsymmetry;

    /// <summary>True when this witness has no polarity content to balance, so the asymmetry
    /// is 0 − 0 and <see cref="ActualVerdict"/> used to read <c>"BALANCED"</c> vacuously;
    /// since 2026-08-07 it reads <c>"DEGENERATE"</c>. This is the DECIDER, and since
    /// 2026-08-07 it is no longer a float-only reading:
    /// <see cref="PolarityCoordinatesResult.IsPolarityDegenerate"/> now takes the exact
    /// structural verdict carried in from the term list, with the float ‖M_anti‖² == 0.0 test
    /// beneath it in the disjunction. Forces lazy evaluation.</summary>
    public bool IsDegenerate => _polarity.Value.IsPolarityDegenerate;

    /// <para><b>DEGENERATE, added 2026-08-07.</b> A witness whose H is Π²-even on this axis
    /// AND whose collapse operators are all Π²-homogeneous on it has NO polarity content:
    /// both halves of M_anti vanish as a theorem, the asymmetry is 0 − 0, and the ratio is a
    /// quotient of two zeros. It is not a balance and must not be counted as one. Decided by
    /// <see cref="IsStructurallyDegenerate"/> from the Pauli letters alone and BEFORE the
    /// Liouvillian is built, with <see cref="IsDegenerate"/> as a second arm beneath it for
    /// the structural test's known false negatives (cancelling +c/−c pairs on the c side,
    /// reachable through the k-body constructor but not through the bond one) where the
    /// float ‖M_anti‖² lands on exact 0.0; routing those through the threshold instead would
    /// compare NaN and report BROKEN.</para>
    /// <summary>Actual verdict: <c>"BALANCED"</c> if
    /// <see cref="ActualRelativeAsymmetry"/> ≤ <see cref="Tolerance"/>, else
    /// <c>"BROKEN"</c>. Forces lazy evaluation.</summary>
    public string ActualVerdict =>
        // The EXACT test first, and it must stay first: it reads Pauli letters, so a silent
        // witness reaches its verdict without ever assembling the 4^N Liouvillian, and
        // Degenerate_Witness_Reaches_Its_Verdict_Without_Building_L pins that. `IsDegenerate`
        // beneath it forces the decomposition. Read it as the FLOAT arm, not as a second
        // structural one: IsPolarityDegenerate is (IsStructurallySilent ?? false) || ‖M_anti‖²
        // == 0.0, and on every witness its first disjunct is computed from the same inputs as
        // arm 1 above, so reaching here at all means arm 1 said false and only the float test
        // can still fire. It exists for the structural test's false negatives (cancelling
        // +c/-c pairs on the c side); without it those rows would reach the threshold, compare
        // NaN, and report BROKEN off a 0/0. NOTE it is dead in the current suite and reachable
        // only through the k-body constructor; J = 0 does NOT reach it, because the bond
        // constructor expands to PauliTerms at coefficient J and the PauliTerm overload skips
        // zero-coefficient terms, so the FIRST arm fires there.
        IsStructurallyDegenerate ? "DEGENERATE"
            : IsDegenerate ? "DEGENERATE"
                : ActualRelativeAsymmetry <= Tolerance ? "BALANCED" : "BROKEN";

    /// <summary>EXACT structural verdict, decided from Pauli letters alone: this
    /// configuration has no polarity content, so it is SILENT rather than balanced.
    /// Outranks <see cref="ActualRelativeAsymmetry"/> in <see cref="ActualVerdict"/>,
    /// because on the Π_X axis the ratio is then a quotient of two zeros and reports
    /// nothing in either direction. Unlike <see cref="IsDegenerate"/>, which reads the
    /// float M_anti and only reaches exact 0.0 at N = 2, this is N-independent.
    /// See <see cref="PolarityCoordinates.IsStructurallyDegenerate"/>.</summary>
    public bool IsStructurallyDegenerate =>
        PolarityCoordinates.IsStructurallyDegenerate(
            Terms, GammaT1.HasValue && GammaT1.Value != 0.0, PauliLetter.X);

    /// <summary>True iff <see cref="ActualVerdict"/> matches
    /// <see cref="ExpectedVerdict"/> bit-exactly.</summary>
    public bool Matches => ActualVerdict == ExpectedVerdict;

    /// <summary>Primary k-body constructor. Threads <c>dephaseLetter = PauliLetter.X</c>
    /// into <see cref="PolarityCoordinates.Decompose(ChainSystem, IReadOnlyList{PauliTerm}, double?, PauliLetter)"/>.</summary>
    public LindbladBitAPiBalanceWitness(
        string witnessName,
        ChainSystem chain,
        IReadOnlyList<PauliTerm> terms,
        double? gammaT1,
        string expectedVerdict,
        double tolerance = 1e-10)
        : base($"F112-X polarity-balance witness: {witnessName}",
               Tier.Tier1Derived,
               "docs/ANALYTICAL_FORMULAS.md F112 + " +
               "docs/proofs/PROOF_F112_CROSS_DEPHASE_VIA_KLEIN_V4.md + " +
               "compute/RCPsiSquared.Core/Symmetry/LindbladBitAPiBalance.cs + " +
               "compute/RCPsiSquared.Diagnostics/Polarity/PolarityCoordinates.cs + " +
               "simulations/f112_klein_v4_cross_dephase_verify.py")
    {
        if (string.IsNullOrWhiteSpace(witnessName))
            throw new ArgumentException("witnessName must be non-empty", nameof(witnessName));
        if (expectedVerdict != "BALANCED" && expectedVerdict != "BROKEN"
            && expectedVerdict != "DEGENERATE")
            throw new ArgumentException(
                $"expectedVerdict must be 'BALANCED', 'BROKEN' or 'DEGENERATE'; " +
                $"got '{expectedVerdict}'",
                nameof(expectedVerdict));
        if (tolerance < 0)
            throw new ArgumentOutOfRangeException(nameof(tolerance),
                $"tolerance must be ≥ 0; got {tolerance}");

        WitnessName = witnessName;
        Chain = chain ?? throw new ArgumentNullException(nameof(chain));
        Terms = terms ?? throw new ArgumentNullException(nameof(terms));
        GammaT1 = gammaT1;
        ExpectedVerdict = expectedVerdict;
        Tolerance = tolerance;

        _polarity = new Lazy<PolarityCoordinatesResult>(
            () => PolarityCoordinates.Decompose(chain, terms, gammaT1, dephaseLetter: PauliLetter.X));
    }

    /// <summary>Convenience ctor for bilinear-bond witnesses: expands
    /// <see cref="PauliPairBondTerm"/> entries across <see cref="ChainSystem.Bonds"/>
    /// into a k-body <see cref="PauliTerm"/> list using the chain's J coupling.</summary>
    public LindbladBitAPiBalanceWitness(
        string witnessName,
        ChainSystem chain,
        IReadOnlyList<PauliPairBondTerm> bondTerms,
        double? gammaT1,
        string expectedVerdict,
        double tolerance = 1e-10)
        : this(witnessName, chain,
               ExpandBondTermsToKBody(chain, bondTerms),
               gammaT1, expectedVerdict, tolerance)
    {
    }

    private static IReadOnlyList<PauliTerm> ExpandBondTermsToKBody(
        ChainSystem chain, IReadOnlyList<PauliPairBondTerm> bondTerms)
    {
        if (chain is null) throw new ArgumentNullException(nameof(chain));
        if (bondTerms is null) throw new ArgumentNullException(nameof(bondTerms));
        var result = new List<PauliTerm>(chain.Bonds.Count * bondTerms.Count);
        foreach (var bond in chain.Bonds)
        {
            foreach (var t in bondTerms)
            {
                result.Add(PauliTerm.TwoSite(
                    chain.N, bond.Site1, t.LetterA, bond.Site2, t.LetterB,
                    coefficient: (Complex)chain.J));
            }
        }
        return result;
    }

    public override string DisplayName =>
        $"F112-X witness '{WitnessName}': expected {ExpectedVerdict}";

    public override string Summary =>
        $"chain N={Chain.N}, {Terms.Count} term(s), γ_T1={FormatGamma(GammaT1)}, dephase=X; " +
        $"expected {ExpectedVerdict}, actual rel asym = {ActualRelativeAsymmetry:E3} → " +
        $"{ActualVerdict} ({(Matches ? "PASS" : "FAIL")})";

    private static string FormatGamma(double? g) =>
        g is null ? "none" : g.Value.ToString("G6");

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("witness name", summary: WitnessName);
            yield return new InspectableNode("chain",
                summary: $"N={Chain.N}, {Chain.HType}, {Chain.Topology}, γ₀={Chain.GammaZero:G6}");
            yield return new InspectableNode("dephase letter", summary: "X (F112-X scope)");
            yield return new InspectableNode("expected verdict", summary: ExpectedVerdict);
            yield return new InspectableNode("actual verdict", summary: ActualVerdict);
            yield return new InspectableNode("agreement", summary: Matches ? "PASS" : "FAIL");
            yield return InspectableNode.RealScalar(
                "relative asymmetry", ActualRelativeAsymmetry, format: "E6");
            yield return InspectableNode.RealScalar(
                "polarity content ‖M_anti‖²", Polarity.Value.MAntiNormSquared, format: "E6");
            yield return new InspectableNode("polarity-degenerate",
                summary: IsStructurallyDegenerate
                    ? "yes, STRUCTURALLY: H is Π²-even on this axis and every c is "
                      + "Π²-homogeneous on it, so both halves of M_anti vanish as a "
                      + "theorem. The verdict is DEGENERATE, not BALANCED, and the tolerance "
                      + $"decides nothing. The float reading agrees: {IsDegenerate}."
                    : "no: this configuration carries polarity content, so the relative "
                      + $"asymmetry is a real measurement. Float ‖M_anti‖² == 0.0: {IsDegenerate} "
                      + "(that float test alone reaches exact 0.0 only at N = 2, which is why "
                      + "the structural one decides).");
            yield return InspectableNode.RealScalar("tolerance", Tolerance, format: "E2");
            yield return new InspectableNode("γ_T1", summary: FormatGamma(GammaT1));
            yield return new InspectableNode("terms",
                summary: $"{Terms.Count} term(s); labels = [{string.Join(", ", Terms.Select(t => t.Label))}]");
            yield return InspectableNode.RealScalar(
                "F81 violation (numerical, X-deph; no closed-form prediction on cross-dephase path)",
                Polarity.Value.F81Violation, format: "E3");
        }
    }

    /// <summary>Five named witnesses spanning the F112-X BALANCED axis under X-dephase:
    /// <list type="number">
    ///   <item><b>Heisenberg_pure_X_DEGENERATE</b>: H = XX+YY+ZZ Heisenberg + pure X-dephasing.
    ///         All three bilinears are bit_a-homogeneous (XX:0, YY:0, ZZ:0 → all even).
    ///         Inside F112-X typed scope, and SILENT within it: bit_a-homogeneous H against a
    ///         bit_a-homogeneous c leaves no polarity content. Expected DEGENERATE (this line
    ///         said BALANCED while the constructor already said DEGENERATE).</item>
    ///   <item><b>YZ_ZY_bit_a_odd_balanced</b>: H = YZ+ZY. Both bilinears are bit_a-odd
    ///         (YZ: 1+0=1; ZY: 0+1=1). Bit_a-homogeneous (both odd) and inside F112-X scope.
    ///         Expected BALANCED.</item>
    ///   <item><b>XY_bit_a_even_DEGENERATE</b>: H = XY. bit_a-even (1+1=0); the canonical
    ///         Π²_X-even bilinear handled by F108 Part 2 + bit_a-homogeneous. Inside F112-X scope
    ///         and silent within it. Expected DEGENERATE.</item>
    ///   <item><b>Heisenberg_with_T1_envelope_DEGENERATE</b>: H = Heisenberg + σ⁻ T1 (γ_T1=0.1).
    ///         Under X-dephase the dissipator c includes c = X (bit_a=1) per site plus σ⁻
    ///         (bit_a=1; σ⁻ = (X − i Y)/2, both X and Y have bit_a=1). c is bit_a-homogeneous;
    ///         F112-X covers this case Tier1Derived, and it is silent within it. Expected
    ///         DEGENERATE.</item>
    ///   <item><b>Xdrive_with_T1_envelope_balanced</b>: H = ω·(X₀+X₁)/2 single-site X-drive
    ///         (ω = 0.13) + σ⁻ T1 (γ_T1 = 0.001). The X-axis analog of F112-Z's f95
    ///         counterexample setup; under X-dephase, X-drive aligns with the dephase axis
    ///         (both bit_a=1) and σ⁻ is also bit_a=1, so c is bit_a-homogeneous. F112-X holds
    ///         and the asymmetry is bit-exact 0. Demonstrates that the F113-style break
    ///         mechanism does NOT trigger on the bit_a axis under X-deph (the analog of
    ///         the BitB Zdrive+T1 BROKEN witness is BALANCED here). Expected BALANCED.</item>
    /// </list>
    /// All five use <see cref="Tolerance"/> = 1e-10. Three of them never reach it (DEGENERATE
    /// outranks the ratio) and the other two measure exactly 0.0, so on this axis the number
    /// is an unexercised guard rather than a cut between measured values; see
    /// <see cref="Tolerance"/>. It is NOT a bit-exact scale, and an earlier version of this
    /// sentence called it one.</summary>
    public static IReadOnlyList<LindbladBitAPiBalanceWitness> StandardSet(ChainSystem chain)
    {
        if (chain is null) throw new ArgumentNullException(nameof(chain));
        if (chain.N != 2)
            throw new ArgumentException(
                $"StandardSet requires N=2 chain (witness 5 is a fixed-ω two-site X-drive); got N={chain.N}",
                nameof(chain));

        var heisenbergTerms = new[]
        {
            new PauliPairBondTerm(PauliLetter.X, PauliLetter.X),
            new PauliPairBondTerm(PauliLetter.Y, PauliLetter.Y),
            new PauliPairBondTerm(PauliLetter.Z, PauliLetter.Z),
        };
        var yzZyTerms = new[]
        {
            new PauliPairBondTerm(PauliLetter.Y, PauliLetter.Z),
            new PauliPairBondTerm(PauliLetter.Z, PauliLetter.Y),
        };
        var xyTerms = new[]
        {
            new PauliPairBondTerm(PauliLetter.X, PauliLetter.Y),
        };

        // Witness 5 ingredients: single-site X-drive at ω/2 on each site (bit_a-axis analog
        // of the BitB Zdrive+T1 witness; under X-deph everything is bit_a=1 homogeneous and
        // F112-X holds, in contrast to F112-Z where Zdrive+T1 BREAKS).
        const double omega = 0.13;
        var xDriveTerms = new List<PauliTerm>
        {
            PauliTerm.SingleSite(chain.N, 0, PauliLetter.X, coefficient: (Complex)(omega / 2.0)),
            PauliTerm.SingleSite(chain.N, 1, PauliLetter.X, coefficient: (Complex)(omega / 2.0)),
        };

        return new[]
        {
            new LindbladBitAPiBalanceWitness(
                witnessName: "Heisenberg_pure_X_DEGENERATE",
                chain: chain,
                bondTerms: heisenbergTerms,
                gammaT1: null,
                expectedVerdict: "DEGENERATE"),
            new LindbladBitAPiBalanceWitness(
                witnessName: "YZ_ZY_bit_a_odd_balanced",
                chain: chain,
                bondTerms: yzZyTerms,
                gammaT1: null,
                expectedVerdict: "BALANCED"),
            new LindbladBitAPiBalanceWitness(
                witnessName: "XY_bit_a_even_DEGENERATE",
                chain: chain,
                bondTerms: xyTerms,
                gammaT1: null,
                expectedVerdict: "DEGENERATE"),
            new LindbladBitAPiBalanceWitness(
                witnessName: "Heisenberg_with_T1_envelope_DEGENERATE",
                chain: chain,
                bondTerms: heisenbergTerms,
                gammaT1: 0.1,
                expectedVerdict: "DEGENERATE"),
            new LindbladBitAPiBalanceWitness(
                witnessName: "Xdrive_with_T1_envelope_balanced",
                chain: chain,
                terms: xDriveTerms,
                gammaT1: 0.001,
                expectedVerdict: "BALANCED"),
        };
    }
}
