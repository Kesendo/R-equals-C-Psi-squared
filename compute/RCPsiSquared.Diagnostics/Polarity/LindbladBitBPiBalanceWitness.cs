using System.Numerics;
using RCPsiSquared.Core.ChainSystems;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Pauli;

namespace RCPsiSquared.Diagnostics.Polarity;

/// <summary>An F112 polarity-balance witness: a named (chain + H + c_ops + γ_T1) triple
/// paired with the expected verdict (<c>"BALANCED"</c> or <c>"BROKEN"</c>) that the
/// computed F112 polarity asymmetry should produce. Self-computes the actual asymmetry
/// lazily via <see cref="PolarityCoordinates.Decompose"/>, so listing many witnesses
/// does not trigger N L-builds upfront.
///
/// <para>Mirrors the <see cref="F87.F87CanonicalWitness"/> shape: expected / actual
/// comparison with a <see cref="Matches"/> boolean. Disagreement signals either a chain
/// configuration that breaks the witness (e.g. wrong N, mismatched dephase letter, missing
/// T1) or a real regression in the <see cref="PolarityCoordinates"/> primitive.</para>
///
/// <para><b>F112 typed scope (Tier1Derived)</b>: for any Lindblad-form Liouvillian with
/// Hermitian H and each c_k bit_b-homogeneous, the asymmetry is exactly 0 bit-exact. See
/// <see cref="LindbladBitBPiBalance"/>. Witnesses with <c>ExpectedVerdict = "BALANCED"</c>
/// either sit strictly inside that scope (witnesses 1-3) or sit in the broader
/// <i>empirical envelope</i> observed across probes 1-14 (witness 4: bit_b-mixed σ⁻ T1 with
/// idle / XY-bond / Heisenberg-style Hermitian H still balances bit-exact). Witnesses with
/// <c>ExpectedVerdict = "BROKEN"</c> are structural counterexamples that fall outside both
/// (witness 5: Z-drive Hermitian H combined with σ⁻ T1, the configuration discovered
/// Welle 2 in <c>simulations/f112_hardware_lens_multi.py</c>).</para>
///
/// <para>Anchor: <c>docs/ANALYTICAL_FORMULAS.md</c> F112 +
/// <c>docs/proofs/PROOF_F112_LINDBLAD_BIT_B_PI_BALANCE.md</c> +
/// <see cref="LindbladBitBPiBalance"/> +
/// <see cref="PolarityCoordinates"/>.</para>
/// </summary>
public sealed class LindbladBitBPiBalanceWitness : Claim
{
    private readonly Lazy<PolarityCoordinatesResult> _polarity;

    /// <summary>Human-readable witness name.</summary>
    public string WitnessName { get; }

    /// <summary>Chain providing N, γ₀, bond layout.</summary>
    public ChainSystem Chain { get; }

    /// <summary>k-body Pauli term list representing the Hamiltonian (already expanded
    /// across bonds for the bilinear witnesses).</summary>
    public IReadOnlyList<PauliTerm> Terms { get; }

    /// <summary>Optional uniform per-site σ⁻ T1 rate. Null or zero means pure Z-dephasing.</summary>
    public double? GammaT1 { get; }

    /// <summary>Expected verdict: <c>"BALANCED"</c> or <c>"BROKEN"</c>.</summary>
    public string ExpectedVerdict { get; }

    /// <summary>True iff <see cref="ExpectedVerdict"/> = <c>"BALANCED"</c>.</summary>
    public bool ExpectedBalanced => ExpectedVerdict == "BALANCED";

    /// <summary>Relative-asymmetry threshold below which the witness reads as
    /// <c>"BALANCED"</c>, applied to
    /// <see cref="PolarityCoordinatesResult.RelativeAsymmetry"/>, the contrast ratio
    /// |asymmetry| / ‖M_anti‖². Default 1e-10 is not a bit-exact scale and must not be called one: it is a
    /// SEPARATION. The in-scope witnesses that carry polarity content measure rel asym
    /// exactly 0.0 (the DEGENERATE ones measure no ratio at all, and read NaN) while the
    /// counterexample measures 7.692080e-3, so 1e-10 sits about seven orders below the
    /// smallest break and eight above the float noise; anything above it is a break, not
    /// rounding. The theorem it serves (asymmetry exactly 0 in F112's typed scope) IS exact,
    /// but a threshold on a float difference of two Frobenius norms cannot witness that;
    /// docs/GLOSSARY.md:304 is the governing rule.
    ///
    /// <para>Where the witness is polarity-degenerate
    /// (<see cref="PolarityCoordinatesResult.IsPolarityDegenerate"/>) the threshold decides
    /// nothing and is not consulted: the ratio is 0/0 and
    /// <see cref="ActualRelativeAsymmetry"/> is NaN. Witnesses 1 and 2 below are in that case.
    /// Read <see cref="ActualVerdict"/> together with <see cref="IsDegenerate"/>.</para></summary>
    public double Tolerance { get; }

    /// <summary>Lazily computed polarity decomposition result. First access triggers the
    /// L-build and Π-decomposition via <see cref="PolarityCoordinates.Decompose"/>.</summary>
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
    /// because on the Π_Z axis the ratio is then a quotient of two zeros and reports
    /// nothing in either direction. Unlike <see cref="IsDegenerate"/>, which reads the
    /// float M_anti and only reaches exact 0.0 at N = 2, this is N-independent.
    /// See <see cref="PolarityCoordinates.IsStructurallyDegenerate"/>.</summary>
    public bool IsStructurallyDegenerate =>
        PolarityCoordinates.IsStructurallyDegenerate(
            Terms, GammaT1.HasValue && GammaT1.Value != 0.0, PauliLetter.Z);

    /// <summary>True iff <see cref="ActualVerdict"/> matches
    /// <see cref="ExpectedVerdict"/> bit-exactly.</summary>
    public bool Matches => ActualVerdict == ExpectedVerdict;

    /// <summary>Primary k-body constructor. Takes the H term list as raw
    /// <see cref="PauliTerm"/> entries (handles single-site, bilinear, and k-body terms
    /// uniformly via <see cref="PolarityCoordinates.Decompose(ChainSystem, IReadOnlyList{PauliTerm}, double?, PauliLetter)"/>).</summary>
    public LindbladBitBPiBalanceWitness(
        string witnessName,
        ChainSystem chain,
        IReadOnlyList<PauliTerm> terms,
        double? gammaT1,
        string expectedVerdict,
        double tolerance = 1e-10)
        : base($"F112 polarity-balance witness: {witnessName}",
               Tier.Tier1Derived,
               "docs/ANALYTICAL_FORMULAS.md F112 + " +
               "docs/proofs/PROOF_F112_LINDBLAD_BIT_B_PI_BALANCE.md + " +
               "compute/RCPsiSquared.Core/Symmetry/LindbladBitBPiBalance.cs + " +
               "compute/RCPsiSquared.Diagnostics/Polarity/PolarityCoordinates.cs + " +
               "experiments/F112_HARDWARE_LENS_KINGSTON.md + " +
               "simulations/f112_hardware_lens_multi.py")
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
            () => PolarityCoordinates.Decompose(chain, terms, gammaT1));
    }

    /// <summary>Convenience ctor for bilinear-bond witnesses: takes
    /// <see cref="PauliPairBondTerm"/> entries and expands them across
    /// <see cref="ChainSystem.Bonds"/> into a k-body <see cref="PauliTerm"/> list using
    /// the chain's J coupling. Equivalent to running the bilinear overload of
    /// <see cref="PolarityCoordinates.Decompose(ChainSystem, IReadOnlyList{PauliPairBondTerm}, double?, PauliLetter)"/>
    /// internally; routes through the k-body path so a single witness type covers both
    /// shapes.</summary>
    public LindbladBitBPiBalanceWitness(
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
        $"F112 witness '{WitnessName}': expected {ExpectedVerdict}";

    public override string Summary =>
        $"chain N={Chain.N}, {Terms.Count} term(s), γ_T1={FormatGamma(GammaT1)}; " +
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
                "F81 violation (inner PiDecomposition sanity)",
                Polarity.Value.F81Violation, format: "E3");
        }
    }

    /// <summary>The five named witnesses spanning the F112 BALANCED / BROKEN axes:
    /// <list type="number">
    ///   <item><b>Heisenberg_pure_Z_DEGENERATE</b>: H = XX+YY+ZZ Heisenberg + pure Z-dephasing.
    ///         Inside F112 typed scope (Hermitian Heisenberg, single-Pauli Z dephase). Expected BALANCED,
    ///         but <b>polarity-degenerate</b>: Heisenberg is <i>truly</i>, so M itself vanishes by the
    ///         F83 closed form and both halves are the zero matrix. The verdict is exactly 0 = 0,
    ///         carrying no information about the threshold.</item>
    ///   <item><b>YZ_ZY_pi2even_DEGENERATE</b>: H = YZ+ZY (F108 non-truly Π²-even).
    ///         Inside F112 typed scope. Expected BALANCED, and <b>also polarity-degenerate</b>, for a
    ///         different reason worth keeping apart: ‖M‖² is large here (256 at N=2, J=1) but Π²-even H
    ///         has no Π²-odd part, so M_anti = L_{H_odd} = 0 by F81 and the halves vanish anyway. This
    ///         is the witness that shows ‖M‖² is the wrong scale: it is nonzero while the object the
    ///         asymmetry measures is empty.</item>
    ///   <item><b>XY_pi2odd_balanced</b>: H = XY (Π²-odd). Inside F112 typed scope.
    ///         Expected BALANCED.</item>
    ///   <item><b>Heisenberg_with_T1_envelope_balanced</b>: H = Heisenberg + σ⁻ T1 (γ_T1=0.1).
    ///         Outside F112 typed scope (σ⁻ is bit_b-mixed). Broader empirical envelope
    ///         observed across probes 1-14 (and re-confirmed Welle 2 across cusp_slowing,
    ///         block_cpsi, chain_gamma0 Tier-A datasets): idle / XY-bond / Heisenberg-style
    ///         Hermitian H combined with σ⁻ T1 still balances bit-exact. Expected BALANCED.</item>
    ///   <item><b>Zdrive_with_T1_envelope_BROKEN</b>: H = ω·(Z₀+Z₁)/2 single-site Z-drive
    ///         (ω = 0.13) + σ⁻ T1 (γ_T1 = 0.001). The structural counterexample, first
    ///         noticed while fitting the f95_angle_steering Tier-A dataset (Kingston,
    ///         2026-05-16), whose protocol supplies a drive of this kind, though on ONE
    ///         qubit of the pair rather than both, so its own value is half of this one. The
    ///         quantity F113 predicts is the <i>absolute</i> asymmetry, not the ratio:
    ///         asymmetry = (4^N/2)·Σ_l ω_l·(γ_pump,l − γ_T1,l), reproduced by
    ///         <c>simulations/f113_break_formula_derivation.py</c>. <b>Mind the factor 2 in
    ///         ω_l</b>: it is the Larmor rate, and the two terms carry ω/2 each as their Pauli
    ///         coefficient, so ω_l = ω on each site and Σ_l ω_l = 2ω = 0.26, not 0.13. That
    ///         gives (16/2)·0.26·(−0.001) = −2.08e-3, met to 2.7e-14 relative
    ///         (−2.079999999999943e-3). Dividing by the polarity content
    ///         ‖M_anti‖² = 2.704080e-1 gives rel asym = 7.692080e-3. (Under the retired ‖M‖²
    ///         denominator the same witness read 3.845528e-3, a factor 2.000266 lower
    ///         (‖M‖² = 5.408880e-1 against ‖M_anti‖² = 2.704080e-1). NOT exactly half, and the
    ///         two quoted values refute that on their own: 3.845528e-3 × 2 = 7.691056e-3. Both are far above any threshold, so the verdict is
    ///         unchanged; this one is a real contrast rather than a ratio against a part of M
    ///         the asymmetry never touches.) Not degenerate: a single-site Z-drive is Π²-odd,
    ///         so the halves are 1.341640e-1 and 1.362440e-1. Expected BROKEN.</item>
    /// </list>
    /// All five use <see cref="Tolerance"/> = 1e-10. That number separates the measured
    /// values rather than certifying exactness: witnesses 3 and 4 read exactly 0.0 and witness
    /// 5 reads 7.692080e-3, so the cut has about seven orders of room on either side and any
    /// rel asym above it counts as BROKEN. Witnesses 1 and 2 are DEGENERATE and read NaN, and
    /// the tolerance never sees them; they used to be counted into the "1-4" above, which was
    /// the vacuous reading in miniature: two of the four values the separation rested on were
    /// zeros with nothing behind them.</summary>
    public static IReadOnlyList<LindbladBitBPiBalanceWitness> StandardSet(ChainSystem chain)
    {
        if (chain is null) throw new ArgumentNullException(nameof(chain));
        if (chain.N != 2)
            throw new ArgumentException(
                $"StandardSet requires N=2 chain (witness 5 is a fixed-ω two-site Z-drive); got N={chain.N}",
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

        // Witness 5 ingredients: single-site Z-drive at ω/2 on each site (Larmor precession).
        const double omega = 0.13;
        var zDriveTerms = new List<PauliTerm>
        {
            PauliTerm.SingleSite(chain.N, 0, PauliLetter.Z, coefficient: (Complex)(omega / 2.0)),
            PauliTerm.SingleSite(chain.N, 1, PauliLetter.Z, coefficient: (Complex)(omega / 2.0)),
        };

        return new[]
        {
            new LindbladBitBPiBalanceWitness(
                witnessName: "Heisenberg_pure_Z_DEGENERATE",
                chain: chain,
                bondTerms: heisenbergTerms,
                gammaT1: null,
                expectedVerdict: "DEGENERATE"),
            new LindbladBitBPiBalanceWitness(
                witnessName: "YZ_ZY_pi2even_DEGENERATE",
                chain: chain,
                bondTerms: yzZyTerms,
                gammaT1: null,
                expectedVerdict: "DEGENERATE"),
            new LindbladBitBPiBalanceWitness(
                witnessName: "XY_pi2odd_balanced",
                chain: chain,
                bondTerms: xyTerms,
                gammaT1: null,
                expectedVerdict: "BALANCED"),
            new LindbladBitBPiBalanceWitness(
                witnessName: "Heisenberg_with_T1_envelope_balanced",
                chain: chain,
                bondTerms: heisenbergTerms,
                gammaT1: 0.1,
                expectedVerdict: "BALANCED"),
            new LindbladBitBPiBalanceWitness(
                witnessName: "Zdrive_with_T1_envelope_BROKEN",
                chain: chain,
                terms: zDriveTerms,
                gammaT1: 0.001,
                expectedVerdict: "BROKEN"),
        };
    }
}
