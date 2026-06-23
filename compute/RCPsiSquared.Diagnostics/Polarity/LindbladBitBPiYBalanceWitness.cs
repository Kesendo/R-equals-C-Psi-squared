using System.Numerics;
using RCPsiSquared.Core.ChainSystems;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Pauli;

namespace RCPsiSquared.Diagnostics.Polarity;

/// <summary>An F112-Y polarity-balance witness: a named (chain + H + c_ops + γ_T1) triple
/// paired with the expected verdict (<c>"BALANCED"</c> or <c>"BROKEN"</c>) that the
/// computed F112-Y polarity asymmetry should produce under <b>Y-dephase Π_Y polarity</b>
/// decomposition. Self-computes the actual asymmetry lazily via
/// <see cref="PolarityCoordinates.Decompose(ChainSystem, IReadOnlyList{PauliTerm}, double?, PauliLetter)"/>
/// with <c>dephaseLetter = PauliLetter.Y</c>, so listing many witnesses does not trigger
/// N L-builds upfront.
///
/// <para>Mirrors <see cref="LindbladBitBPiBalanceWitness"/> (F112-Z) exactly except for
/// the Y-dephase routing. The 5-witness <see cref="StandardSet"/> REUSES the same
/// term-lists as the F112-Z witness set because the F112-Y hypothesis on c (bit_b-
/// homogeneous) is identical to F112-Z (same bit_b axis per F38); only the polarity
/// decomposition axis differs (Π_Y vs Π_Z).</para>
///
/// <para><b>F112-Y typed scope (Tier1Derived)</b>: for any Lindblad-form Liouvillian
/// with any H (Hermitian or non-Hermitian) and each c_k bit_b-homogeneous, the asymmetry
/// is exactly 0 bit-exact. See <see cref="LindbladBitBPiYBalance"/>.</para>
///
/// <para><b>Welle 13 derivation route</b>: F112-Y is derived via Route 1 (axis-direct
/// re-run of Welle-11 lemmas with d = Y substituted for d = Z); D-conjugation from
/// F112-Z is NOT available because D-conjugation lacks a Hilbert-space unitary lift
/// (per <c>PROOF_F112_CROSS_DEPHASE_VIA_KLEIN_V4.md</c> section (d) Remark).</para>
///
/// <para>Anchor: <c>docs/ANALYTICAL_FORMULAS.md</c> F112 +
/// <c>docs/proofs/PROOF_F112_CROSS_DEPHASE_VIA_KLEIN_V4.md</c> +
/// <c>docs/proofs/PROOF_F112_LINDBLAD_BIT_B_PI_BALANCE.md</c> +
/// <see cref="LindbladBitBPiYBalance"/> + <see cref="LindbladBitBPiBalance"/> +
/// <see cref="PolarityCoordinates"/>.</para>
/// </summary>
public sealed class LindbladBitBPiYBalanceWitness : Claim
{
    private readonly Lazy<PolarityCoordinatesResult> _polarity;

    /// <summary>Human-readable witness name.</summary>
    public string WitnessName { get; }

    /// <summary>Chain providing N, γ₀, bond layout.</summary>
    public ChainSystem Chain { get; }

    /// <summary>k-body Pauli term list representing the Hamiltonian.</summary>
    public IReadOnlyList<PauliTerm> Terms { get; }

    /// <summary>Optional uniform per-site σ⁻ T1 rate. Null or zero means pure Y-dephasing.</summary>
    public double? GammaT1 { get; }

    /// <summary>Expected verdict: <c>"BALANCED"</c> or <c>"BROKEN"</c>.</summary>
    public string ExpectedVerdict { get; }

    /// <summary>True iff <see cref="ExpectedVerdict"/> = <c>"BALANCED"</c>.</summary>
    public bool ExpectedBalanced => ExpectedVerdict == "BALANCED";

    /// <summary>Relative-asymmetry threshold below which the witness reads as
    /// <c>"BALANCED"</c>. Default 1e-10 matches the bit-exact balance scale of F112-Y's
    /// typed scope.</summary>
    public double Tolerance { get; }

    /// <summary>Lazily computed polarity decomposition result. First access triggers the
    /// L-build and Π_Y-decomposition via <see cref="PolarityCoordinates.Decompose"/>.</summary>
    public Lazy<PolarityCoordinatesResult> Polarity => _polarity;

    /// <summary>Actual relative asymmetry: |Asymmetry| / max(‖M‖², 1e-15). Forces lazy
    /// evaluation of <see cref="Polarity"/>.</summary>
    public double ActualRelativeAsymmetry
    {
        get
        {
            var pol = _polarity.Value;
            return Math.Abs(pol.Asymmetry) / Math.Max(pol.MNormSquared, 1e-15);
        }
    }

    /// <summary>Actual verdict: <c>"BALANCED"</c> if
    /// <see cref="ActualRelativeAsymmetry"/> ≤ <see cref="Tolerance"/>, else
    /// <c>"BROKEN"</c>. Forces lazy evaluation.</summary>
    public string ActualVerdict =>
        ActualRelativeAsymmetry <= Tolerance ? "BALANCED" : "BROKEN";

    /// <summary>True iff <see cref="ActualVerdict"/> matches
    /// <see cref="ExpectedVerdict"/> bit-exactly.</summary>
    public bool Matches => ActualVerdict == ExpectedVerdict;

    /// <summary>Primary k-body constructor. Threads <c>dephaseLetter = PauliLetter.Y</c>
    /// into <see cref="PolarityCoordinates.Decompose(ChainSystem, IReadOnlyList{PauliTerm}, double?, PauliLetter)"/>.</summary>
    public LindbladBitBPiYBalanceWitness(
        string witnessName,
        ChainSystem chain,
        IReadOnlyList<PauliTerm> terms,
        double? gammaT1,
        string expectedVerdict,
        double tolerance = 1e-10)
        : base($"F112-Y polarity-balance witness: {witnessName}",
               Tier.Tier1Derived,
               "docs/ANALYTICAL_FORMULAS.md F112 + " +
               "docs/proofs/PROOF_F112_CROSS_DEPHASE_VIA_KLEIN_V4.md + " +
               "docs/proofs/PROOF_F112_LINDBLAD_BIT_B_PI_BALANCE.md + " +
               "compute/RCPsiSquared.Core/Symmetry/LindbladBitBPiYBalance.cs + " +
               "compute/RCPsiSquared.Core/Symmetry/LindbladBitBPiBalance.cs + " +
               "compute/RCPsiSquared.Diagnostics/Polarity/PolarityCoordinates.cs + " +
               "simulations/f112_klein_v4_cross_dephase_verify.py")
    {
        if (string.IsNullOrWhiteSpace(witnessName))
            throw new ArgumentException("witnessName must be non-empty", nameof(witnessName));
        if (expectedVerdict != "BALANCED" && expectedVerdict != "BROKEN")
            throw new ArgumentException(
                $"expectedVerdict must be 'BALANCED' or 'BROKEN'; got '{expectedVerdict}'",
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
            () => PolarityCoordinates.Decompose(chain, terms, gammaT1, dephaseLetter: PauliLetter.Y));
    }

    /// <summary>Convenience ctor for bilinear-bond witnesses.</summary>
    public LindbladBitBPiYBalanceWitness(
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
        $"F112-Y witness '{WitnessName}': expected {ExpectedVerdict}";

    public override string Summary =>
        $"chain N={Chain.N}, {Terms.Count} term(s), γ_T1={FormatGamma(GammaT1)}, dephase=Y; " +
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
            yield return new InspectableNode("dephase letter", summary: "Y (F112-Y scope)");
            yield return new InspectableNode("expected verdict", summary: ExpectedVerdict);
            yield return new InspectableNode("actual verdict", summary: ActualVerdict);
            yield return new InspectableNode("agreement", summary: Matches ? "PASS" : "FAIL");
            yield return InspectableNode.RealScalar(
                "relative asymmetry", ActualRelativeAsymmetry, format: "E6");
            yield return InspectableNode.RealScalar("tolerance", Tolerance, format: "E2");
            yield return new InspectableNode("γ_T1", summary: FormatGamma(GammaT1));
            yield return new InspectableNode("terms",
                summary: $"{Terms.Count} term(s); labels = [{string.Join(", ", Terms.Select(t => t.Label))}]");
            yield return InspectableNode.RealScalar(
                "F81 violation (numerical, Y-deph; no closed-form prediction on cross-dephase path)",
                Polarity.Value.F81Violation, format: "E3");
        }
    }

    /// <summary>Five named witnesses spanning the F112-Y BALANCED / BROKEN axes under
    /// Y-dephase. Term-lists are identical to <see cref="LindbladBitBPiBalanceWitness.StandardSet"/>
    /// (F112-Y shares the bit_b-homogeneity hypothesis with F112-Z); only the dephase letter
    /// changes:
    /// <list type="number">
    ///   <item><b>Heisenberg_pure_Y_balanced</b>: H = XX+YY+ZZ + pure Y-dephasing.
    ///         Inside F112-Y typed scope. Expected BALANCED.</item>
    ///   <item><b>YZ_ZY_pi2even_Y_balanced</b>: H = YZ+ZY (Π²_Y-even).
    ///         Inside F112-Y typed scope. Expected BALANCED.</item>
    ///   <item><b>XY_pi2odd_Y_balanced</b>: H = XY (Π²_Y-odd). Inside F112-Y typed scope.
    ///         Expected BALANCED.</item>
    ///   <item><b>Heisenberg_with_T1_envelope_Y_balanced</b>: H = Heisenberg + σ⁻ T1 (γ_T1=0.1)
    ///         under Y-dephase. Under Y-deph the dephasing c is Y (bit_b = 1) per site; σ⁻
    ///         is bit_b-mixed (σ⁻ = (X − i Y)/2, X has bit_b=0 and Y has bit_b=1).
    ///         Outside the F112-Y typed Tier1Derived scope (bit_b-homogeneous c hypothesis
    ///         broken by σ⁻), but the broader empirical envelope (mirroring the F112-Z
    ///         envelope across probes 1-14) still gives BALANCED for idle / XY-bond /
    ///         Heisenberg-style Hermitian H. Expected BALANCED.</item>
    ///   <item><b>Zdrive_with_T1_envelope_Y_BROKEN</b>: H = ω·(Z₀+Z₁)/2 single-site Z-drive
    ///         (ω = 0.13) + σ⁻ T1 (γ_T1 = 0.001) under Y-dephase. The Y-dephase analog of the
    ///         F112-Z f95 BROKEN counterexample; both Z-drive and σ⁻ produce non-Hermitian
    ///         Π_Y-eigenspace coupling identically to the Z-dephase case (the F113 [Z, σ⁻]
    ///         = −2σ⁻ commutator structure is dephase-letter independent at the H × c
    ///         interaction level, with the Y-deph Π_Y producing the same +i/−i imbalance
    ///         as the Z-deph Π_Z up to the Y/Z 2-cycle phase). Expected BROKEN.</item>
    /// </list>
    /// All five use <see cref="Tolerance"/> = 1e-10.</summary>
    public static IReadOnlyList<LindbladBitBPiYBalanceWitness> StandardSet(ChainSystem chain)
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

        const double omega = 0.13;
        var zDriveTerms = new List<PauliTerm>
        {
            PauliTerm.SingleSite(chain.N, 0, PauliLetter.Z, coefficient: (Complex)(omega / 2.0)),
            PauliTerm.SingleSite(chain.N, 1, PauliLetter.Z, coefficient: (Complex)(omega / 2.0)),
        };

        return new[]
        {
            new LindbladBitBPiYBalanceWitness(
                witnessName: "Heisenberg_pure_Y_balanced",
                chain: chain,
                bondTerms: heisenbergTerms,
                gammaT1: null,
                expectedVerdict: "BALANCED"),
            new LindbladBitBPiYBalanceWitness(
                witnessName: "YZ_ZY_pi2even_Y_balanced",
                chain: chain,
                bondTerms: yzZyTerms,
                gammaT1: null,
                expectedVerdict: "BALANCED"),
            new LindbladBitBPiYBalanceWitness(
                witnessName: "XY_pi2odd_Y_balanced",
                chain: chain,
                bondTerms: xyTerms,
                gammaT1: null,
                expectedVerdict: "BALANCED"),
            new LindbladBitBPiYBalanceWitness(
                witnessName: "Heisenberg_with_T1_envelope_Y_balanced",
                chain: chain,
                bondTerms: heisenbergTerms,
                gammaT1: 0.1,
                expectedVerdict: "BALANCED"),
            new LindbladBitBPiYBalanceWitness(
                witnessName: "Zdrive_with_T1_envelope_Y_BROKEN",
                chain: chain,
                terms: zDriveTerms,
                gammaT1: 0.001,
                expectedVerdict: "BROKEN"),
        };
    }
}
