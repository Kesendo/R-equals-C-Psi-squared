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
/// Welle 2 in <c>simulations/_f112_hardware_lens_multi.py</c>).</para>
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
    /// <c>"BALANCED"</c>. Computed as |asymmetry| / max(‖M‖², 1e-15) per
    /// <c>_f112_hardware_lens_multi.py</c>. Default 1e-10 matches the bit-exact balance
    /// scale of F112's typed scope.</summary>
    public double Tolerance { get; }

    /// <summary>Lazily computed polarity decomposition result. First access triggers the
    /// L-build and Π-decomposition via <see cref="PolarityCoordinates.Decompose"/>.</summary>
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
               "simulations/_f112_hardware_lens_multi.py")
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
    ///   <item><b>Heisenberg_pure_Z_balanced</b>: H = XX+YY+ZZ Heisenberg + pure Z-dephasing.
    ///         Inside F112 typed scope (Hermitian Heisenberg, single-Pauli Z dephase). Expected BALANCED.</item>
    ///   <item><b>YZ_ZY_pi2even_balanced</b>: H = YZ+ZY (F108 non-truly Π²-even).
    ///         Inside F112 typed scope. Expected BALANCED.</item>
    ///   <item><b>XY_pi2odd_balanced</b>: H = XY (Π²-odd). Inside F112 typed scope.
    ///         Expected BALANCED.</item>
    ///   <item><b>Heisenberg_with_T1_envelope_balanced</b>: H = Heisenberg + σ⁻ T1 (γ_T1=0.1).
    ///         Outside F112 typed scope (σ⁻ is bit_b-mixed). Broader empirical envelope
    ///         observed across probes 1-14 (and re-confirmed Welle 2 across cusp_slowing,
    ///         block_cpsi, chain_gamma0 Tier-A datasets): idle / XY-bond / Heisenberg-style
    ///         Hermitian H combined with σ⁻ T1 still balances bit-exact. Expected BALANCED.</item>
    ///   <item><b>Zdrive_with_T1_envelope_BROKEN</b>: H = ω·(Z₀+Z₁)/2 single-site Z-drive
    ///         (ω = 0.13) + σ⁻ T1 (γ_T1 = 0.001). The structural counterexample discovered
    ///         Welle 2, 2026-05-26, via the f95_angle_steering Tier-A dataset (Kingston,
    ///         2026-05-16). Synthetic isolation in <c>_f112_hardware_lens_multi.py</c>
    ///         confirms rel asym ≈ 3.85e-3 bit-exact. Expected BROKEN.</item>
    /// </list>
    /// All five use <see cref="Tolerance"/> = 1e-10 to keep the BALANCED/BROKEN cut at
    /// the bit-exact scale of F112 (any rel asym > 1e-10 counts as BROKEN).</summary>
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
                witnessName: "Heisenberg_pure_Z_balanced",
                chain: chain,
                bondTerms: heisenbergTerms,
                gammaT1: null,
                expectedVerdict: "BALANCED"),
            new LindbladBitBPiBalanceWitness(
                witnessName: "YZ_ZY_pi2even_balanced",
                chain: chain,
                bondTerms: yzZyTerms,
                gammaT1: null,
                expectedVerdict: "BALANCED"),
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
