using System;
using System.Collections.Generic;
using System.Linq;
using System.Numerics;
using RCPsiSquared.Core.ChainSystems;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Pauli;

namespace RCPsiSquared.Diagnostics.F87;

/// <summary>An F103 §7 bipartite-chirality witness: a named diagonal-cell pair (k-body
/// templates + dephase letter) paired with its expected F87 class. Self-computes, lazily, both
/// the bipartite criterion's verdict and the actual F87 trichotomy via
/// <see cref="BipartiteChirality.Classify"/>, and checks they agree.
///
/// <para>Mirrors <see cref="F87CanonicalWitness"/> and
/// <see cref="Polarity.LindbladBitBPiBalanceWitness"/>: expected / actual / predicted with a
/// <see cref="Matches"/> boolean. A failure signals either a broken witness (wrong N, dephase
/// letter, or term set) or a regression in the criterion vs the trichotomy.</para>
///
/// <para><b>Claim (Tier1Candidate):</b> a diagonal-cell pair is soft iff its hopping graph is
/// bipartite in the dephasing basis. <c>bipartite ⟹ soft</c> is derived from the palindrome
/// plus a chiral sublattice K (<see cref="Core.Symmetry.ChiralKClaim"/>); the converse
/// non-bipartite ⟹ hard is now derived (PROOF_F103 §7.5, 2026-06-04) modulo the first-order-block
/// premise: the K3 triangle obstructs the chiral functional that would supply the gain channel's
/// −N reflection mode and pair its +N population Perron mode, and any palindromizer forces a
/// spectral palindrome (no non-chiral escape). The un-derived first-order reduction keeps this
/// Tier1Candidate.</para>
///
/// <para>Anchor: <c>docs/proofs/PROOF_F103_F87_Z2_CUBED_REFINEMENT.md</c> §7 +
/// <see cref="Core.Symmetry.ChiralKClaim"/> + <see cref="PauliPairTrichotomy"/> +
/// <c>simulations/f87_42_8_bipartite_fullcell.py</c> +
/// <c>simulations/f87_bipartite_chiral_witness.py</c>.</para></summary>
public sealed class F87DiagonalCellBipartiteWitness : Claim
{
    private readonly Lazy<BipartiteChiralityResult> _criterion;

    /// <summary>Human-readable witness name.</summary>
    public string WitnessName { get; }

    /// <summary>Chain providing N, γ₀, bond layout.</summary>
    public ChainSystem Chain { get; }

    /// <summary>k-body Pauli templates summed across the chain by the sliding-window builder.</summary>
    public IReadOnlyList<PauliTerm> Terms { get; }

    /// <summary>Single-letter dephasing axis.</summary>
    public PauliLetter DephaseLetter { get; }

    /// <summary>Expected F87 class (the documented anchor).</summary>
    public TrichotomyClass ExpectedClass { get; }

    /// <summary>Lazily computed criterion result; first access triggers the H/L builds.</summary>
    public Lazy<BipartiteChiralityResult> Criterion => _criterion;

    /// <summary>Whether the hopping graph is bipartite in the dephasing basis.</summary>
    public bool IsBipartite => _criterion.Value.IsBipartite;

    /// <summary>The bipartite criterion's predicted class (Soft if bipartite, else Hard).</summary>
    public TrichotomyClass PredictedClass => _criterion.Value.PredictedClass;

    /// <summary>The actual F87 trichotomy verdict.</summary>
    public TrichotomyClass ActualClass => _criterion.Value.ActualClass;

    /// <summary>True iff the criterion's prediction equals the actual verdict (the core claim).</summary>
    public bool CriterionAgrees => _criterion.Value.Agrees;

    /// <summary>True iff the criterion agrees with the actual verdict and that verdict is the
    /// expected one.</summary>
    public bool Matches => CriterionAgrees && ActualClass == ExpectedClass;

    public F87DiagonalCellBipartiteWitness(
        string witnessName, ChainSystem chain, IReadOnlyList<PauliTerm> terms,
        PauliLetter dephaseLetter, TrichotomyClass expectedClass)
        : base($"F103 §7 bipartite-chirality witness: {witnessName}",
               Tier.Tier1Candidate,
               "docs/proofs/PROOF_F103_F87_Z2_CUBED_REFINEMENT.md §7 + " +
               "compute/RCPsiSquared.Core/Symmetry/ChiralKClaim.cs + " +
               "compute/RCPsiSquared.Diagnostics/F87/PauliPairTrichotomy.cs + " +
               "simulations/f87_42_8_bipartite_fullcell.py + " +
               "simulations/f87_bipartite_chiral_witness.py")
    {
        if (string.IsNullOrWhiteSpace(witnessName))
            throw new ArgumentException("witnessName must be non-empty", nameof(witnessName));
        Chain = chain ?? throw new ArgumentNullException(nameof(chain));
        Terms = terms ?? throw new ArgumentNullException(nameof(terms));
        if (Terms.Count == 0)
            throw new ArgumentException("terms must be non-empty", nameof(terms));

        WitnessName = witnessName;
        DephaseLetter = dephaseLetter;
        ExpectedClass = expectedClass;
        _criterion = new Lazy<BipartiteChiralityResult>(
            () => BipartiteChirality.Classify(chain, terms, dephaseLetter));
    }

    public override string DisplayName =>
        $"F103 §7 witness '{WitnessName}': {DephaseLetter}-deph, expect {ExpectedClass}";

    public override string Summary =>
        $"chain N={Chain.N}, [{TermsLabel()}], {DephaseLetter}-deph; bipartite={IsBipartite} → " +
        $"predicted {PredictedClass}, actual {ActualClass}, expected {ExpectedClass} " +
        $"({(Matches ? "PASS" : "FAIL")})";

    private string TermsLabel() => string.Join(" + ", Terms.Select(t => t.Label));

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("witness name", summary: WitnessName);
            yield return new InspectableNode("chain",
                summary: $"N={Chain.N}, {Chain.Topology}, γ₀={Chain.GammaZero:G6}");
            yield return new InspectableNode("dephase letter", summary: DephaseLetter.ToString());
            yield return new InspectableNode("terms", summary: TermsLabel());
            yield return new InspectableNode("hopping graph bipartite (in dephase basis)",
                summary: IsBipartite ? "yes" : "no");
            yield return new InspectableNode("chiral K verified (KHK = −H)",
                summary: _criterion.Value.ChiralKVerified ? "yes" : "no");
            yield return new InspectableNode("predicted class (bipartite criterion)",
                summary: PredictedClass.ToString());
            yield return new InspectableNode("actual class (F87 trichotomy)",
                summary: ActualClass.ToString());
            yield return new InspectableNode("expected class", summary: ExpectedClass.ToString());
            yield return new InspectableNode("criterion agrees with actual",
                summary: CriterionAgrees ? "yes" : "no");
            yield return new InspectableNode("verdict", summary: Matches ? "PASS" : "FAIL");
        }
    }

    /// <summary>Four canonical witnesses spanning the diagonal-cell mechanism:
    /// <list type="number">
    ///   <item><b>XXZ_ZXX_Z_soft</b>: Z-dephasing, diagonal letters at the window ends (even
    ///         hopping cycles), bipartite, soft.</item>
    ///   <item><b>XXZ_XZX_Z_hard_oddcycle</b>: Z-dephasing, opposite position-parity closes an
    ///         odd cycle, non-bipartite, hard (rule (b)).</item>
    ///   <item><b>ZZZ_XXZ_Z_hard_template</b>: Z-dephasing, ZZZ is a pure-D template (diagonal),
    ///         it lifts H's diagonal so no chiral K exists, non-bipartite, hard (rule (a)).</item>
    ///   <item><b>ZZX_XZZ_X_soft</b>: X-dephasing, the Z↔X mirror of witness 1; bipartite only
    ///         in the X eigenbasis, so it exercises the dephase-letter-relative rotation.</item>
    /// </list></summary>
    public static IReadOnlyList<F87DiagonalCellBipartiteWitness> StandardSet(ChainSystem chain)
    {
        if (chain is null) throw new ArgumentNullException(nameof(chain));
        static PauliTerm T(params PauliLetter[] letters) => new(letters, Complex.One);
        const PauliLetter X = PauliLetter.X, Z = PauliLetter.Z;
        return new[]
        {
            new F87DiagonalCellBipartiteWitness("XXZ_ZXX_Z_soft", chain,
                new[] { T(X, X, Z), T(Z, X, X) }, Z, TrichotomyClass.Soft),
            new F87DiagonalCellBipartiteWitness("XXZ_XZX_Z_hard_oddcycle", chain,
                new[] { T(X, X, Z), T(X, Z, X) }, Z, TrichotomyClass.Hard),
            new F87DiagonalCellBipartiteWitness("ZZZ_XXZ_Z_hard_template", chain,
                new[] { T(Z, Z, Z), T(X, X, Z) }, Z, TrichotomyClass.Hard),
            new F87DiagonalCellBipartiteWitness("ZZX_XZZ_X_soft", chain,
                new[] { T(Z, Z, X), T(X, Z, Z) }, X, TrichotomyClass.Soft),
        };
    }
}
