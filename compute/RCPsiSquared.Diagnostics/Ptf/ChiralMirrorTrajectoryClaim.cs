using System;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Pauli;
using RCPsiSquared.Core.States;
using ComplexMatrix = MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex>;

namespace RCPsiSquared.Diagnostics.Ptf;

/// <summary>The K₁ chiral mirror rate law, PTF's only surviving exact law (EQ-014), derived
/// 2026-06-10 as a site-wise trajectory identity. K₁ = Π_{l odd} Z_l (the odd-sublattice Z
/// product, 0-indexed) anticommutes with every XY bond term and with the J-bond defect
/// V = δJ·½(X_bX_{b+1} + Y_bY_{b+1}), while commuting with every dephasing Z_l. Conjugation by
/// K₁ (flips H+V to −(H+V), dissipator invariant) chained with complex conjugation (H, V real;
/// flips −(H+V) back) gives P_i(t; H+V, K₁ψ) = P_i(t; H+V, ψ) exactly for every site i and
/// time t, for any real initial state ψ. The single-excitation sine modes satisfy
/// K₁ψ_k = ψ_{N+1−k} exactly, and any sublattice sign is absorbed by the U(1) phase
/// e^{iπN̂} = Π_l Z_l, so the PTF pair states obey P_i(t; φ_k) = P_i(t; φ_{N+1−k}). The
/// EQ-014 Σ-mirror law Σ_i f_i(ψ_k) = Σ_i f_i(ψ_{N+1−k}) (machine-exact at N = 5, 7, 8) is
/// the corollary obtained by summing the site-wise f_i equality over sites.
///
/// <para>The eigenvector/dynamics side of the same sublattice chirality whose eigenvalue side
/// is <see cref="Core.Symmetry.ChiralKClaim"/> (spectrum inversion E_{N+1−k} = −E_k;
/// bipartite ⟹ soft). Anchor: <c>docs/proofs/PROOF_PTF_CHIRAL_MIRROR_RATE_LAW.md</c> +
/// <c>simulations/ptf_chiral_mirror_trajectory.py</c> + <c>review/EQ014_FINDINGS.md</c>.</para></summary>
public sealed class ChiralMirrorTrajectoryClaim : Claim
{
    /// <summary>One pure-algebra self-check of the K₁ machinery on small dense matrices.</summary>
    public readonly record struct BatteryCase(string Name, string Detail, string Expected, string Actual)
    {
        public bool Passes => string.Equals(Expected, Actual, StringComparison.Ordinal);
    }

    public IReadOnlyList<BatteryCase> Cases { get; }
    public int PassCount => Cases.Count(c => c.Passes);

    public ChiralMirrorTrajectoryClaim()
        : base("PTF K₁ chiral mirror: P_i(t; φ_k) = P_i(t; φ_{N+1−k}) site-wise trajectory identity (K₁-conjugation + complex conjugation + U(1) sign absorption); the EQ-014 Σ f_i mirror law is the summed corollary",
               Tier.Tier1Derived,
               "docs/proofs/PROOF_PTF_CHIRAL_MIRROR_RATE_LAW.md + " +
               "simulations/ptf_chiral_mirror_trajectory.py + " +
               "review/EQ014_FINDINGS.md (the Σ-law, machine-exact at N = 5, 7, 8)")
    {
        Cases = BuildBattery();
    }

    public string Step1Algebra =>
        "K₁ = Π_{l odd} Z_l: every XY bond term has exactly one odd endpoint, so K₁HK₁ = −H and " +
        "K₁VK₁ = −V (the defect bond included); K₁ is a Z string, so [K₁, Z_l] = 0 and the " +
        "Z-dephasing dissipator is K₁-invariant (any site rates).";

    public string Step2Conjugation =>
        "ρ̃(t) = K₁ρ(t)K₁ solves the Lindblad equation with H+V replaced by −(H+V) and the same " +
        "dissipator; K₁ acts per site as I or Z, so single-site reduced states transform unitarily " +
        "and site purities are K₁-invariant: P_i(t; −(H+V), K₁ψ) = P_i(t; H+V, ψ).";

    public string Step3Reality =>
        "H and V are real in the computational basis and the Z-dephasing action is real, so complex " +
        "conjugation maps −H' dynamics to +H' dynamics on the conjugated state; for real initial " +
        "states P_i(t; −H', ψ) = P_i(t; H', ψ). The pair states φ_k are real.";

    public string Step4Modes =>
        "K₁ψ_k = ψ_{N+1−k} exactly (sine identity, 0-indexed odd sublattice; no sign); the " +
        "even-sublattice complement picks a minus sign, absorbed by e^{iπN̂} = Π_l Z_l, which " +
        "commutes with H, V, and the dissipator and preserves site purities. Hence " +
        "P_i(t; φ_k) = P_i(t; φ_{N+1−k}) for every site and time, exactly.";

    public string Scope =>
        "Needs: real H and V, K₁-odd H and V (one chirality for chain AND defect), K₁-invariant " +
        "dissipator, real initial state modulo a U(1) phase. Does NOT need: uniformity, the specific " +
        "XY form, single-excitation states, or perturbation theory (exact at every δJ). Breaks under: " +
        "Z-fields (K₁-even component), complex hopping phases, T1 amplitude damping.";

    public string History =>
        "Found 2026-04 in EQ-014 as a Σ-law (the survivor of the closure-law retraction; " +
        "Σ ln α_i = 0 fell, the k ↔ N+1−k pairing stayed machine-exact at N = 5, 7, 8). Lived " +
        "untyped in two docs and a comment in PerturbationMatrixElements.cs. Derived 2026-06-10 as " +
        "a site-wise trajectory identity via the involution + sign-table idiom of the " +
        "windowed-converse wave, strengthening it twice: Σ-level → site-wise, fitted-rate → trajectory.";

    public string Sibling =>
        "ChiralKClaim rides the same sublattice chirality on the eigenvalue side (KHK = −H ⟹ " +
        "E_{N+1−k} = −E_k spectrum inversion; bipartite ⟹ soft). This claim is the " +
        "eigenvector/dynamics side: the same K₁ forces whole purity trajectories to coincide pairwise.";

    public override string DisplayName =>
        "PTF K₁ chiral mirror (site-wise trajectory identity, Tier1Derived)";

    public override string Summary =>
        "K₁-conjugation + reality + U(1) absorption ⟹ P_i(t; φ_k) = P_i(t; φ_{N+1−k}) exactly for " +
        "every site and time; the EQ-014 Σ f_i mirror law (machine-exact at N = 5, 7, 8) is the " +
        $"summed corollary; {PassCount}/{Cases.Count} PASS ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("Step 1: algebra (K₁HK₁ = −H, K₁VK₁ = −V, [K₁, Z_l] = 0)", summary: Step1Algebra);
            yield return new InspectableNode("Step 2: unitary conjugation (purities K₁-invariant)", summary: Step2Conjugation);
            yield return new InspectableNode("Step 3: complex conjugation (real H, V flip back)", summary: Step3Reality);
            yield return new InspectableNode("Step 4: mode mapping + U(1) sign absorption", summary: Step4Modes);
            yield return new InspectableNode("Scope (what it needs, what it does not)", summary: Scope);
            yield return new InspectableNode("History (EQ-014 Σ-law → site-wise identity)", summary: History);
            yield return new InspectableNode("Sibling: ChiralKClaim (eigenvalue side of the same K₁)", summary: Sibling);
            foreach (var c in Cases)
                yield return new InspectableNode(c.Name,
                    summary: $"{c.Detail}; expected {c.Expected}, got {c.Actual}, " + (c.Passes ? "PASS" : "FAIL"));
        }
    }

    /// <summary>Pure-algebra battery on 16×16 dense matrices at N = 4 (no Liouvillian, no
    /// propagation): K₁-parity of the XY chain and the defect bond, dissipator invariance,
    /// the exact sine-mode mapping K₁ψ_k = ψ_{N+1−k}, and the U(1) phase's commutation.</summary>
    private static IReadOnlyList<BatteryCase> BuildBattery()
    {
        const int N = 4;
        const double tol = 1e-12;
        string Zero(double v) => v < tol ? "0" : v.ToString("E3", CultureInfo.InvariantCulture);

        ComplexMatrix BondTerm(int l) =>
            0.5 * (PauliString.SiteOp(N, l, PauliLetter.X) * PauliString.SiteOp(N, l + 1, PauliLetter.X)
                 + PauliString.SiteOp(N, l, PauliLetter.Y) * PauliString.SiteOp(N, l + 1, PauliLetter.Y));

        var H = BondTerm(0) + BondTerm(1) + BondTerm(2);          // uniform XY chain, J = 1
        var V = 0.1 * BondTerm(0);                                 // δJ = 0.1 defect on bond (0,1)
        var K1 = PauliString.SiteOp(N, 1, PauliLetter.Z) * PauliString.SiteOp(N, 3, PauliLetter.Z);
        var uOne = PauliString.SiteOp(N, 0, PauliLetter.Z) * PauliString.SiteOp(N, 1, PauliLetter.Z)
                 * PauliString.SiteOp(N, 2, PauliLetter.Z) * PauliString.SiteOp(N, 3, PauliLetter.Z);

        var cases = new List<BatteryCase>
        {
            new(Name: "K₁HK₁ = −H (chain K-parity)",
                Detail: "‖K₁HK₁ + H‖_F at N=4, uniform XY",
                Expected: "0",
                Actual: Zero((K1 * H * K1 + H).FrobeniusNorm())),
            new(Name: "K₁VK₁ = −V (defect K-parity)",
                Detail: "‖K₁VK₁ + V‖_F, δJ=0.1 on bond (0,1)",
                Expected: "0",
                Actual: Zero((K1 * V * K1 + V).FrobeniusNorm())),
            new(Name: "[K₁, Z_l] = 0 (dissipator invariance)",
                Detail: "max_l ‖K₁Z_l − Z_lK₁‖_F",
                Expected: "0",
                Actual: Zero(Enumerable.Range(0, N).Max(l =>
                {
                    var zl = PauliString.SiteOp(N, l, PauliLetter.Z);
                    return (K1 * zl - zl * K1).FrobeniusNorm();
                }))),
            new(Name: "U(1) phase commutes with H",
                Detail: "‖[Π_l Z_l, H]‖_F (e^{iπN̂} sign absorber)",
                Expected: "0",
                Actual: Zero((uOne * H - H * uOne).FrobeniusNorm())),
        };

        double worstMode = 0.0;
        for (int k = 1; k <= N; k++)
        {
            var mapped = K1 * BondingMode.Build(N, k);
            var partner = BondingMode.Build(N, N + 1 - k);
            worstMode = Math.Max(worstMode, (mapped - partner).L2Norm());
        }
        cases.Add(new BatteryCase(
            Name: "K₁ψ_k = ψ_{N+1−k} exactly (no sign)",
            Detail: "max_k ‖K₁ψ_k − ψ_{N+1−k}‖₂ at N=4",
            Expected: "0",
            Actual: Zero(worstMode)));

        return cases;
    }
}
