using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.F1;

/// <summary>F49 cross-term closed form under non-uniform γ (Tier 1 derived; bit-exact at
/// N = 3, 4, 5 across Heisenberg, Ising, XY, soft XY+YX H-classes).
///
/// <para>Extends the F49 uniform-γ identity <c>‖{L_H, L_Dc}‖²_F = 4γ²·(N−2)·‖L_H‖²_F</c>
/// to arbitrary site-dependent {γ_l}, with the F1-centered dissipator L_Dc := L_D + σ·I
/// (σ = Σ_l γ_l). The closed form splits the cross-term into a spectator part (per-bond,
/// depending on Σ_{m ∉ bond} γ_m²) plus a bond-asymmetry part (per-bond, depending on
/// (γ_i − γ_j)² with a per-Pauli-class coefficient G(bond, H) = <c>4·‖L_{ZZ-part of H}^bond‖²_F</c>;
/// only the ZZ-fraction of each bond Hamiltonian carries (γ_i − γ_j)² sensitivity):</para>
///
/// <code>
///     ‖{L_H, L_Dc}‖²_F  =  4 · Σ_b ‖L_H^bond_b‖²_F · Σ_{m ∉ bond_b} γ_m²        (spectator)
///                        +     Σ_b G(bond_b, H) · (γ_{i_b} − γ_{j_b})²            (bond-asymmetry)
/// </code>
///
/// <para>Per-class G fractions (G/‖L_H^bond‖²_F): see the <see cref="GHeisenbergFraction"/>,
/// <see cref="GIsingFraction"/>, <see cref="GXyFraction"/>, <see cref="GSoftXyYxFraction"/>
/// constants.</para>
///
/// <para><b>Convention on <c>bondNormSquared</c>:</b> the per-bond Frobenius norm squared
/// <c>‖L_H^bond_b‖²_F</c> is computed on the FULL N-qubit Pauli-string operator space (spectator
/// I-tensors included). For a single Heisenberg J = 1 bond this equals 384 at N = 3, 1536 at N = 4,
/// 6144 at N = 5; each step multiplies by 4 (one additional spectator contributes <c>tr(I_4) = 4</c>
/// to the Frobenius-norm tensor calculation). Intrinsic local 2-qubit norms are 96 / 32 / 64 / 64
/// for Heisenberg / Ising / XY / XY+YX respectively. Matches <c>_bond_LH_norm_sq</c> in
/// <c>simulations/_f49_nonuniform_gamma_crossterm_verify.py</c>; the Predict methods are linear in
/// <c>bondNormSquared</c>, so any consistent convention works as long as caller and theorem agree.</para>
///
/// <para>Uniform γ_l ≡ γ recovers F49's <c>4γ²·(N−2)·‖L_H‖²_F</c>
/// ([F49 / PROOF_CROSS_TERM_FORMULA]): bond-asymmetry vanishes (γ_i = γ_j); the spectator part
/// collapses via the disjoint-bond-supports lemma.</para>
///
/// <para>Anchor: <c>docs/proofs/PROOF_F49_NONUNIFORM_GAMMA_EXTENSION.md</c> (Steps 1-6).
/// Closes the "Open follow-ups" item of <c>docs/proofs/PROOF_F1_NONUNIFORM_GAMMA.md</c>.
/// Verification: <c>simulations/_f49_nonuniform_gamma_crossterm_verify.py</c> (Phase 1
/// exploratory at commit <c>1c6701c</c>, Phase 2 assertion block; bit-exact at N = 3, 4, 5
/// across 4 H-classes).</para>
/// </summary>
public sealed class F49NonUniformCrossTermClaim : Claim
{
    /// <summary>Prefactor on the spectator part: <c>4 · Σ_b ‖L_H^bond_b‖² · Σ_{m ∉ bond_b} γ_m²</c>.
    /// The 4 comes from <c>(d_α^c + d_β^c)² = 4·(spectator-sum)² + A²·(γ_i−γ_j)²</c> after the
    /// bond-sum rule eliminates the cross-term (Step 3 of the proof).</summary>
    public const double SpectatorPrefactor = 4.0;

    /// <summary>G(bond, Heisenberg) / ‖L_H^bond‖²_F = 4/3. Heisenberg J·(XX+YY+ZZ) decomposes
    /// 2/3 XY-class (XX + YY, A = 0) plus 1/3 ZZ-class (ZZ, A = ±2); the ZZ third inflates
    /// the bond-asymmetry coefficient to <c>4·(1/3) = 4/3</c>.</summary>
    public const double GHeisenbergFraction = 4.0 / 3.0;

    /// <summary>G(bond, Ising) / ‖L_H^bond‖²_F = 4. Ising J·ZZ is entirely ZZ-class, so
    /// every transition has A = ±2 ⟹ <c>G = 4·‖L_H^bond‖²</c>.</summary>
    public const double GIsingFraction = 4.0;

    /// <summary>G(bond, XY) / ‖L_H^bond‖²_F = 0. XY model J·(XX+YY) has no ZZ component; both
    /// XX and YY are XY-class (both bond Paulis in {X, Y}), giving A = 0 on every transition.</summary>
    public const double GXyFraction = 0.0;

    /// <summary>G(bond, soft XY+YX) / ‖L_H^bond‖²_F = 0. Soft Π²-odd J·(XY+YX) has no ZZ component;
    /// XY and YX are both XY-class (both bond Paulis in {X, Y}), giving A = 0 on every transition.</summary>
    public const double GSoftXyYxFraction = 0.0;

    public F49NonUniformCrossTermClaim()
        : base("F49 non-uniform γ cross-term: ‖{L_H, L_Dc}‖² = 4·Σ_b ‖L_H^bond‖²·Σ_{m∉bond}γ_m² + Σ_b G(bond,H)·(γ_i−γ_j)²",
               Tier.Tier1Derived,
               "docs/proofs/PROOF_F49_NONUNIFORM_GAMMA_EXTENSION.md + " +
               "docs/proofs/PROOF_CROSS_TERM_FORMULA.md (uniform-γ parent) + " +
               "docs/proofs/PROOF_F1_NONUNIFORM_GAMMA.md (sibling H-block γ-independence) + " +
               "compute/RCPsiSquared.Core/F1/F1T1ResidualClosedForm.cs + " +
               "compute/RCPsiSquared.Core/F1/F1DepolResidualClosedForm.cs (sibling dissipator-block closed forms)")
    { }

    /// <summary>Predicted <c>‖{L_H, L_Dc}‖²_F</c> for a Heisenberg chain (XX+YY+ZZ on every nearest-neighbour
    /// bond, uniform J) at chain length N with per-site γ rates {γ_l}. <paramref name="bondNormSquared"/> is
    /// the caller-supplied <c>‖L_H^bond‖²_F</c> on the full N-qubit operator space (e.g., 384 at N=3, J=1;
    /// 1536 at N=4, J=1; 6144 at N=5, J=1; see XML <c>summary</c>).</summary>
    public static double PredictHeisenbergChain(int N, IReadOnlyList<double> gamma, double bondNormSquared) =>
        PredictChain(N, gamma, bondNormSquared, GHeisenbergFraction);

    /// <summary>Predicted <c>‖{L_H, L_Dc}‖²_F</c> for an Ising chain (ZZ-only on every nearest-neighbour
    /// bond, uniform J) at chain length N with per-site γ rates {γ_l}.</summary>
    public static double PredictIsingChain(int N, IReadOnlyList<double> gamma, double bondNormSquared) =>
        PredictChain(N, gamma, bondNormSquared, GIsingFraction);

    /// <summary>Predicted <c>‖{L_H, L_Dc}‖²_F</c> for an XY chain (XX+YY on every nearest-neighbour bond,
    /// uniform J) at chain length N with per-site γ rates {γ_l}. G-fraction is 0, so the formula reduces
    /// to the spectator-only part.</summary>
    public static double PredictXyChain(int N, IReadOnlyList<double> gamma, double bondNormSquared) =>
        PredictChain(N, gamma, bondNormSquared, GXyFraction);

    /// <summary>Predicted <c>‖{L_H, L_Dc}‖²_F</c> for a soft Π²-odd chain (XY+YX on every nearest-neighbour
    /// bond, uniform J) at chain length N with per-site γ rates {γ_l}. G-fraction is 0, so the formula
    /// reduces to the spectator-only part.</summary>
    public static double PredictSoftXyYxChain(int N, IReadOnlyList<double> gamma, double bondNormSquared) =>
        PredictChain(N, gamma, bondNormSquared, GSoftXyYxFraction);

    /// <summary>Predicted <c>‖{L_H, L_Dc}‖²_F</c> for an arbitrary graph topology with arbitrary per-bond
    /// Hamiltonian. <paramref name="bondEdges"/> is the list of bonds (each (i, j) with 0 ≤ i &lt; j &lt; N);
    /// <paramref name="bondNormSquaredPerBond"/> is the per-bond <c>‖L_H^bond_b‖²_F</c> on the full N-qubit
    /// operator space; <paramref name="gFractionPerBond"/> is the per-bond G(bond_b, H) / ‖L_H^bond_b‖²_F
    /// (4/3 for Heisenberg, 4 for Ising, 0 for XY-class; see the per-class constants).</summary>
    public static double Predict(
        int N,
        IReadOnlyList<double> gamma,
        IReadOnlyList<(int i, int j)> bondEdges,
        IReadOnlyList<double> bondNormSquaredPerBond,
        IReadOnlyList<double> gFractionPerBond)
    {
        ValidateCommon(N, gamma);
        if (bondNormSquaredPerBond.Count != bondEdges.Count)
            throw new ArgumentException(
                $"bondNormSquaredPerBond length ({bondNormSquaredPerBond.Count}) must equal bondEdges length ({bondEdges.Count}).",
                nameof(bondNormSquaredPerBond));
        if (gFractionPerBond.Count != bondEdges.Count)
            throw new ArgumentException(
                $"gFractionPerBond length ({gFractionPerBond.Count}) must equal bondEdges length ({bondEdges.Count}).",
                nameof(gFractionPerBond));
        for (int b = 0; b < bondEdges.Count; b++)
        {
            (int i, int j) = bondEdges[b];
            if (i < 0 || j < 0 || i >= N || j >= N)
                throw new ArgumentException(
                    $"bondEdges[{b}] = ({i}, {j}) is out of range; both indices must lie in [0, {N}).",
                    nameof(bondEdges));
            if (i >= j)
                throw new ArgumentException(
                    $"bondEdges[{b}] = ({i}, {j}) violates i < j ordering.",
                    nameof(bondEdges));
        }

        double total = 0.0;
        for (int b = 0; b < bondEdges.Count; b++)
        {
            (int i, int j) = bondEdges[b];
            total += BondContribution(N, gamma, i, j, bondNormSquaredPerBond[b], gFractionPerBond[b]);
        }
        return total;
    }

    /// <summary>Shared chain specialisation: nearest-neighbour bonds (0,1), (1,2), …, (N−2, N−1)
    /// with uniform per-bond <c>‖L_H^bond‖²_F</c> and uniform G-fraction.</summary>
    private static double PredictChain(int N, IReadOnlyList<double> gamma, double bondNormSquared, double gFraction)
    {
        ValidateCommon(N, gamma);
        if (bondNormSquared < 0)
            throw new ArgumentOutOfRangeException(nameof(bondNormSquared),
                $"bondNormSquared must be ≥ 0; got {bondNormSquared}");
        double total = 0.0;
        for (int b = 0; b < N - 1; b++)
            total += BondContribution(N, gamma, i: b, j: b + 1, bondNormSquared, gFraction);
        return total;
    }

    /// <summary>Single-bond contribution to <c>‖{L_H, L_Dc}‖²_F</c>:
    /// <c>4·bondNormSq·Σ_{m ∉ {i,j}} γ_m² + gFraction·bondNormSq·(γ_i − γ_j)²</c>.</summary>
    private static double BondContribution(
        int N, IReadOnlyList<double> gamma, int i, int j, double bondNormSq, double gFraction)
    {
        double spectatorSum = 0.0;
        for (int m = 0; m < N; m++)
        {
            if (m == i || m == j) continue;
            spectatorSum += gamma[m] * gamma[m];
        }
        double deltaGamma = gamma[i] - gamma[j];
        return SpectatorPrefactor * bondNormSq * spectatorSum
             + gFraction * bondNormSq * deltaGamma * deltaGamma;
    }

    private static void ValidateCommon(int N, IReadOnlyList<double> gamma)
    {
        if (N < 2) throw new ArgumentOutOfRangeException(nameof(N), $"N must be ≥ 2; got {N}");
        if (gamma.Count != N)
            throw new ArgumentException(
                $"gamma length ({gamma.Count}) must equal N ({N}).", nameof(gamma));
    }

    public override string DisplayName =>
        "F49 non-uniform γ cross-term: ‖{L_H, L_Dc}‖² = 4·Σ_b ‖L_H^bond‖²·Σ_{m∉bond}γ_m² + Σ_b G(bond,H)·(γ_i−γ_j)²";

    public override string Summary =>
        "F49 uniform-γ extension under site-dependent γ; spectator part 4·Σ_b ‖L_H^bond‖²·Σ_{m∉bond}γ_m² + " +
        "bond-asymmetry part Σ_b G(bond,H)·(γ_i−γ_j)²; G/‖L_H^bond‖² = 4/3 (Heisenberg) / 4 (Ising) / 0 (XY, XY+YX); " +
        "uniform γ recovers F49's 4γ²·(N−2)·‖L_H‖²";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("statement",
                summary: "‖{L_H, L_Dc}‖²_F = 4·Σ_b ‖L_H^bond_b‖²·Σ_{m∉bond_b}γ_m² + Σ_b G(bond_b, H)·(γ_{i_b}−γ_{j_b})²");
            yield return InspectableNode.RealScalar("spectator prefactor", SpectatorPrefactor);
            yield return InspectableNode.RealScalar("G(bond, Heisenberg) / ‖L_H^bond‖² = 4/3 (ZZ is 1/3 of bond norm)",
                GHeisenbergFraction);
            yield return InspectableNode.RealScalar("G(bond, Ising) / ‖L_H^bond‖² = 4 (ZZ is 100% of bond norm)",
                GIsingFraction);
            yield return InspectableNode.RealScalar("G(bond, XY) / ‖L_H^bond‖² = 0 (no ZZ content)",
                GXyFraction);
            yield return InspectableNode.RealScalar("G(bond, soft XY+YX) / ‖L_H^bond‖² = 0 (no ZZ content)",
                GSoftXyYxFraction);
            yield return new InspectableNode("spectator-vs-asymmetry split",
                summary: "Step 3 of the proof: cross-term (d_α^c + d_β^c)² decomposes into " +
                         "4·(spectator-sum)² + A²·(γ_i−γ_j)² after the bond-sum rule kills the spectator-bond " +
                         "cross-term; spectator-sum is per-bond Σ_{m ∉ bond} γ_m·ε_m(α) and squares to " +
                         "‖L_H^bond‖²·Σ_{m ∉ bond} γ_m² after averaging over spectator letters (ε² = 1)");
            yield return new InspectableNode("A-classification per bond Pauli class",
                summary: "Step 5: only ZZ-class bond terms (both Paulis in {I, Z}) generate A = ±2 transitions; " +
                         "XY-class bond terms (both Paulis in {X, Y}: XX, YY, XY, YX) generate A = 0; hence " +
                         "G(bond, H) = 4·‖L_{ZZ-part of H}^bond‖²");
            yield return new InspectableNode("bondNormSquared convention",
                summary: "‖L_H^bond_b‖²_F is the per-bond Frobenius norm on the full N-qubit operator space " +
                         "(spectator I-tensors included). For Heisenberg J=1: 384 at N=3, 1536 at N=4, 6144 at N=5; " +
                         "Ising J=1: 128 / 512 / 2048; XY J=1: 256 / 1024 / 4096; XY+YX J=1: 256 / 1024 / 4096. " +
                         "Predict methods are linear in bondNormSquared; matches the convention of " +
                         "simulations/_f49_nonuniform_gamma_crossterm_verify.py (Phase 1 verification anchor)");
            yield return new InspectableNode("uniform-γ recovery",
                summary: "γ_l ≡ γ ⟹ bond-asymmetry vanishes (γ_i − γ_j = 0), spectator part collapses to " +
                         "4γ²·(N−2)·Σ_b ‖L_H^bond‖² = 4γ²·(N−2)·‖L_H‖² (disjoint-bond-supports lemma); " +
                         "this is the F49 uniform-γ identity PROOF_CROSS_TERM_FORMULA");
            yield return new InspectableNode("verification",
                summary: "bit-exact at N = 3, 4, 5 across Heisenberg / Ising / XY / soft XY+YX " +
                         "(simulations/_f49_nonuniform_gamma_crossterm_verify.py, Phase 1 commit 1c6701c + Phase 2 assertions)");
            yield return new InspectableNode("closes follow-up",
                summary: "the 'Open follow-ups' item in PROOF_F1_NONUNIFORM_GAMMA.md (the F49 non-uniform γ gap " +
                         "at N=3 Heisenberg γ=[0.1, 0.2, 0.3]) is closed by this proof + typed claim");
        }
    }
}
