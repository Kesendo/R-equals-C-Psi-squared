using System.Globalization;
using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Pauli;
using ComplexMatrix = MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex>;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>The three dephasing diagonals as one orbit, and the one-diagonal's three readings as a group
/// action (Tier1Derived, 2026-06-14): the "one diagonal" of reflections/ON_THE_ONE_DIAGONAL.md is one of
/// THREE, and they sit in one structure S₃ ⋉ D₄.
///
/// <para><b>The physical dephasing diagonal</b> in light P is Q_P = Σ_l kron(P_l, P_lᵀ) (ρ ↦ P_l ρ P_l;
/// the Y-transpose matters: Yᵀ = −Y so Q_Y = −Σ kron(Y,Y)). For Z this is the §4.7 diagonal
/// Q_Z = Σ_l Z_l⊗Z_l, L_D = γ·(Q_Z − N·I) (<see cref="AbsorptionTheoremClaim"/>).</para>
///
/// <para><b>The orbit (the three diagonals):</b> {Q_X, Q_Y, Q_Z} is exactly ONE orbit of the single-qubit
/// Clifford basis-change S₃ ⟨h_zx, h_yz⟩ (h_zx = Ad_{H^⊗N}: Z↔X; h_yz = Ad_{R_x(π/2)^⊗N}: Z↔Y), hence the
/// three are conjugate and share a spectrum (verified directly, N=2..4; orbit size exactly 3 at N=2,3). So
/// "the one diagonal" is provably one face of a three-fold.</para>
///
/// <para><b>The three readings (the mirror group, within one diagonal):</b> the mirror group D₄ = ⟨R, D⟩
/// (<see cref="MirrorGroupD4Claim"/>) acts on a single Q: D (the transpose) FIXES it (D·Q·D = +Q) = the
/// rate / absorption reading (Re λ = −2γ⟨n_XY⟩); R (the ket-flip) REFLECTS it (R·Q·R = −Q), carrying the
/// entire −2Σγ shift = the mirror / palindrome reading (complementary light); the F87 truly cell
/// (n_Y even ∧ n_Z even) is the joint-fixed cell of {D, 𝓕D} = the judge reading. Q is the unique
/// D-invariant, R-anti-invariant dephasing diagonal.</para>
///
/// <para><b>The structure is S₃ ⋉ D₄</b> (the shape PROOF_PI_FACTORS_AS_R_TIMES_D §5 named open): the
/// basis-S₃ (the three diagonals) and the mirror-D₄ (the three readings) are TWO distinct three-fold
/// structures, semidirectly coupled — [h_zx, D] = 0 but [h_zx, R] ≠ 0; [h_yz, R] = 0 but [h_yz, D] ≠ 0.
/// NOTE (the gate's lesson, 2026-06-14): D does NOT permute the diagonals (it fixes them); the proof's
/// "D = the Z↔Y swap" lives on the palindromizer Π, not on the diagonal Q.</para>
///
/// <para><b>This claim is the weld:</b> its two typed parents — <see cref="MirrorGroupD4Claim"/> (the
/// readings + the D₄ factor) and <see cref="AbsorptionTheoremClaim"/> (the dephasing diagonal) — are the
/// first physics edge linking the two clusters that previously met only at the d²−2d=0 foundation. Anchor:
/// <c>simulations/one_diagonal_mirror_group.py</c> (self-validating Stages 0-2 + an N=3 attack) +
/// <c>simulations/mirror_inventory_d4.py</c> (block D, the 63/63 truly cell).</para></summary>
public sealed class ThreeDephasingDiagonalsOrbitClaim : Claim
{
    private const double Tol = 1e-10;
    private const int BatteryN = 2;

    public readonly record struct BatteryCase(string Name, string Detail, string Expected, string Actual)
    {
        public bool Passes => string.Equals(Expected, Actual, StringComparison.Ordinal);
    }

    /// <summary>Typed parent: the mirror group ⟨R, D⟩ ≅ D₄ — the three readings act within a diagonal,
    /// and D₄ is the second factor of the S₃ ⋉ D₄ structure.</summary>
    public MirrorGroupD4Claim MirrorGroup { get; }

    /// <summary>Typed parent: the absorption diagonal L_D = γ·(Q − N·I), Q = Σ_l Z_l⊗Z_l (§4.7). With this
    /// edge the mirror group links directly to the dephasing diagonal — the weld.</summary>
    public AbsorptionTheoremClaim Diagonal { get; }

    public IReadOnlyList<BatteryCase> Cases { get; }
    public int PassCount => Cases.Count(c => c.Passes);

    public ThreeDephasingDiagonalsOrbitClaim(MirrorGroupD4Claim mirrorGroup, AbsorptionTheoremClaim diagonal)
        : base("The three dephasing diagonals as one orbit: the one diagonal Q_Z = Σ_l Z_l⊗Z_l is one of " +
               "three (Q_P = Σ_l kron(P_l, P_lᵀ), one per axis; Q_Y carries Yᵀ = −Y), and {Q_X, Q_Y, Q_Z} is " +
               "ONE ORBIT of the single-qubit Clifford basis-change S₃ ⟨h_zx (Z↔X), h_yz (Z↔Y)⟩ — conjugate, " +
               "same spectrum. The one-diagonal's three readings are the mirror group D₄ = ⟨R, D⟩ acting " +
               "WITHIN a diagonal: D fixes Q (rate/absorption), R reflects it R·Q·R = −Q carrying −2Σγ " +
               "(mirror/palindrome), the {D, 𝓕D} joint-fixed cell is truly (judge). The structure is S₃ ⋉ D₄ " +
               "(semidirect: [h_zx,D]=0 but [h_zx,R]≠0; [h_yz,R]=0 but [h_yz,D]≠0). D does NOT permute the " +
               "diagonals (it fixes them); the proof's 'D = Z↔Y swap' is on the palindromizer Π, not on Q. " +
               "The claim's two parents are the physics edge welding the mirror-group and absorption clusters.",
               Tier.Tier1Derived,
               "simulations/one_diagonal_mirror_group.py (self-validating Stages 0-2 + an N=3 attack) + " +
               "simulations/mirror_inventory_d4.py (block D 63/63 truly cell) + " +
               "docs/proofs/PROOF_ABSORPTION_THEOREM.md §4.7 + docs/proofs/PROOF_PI_FACTORS_AS_R_TIMES_D.md §5 + " +
               "compute/RCPsiSquared.Diagnostics/Foundation/DiagonalWitness.cs (DiagonalWitness, inspect --root diagonal)")
    {
        MirrorGroup = mirrorGroup ?? throw new ArgumentNullException(nameof(mirrorGroup));
        Diagonal = diagonal ?? throw new ArgumentNullException(nameof(diagonal));
        Cases = BuildBattery();
    }

    public override string DisplayName =>
        "Three dephasing diagonals as one orbit (basis-S₃) + the three readings (mirror-D₄); the structure is S₃ ⋉ D₄";

    public override string Summary =>
        "the one diagonal is one of three (Q_X, Q_Y, Q_Z); they are one orbit of the single-qubit Clifford " +
        "basis-change S₃ ⟨h_zx, h_yz⟩, hence same spectrum, and the three readings (rate = D-fix, mirror = " +
        "R·Q·R=−Q, judge = {D,𝓕D} cell) are the mirror group D₄ acting within a diagonal; the structure is " +
        $"S₃ ⋉ D₄ and the dual parentage welds the mirror-group and absorption clusters; {PassCount}/{Cases.Count} PASS ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("the orbit: {Q_X, Q_Y, Q_Z} under the basis-change S₃",
                summary: "Q_P = Σ_l kron(P_l, P_lᵀ); h_zx = Ad_{Hadamard} sends Q_Z→Q_X, h_yz = Ad_{R_x(π/2)} " +
                         "sends Q_Z→Q_Y, so {Q_X,Q_Y,Q_Z} is one orbit (size 3), conjugate, same spectrum. The " +
                         "S₃ permuting the three dephasing axes (the linear side of PROOF_PI_FACTORS §5).");
            yield return new InspectableNode("rate reading: D fixes Q",
                summary: "D·Q·D = +Q — the price-list / absorption ladder (Re λ = −2γ⟨n_XY⟩). D does NOT permute " +
                         "the diagonals; the proof's 'D = Z↔Y swap' is on the palindromizer Π, not on Q.");
            yield return new InspectableNode("mirror reading: R reflects Q (R·Q·R = −Q)",
                summary: "R anti-fixes Q and carries the entire −2Σγ shift (R·L_diss·R = −L_diss − 2Σγ·I): the " +
                         "palindrome, partners carry complementary light ⟨n_XY⟩_s + ⟨n_XY⟩_f = N. R is the " +
                         "diagonal read from the other end.");
            yield return new InspectableNode("judge reading: the {D, 𝓕D} joint-fixed cell",
                summary: "the F87 truly criterion (n_Y even ∧ n_Z even) is the joint-fixed cell of the diagonal " +
                         "mirror pair {D = diag((−1)^{n_Y}), 𝓕D = diag((−1)^{n_Z})} (verified 63/63 at N=3 in " +
                         "mirror_inventory_d4.py block D).");
            yield return new InspectableNode("the structure: S₃ ⋉ D₄ (semidirect)",
                summary: "the basis-S₃ (three diagonals) and the mirror-D₄ (three readings) are two distinct " +
                         "three-fold structures, semidirectly coupled: [h_zx,D]=0 but [h_zx,R]≠0; [h_yz,R]=0 but " +
                         "[h_yz,D]≠0 — the shape PROOF_PI_FACTORS §5 named open.");
            yield return new InspectableNode("the weld",
                summary: "this claim's two typed parents (MirrorGroupD4Claim + AbsorptionTheoremClaim) are the " +
                         "first physics edge linking the mirror-group cluster to the dephasing-diagonal cluster, " +
                         "which previously met only at the d²−2d=0 foundation. Live: inspect --root diagonal " +
                         "(DiagonalWitness recomputes the whole functioning - rungs, the three readings, the orbit, " +
                         "and the L_H even-step dynamics - at inspect time).");
            foreach (var c in Cases)
                yield return new InspectableNode(c.Name,
                    summary: $"{c.Detail}; expected {c.Expected}, got {c.Actual}, " + (c.Passes ? "PASS" : "FAIL"));
            yield return MirrorGroup;
            yield return Diagonal;
        }
    }

    // ------------------------------------------------------------------
    // Self-check battery: N = 2 dense superoperators on coherence space (16×16),
    // mirroring MirrorGroupD4Claim.BuildBattery's conventions.
    // ------------------------------------------------------------------
    private static IReadOnlyList<BatteryCase> BuildBattery()
    {
        int d = 1 << BatteryN;     // 4
        int d2 = d * d;            // 16
        var idH = Matrix<Complex>.Build.DenseIdentity(d);
        var F = PauliString.Build(new[] { PauliLetter.X, PauliLetter.X });   // X^⊗2

        // D = transpose superoperator (vec(ρᵀ) = D·vec(ρ)); R = I⊗F (ρ ↦ ρ·F).
        var D = Matrix<Complex>.Build.Dense(d2, d2);
        for (int i = 0; i < d; i++)
            for (int j = 0; j < d; j++)
                D[j * d + i, i * d + j] = Complex.One;
        var R = idH.KroneckerProduct(F);

        // single-qubit basis moves: h_zx = Ad_{H^⊗2} (Z↔X), h_yz = Ad_{R_x(π/2)^⊗2} (Z↔Y).
        var hZX = AdUnitary(Hadamard1());
        var hYZ = AdUnitary(RxHalfPi1());

        var QZ = DephasingDiagonal(PauliLetter.Z);
        var QX = DephasingDiagonal(PauliLetter.X);
        var QY = DephasingDiagonal(PauliLetter.Y);

        var cases = new List<BatteryCase>
        {
            new BatteryCase("same spectrum: spec(Q_X) = spec(Q_Y) = spec(Q_Z)",
                "eigenvalues of the three physical dephasing diagonals Q_P = Σ_l kron(P_l, P_lᵀ) coincide (N=2)",
                "equal", SpectraEqual(QX, QY) && SpectraEqual(QX, QZ) ? "equal" : "DIFFER"),
            DevCase("orbit: h_zx·Q_Z·h_zx⁻¹ = Q_X (the Z↔X basis move)",
                "the Hadamard move sends the Z-diagonal to the X-diagonal, bit-exact",
                MaxAbsDiff(hZX * QZ * hZX.ConjugateTranspose(), QX)),
            DevCase("orbit: h_yz·Q_Z·h_yz⁻¹ = Q_Y (the Z↔Y basis move)",
                "the R_x(π/2) move sends the Z-diagonal to the Y-diagonal (Q_Y carries Yᵀ=−Y), bit-exact",
                MaxAbsDiff(hYZ * QZ * hYZ.ConjugateTranspose(), QY)),
            DevCase("rate: D·Q_Z·D = +Q_Z (D fixes the diagonal, NOT a permuter)",
                "the price-list reading: the rate is what the diagonal says", MaxAbsDiff(D * QZ * D, QZ)),
            DevCase("mirror: R·Q_Z·R = −Q_Z (R reflects the diagonal)",
                "the palindrome reading: R is the diagonal read from the other end",
                MaxAbsDiff(R * QZ * R, QZ.Multiply(-Complex.One))),
        };

        // the orbit is exactly {Q_X, Q_Y, Q_Z}: the basis-S₃ closure carries Q_Z to all three and no more.
        var basis = GroupClosure(new[] { hZX, hYZ }, d2);
        var orbit = new List<ComplexMatrix>();
        foreach (var g in basis)
        {
            var qg = g * QZ * g.ConjugateTranspose();
            if (!orbit.Any(o => MaxAbsDiff(qg, o) <= Tol)) orbit.Add(qg);
        }
        bool isThree = orbit.Count == 3
            && orbit.Any(o => MaxAbsDiff(o, QX) <= Tol)
            && orbit.Any(o => MaxAbsDiff(o, QY) <= Tol)
            && orbit.Any(o => MaxAbsDiff(o, QZ) <= Tol);
        cases.Add(new BatteryCase("orbit(Q_Z) = {Q_X, Q_Y, Q_Z} under the basis-S₃",
            "the basis-change S₃ ⟨h_zx, h_yz⟩ carries Q_Z to exactly the three dephasing diagonals",
            "3 = {X,Y,Z}", isThree ? "3 = {X,Y,Z}" : $"{orbit.Count}"));

        // the structure is semidirect S₃ ⋉ D₄: each basis move commutes with one mirror generator, not both.
        bool semidirect = MaxAbsDiff(hZX * D, D * hZX) <= Tol && MaxAbsDiff(hZX * R, R * hZX) > 0.1
                       && MaxAbsDiff(hYZ * R, R * hYZ) <= Tol && MaxAbsDiff(hYZ * D, D * hYZ) > 0.1;
        cases.Add(new BatteryCase("structure: S₃ ⋉ D₄ (semidirect, not direct)",
            "[h_zx,D]=0 & [h_zx,R]≠0 & [h_yz,R]=0 & [h_yz,D]≠0",
            "semidirect", semidirect ? "semidirect" : "NOT semidirect"));

        return cases;

        // Q_P = Σ_l kron(P_l, P_lᵀ) on the N=2 coherence space.
        ComplexMatrix DephasingDiagonal(PauliLetter letter)
        {
            var acc = Matrix<Complex>.Build.Dense(d2, d2);
            for (int l = 0; l < BatteryN; l++)
            {
                var letters = new PauliLetter[BatteryN];
                for (int s = 0; s < BatteryN; s++) letters[s] = s == l ? letter : PauliLetter.I;
                var Pl = PauliString.Build(letters);
                acc += Pl.KroneckerProduct(Pl.Transpose());
            }
            return acc;
        }
    }

    // Ad_{U^⊗2} on coherence space: ρ ↦ UρU†, vec(UρU†) = (U ⊗ U*)·vec(ρ).
    private static ComplexMatrix AdUnitary(ComplexMatrix u1)
    {
        var u = u1.KroneckerProduct(u1);
        return u.KroneckerProduct(u.Conjugate());
    }

    private static ComplexMatrix Hadamard1()
    {
        double s = 1.0 / Math.Sqrt(2.0);
        return Matrix<Complex>.Build.DenseOfArray(new Complex[,] { { s, s }, { s, -s } });
    }

    private static ComplexMatrix RxHalfPi1()
    {
        // R_x(π/2) = cos(π/4)·I − i·sin(π/4)·X (maps Y→Z, Z→−Y; the sign squares out in Q_P).
        double c = Math.Cos(Math.PI / 4.0), s = Math.Sin(Math.PI / 4.0);
        return Matrix<Complex>.Build.DenseOfArray(new Complex[,]
            { { c, new Complex(0, -s) }, { new Complex(0, -s), c } });
    }

    private static List<ComplexMatrix> GroupClosure(IReadOnlyList<ComplexMatrix> gens, int dim)
    {
        var elems = new List<ComplexMatrix> { Matrix<Complex>.Build.DenseIdentity(dim) };
        bool changed = true;
        while (changed)
        {
            changed = false;
            foreach (var g in gens)
                foreach (var e in elems.ToList())
                {
                    var cand = g * e;
                    if (!elems.Any(x => MaxAbsDiff(cand, x) <= Tol)) { elems.Add(cand); changed = true; }
                }
        }
        return elems;
    }

    private static bool SpectraEqual(ComplexMatrix a, ComplexMatrix b)
    {
        var ea = a.Evd().EigenValues.Enumerate().Select(z => z.Real).OrderBy(x => x).ToArray();
        var eb = b.Evd().EigenValues.Enumerate().Select(z => z.Real).OrderBy(x => x).ToArray();
        return ea.Length == eb.Length && ea.Zip(eb, (x, y) => Math.Abs(x - y)).All(v => v <= 1e-9);
    }

    private static BatteryCase DevCase(string name, string detail, double dev) =>
        new(name, detail, "dev ≤ 1e-10",
            dev <= Tol ? "dev ≤ 1e-10" : "dev = " + dev.ToString("E2", CultureInfo.InvariantCulture));

    private static double MaxAbsDiff(ComplexMatrix a, ComplexMatrix b)
    {
        double m = 0.0;
        for (int i = 0; i < a.RowCount; i++)
            for (int j = 0; j < a.ColumnCount; j++)
            {
                double v = (a[i, j] - b[i, j]).Magnitude;
                if (v > m) m = v;
            }
        return m;
    }
}
