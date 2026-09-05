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
/// THREE, permuted by the letter moves, while a second three-fold (the mirror D₄) acts inside a single
/// diagonal. The two do NOT assemble into a semidirect product (corrected 2026-08-01, see below).
///
/// <para><b>The physical dephasing diagonal</b> in light P is Q_P = Σ_l kron(P_l, P_lᵀ) (ρ ↦ P_l ρ P_l;
/// the Y-transpose matters: Yᵀ = −Y so Q_Y = −Σ kron(Y,Y)). For Z this is the §4.7 diagonal
/// Q_Z = Σ_l Z_l⊗Z_l, L_D = γ·(Q_Z − N·I) (<see cref="AbsorptionTheoremClaim"/>).</para>
///
/// <para><b>The orbit (the three diagonals):</b> {Q_X, Q_Y, Q_Z} is exactly ONE orbit of the single-qubit
/// letter moves (h_zx = Ad_{H^⊗N}: Z↔X; h_yz = Ad_{R_x(π/2)^⊗N}: Z↔Y), hence the three are conjugate and
/// share a spectrum (verified directly, N=2..4; orbit size exactly 3 at N=2,3). So "the one diagonal" is
/// provably one face of a three-fold. <b>Naming, corrected 2026-08-01:</b> ⟨h_zx, h_yz⟩ is NOT S₃ — R_x(π/2)
/// is a quarter-turn on the letters, so h_yz has order 4 and the closure has order 24 (the single-qubit
/// Clifford group mod phase). The genuine letter-S₃ has order 6 and uses the INVOLUTIVE transposition
/// t_yz = Ad of (Y+Z)/√2. Both give the same orbit of 3: the order-6 group acts faithfully, the order-24
/// one through a kernel of order 4, because Q_P is quadratic in P and so cannot see a sign-only move.</para>
///
/// <para><b>The three readings (the mirror group, within one diagonal):</b> the mirror group D₄ = ⟨R, D⟩
/// (<see cref="MirrorGroupD4Claim"/>) acts on a single Q: D (the transpose) FIXES it (D·Q·D = +Q) = the
/// rate / absorption reading (Re λ = −2γ⟨n_XY⟩); R (the ket-flip) REFLECTS it (R·Q·R = −Q), carrying the
/// entire −2Σγ shift = the mirror / palindrome reading (complementary light); the F87 truly cell
/// (n_Y even ∧ n_Z even) is the joint-fixed cell of {D, 𝓕D} = the judge reading. Q is the unique
/// D-invariant, R-anti-invariant dephasing diagonal.</para>
///
/// <para><b>The structure is NOT S₃ ⋉ D₄</b> (corrected 2026-08-01; the order-48 shape
/// PROOF_PI_FACTORS_AS_R_TIMES_D §5 expected was disproved in that same §5 on 2026-06-15, and this claim
/// did not follow until review caught it). The letter three-fold (the three diagonals) and the mirror-D₄
/// (the three readings) are TWO distinct three-fold structures, and each letter move commutes with one
/// mirror generator but not the other — [h_zx, D] = 0 but [h_zx, R] ≠ 0; [h_yz, R] = 0 but [h_yz, D] ≠ 0.
/// That pattern is NOT evidence of a semidirect product; relations of that shape hold in any group
/// containing both factors. A semidirect product needs the letter moves to NORMALIZE D₄, and they do not:
/// h_zx·R·h_zx⁻¹ is the one-sided multiplication by Z^⊗N, outside the eight elements of ⟨R, D⟩, and the
/// coherence-space closure ⟨R, D, h_zx, t_yz⟩ has order 96·2^N (384 at N=2, 768 at N=3), not 48.
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

    /// <summary>Typed parent: the mirror group ⟨R, D⟩ ≅ D₄ — the three readings act within a diagonal.
    /// It is the second three-fold, not the second factor of a product: the letter moves do not
    /// normalize it.</summary>
    public MirrorGroupD4Claim MirrorGroup { get; }

    /// <summary>Typed parent: the absorption diagonal L_D = γ·(Q − N·I), Q = Σ_l Z_l⊗Z_l (§4.7). With this
    /// edge the mirror group links directly to the dephasing diagonal — the weld.</summary>
    public AbsorptionTheoremClaim Diagonal { get; }

    public IReadOnlyList<BatteryCase> Cases { get; }
    public int PassCount => Cases.Count(c => c.Passes);

    public ThreeDephasingDiagonalsOrbitClaim(MirrorGroupD4Claim mirrorGroup, AbsorptionTheoremClaim diagonal)
        : base("The three dephasing diagonals as one orbit: the one diagonal Q_Z = Σ_l Z_l⊗Z_l is one of " +
               "three (Q_P = Σ_l kron(P_l, P_lᵀ), one per axis; Q_Y carries Yᵀ = −Y), and {Q_X, Q_Y, Q_Z} is " +
               "ONE ORBIT of the single-qubit letter moves ⟨h_zx (Z↔X), h_yz (Z↔Y)⟩ — conjugate, " +
               "same spectrum. That group has order 24, not 6; the genuine letter-S₃ (order 6) uses the " +
               "involutive t_yz = Ad of (Y+Z)/√2, and both induce the same orbit of three. The " +
               "one-diagonal's three readings are the mirror group D₄ = ⟨R, D⟩ acting " +
               "WITHIN a diagonal: D fixes Q (rate/absorption), R reflects it R·Q·R = −Q carrying −2Σγ " +
               "(mirror/palindrome), the {D, 𝓕D} joint-fixed cell is truly (judge). The two three-folds do " +
               "NOT form a semidirect product S₃ ⋉ D₄: the letter moves do not normalize D₄ " +
               "(h_zx·R·h_zx⁻¹ = one-sided Z^⊗N, outside ⟨R,D⟩; closure order 96·2^N, not 48). D does NOT permute the " +
               "diagonals (it fixes them); the proof's 'D = Z↔Y swap' is on the palindromizer Π, not on Q. " +
               "The claim's two parents are the physics edge welding the mirror-group and absorption clusters.",
               Tier.Tier1Derived,
               "simulations/one_diagonal_mirror_group.py (self-validating Stages 0-2 + an N=3 attack) + " +
               "simulations/mirror_inventory_d4.py (block D 63/63 truly cell) + " +
               "docs/proofs/PROOF_ABSORPTION_THEOREM.md §4.7 + docs/proofs/PROOF_PI_FACTORS_AS_R_TIMES_D.md §5 + " +
               "compute/RCPsiSquared.Diagnostics/Foundation/DiagonalWitness.cs (DiagonalWitness, inspect --root diagonal) + " +
               "docs/THE_THREE_DIAGONALS.md (synthesis: the one diagonal as one of three, the basis-S₃ orbit)")
    {
        MirrorGroup = mirrorGroup ?? throw new ArgumentNullException(nameof(mirrorGroup));
        Diagonal = diagonal ?? throw new ArgumentNullException(nameof(diagonal));
        Cases = BuildBattery();
    }

    public override string DisplayName =>
        "Three dephasing diagonals as one letter orbit + the three readings (mirror-D₄); the two three-folds do not form S₃ ⋉ D₄";

    public override string Summary =>
        "the one diagonal is one of three (Q_X, Q_Y, Q_Z); they are one orbit of the single-qubit letter " +
        "moves ⟨h_zx, h_yz⟩ (order 24; the genuine letter-S₃ is the order-6 ⟨h_zx, t_yz⟩), hence same " +
        "spectrum, and the three readings (rate = D-fix, mirror = R·Q·R=−Q, judge = {D,𝓕D} cell) are the " +
        "mirror group D₄ acting within a diagonal; the two three-folds do NOT form S₃ ⋉ D₄ (the letter " +
        $"moves do not normalize D₄) and the dual parentage welds the mirror-group and absorption clusters; {PassCount}/{Cases.Count} PASS ({Tier.Label()})";

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
            yield return new InspectableNode("the structure: two three-folds, NOT a semidirect product",
                summary: "the letter three-fold (three diagonals) and the mirror-D₄ (three readings) are two " +
                         "distinct structures. Each letter move commutes with one mirror generator and not the " +
                         "other ([h_zx,D]=0 but [h_zx,R]≠0; [h_yz,R]=0 but [h_yz,D]≠0), which is NOT a semidirect " +
                         "product: h_zx·R·h_zx⁻¹ = one-sided Z^⊗N leaves ⟨R,D⟩, so the letter moves do not " +
                         "normalize D₄, and the closure has order 96·2^N, not the 48 of S₃ ⋉ D₄. The order-48 " +
                         "shape PROOF_PI_FACTORS §5 once named open was disproved there on 2026-06-15.");
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
                "the Hadamard move sends the Z-diagonal to the X-diagonal, to 1e-10 on the float route (the conjugation itself is dyadic)",
                MaxAbsDiff(hZX * QZ * hZX.ConjugateTranspose(), QX)),
            DevCase("orbit: h_yz·Q_Z·h_yz⁻¹ = Q_Y (the Z↔Y basis move)",
                "the R_x(π/2) move sends the Z-diagonal to the Y-diagonal (Q_Y carries Yᵀ=−Y), to 1e-10 on the float route",
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
        cases.Add(new BatteryCase("orbit(Q_Z) = {Q_X, Q_Y, Q_Z} under the letter moves",
            "⟨h_zx, h_yz⟩ carries Q_Z to exactly the three dephasing diagonals",
            "3 = {X,Y,Z}", isThree ? "3 = {X,Y,Z}" : $"{orbit.Count}"));

        // the group NAME: ⟨h_zx, h_yz⟩ has order 24 (R_x(π/2) is a quarter-turn), NOT 6. The genuine
        // letter-S₃ is ⟨h_zx, t_yz⟩ with the involutive t_yz, order 6, and it induces the same orbit.
        var tYZ = AdUnitary(TransposeYZ1());
        var letterS3 = GroupClosure(new[] { hZX, tYZ }, d2);
        cases.Add(new BatteryCase("group order: ⟨h_zx, h_yz⟩ = 24, ⟨h_zx, t_yz⟩ = 6",
            "R_x(π/2) is order 4, so ⟨h_zx, h_yz⟩ is the Clifford group mod phase; only the involutive "
            + "t_yz = (Y+Z)/√2 generates the order-6 letter-S₃",
            "24 and 6", $"{basis.Count} and {letterS3.Count}"));

        // the structure is NOT a semidirect product: a semidirect product would need the letter moves to
        // NORMALIZE D₄, and h_zx·R·h_zx⁻¹ = the one-sided multiplication by Z^⊗N, which leaves ⟨R, D⟩.
        // (Until 2026-08-01 this battery asserted "semidirect" from the commutator pattern below, which
        // does not imply it; the commutators are kept as their own case.)
        var d4 = GroupClosure(new[] { R, D }, d2);
        var conjR = hZX * R * hZX.ConjugateTranspose();
        bool normalizes = d4.Any(g => MaxAbsDiff(conjR, g) <= Tol);
        cases.Add(new BatteryCase("structure: the letter moves do NOT normalize D₄ (so no S₃ ⋉ D₄)",
            "h_zx·R·h_zx⁻¹ must leave the 8-element ⟨R, D⟩ for the semidirect product to fail",
            "outside ⟨R,D⟩", normalizes ? "INSIDE ⟨R,D⟩" : "outside ⟨R,D⟩"));

        bool commutatorPattern = MaxAbsDiff(hZX * D, D * hZX) <= Tol && MaxAbsDiff(hZX * R, R * hZX) > 0.1
                       && MaxAbsDiff(hYZ * R, R * hYZ) <= Tol && MaxAbsDiff(hYZ * D, D * hYZ) > 0.1;
        cases.Add(new BatteryCase("commutators: each letter move commutes with exactly one mirror generator",
            "[h_zx,D]=0 & [h_zx,R]≠0 & [h_yz,R]=0 & [h_yz,D]≠0 (a fact about the factors, not a product)",
            "as expected", commutatorPattern ? "as expected" : "DIFFERENT"));

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

    /// <summary>The INVOLUTIVE Y↔Z transposition Clifford, (Y + Z)/√2. Order 2, unlike R_x(π/2), so
    /// ⟨h_zx, t_yz⟩ is the genuine order-6 letter-S₃ rather than the order-24 Clifford group.</summary>
    private static ComplexMatrix TransposeYZ1()
    {
        double s = 1.0 / Math.Sqrt(2.0);
        // (Y + Z)/√2 with Y = [[0,-i],[i,0]], Z = [[1,0],[0,-1]].
        return Matrix<Complex>.Build.DenseOfArray(new Complex[,]
            { { s, new Complex(0, -s) }, { new Complex(0, s), -s } });
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
