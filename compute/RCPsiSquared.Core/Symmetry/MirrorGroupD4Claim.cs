using System.Globalization;
using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Pauli;
using ComplexMatrix = MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex>;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>The mirror group D₄ (Tier1Derived, 2026-06-10): the canonical palindromizer
/// factors as Π_Z = R·D, and the repository's mirror inventory closes into one dihedral
/// group of order 8.
///
/// <para><b>The two generators</b>, on the coherence space of an N-qubit chain in the
/// row-stacking (C-order) vec convention (|i⟩⟨j| ↦ e_i ⊗ e_j, kron(A, B): ρ ↦ A·ρ·Bᵀ),
/// with F = X^⊗N (F² = I, Fᵀ = F):</para>
/// <list type="bullet">
///   <item><b>R</b> = I ⊗ F, the ket reflection R(ρ) = ρ·F: the windowed-converse spine's
///         one-sided involution (PROOF_F87_WINDOWED_MONOMIAL_CONVERSE §2). It flips the
///         ket index of every coherence, j ↦ j̄ (all bits flipped).</item>
///   <item><b>D</b>, the transpose superoperator D(ρ) = ρᵀ (the SWAP permutation on
///         coherence space). On the Pauli basis D = diag((−1)^{n_Y}): F114's / Welle 12's
///         diagonal mirror (<see cref="Pi2KleinV4DephaseSwapGroup.BuildD"/>), because
///         transposition flips the sign of each Y letter and nothing else.</item>
/// </list>
///
/// <para><b>Theorem (factorization):</b> Π_Z = R ∘ D, D applied first: Π_Z(ρ) = ρᵀ·X^⊗N.
/// Per site σ ↦ σᵀ·X reproduces the April rule I → X, X → I, Y → iZ, Z → iY with no extra
/// phase (the hard-won factors i fall out of Yᵀ = −Y meeting YX = −iZ). The opposite order
/// is the inverse: Π_Y = D ∘ R = Π_Z⁻¹.</para>
///
/// <para><b>Theorem (the group):</b> ⟨R, D⟩ ≅ D₄, order 8. Rotations
/// {I, Π_Z, 𝓕 = Π_Z², Π_Y = Π_Z³} with 𝓕 = F⊗F the charge conjugation (F1², the center);
/// reflections {D, 𝓕D = diag((−1)^{n_Z}), R, 𝓕R = F⊗I}. The two reflection classes have
/// Pauli-basis meaning: the square's diagonal mirrors {D, 𝓕D} are literally the diagonal
/// matrices (pure signs), the edge mirrors {R, 𝓕R} are the one-sided F multiplications
/// that move strings. The windowed-converse spine V₄ {I, 𝓕, R, 𝓕R} is the Klein subgroup
/// of one-sided and two-sided F multiplications; the diagonal Klein subgroup {I, 𝓕, D, 𝓕D}
/// is the sign-grading square; the rotation subgroup is the palindromizer family. All three
/// intersect in the center {I, 𝓕}. Welle 12's D·Π_Z·D = Π_Y is the dihedral inversion
/// relation s·r·s = r⁻¹ in disguise.</para>
///
/// <para><b>The palindrome splits along the generators:</b> D·L_H·D = −L_H (transposition
/// reverses commutators; n_Y-even H, F114's ε = −1) and D fixes the Z-dephasing dissipator
/// (the rate depends only on i ⊕ j); R·L_H·R = +L_H (XXZ terms are F-conjugation-even) and
/// R reflects the dissipator, R·L_diss·R = −L_diss − 2Σγ·I, carrying the entire shift
/// (flipping j complements the set of lit sites). The product is the April palindrome
/// Π·L·Π⁻¹ = −L − 2Σγ·I: April multiplied the generators, June kept them separate.</para>
///
/// <para><b>The cube filled (§7):</b> the polarity cube's three axes
/// (<see cref="KleinEightCellClaim"/>) are characters of three specific mirrors of the
/// Pauli algebra: bit_a = char(Ad_{Z^⊗N}), bit_b = char(Ad_{X^⊗N}), y_par = char(θ) with
/// θ the transpose. Unitary conjugations can never read a single-letter parity (conjugating
/// by a Pauli string flips exactly the two letters that anticommute with it), so the
/// conjugation mirrors span only the even Klein square; <b>y_par is the antiautomorphism
/// axis</b>: θ alone reads (−1)^{n_Y} (Y is the only antisymmetric Pauli), and its
/// composites deliver (−1)^{n_X} and (−1)^{n_Z} (the latter = 𝓕D on coherence space).
/// The quadratic-to-cubic step of the polarity cube and the Π = R·D factorization are one
/// move seen twice. Corollary: the F85/F87 truly criterion (n_Y even AND n_Z even) is the
/// joint-fixed cell of the diagonal Klein subgroup.</para>
///
/// <para><b>Deliberately outside</b> (PROOF_PI_FACTORS_AS_R_TIMES_D §5): the sublattice
/// chirality K₁ (<see cref="ChiralKClaim"/>, grades by site, not letters), the F116 golden
/// router W (two-sided P ≠ Q, non-involutive; it covers the n_Z-odd ceiling scope that D₄
/// provably cannot), the crossover mirror M (continuous R_z(π/4) object), F71's spatial
/// bond mirror, and the dephase-letter swaps Q_zx / Q_yx of
/// <see cref="Pi2KleinV4DephaseSwapGroup"/> (the Z↔Y transposition is D itself and is
/// inside; the other two need the X↔Z basis permutation and are outside; adjoining them is
/// the named-open S₃ ⋉ D₄ completion).</para>
///
/// <para><b>Layer note:</b> the windowed-converse spine claim
/// (WindowedConverseThresholdClaim) lives in RCPsiSquared.Diagnostics, which Core cannot
/// reference; the spine edge is carried in prose and through the
/// PROOF_F87_WINDOWED_MONOMIAL_CONVERSE §2 anchor, not as a typed parent.</para>
///
/// <para><b>Self-check battery:</b> built in the constructor at N = 2 (16×16 dense
/// superoperators on coherence space, &lt; 100 ms): Π_Z = R·D, Π_Z² = 𝓕, Π_Z⁴ = I,
/// D·Π_Z·D = Π_Z⁻¹ = Π_Y, closure |⟨R, D⟩| = 8, spine V₄ membership, the diagonal-mirror
/// identifications D = diag((−1)^{n_Y}) and 𝓕D = diag((−1)^{n_Z}), and the three cube
/// characters on all 16 Pauli strings. Mirrors blocks A/B/G of
/// <c>simulations/mirror_inventory_d4.py</c> (which additionally verifies N = 1..3, the
/// palindrome factorization rows on XXZ, the truly-cell equivalence 63/63, and an F114
/// N = 5 spot check, all dev 0.00e+00).</para></summary>
public sealed class MirrorGroupD4Claim : Claim
{
    private const double Tol = 1e-12;
    private const int BatteryN = 2;

    /// <summary>One self-check tying the claim to the D₄ identities at N = 2.</summary>
    public readonly record struct BatteryCase(string Name, string Detail, string Expected, string Actual)
    {
        public bool Passes => string.Equals(Expected, Actual, StringComparison.Ordinal);
    }

    /// <summary>Typed parent: the polarity cube (Z₂³ 8-cell decomposition) whose three
    /// axes this claim characterizes as characters of Ad_{Z^⊗N}, Ad_{X^⊗N}, and θ.</summary>
    public KleinEightCellClaim Cube { get; }

    /// <summary>Typed parent: F114, D's sign law D·L_σ·D = ε(σ)·L_σ with
    /// ε(σ) = (−1)^{n_Y(σ)+1}. Once D is recognized as the transpose superoperator the law
    /// is transparent (one sign from commutator reversal, one from σᵀ = (−1)^{n_Y}·σ);
    /// it is D's row of the palindrome split.</summary>
    public CommutatorDConjugationSign F114 { get; }

    /// <summary>Typed parent: Welle 12's Klein-V₄ dephase-swap group, owner of the
    /// diagonal D (<see cref="Pi2KleinV4DephaseSwapGroup.BuildD"/>). Its flagship identity
    /// D·Π_Z·D = Π_Y is re-read here as the dihedral inversion relation s·r·s = r⁻¹.</summary>
    public Pi2KleinV4DephaseSwapGroup KleinV4 { get; }

    public IReadOnlyList<BatteryCase> Cases { get; }
    public int PassCount => Cases.Count(c => c.Passes);

    public MirrorGroupD4Claim(
        KleinEightCellClaim cube,
        CommutatorDConjugationSign f114,
        Pi2KleinV4DephaseSwapGroup kleinV4)
        : base("Mirror group D₄: Π_Z = R·D (Π_Z(ρ) = ρᵀ·X^⊗N, D the transpose first, R = I⊗F the " +
               "ket reflection second; the opposite order is Π_Y = Π_Z⁻¹); ⟨R, D⟩ ≅ D₄ of order 8 " +
               "(rotations {I, Π_Z, 𝓕, Π_Y}; reflections {D, 𝓕D, R, 𝓕R}) with the windowed-converse " +
               "spine V₄ {I, 𝓕, R, 𝓕R} as Klein subgroup; the palindrome splits along the generators " +
               "(D flips L_H, R reflects the dissipator and carries −2Σγ); the polarity cube axes are " +
               "characters: bit_a = char(Ad_{Z^⊗N}), bit_b = char(Ad_{X^⊗N}), y_par = char(θ), the " +
               "antiautomorphism axis. Tier1Derived (signed-permutation identities, exact)",
               Tier.Tier1Derived,
               "docs/proofs/PROOF_PI_FACTORS_AS_R_TIMES_D.md + " +
               "simulations/mirror_inventory_d4.py + " +
               "docs/proofs/PROOF_F87_WINDOWED_MONOMIAL_CONVERSE.md (§2 spine involutions 𝓕 and R) + " +
               "docs/ANALYTICAL_FORMULAS.md F114 + " +
               "compute/RCPsiSquared.Core/Symmetry/PiOperator.cs")
    {
        Cube = cube ?? throw new ArgumentNullException(nameof(cube));
        F114 = f114 ?? throw new ArgumentNullException(nameof(f114));
        KleinV4 = kleinV4 ?? throw new ArgumentNullException(nameof(kleinV4));
        Cases = BuildBattery();
    }

    /// <summary>The factorization theorem in one line.</summary>
    public string Factorization =>
        "Π_Z = R ∘ D, D first: Π_Z(ρ) = ρᵀ·X^⊗N. Per site σ ↦ σᵀ·X reproduces the April rule " +
        "I → X, X → I, Y → iZ, Z → iY with no extra phase; the opposite order is Π_Y = D ∘ R = Π_Z⁻¹.";

    /// <summary>The group theorem in one line.</summary>
    public string GroupClosure =>
        "⟨R, D⟩ ≅ D₄ (order 8): rotations {I, Π_Z, 𝓕 = Π_Z² = F⊗F, Π_Y = Π_Z³}; reflections " +
        "{D = diag((−1)^{n_Y}), 𝓕D = diag((−1)^{n_Z}), R = I⊗F, 𝓕R = F⊗I}. Defining relations " +
        "R² = D² = I, Π_Z⁴ = I, D·Π_Z·D = Π_Z⁻¹ (= Welle 12's D·Π_Z·D = Π_Y).";

    /// <summary>The palindrome split in one line.</summary>
    public string PalindromeSplit =>
        "D·L_H·D = −L_H and D·L_diss·D = +L_diss; R·L_H·R = +L_H and R·L_diss·R = −L_diss − 2Σγ·I. " +
        "Each generator does one job; the product is the April palindrome Π·L·Π⁻¹ = −L − 2Σγ·I.";

    /// <summary>The cube-character theorem in one line.</summary>
    public string CubeCharacters =>
        "bit_a = char(Ad_{Z^⊗N}), bit_b = char(Ad_{X^⊗N}), y_par = char(θ). Conjugations span only " +
        "the even Klein square (they flip two letters at once); the transpose θ is the antiautomorphism " +
        "that reads the single-letter parity (−1)^{n_Y} and lifts the square to the cube.";

    public override string DisplayName =>
        "Mirror group D₄ (Π_Z = R·D; ⟨R, D⟩ ≅ D₄; spine V₄ ⊂ D₄; cube axes = characters)";

    public override string Summary =>
        "the canonical palindromizer factors as Π_Z = R·D (transpose, then ket reflection by X^⊗N) and " +
        "⟨R, D⟩ closes into one dihedral D₄ holding every named mirror of the palindrome story: the April " +
        "Π_Z/Π_Y rotations, the spine V₄ {I, 𝓕, R, 𝓕R} as Klein subgroup, F114's D and the fourth mirror " +
        "𝓕D = diag((−1)^{n_Z}). D flips L_H, R carries the −2Σγ shift; the polarity cube axes are the " +
        "characters bit_a/bit_b/y_par of Ad_{Z^⊗N}/Ad_{X^⊗N}/θ, with y_par the antiautomorphism axis no " +
        $"unitary conjugation can read; {PassCount}/{Cases.Count} battery PASS ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("Factorization Π_Z = R·D", summary: Factorization);
            yield return new InspectableNode("Group closure ⟨R, D⟩ ≅ D₄", summary: GroupClosure);
            yield return new InspectableNode("Palindrome split along the generators", summary: PalindromeSplit);
            yield return new InspectableNode("The cube filled (three axes = three characters)", summary: CubeCharacters);
            yield return new InspectableNode("Welle 12 re-read",
                summary: "D·Π_Z·D = Π_Y (PROOF_D_PI_Z_EQUALS_PI_Y_UNIVERSAL_N) is the dihedral inversion " +
                         "relation s·r·s = r⁻¹ with r = Π_Z, s = D: the Z↔Y dephase-letter swap and " +
                         "running the palindromizer backwards are the same operation.");
            yield return new InspectableNode("Truly criterion as character cell",
                summary: "The diagonal Klein subgroup {I, 𝓕, D, 𝓕D} is jointly diagonal on Pauli strings " +
                         "with eigenvalues ((−1)^{n_Y}, (−1)^{n_Z}); the F85/F87 truly criterion (n_Y even " +
                         "AND n_Z even) is exactly the joint-fixed cell (verified 63/63 at N = 3 in the " +
                         "committed verifier).");
            yield return new InspectableNode("Deliberately outside",
                summary: "K₁ sublattice chirality (grades by site), the F116 golden router W (two-sided, " +
                         "non-involutive; covers the n_Z-odd ceiling scope D₄ provably cannot enter), the " +
                         "crossover mirror M (continuous R_z(π/4)), F71's spatial bond mirror, and the " +
                         "Q_zx/Q_yx dephase swaps (S₃ ⋉ D₄ completion named open in the proof §5).");
            yield return new InspectableNode("Spine claim (prose edge only)",
                summary: "The windowed-converse spine's typed claim (WindowedConverseThresholdClaim) lives " +
                         "in RCPsiSquared.Diagnostics; Core cannot reference it, so the edge is carried by " +
                         "the PROOF_F87_WINDOWED_MONOMIAL_CONVERSE §2 anchor instead of a ctor parent.");
            yield return new InspectableNode("Typed parents",
                summary: $"KleinEightCellClaim ({Cube.Tier.Label()}): the cube whose axes become characters; " +
                         $"CommutatorDConjugationSign ({F114.Tier.Label()}): F114, D's sign law = D's row of " +
                         $"the palindrome split; Pi2KleinV4DephaseSwapGroup ({KleinV4.Tier.Label()}): owner " +
                         "of D, its D·Π_Z·D = Π_Y re-read as the dihedral inversion.");
            foreach (var c in Cases)
                yield return new InspectableNode(c.Name,
                    summary: $"{c.Detail}; expected {c.Expected}, got {c.Actual}, " + (c.Passes ? "PASS" : "FAIL"));
        }
    }

    // ------------------------------------------------------------------
    // Self-check battery: N = 2 dense superoperators on coherence space.
    // ------------------------------------------------------------------

    private static IReadOnlyList<BatteryCase> BuildBattery()
    {
        int d = 1 << BatteryN;     // 4
        int d2 = d * d;            // 16 coherence dimensions = 4^N Pauli strings

        var idH = Matrix<Complex>.Build.DenseIdentity(d);
        var xx = new[] { PauliLetter.X, PauliLetter.X };
        var zz = new[] { PauliLetter.Z, PauliLetter.Z };
        var F = PauliString.Build(xx);          // X^⊗2
        var ZN = PauliString.Build(zz);         // Z^⊗2

        // Row-stacking vec convention: |i⟩⟨j| ↦ e_i ⊗ e_j; kron(A, B): ρ ↦ A·ρ·Bᵀ.
        var D = Matrix<Complex>.Build.Dense(d2, d2);
        for (int i = 0; i < d; i++)
            for (int j = 0; j < d; j++)
                D[j * d + i, i * d + j] = Complex.One;     // vec(ρᵀ) = D·vec(ρ)
        var R = idH.KroneckerProduct(F);                    // ρ ↦ ρ·F (Fᵀ = F)
        var FF = F.KroneckerProduct(F);                     // 𝓕: ρ ↦ F·ρ·F
        var FR = F.KroneckerProduct(idH);                   // 𝓕R: ρ ↦ F·ρ

        // Pauli-string basis → coherence basis: M[:, k] = vec(σ_k); M⁻¹ = M†/2^N.
        var M = Matrix<Complex>.Build.Dense(d2, d2);
        for (int k = 0; k < d2; k++)
        {
            var sigma = PauliString.Build(PauliIndex.FromFlat(k, BatteryN));
            for (int i = 0; i < d; i++)
                for (int j = 0; j < d; j++)
                    M[i * d + j, k] = sigma[i, j];
        }
        var mInv = M.ConjugateTranspose().Divide(d);

        var pi = M * PiOperator.BuildFull(BatteryN, PauliLetter.Z) * mInv;
        var piY = M * PiOperator.BuildFull(BatteryN, PauliLetter.Y) * mInv;
        var piInv = pi * pi * pi;
        var idS = Matrix<Complex>.Build.DenseIdentity(d2);

        var cases = new List<BatteryCase>
        {
            DevCase("Π_Z = R·D (factorization)",
                "max|Π_coh − R·D| on the 16×16 coherence space; Π_coh = M·PiOperator.BuildFull(2, Z)·M⁻¹, D applied first",
                MaxAbsDiff(pi, R * D)),
            DevCase("Π_Z² = 𝓕 = F⊗F (charge conjugation, F1²)",
                "max|Π² − F⊗F|; the center of D₄",
                MaxAbsDiff(pi * pi, FF)),
            DevCase("Π_Z⁴ = I (order 4)",
                "max|Π·Π³ − I|; R, D involutions with order-4 product is the dihedral presentation",
                MaxAbsDiff(pi * piInv, idS)),
            DevCase("D·Π_Z·D = Π_Z⁻¹ = Π_Y (dihedral inversion = Welle 12)",
                "max over |D·Π·D − Π³| and |D·Π·D − Π_Y(PiOperator)|",
                Math.Max(MaxAbsDiff(D * pi * D, piInv), MaxAbsDiff(D * pi * D, piY))),
        };

        // Group closure from the two generators.
        var elems = new List<ComplexMatrix> { idS };
        bool changed = true;
        while (changed)
        {
            changed = false;
            foreach (var g in new[] { R, D })
            {
                foreach (var e in elems.ToList())
                {
                    var cand = g * e;
                    if (!elems.Any(x => MaxAbsDiff(cand, x) <= Tol))
                    {
                        elems.Add(cand);
                        changed = true;
                    }
                }
            }
        }
        cases.Add(new BatteryCase(
            Name: "|⟨R, D⟩| = 8 (D₄ closure)",
            Detail: "closure of {R, D} under multiplication on the 16×16 coherence space",
            Expected: "8",
            Actual: elems.Count.ToString(CultureInfo.InvariantCulture)));

        var spine = new[] { idS, FF, R, FR };
        int spineHits = spine.Count(s => elems.Any(x => MaxAbsDiff(s, x) <= Tol));
        cases.Add(new BatteryCase(
            Name: "spine V₄ {I, 𝓕, R, 𝓕R} ⊂ ⟨R, D⟩ (Klein subgroup)",
            Detail: "membership of the windowed-converse spine involution set in the closure",
            Expected: "4/4",
            Actual: spineHits.ToString(CultureInfo.InvariantCulture) + "/4"));

        // Diagonal mirrors in the Pauli basis: D = diag((−1)^{n_Y}), 𝓕D = diag((−1)^{n_Z}).
        var dPauli = mInv * D * M;
        var fdPauli = mInv * (FF * D) * M;
        double devDiag = 0.0;
        for (int k = 0; k < d2; k++)
        {
            var letters = PauliIndex.FromFlat(k, BatteryN);
            int nY = letters.Count(l => l == PauliLetter.Y);
            int nZ = letters.Count(l => l == PauliLetter.Z);
            for (int r = 0; r < d2; r++)
            {
                Complex expectD = r == k ? SignOf(nY) : Complex.Zero;
                Complex expectFd = r == k ? SignOf(nZ) : Complex.Zero;
                devDiag = Math.Max(devDiag, (dPauli[r, k] - expectD).Magnitude);
                devDiag = Math.Max(devDiag, (fdPauli[r, k] - expectFd).Magnitude);
            }
        }
        cases.Add(DevCase("D = diag((−1)^{n_Y}), 𝓕D = diag((−1)^{n_Z}) in the Pauli basis",
            "the two diagonal mirrors of the square are literally the diagonal sign matrices (F114's D + the fourth mirror)",
            devDiag));

        // Cube characters on all 16 strings: bit_a / bit_b / y_par.
        int okA = 0, okB = 0, okY = 0;
        for (int k = 0; k < d2; k++)
        {
            var letters = PauliIndex.FromFlat(k, BatteryN);
            var sigma = PauliString.Build(letters);
            int bitA = 0, bitB = 0, nY = 0;
            foreach (var l in letters)
            {
                bitA += l.BitA();
                bitB += l.BitB();
                if (l == PauliLetter.Y) nY++;
            }
            if (MaxAbsDiff(ZN * sigma * ZN, sigma.Multiply(SignOf(bitA))) <= Tol) okA++;
            if (MaxAbsDiff(F * sigma * F, sigma.Multiply(SignOf(bitB))) <= Tol) okB++;
            if (MaxAbsDiff(sigma.Transpose(), sigma.Multiply(SignOf(nY))) <= Tol) okY++;
        }
        cases.Add(new BatteryCase(
            Name: "bit_a = char(Ad_{Z^⊗N})",
            Detail: "Z^⊗2·σ·Z^⊗2 = (−1)^{bit_a}·σ on all 16 strings",
            Expected: "16/16",
            Actual: okA.ToString(CultureInfo.InvariantCulture) + "/16"));
        cases.Add(new BatteryCase(
            Name: "bit_b = char(Ad_{X^⊗N})",
            Detail: "X^⊗2·σ·X^⊗2 = (−1)^{bit_b}·σ on all 16 strings",
            Expected: "16/16",
            Actual: okB.ToString(CultureInfo.InvariantCulture) + "/16"));
        cases.Add(new BatteryCase(
            Name: "y_par = char(θ) (the antiautomorphism axis)",
            Detail: "σᵀ = (−1)^{n_Y}·σ on all 16 strings; the parity no unitary conjugation reads",
            Expected: "16/16",
            Actual: okY.ToString(CultureInfo.InvariantCulture) + "/16"));

        return cases;
    }

    private static Complex SignOf(int parity) => (parity & 1) == 0 ? Complex.One : -Complex.One;

    private static BatteryCase DevCase(string name, string detail, double dev) =>
        new(Name: name,
            Detail: detail,
            Expected: "dev ≤ 1e-12",
            Actual: dev <= Tol
                ? "dev ≤ 1e-12"
                : "dev = " + dev.ToString("E2", CultureInfo.InvariantCulture));

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
