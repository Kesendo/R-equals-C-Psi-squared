using System.Globalization;
using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Pauli;
using RCPsiSquared.Core.Symmetry;
using ComplexMatrix = MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex>;

namespace RCPsiSquared.Diagnostics.Foundation;

/// <summary>The live mirror group D₄ (<c>inspect --root mirrorgroup</c>): the operator-algebra side of
/// the palindrome. Where <see cref="DiagonalWitness"/> lives on the dephasing diagonal Q (the shadow)
/// and <see cref="MirrorSystem"/> (<c>--root mirror</c>) reads the palindrome off the SPECTRUM
/// (each decay rate r pairs with 2σ − r), this one recomputes, live at inspect time, the operator
/// factorization the palindrome is made of: Π_Z = R·D and the whole dihedral group it generates.
///
/// <para>The typed claim <see cref="RCPsiSquared.Core.Symmetry.MirrorGroupD4Claim"/> records this
/// structure with a fixed N = 2 self-check battery; this witness ports it into a live, N-parameterized
/// telescope and computes the one thing neither layer surfaces live: the §3 palindrome split
/// generator-by-generator, WITH the −2Σγ shift on the dissipator rows that the diagonal witness
/// recenters away (showing R·Q·R = −Q). The palindrome product Π_Z·L·Π_Z⁻¹ = −L − 2Σγ is shown as the
/// closure of R·D — the SAME palindrome <see cref="MirrorSystem"/> reads spectrally, here from the two
/// generators. Origin: <c>reflections/D_PI_Z_EQUALS_PI_Y.md</c> (D·Π_Z·D = Π_Y) +
/// <c>docs/proofs/PROOF_PI_FACTORS_AS_R_TIMES_D.md</c> (Π_Z = R·D, the D₄ closure, §3, §5) +
/// <c>docs/THE_THREE_DIAGONALS.md</c> (the one diagonal read three ways by D₄, and one of three).</para>
///
/// <para>Five live nodes: (1) the factorization Π_Z = R·D; (2) the group ⟨R, D⟩ ≅ D₄ (order 8, the
/// dihedral inversion D·Π_Z·D = Π_Y, the center 𝓕 = Π_Z²); (3) the §3 palindrome split along the
/// generators, with the shift; (4) the palindrome product −L − 2Σγ; (5) §5's deliberately-outside
/// boundary. The first two devs are the GATE: Π_coh = M·PiOperator·M⁻¹ and R·D agree only when both are
/// built in one consistent vec convention (the exact (−1)^{n_Y} twist the origin reflection is about).</para>
///
/// <para>Guard: operator dimension 4^N built only for 4^N ≤ <see cref="MaxDim"/> (= 1024), i.e. N ≤ 5.
/// Convention: row-stacking vec, kron(A, B): ρ ↦ AρBᵀ, matching the diagonal witness, the typed battery,
/// and <see cref="PiOperator"/> bridged through the Pauli↔coherence transform M.</para></summary>
public sealed class MirrorGroupWitness : IInspectable
{
    private static readonly CultureInfo Inv = CultureInfo.InvariantCulture;
    private const double Tol = 1e-9;

    /// <summary>Largest coherence-space dimension 4^N the live build will materialise.</summary>
    public const int MaxDim = 1024;

    public int N { get; }
    public double Gamma { get; }
    public double Delta { get; }
    public double J { get; }
    public int D { get; }       // 2^N
    public int Dim => D * D;    // 4^N coherence space

    // built once in the ctor (cheap at the default N=3; on-demand for larger N).
    private readonly ComplexMatrix _R, _Dsuper, _piZ, _piY, _FF, _LH, _Ldiss, _L;

    public MirrorGroupWitness(int n = 3, double gamma = 0.05, double delta = 0.7, double j = 1.0)
    {
        if (n < 2) throw new ArgumentOutOfRangeException(nameof(n), $"N must be ≥ 2; got {n}");
        if (gamma <= 0) throw new ArgumentOutOfRangeException(nameof(gamma), $"gamma must be positive; got {gamma}");
        long dim = 1L << (2 * n);   // 4^N
        if (dim > MaxDim)
            throw new ArgumentOutOfRangeException(nameof(n),
                $"4^N = {dim} exceeds the live-build guard {MaxDim} for N={n}; pick N ≤ 5.");

        N = n; Gamma = gamma; Delta = delta; J = j; D = 1 << n;

        _Dsuper = TransposeSuper(N);                  // D(ρ) = ρᵀ
        _R = KetReflection(N);                        // R(ρ) = ρ·F, F = X^⊗N
        _FF = ChargeConjugation(N);                   // 𝓕 = F⊗F: ρ ↦ F·ρ·F
        var (m, mInv) = PauliCoherenceTransform(N);   // M·σ_k = vec(σ_k), M⁻¹ = M†/2^N
        _piZ = m * PiOperator.BuildFull(N, PauliLetter.Z) * mInv;   // Π_Z on coherence space
        _piY = m * PiOperator.BuildFull(N, PauliLetter.Y) * mInv;   // Π_Y on coherence space
        _LH = HamiltonianSuper(BuildChainH(N, J, Delta), N);        // L_H = −i[H,·]
        _Ldiss = DissipatorSuper(N, Gamma);                        // L_diss = diag(−2γ·popcount(i⊕j))
        _L = _LH + _Ldiss;
    }

    public string DisplayName => $"the mirror group (N={N}, γ={Gamma.ToString("0.###", Inv)})";

    public string Summary
    {
        get
        {
            double fact = FactorizationDev();
            double inv = DihedralInversionDev();
            int order = GroupOrder();
            double prod = PalindromeProductDev();
            return $"N={N}: Π_Z=R·D (dev {fact:0.0e+00}); D·Π_Z·D=Π_Y (dev {inv:0.0e+00}); "
                 + $"|⟨R,D⟩|={order}; palindrome Π·L·Π⁻¹=−L−2Σγ (dev {prod:0.0e+00})";
        }
    }

    public IEnumerable<IInspectable> Children
    {
        get
        {
            yield return TheFactorization();
            yield return TheGroup();
            yield return ThePalindromeSplit();
            yield return ThePalindrome();
            yield return DeliberatelyOutside();
        }
    }

    public InspectablePayload Payload => InspectablePayload.Empty;

    // ============================ public live readouts (the gate + the rest) ============================
    /// <summary>The proof's headline, live: max|Π_Z − R·D|. The first half of the GATE — fires if Π_coh
    /// (from PiOperator) and the R·D product are in mismatched vec conventions.</summary>
    public double FactorizationDev() => MaxAbsDiff(_piZ, _R * _Dsuper);

    /// <summary>reflections/D_PI_Z_EQUALS_PI_Y, live: max|D·Π_Z·D − Π_Y| (the dihedral inversion
    /// s·r·s = r⁻¹, Π_Y = Π_Z⁻¹). The second half of the GATE.</summary>
    public double DihedralInversionDev() => MaxAbsDiff(_Dsuper * _piZ * _Dsuper, _piY);

    /// <summary>|⟨R, D⟩|: the closure of the two generators (= 8, the dihedral group of the square).</summary>
    public int GroupOrder() => ClosureSize(new[] { _R, _Dsuper });

    /// <summary>The §3 palindrome split along the generators, WITH the −2Σγ shift on the dissipator row
    /// (the row the diagonal witness recenters away): devs of D·L_H·D=−L_H, D·L_diss·D=+L_diss,
    /// R·L_H·R=+L_H, R·L_diss·R=−L_diss−2Σγ·I.</summary>
    public (double DLh, double DLdiss, double RLh, double RLdiss) PalindromeSplitDevs()
    {
        var shiftI = Matrix<Complex>.Build.DenseIdentity(Dim).Multiply(new Complex(2.0 * N * Gamma, 0));   // 2Σγ·I (uniform γ)
        double dLh = MaxAbsDiff(_Dsuper * _LH * _Dsuper, _LH.Multiply(-Complex.One));
        double dLdiss = MaxAbsDiff(_Dsuper * _Ldiss * _Dsuper, _Ldiss);
        double rLh = MaxAbsDiff(_R * _LH * _R, _LH);
        double rLdiss = MaxAbsDiff(_R * _Ldiss * _R, _Ldiss.Multiply(-Complex.One) - shiftI);
        return (dLh, dLdiss, rLh, rLdiss);
    }

    /// <summary>The palindrome itself, first live in C#: max|Π_Z·L·Π_Z⁻¹ − (−L − 2Σγ·I)|, L = L_H + L_diss,
    /// Π_Z⁻¹ = Π_Z³ (order 4).</summary>
    public double PalindromeProductDev()
    {
        var piInv = _piZ * _piZ * _piZ;   // Π_Z⁻¹ = Π_Z³
        var rhs = _L.Multiply(-Complex.One)
                - Matrix<Complex>.Build.DenseIdentity(Dim).Multiply(new Complex(2.0 * N * Gamma, 0));
        return MaxAbsDiff(_piZ * _L * piInv, rhs);
    }

    // ============================ nodes ============================
    private IInspectable TheFactorization()
    {
        double dev = FactorizationDev();
        return new InspectableNode("the factorization (Π_Z = R·D — the proof's headline)",
            summary: "Π_Z(ρ) = ρᵀ·X^⊗N: transpose (D) first, then ket-reflect (R). Per site σ ↦ σᵀ·X reproduces "
                   + "the April rule I→X, X→I, Y→iZ, Z→iY (the phase i falls out of Yᵀ=−Y meeting YX=−iZ); the two "
                   + $"hinges of the one operator the diagonal witness only ever sees apart. dev = {dev:0.0e+00}");
    }

    private IInspectable TheGroup()
    {
        int order = GroupOrder();
        double inv = DihedralInversionDev();
        double center = MaxAbsDiff(_piZ * _piZ, _FF);
        return new InspectableNode("the group (⟨R, D⟩ ≅ D₄, order 8)",
            summary: $"two involutions R, D whose product Π_Z has order 4 ⟹ the dihedral group of the square; "
                   + $"|⟨R,D⟩| = {order}. Rotations {{I, Π_Z, 𝓕, Π_Y}}, reflections {{D, 𝓕D, R, 𝓕R}}.",
            children: new IInspectable[]
            {
                new InspectableNode("the dihedral inversion: D·Π_Z·D = Π_Y = Π_Z⁻¹ (reflections/D_PI_Z_EQUALS_PI_Y, live)",
                    summary: "conjugating the rotation by a reflection inverts it (s·r·s = r⁻¹): the Z↔Y dephase-letter "
                           + $"swap and running the palindromizer backwards are one operation. dev = {inv:0.0e+00}"),
                new InspectableNode("𝓕 = Π_Z² (the center, charge conjugation F⊗F)",
                    summary: "the square of the palindromizer is the charge conjugation 𝓕 = X^⊗N·(·)·X^⊗N (F1²), the "
                           + $"center of D₄ that commutes with every mirror. dev = {center:0.0e+00}"),
                new InspectableNode("𝓕D = diag((−1)^{n_Z}) (the fourth mirror, §4e)",
                    summary: "the product of the center and the transpose, ρ ↦ F·ρᵀ·F; it grades by Z-content and is the "
                           + "second leg of the F87 truly cell, the mirror the repository never named until the "
                           + "factorization closed (verified in MirrorGroupD4Claim's Pauli-basis battery)."),
            });
    }

    private IInspectable ThePalindromeSplit()
    {
        var (dLh, dLdiss, rLh, rLdiss) = PalindromeSplitDevs();
        return new InspectableNode("the palindrome split (§3, along the generators — WITH the −2Σγ shift)",
            summary: "each generator does one job; their product is the palindrome. This is the full §3 column the "
                   + "diagonal witness recenters away (it shows R·Q·R = −Q, dropping the shift that is the center).",
            children: new IInspectable[]
            {
                new InspectableNode("D·L_H·D = −L_H (D carries the Hamiltonian flip)",
                    summary: "transposition reverses commutators and every XXZ term has n_Y even (σᵀ=+σ); F114's ε=−1 "
                           + $"in action. dev = {dLh:0.0e+00}"),
                new InspectableNode("D·L_diss·D = +L_diss (D fixes the dissipator)",
                    summary: $"the dephasing rate depends only on i⊕j, symmetric under i↔j. dev = {dLdiss:0.0e+00}"),
                new InspectableNode("R·L_H·R = +L_H (R fixes the Hamiltonian commutator)",
                    summary: "conjugating the ket leg by F leaves every XXZ term invariant (n_Y+n_Z even). "
                           + $"dev = {rLh:0.0e+00}"),
                new InspectableNode("R·L_diss·R = −L_diss − 2Σγ·I (R reflects it, carries the ENTIRE shift)",
                    summary: "flipping the ket index complements the lit sites: the rate −2Σ_lit γ becomes −2Σγ + 2Σ_lit γ, "
                           + $"so the constant −2Σγ appears. The palindrome's center, made visible. dev = {rLdiss:0.0e+00}"),
            });
    }

    private IInspectable ThePalindrome()
    {
        double dev = PalindromeProductDev();
        return new InspectableNode("the palindrome (the product — the closure of the factorization)",
            summary: $"Π_Z·L·Π_Z⁻¹ = −L − 2Σγ·I on the XXZ chain (Δ={Delta.ToString("0.###", Inv)}, "
                   + $"γ={Gamma.ToString("0.###", Inv)}): the four §3 rows multiplied back into one conjugation. The same "
                   + $"palindrome --root mirror reads off the spectrum (PalindromeHolds), here derived from R·D. dev = {dev:0.0e+00}");
    }

    private static IInspectable DeliberatelyOutside()
    {
        return new InspectableNode("deliberately outside (§5 — the honest boundary)",
            summary: "five mirrors kept OUT of D₄, each for a stated reason; keeping them out is part of the result.",
            children: new IInspectable[]
            {
                new InspectableNode("K₁ = Π_{l odd} Z_l (sublattice chirality)",
                    summary: "grades by SITE (even/odd sublattice), not by Pauli letter; its mirror lives on H "
                           + "(K₁HK₁ = −H), not the letter algebra."),
                new InspectableNode("the golden router W (F116)",
                    summary: "two-sided P ≠ Q golden-frame product, non-involutive per site; it covers the n_Z-odd "
                           + "ceiling scope that D₄ provably cannot enter."),
                new InspectableNode("the crossover mirror M (√90°)",
                    summary: "a continuous R_z(π/4) conjugation; D₄ samples that rotation circle only at right angles "
                           + "(the T-gate to D₄'s S-gate)."),
                new InspectableNode("F71's bond mirror",
                    summary: "a spatial reflection of the chain (site k ↔ N+1−k), not an operator-space letter map."),
                new InspectableNode("the dephase-letter swaps Q_zx, Q_yx (the S₃ ⋉ D₄ completion, named open)",
                    summary: "the Z↔Y swap IS D (inside ⟨R,D⟩); the other two need the X↔Z basis permutation and are "
                           + "outside. Adjoining them is the open completion the proof §5 named."),
            });
    }

    // ============================ builders (reused conventions) ============================
    private static ComplexMatrix TransposeSuper(int n)
    {
        int d = 1 << n, d2 = d * d;
        var t = Matrix<Complex>.Build.Dense(d2, d2);
        for (int i = 0; i < d; i++)
            for (int jj = 0; jj < d; jj++)
                t[jj * d + i, i * d + jj] = Complex.One;   // vec(ρᵀ) = D·vec(ρ)
        return t;
    }

    private static ComplexMatrix KetReflection(int n)
    {
        int d = 1 << n;
        var idH = Matrix<Complex>.Build.DenseIdentity(d);
        var letters = new PauliLetter[n];
        for (int s = 0; s < n; s++) letters[s] = PauliLetter.X;
        var f = PauliString.Build(letters);                 // X^⊗N
        return idH.KroneckerProduct(f);                     // R: ρ ↦ ρ·F
    }

    private static ComplexMatrix ChargeConjugation(int n)
    {
        var letters = new PauliLetter[n];
        for (int s = 0; s < n; s++) letters[s] = PauliLetter.X;
        var f = PauliString.Build(letters);                 // X^⊗N (Fᵀ = F)
        return f.KroneckerProduct(f);                       // 𝓕 = F⊗F: ρ ↦ F·ρ·F
    }

    // M·σ_k = vec(σ_k) (row-stacking); M⁻¹ = M†/2^N. Bridges PiOperator's Pauli basis to coherence space.
    private static (ComplexMatrix M, ComplexMatrix MInv) PauliCoherenceTransform(int n)
    {
        int d = 1 << n, d2 = d * d;
        var m = Matrix<Complex>.Build.Dense(d2, d2);
        for (int k = 0; k < d2; k++)
        {
            var sigma = PauliString.Build(PauliIndex.FromFlat(k, n));
            for (int i = 0; i < d; i++)
                for (int jj = 0; jj < d; jj++)
                    m[i * d + jj, k] = sigma[i, jj];
        }
        return (m, m.ConjugateTranspose().Divide(d));
    }

    private static ComplexMatrix BuildChainH(int n, double j, double delta)
    {
        int d = 1 << n;
        var h = Matrix<Complex>.Build.Dense(d, d);
        for (int b = 0; b < n - 1; b++)
        {
            var xx = PauliString.SiteOp(n, b, PauliLetter.X) * PauliString.SiteOp(n, b + 1, PauliLetter.X);
            var yy = PauliString.SiteOp(n, b, PauliLetter.Y) * PauliString.SiteOp(n, b + 1, PauliLetter.Y);
            var zz = PauliString.SiteOp(n, b, PauliLetter.Z) * PauliString.SiteOp(n, b + 1, PauliLetter.Z);
            h += (xx + yy).Multiply(new Complex(j, 0)) + zz.Multiply(new Complex(delta, 0));
        }
        return h;
    }

    private static ComplexMatrix HamiltonianSuper(ComplexMatrix h, int n)
    {
        int d = 1 << n;
        var idH = Matrix<Complex>.Build.DenseIdentity(d);
        var lh = h.KroneckerProduct(idH) - idH.KroneckerProduct(h.Transpose());   // [H,·]
        return lh.Multiply(new Complex(0, -1));                                     // −i[H,·]
    }

    private static ComplexMatrix DissipatorSuper(int n, double gamma)
    {
        int d = 1 << n, d2 = d * d;
        var ld = Matrix<Complex>.Build.Dense(d2, d2);
        for (int i = 0; i < d; i++)
            for (int jj = 0; jj < d; jj++)
            {
                int k = System.Numerics.BitOperations.PopCount((uint)(i ^ jj));
                ld[i * d + jj, i * d + jj] = new Complex(-2.0 * gamma * k, 0);   // Re λ = −2γ·popcount(i⊕j)
            }
        return ld;
    }

    // ============================ helpers ============================
    private static int ClosureSize(ComplexMatrix[] gens)
    {
        int dim = gens[0].RowCount;
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
        return elems.Count;
    }

    private static double MaxAbsDiff(ComplexMatrix a, ComplexMatrix b)
    {
        double m = 0.0;
        for (int i = 0; i < a.RowCount; i++)
            for (int jj = 0; jj < a.ColumnCount; jj++)
            {
                double v = (a[i, jj] - b[i, jj]).Magnitude;
                if (v > m) m = v;
            }
        return m;
    }
}
