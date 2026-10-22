using RCPsiSquared.Core.F86;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Numerics;
using RCPsiSquared.Core.Symmetry;
using Mat = RCPsiSquared.Core.Numerics.GaussianRational[,];

namespace RCPsiSquared.Diagnostics.Foundation;

/// <summary>Live witness for <see cref="TwoByTwoExceptionalPointClaim"/>: the rows of
/// <c>docs/THE_DOUBLE_ROOT.md</c>'s section "One 2 × 2 in many settings", each at a point with integer entries on its exceptional
/// point, read exactly over ℚ(i) through the claim's own law (<see cref="TwoByTwoExceptionalPointClaim.IsJordanBlock"/>):
/// on the point one 2 × 2 Jordan block, one step off it none (the control). The thirteen rows hold nine distinct
/// on-point matrices (several settings meet at the same integer block), counted live. The coalesced vector is
/// self-orthogonal in the complex-symmetric rows and not in the real opposite-sign ones; that reading checks the label
/// against the matrix's symmetry and nothing physical (a non-zero real vector has vᵀv > 0, the kernel vector (b, −a) of a
/// symmetric nilpotent [[a, b], [b, −a]] has vᵀv = a² + b² = 0). The similarity diag(1, i) carries the ring's writing into
/// the complex-symmetric one and its vector (1, −1), vᵀv = 2, into (1, i), vᵀv = 0. F155's dimer: the single-excitation
/// block is a Jordan block at g = J, so H is nilpotent there, and the no-jump generator ρ ↦ −i(Hρ − ρH†) = −i(H⊗1 − 1⊗H̄)
/// carries a Jordan block of size 3, the Kronecker sum of H's nilpotent 2 × 2 block with its conjugate (Jordan type
/// 3 ⊕ 1; the whole G is 3 ⊕ 1 ⊕ 2⁴ ⊕ 1⁴), read from exact ranks of its powers. A
/// companion matrix (the rows whose page writes a polynomial) is a Jordan block at any double root, so its on-point
/// reading cannot fail; those rows are marked, and for the N = 2 pair the physical block that carries the verdict is the
/// monotonicity row's (the N = 3 horizon's defectiveness rests on EpCharacterWitness, not on a block read here).
///
/// <para>Owners: the sodium points come from <see cref="SodiumTwoMirrorsWitness"/> (the bright centre and ω are a scalar
/// shift and change no verdict, so only 2q = 3J is read there); the F86 toy is asked of <see cref="TwoLevelEpModel"/>,
/// whose discriminant, centre and eigenvalues must equal the law's at an integer EP point and one step off; and
/// <see cref="QEpLaw"/>'s value, set as the coupling, must put the toy's block on d = 0 and a quarter step off it not.
/// The horizon's Q* = 2/√κ_N is not read from an owner here (CoherenceHorizonClaim exposes no value): its rows are the
/// companions λ² + 4γλ + κJ² at κ = 4, 2. The toy's Q_EP and the horizon's Q* are critical damping of a 2 × 2 whose
/// diagonal splitting is 4γ, the same shape; they are different modes (docs/Q_REGIME_ANCHORS.md), and no gate here
/// joins them.</para>
///
/// <para>Typed home <see cref="TwoByTwoExceptionalPointClaim"/>; the Python twin is
/// <c>simulations/two_by_two_ep_family.py</c>.</para></summary>
public sealed class TwoByTwoExceptionalPointWitness : IInspectable
{
    /// <summary>One row at its point and one step off it.</summary>
    public sealed record RowReading(string Setting, string Writing, bool Companion, bool JordanAtPoint, bool NoneOffPoint,
        bool? SelfOrthogonal)
    {
        public bool Holds => JordanAtPoint && NoneOffPoint;
    }

    public IReadOnlyList<RowReading> Rows { get; }
    /// <summary>How many distinct on-point matrices the rows hold.</summary>
    public int DistinctOnPointMatrices { get; }
    /// <summary>The dimer's single-excitation block [[ig, J], [J, −ig]] at g = J = 1 is a Jordan block.</summary>
    public bool DimerBlockIsJordan { get; }
    /// <summary>The nullities of the no-jump generator's powers G, G², G³, G⁴ at g = J = 1.</summary>
    public IReadOnlyList<int> DimerNoJumpNullities { get; }
    /// <summary>A Jordan block of size 3 and none larger in the no-jump generator.</summary>
    public bool DimerNoJumpSizeThree => DimerNoJumpNullities[2] > DimerNoJumpNullities[1] && DimerNoJumpNullities[3] == DimerNoJumpNullities[2];
    /// <summary>diag(1, i)⁻¹·[[0, 2], [−2, −4]]·diag(1, i) = [[0, 2i], [2i, −4]], and the ring's (1, −1) with vᵀv = 2
    /// becomes (1, i) with vᵀv = 0.</summary>
    public bool WritingCarriesSelfOrthogonality { get; }
    /// <summary>Every complex-symmetric row self-orthogonal, every real opposite-sign row not: a check of the label against
    /// the matrix's symmetry.</summary>
    public bool SelfOrthogonalityFollowsTheWriting =>
        Rows.Where(r => r.SelfOrthogonal is not null).All(r => r.SelfOrthogonal == (r.Writing == "complex symmetric"));
    /// <summary>TwoLevelEpModel, the F86 toy's owner, agrees exactly with the law at an integer EP point and one step off.</summary>
    public bool OwnerToyAgrees { get; }
    /// <summary>QEpLaw's value 2/g_eff, set as J at γ₀ = 1, puts the toy's block [[−2, iJg], [iJg, −6]] on d = 0 and a
    /// quarter step in J off it not, at g_eff = 1, 2, 4 (dyadic, so the double is read into ℚ exactly).</summary>
    public bool QEpLawSitsOnTheToyPoint { get; }
    public bool Holds => OwnerToyAgrees && QEpLawSitsOnTheToyPoint && Rows.All(r => r.Holds) && DimerBlockIsJordan && DimerNoJumpSizeThree && WritingCarriesSelfOrthogonality;   // SelfOrthogonalityFollowsTheWriting is a reading of the labels, not counted

    private static GaussianRational G(int re, int im = 0) => new(re, im);

    public TwoByTwoExceptionalPointWitness()
    {
        var rows = new List<RowReading>();
        var onPoint = new HashSet<string>();
        void Row(string setting, string writing, bool companion, Mat at, Mat off, bool readVector)
        {
            onPoint.Add($"{at[0, 0]}|{at[0, 1]}|{at[1, 0]}|{at[1, 1]}");
            bool? self = readVector
                ? TwoByTwoExceptionalPointClaim.BilinearSquare(TwoByTwoExceptionalPointClaim.CoalescedVector(at)).IsZero
                : null;
            rows.Add(new RowReading(setting, writing, companion, TwoByTwoExceptionalPointClaim.IsJordanBlock(at),
                !TwoByTwoExceptionalPointClaim.IsJordanBlock(off), self));
        }
        Mat Companion(int a1, int a0) => new[,] { { G(0), G(1) }, { G(-a0), G(-a1) } };

        // closed blocks of a generator, each at an integer point on its exceptional point and one step off it
        Row("two dephased qubits, the monotonicity block [[0, −4J], [4J, −4γ]] at γ = 2J (J = 1, γ = 2; off: γ = 3)", "real opposite sign", false,
            new[,] { { G(0), G(-4) }, { G(4), G(-8) } }, new[,] { { G(0), G(-4) }, { G(4), G(-12) } }, true);
        Row("the horizon at N = 2, λ² + 4γλ + 4J² at J = γ = 1 (off: J = 2)", "companion", true, Companion(4, 4), Companion(4, 16), false);
        Row("the horizon at N = 3, λ² + 4γλ + 2J² at γ = 1, J² = 2 (off: J² = 3)", "companion", true, Companion(4, 4), Companion(4, 6), false);
        Row("γ₀'s oscillator z̈ + 4γż + 4J²z at J = γ = 1 (off: J = 2)", "companion", true, Companion(4, 4), Companion(4, 16), false);
        Row("the ring's {0, 2} block [[0, B], [−B, −4γ]] at B = 2γ = 2 (off: B = 3)", "real opposite sign", false,
            new[,] { { G(0), G(2) }, { G(-2), G(-4) } }, new[,] { { G(0), G(3) }, { G(-3), G(-4) } }, true);
        Row("the proton doublet [[0, 2J], [−2J, −2γ]] at γ = 2J = 2 (off: γ = 3)", "real opposite sign", false,
            new[,] { { G(0), G(2) }, { G(-2), G(-4) } }, new[,] { { G(0), G(2) }, { G(-2), G(-6) } }, true);
        Row("the N = 2 edge block [[−2γ₀, 2iJ], [2iJ, −2γ₁]] at γ₁ = γ₀ + 2J (γ₀ = J = 1, γ₁ = 3; off: γ₁ = 4)", "complex symmetric", false,
            new[,] { { G(-2), G(0, 2) }, { G(0, 2), G(-6) } }, new[,] { { G(-2), G(0, 2) }, { G(0, 2), G(-8) } }, true);
        // the sodium pairs at the owner's points (SodiumTwoMirrorsWitness); the centre (Spin32QuadrupolarMirrorsClaim) and ω are
        // a scalar shift, kept so the block is the page's, and change no verdict
        Mat Bright(long j0, long j1, long j2, long q, long omega, int order)
        {
            long c2 = Spin32QuadrupolarMirrorsClaim.BrightCentre(j0, j1, j2);
            long kappa = order == 1 ? 3 * j2 : -3 * j1;   // the page: +3J₂ for the satellites, −3J₁ for the double quanta
            var shift = new GaussianRational(new BigRational(-c2), new BigRational(-order * omega));
            return new[,] { { shift + G(0, (int)(-2 * q)), G((int)kappa) }, { G((int)kappa), shift + G(0, (int)(2 * q)) } };
        }
        var sat = SodiumTwoMirrorsWitness.SatelliteEpPoint;
        var dq = SodiumTwoMirrorsWitness.DoubleQuantumEpPoint;
        Row($"the sodium satellites −c₂ − iω + [[−2iq, 3J₂], [3J₂, 2iq]] at the owner's point (J₀, J₁, J₂, q, ω) = {sat}, c₂ from " +
            "Spin32QuadrupolarMirrorsClaim.BrightCentre (off: q + 1)", "complex symmetric", false,
            Bright(sat.J0, sat.J1, sat.J2, sat.Q, sat.Omega, 1), Bright(sat.J0, sat.J1, sat.J2, sat.Q + 1, sat.Omega, 1), true);
        Row($"the sodium double quanta −c₂ − 2iω + [[−2iq, −3J₁], [−3J₁, 2iq]] at the owner's point {dq} (off: q + 1)", "complex symmetric", false,
            Bright(dq.J0, dq.J1, dq.J2, dq.Q, dq.Omega, 2), Bright(dq.J0, dq.J1, dq.J2, dq.Q + 1, dq.Omega, 2), true);
        Row("the neural pair −s + [[x, −h], [h, −x]] at x = h = 1, s = 2 (off: x = 2)", "real opposite sign", false,
            new[,] { { G(-1), G(-1) }, { G(1), G(-3) } }, new[,] { { G(0), G(-1) }, { G(1), G(-4) } }, true);
        // reductions, toys and controls
        Row("the F86 toy [[−2γ₀(2k−1), iJg], [iJg, −2γ₀(2k+1)]] at Jg = 2γ₀ (γ₀ = k = 1, Jg = 2; off: Jg = 3)", "complex symmetric", false,
            new[,] { { G(-2), G(0, 2) }, { G(0, 2), G(-6) } }, new[,] { { G(-2), G(0, 3) }, { G(0, 3), G(-6) } }, true);
        Row("the seed toy [[s, 1], [1, −s]] at s = i (off: s = 2i)", "complex symmetric", false,
            new[,] { { G(0, 1), G(1) }, { G(1), G(0, -1) } }, new[,] { { G(0, 2), G(1) }, { G(1), G(0, -2) } }, true);
        Row("the Route-B control −4 + i[[q − 1, ε], [ε, −(q − 1)]] at q − 1 = iε, ε = 1 (off: ε = 2)", "complex symmetric", false,
            new[,] { { G(-5), G(0, 1) }, { G(0, 1), G(-3) } }, new[,] { { G(-5), G(0, 2) }, { G(0, 2), G(-3) } }, true);
        Rows = rows;
        DistinctOnPointMatrices = onPoint.Count;

        // the F86 toy's owner, TwoLevelEpModel, asked at an integer point on its EP (γ₀ = 1, J = 2, g_eff = 1, k = 1) and one
        // step off (J = 3): its discriminant and eigenvalues are exact there (integers, √0), compared with == against the
        // law's block [[−2γ₀(2k−1), iJg], [iJg, −2γ₀(2k+1)]]: d = its Discriminant, −c = its TraceCentre
        var toy = new[,] { { G(-2), G(0, 2) }, { G(0, 2), G(-6) } };
        var onModel = new TwoLevelEpModel(1.0, 2.0, 1.0);
        var offModel = new TwoLevelEpModel(1.0, 3.0, 1.0);
        var toyOff = new[,] { { G(-2), G(0, 3) }, { G(0, 3), G(-6) } };
        OwnerToyAgrees = onModel.Regime == EpRegime.AtEp && onModel.Discriminant == 0.0
            && offModel.Regime != EpRegime.AtEp
            && TwoByTwoExceptionalPointClaim.D(toyOff) == G((int)offModel.Discriminant)
            && TwoByTwoExceptionalPointClaim.HalfTrace(toy) == G((int)onModel.TraceCentre)
            && onModel.Eigenvalues.LamPlus == onModel.TraceCentre && onModel.Eigenvalues.LamMinus == onModel.TraceCentre
            && TwoByTwoExceptionalPointClaim.IsJordanBlock(toy) && !TwoByTwoExceptionalPointClaim.IsJordanBlock(toyOff);

        bool qAtD0 = true;
        foreach (int g in new[] { 1, 2, 4 })
        {
            // the owner's Q_EP read into ℚ exactly (2/g at g = 1, 2, 4 is a quarter-integer double), set as J at γ₀ = 1,
            // and the toy block [[−2, iJg], [iJg, −6]] built from it must sit on d = 0, a quarter step in J off it not
            double qEp = new QEpLaw(g).Value;
            var jRat = new BigRational(new System.Numerics.BigInteger(qEp * 4), 4);
            bool exactRead = (double)(qEp * 4) == Math.Floor(qEp * 4);
            GaussianRational Off(BigRational jj) => new(BigRational.Zero, jj * new BigRational(g));
            var at = new[,] { { G(-2), Off(jRat) }, { Off(jRat), G(-6) } };
            var step = jRat + new BigRational(1, 4);
            var off = new[,] { { G(-2), Off(step) }, { Off(step), G(-6) } };
            qAtD0 &= exactRead && TwoByTwoExceptionalPointClaim.D(at).IsZero && !TwoByTwoExceptionalPointClaim.D(off).IsZero;
        }
        QEpLawSitsOnTheToyPoint = qAtD0;

        // F155's dimer: single-excitation block and the no-jump generator, g = J = 1
        // the single-excitation block read off the built H on |01⟩, |10⟩, at g = J = 1 and one step off (g = 2)
        static Mat SeBlock(Mat hh) => new[,] { { hh[1, 1], hh[1, 2] }, { hh[2, 1], hh[2, 2] } };
        var h = DimerHamiltonian(1, 1);
        DimerBlockIsJordan = TwoByTwoExceptionalPointClaim.IsJordanBlock(SeBlock(h))
                             && !TwoByTwoExceptionalPointClaim.IsJordanBlock(SeBlock(DimerHamiltonian(1, 2)));
        var gen = NoJumpGenerator(h);
        var power = gen;
        var nullities = new List<int>();
        for (int p = 1; p <= 4; p++)
        {
            nullities.Add(16 - Rank(power));
            power = Mul(power, gen);
        }
        DimerNoJumpNullities = nullities;

        // the writing: diag(1, i) carries the ring's block to the complex-symmetric one
        var ring = new[,] { { G(0), G(2) }, { G(-2), G(-4) } };
        var d = new[,] { { G(1), G(0) }, { G(0), G(0, 1) } };
        var dInv = new[,] { { G(1), G(0) }, { G(0), G(0, -1) } };
        var turned = Mul(Mul(dInv, ring), d);
        bool form = Equal(turned, new[,] { { G(0), G(0, 2) }, { G(0, 2), G(-4) } });
        var v = (G(1), G(-1));
        var vTurned = (dInv[0, 0] * v.Item1, dInv[1, 1] * v.Item2);
        WritingCarriesSelfOrthogonality = form
            && TwoByTwoExceptionalPointClaim.BilinearSquare(v) == G(2)
            && TwoByTwoExceptionalPointClaim.BilinearSquare(vTurned).IsZero
            && vTurned == (G(1), G(0, 1));
    }

    /// <summary>H = J(XX + YY)/2 + ig(Z₀ − Z₁)/2 on |00⟩, |01⟩, |10⟩, |11⟩: J on the |01⟩, |10⟩ flip-flop,
    /// ig·(z₀ − z₁)/2 on the diagonal.</summary>
    private static Mat DimerHamiltonian(int coupling, int gain)
    {
        var h = new GaussianRational[4, 4];
        for (int i = 0; i < 4; i++) for (int j = 0; j < 4; j++) h[i, j] = GaussianRational.Zero;
        h[1, 2] = G(coupling);
        h[2, 1] = G(coupling);
        for (int b = 0; b < 4; b++)
        {
            int z0 = (b & 2) == 0 ? 1 : -1, z1 = (b & 1) == 0 ? 1 : -1;
            h[b, b] = new GaussianRational(BigRational.Zero, new BigRational(gain * (z0 - z1), 2));
        }
        return h;
    }

    /// <summary>ρ ↦ −i(Hρ − ρH†) as a 16 × 16 matrix on the row-major vec of ρ: −i(H⊗1 − 1⊗H̄).</summary>
    private static Mat NoJumpGenerator(Mat h)
    {
        var g = new GaussianRational[16, 16];
        for (int r = 0; r < 16; r++) for (int c = 0; c < 16; c++) g[r, c] = GaussianRational.Zero;
        var minusI = new GaussianRational(BigRational.Zero, -BigRational.One);
        for (int i = 0; i < 4; i++)
            for (int j = 0; j < 4; j++)
                for (int k = 0; k < 4; k++)
                {
                    g[4 * i + j, 4 * k + j] += minusI * h[i, k];                 // −i·H·ρ
                    g[4 * i + j, 4 * i + k] -= minusI * h[j, k].Conjugate;      // +i·ρ·H†: (ρH†)_{ij} = Σ_k ρ_{ik} conj(H_{jk})
                }
        return g;
    }

    private static Mat Mul(Mat a, Mat b)
    {
        int n = a.GetLength(0), m = b.GetLength(1), l = a.GetLength(1);
        var r = new GaussianRational[n, m];
        for (int i = 0; i < n; i++)
            for (int j = 0; j < m; j++)
            {
                var s = GaussianRational.Zero;
                for (int k = 0; k < l; k++)
                    if (!a[i, k].IsZero && !b[k, j].IsZero) s += a[i, k] * b[k, j];
                r[i, j] = s;
            }
        return r;
    }

    private static bool Equal(Mat a, Mat b)
    {
        for (int i = 0; i < a.GetLength(0); i++)
            for (int j = 0; j < a.GetLength(1); j++)
                if (a[i, j] != b[i, j]) return false;
        return true;
    }

    private static int Rank(Mat source)
    {
        int n = source.GetLength(0), m = source.GetLength(1);
        var a = (GaussianRational[,])source.Clone();
        int rank = 0;
        for (int col = 0; col < m && rank < n; col++)
        {
            int pivot = -1;
            for (int r = rank; r < n; r++)
                if (!a[r, col].IsZero) { pivot = r; break; }
            if (pivot < 0) continue;
            for (int c = 0; c < m; c++) (a[rank, c], a[pivot, c]) = (a[pivot, c], a[rank, c]);
            for (int r = 0; r < n; r++)
            {
                if (r == rank || a[r, col].IsZero) continue;
                var f = a[r, col] / a[rank, col];
                for (int c = col; c < m; c++) a[r, c] -= f * a[rank, c];
            }
            rank++;
        }
        return rank;
    }

    private static string Verdict(bool ok, string yes, string no) => ok ? yes : no;

    public string DisplayName => "One 2 × 2 in many settings, read live";

    public string Summary =>
        $"{Rows.Count} settings ({DistinctOnPointMatrices} distinct on-point matrices) of the one closed 2 × 2, each at an integer point on its exceptional point and one step off it, exact over ℚ(i): " +
        $"{Verdict(Rows.All(r => r.Holds), "every one a Jordan block on its point and none off it", "a row FAILS")} " +
        $"({Rows.Count(r => r.Companion)} of them companion matrices, whose on-point reading cannot fail); the coalesced vector " +
        $"{Verdict(SelfOrthogonalityFollowsTheWriting, "self-orthogonal exactly in the complex-symmetric rows", "does NOT follow the writing")}, " +
        $"and diag(1, i) {Verdict(WritingCarriesSelfOrthogonality, "turns the ring's (1, −1), vᵀv = 2, into (1, i), vᵀv = 0", "does NOT act as stated")}; " +
        $"F155's dimer {Verdict(DimerBlockIsJordan, "a Jordan block at g = J", "NOT a Jordan block")}, its no-jump generator " +
        $"{Verdict(DimerNoJumpSizeThree, "with a block of size 3 there", "WITHOUT a block of size 3")} (nullities {string.Join(", ", DimerNoJumpNullities)}). " +
        $"Owners: TwoLevelEpModel {Verdict(OwnerToyAgrees, "agrees exactly on and off its point", "DISAGREES")}; QEpLaw's value " +
        $"{Verdict(QEpLawSitsOnTheToyPoint, "puts the toy's block on d = 0", "does NOT put the toy's block on d = 0")}.";

    public IEnumerable<IInspectable> Children
    {
        get
        {
            foreach (var r in Rows)
                yield return new InspectableNode(
                    displayName: r.Setting,
                    summary: $"{r.Writing}{(r.Companion ? " (companion: the on-point reading cannot fail)" : "")}: " +
                             $"{Verdict(r.JordanAtPoint, "a Jordan block on its point", "NOT a Jordan block on its point")}, " +
                             $"{Verdict(r.NoneOffPoint, "none one step off", "STILL one off the point")}" +
                             (r.SelfOrthogonal is { } so ? $"; the coalesced vector {(so ? "self-orthogonal" : "not self-orthogonal")}" : ""),
                    provenance: NodeProvenance.Live);
            yield return new InspectableNode(
                displayName: "the owners asked",
                summary: $"TwoLevelEpModel at γ₀ = 1, J = 2, g_eff = 1 (and J = 3): " +
                         $"{Verdict(OwnerToyAgrees, "discriminant 0.0, regime AtEp, centre −4 = the law's half trace, off the point its discriminant = d", "DISAGREES with the law")}; " +
                         $"QEpLaw at g_eff = 1, 2, 4, its value set as J: " +
                         $"{Verdict(QEpLawSitsOnTheToyPoint, "the toy's block on d = 0, a quarter step off it not", "NOT on d = 0")}",
                provenance: NodeProvenance.Live);
            yield return new InspectableNode(
                displayName: "F155's PT dimer, a Hamiltonian",
                summary: $"single-excitation block [[ig, J], [J, −ig]] at g = J {Verdict(DimerBlockIsJordan, "a Jordan block", "NOT a Jordan block")}; " +
                         $"the no-jump generator ρ ↦ −i(Hρ − ρH†): nullities of G, G², G³, G⁴ {string.Join(", ", DimerNoJumpNullities)}, " +
                         $"{Verdict(DimerNoJumpSizeThree, "a Jordan block of size 3 and none larger", "NOT size 3")} (H nilpotent there, " +
                         "the Kronecker sum of its 2 × 2 block with the conjugate, 3 ⊕ 1)",
                provenance: NodeProvenance.Live);
        }
    }

    public InspectablePayload Payload => InspectablePayload.Empty;
}
