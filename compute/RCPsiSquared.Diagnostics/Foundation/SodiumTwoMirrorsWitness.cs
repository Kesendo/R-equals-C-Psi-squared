using System.Numerics;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Numerics;
using RCPsiSquared.Core.Symmetry;
using Mat = RCPsiSquared.Core.Numerics.GaussianRational[,];

namespace RCPsiSquared.Diagnostics.Foundation;

/// <summary>The live lab for F167 (<c>experiments/THE_SODIUM_NUCLEUS_TWO_MIRRORS.md</c>): one spin
/// 3/2, levels m = 3/2, 1/2, −1/2, −3/2, under secular quadrupolar relaxation at infinite
/// temperature, dissipator Σ_m J_|m|·(T₂ₘρT₂ₘ† − ½{T₂ₘ†T₂ₘ, ρ}), Hamiltonian ω_Q·T₂₀ + ω·I_z.
///
/// <para>The page's verifier works with sympy symbols, a Gröbner basis and F158 through Pauli
/// strings. This witness takes the other road: it builds the 16 × 16 generator from the SPIN
/// MATRICES, entry by entry on matrix units, and reads everything exactly over ℚ(i). The
/// irrationals are handled, not rounded: I_± has the entries √3, 2, √3, so the spin matrices are
/// held as A + √3·B with A, B rational, and the witness checks that every T₂,ₘ with m ≠ 0 has a
/// vanishing rational part (the middle entry 2 is killed in T₂,±₁ and absent in T₂,±₂), so
/// T₂,ₘ = √3·R₂,ₘ with R rational and the quadratic dissipator carries 3·R ρ R† exactly;
/// T₂₀ = D/√6 with D = 3I_z² − 15/4 = 3·diag(1, −1, −1, 1) (checked), so its dissipator carries
/// D ρ D/6, and the splitting enters as the integer q, the coefficient of diag(1, −1, −1, 1) in H,
/// q = ω_Q·√6/2 (ω = √6·ω_Q is then ω = 2q). Integer inputs J₀, J₁, J₂, q, ω; no double anywhere
/// in a gate.</para>
///
/// <para>Children, the gate first: the five Dirac matrices (anticommuting, squaring to 3/2, 3, 3,
/// 3, 3, the strings ZZ, ZX, ZY, XI, YI), read once from the √3 representation and once, at
/// spins 1, 3/2, 2, 5/2, from a rational diagonal similarity (I_+ ↦ X − m(m + 1) on the
/// superdiagonal, I_− ↦ 1), the control spins breaking it; the one map S(X) = Γ̂₀·Xᵀ with
/// S·L·S⁻¹ = −L + 2·L_{J₀} − 6(J₁ + J₂)·1 at every point, beside two controls (Γ̂₀·X without the
/// transpose, Xᵀ without Γ̂₀) that must break it; the bright projector P commuting with L and S,
/// L_{J₀} = −3J₀·P, the dark and bright parts paired about −c₁ and −c₂; and the global-centre
/// theorem at a fixed set of exact points covering both outcomes, each read by route A (F158's two
/// ends, dim ker L against dim ker(L + 2s) as exact ranks, only for J ≥ 0) and route B (the
/// characteristic polynomial χ with χ(−x − 2s) = χ(x) identically, any real J), the two routes
/// against each other and against the claim's <see cref="Spin32QuadrupolarMirrorsClaim.HasGlobalCentre"/>;
/// and the handshake with NMR: the textbook biexponential relaxation (Hubbard's weights 1/5 and 4/5
/// for I_z, 3/5 and 2/5 for I_+ at ω_Q = 0), the 3:4:3 triplet, the central transition and the
/// triple quantum blind to J₀ and to the splitting, extreme narrowing, the rank-2 coherence T₂,₁ as
/// the satellites' mirror image, each as an exact eigenoperator equation on string parts, with a
/// splitting and swapped rate assignments as the controls.</para>
///
/// <para>Typed home <see cref="Spin32QuadrupolarMirrorsClaim"/>.</para></summary>
public sealed class SodiumTwoMirrorsWitness : IInspectable
{
    /// <summary>The largest |input| accepted; every comparison is exact at any size, this only
    /// keeps the inspect-time arithmetic small.</summary>
    public const long MaxMagnitude = 1000;

    private const int Levels = 4, Dim = 16;

    /// <summary>One exact point of the global-centre theorem. <see cref="NearKernel"/> and
    /// <see cref="FarKernel"/> are route A (null where a density is negative, outside F158's class,
    /// and where 2s = 0, every density zero, where the two ends are one matrix and the comparison says nothing);
    /// <see cref="CharpolySymmetric"/> is route B.</summary>
    public sealed record CentrePoint(long J0, long J1, long J2, long Q, long Omega, long TwoS,
        bool Rule, int? NearKernel, int? FarKernel, bool CharpolySymmetric,
        bool MirrorIdentityHolds, bool SectorMirrorsHold)
    {
        public bool RouteAApplies => NearKernel.HasValue;
        public bool? RouteAPairs => NearKernel.HasValue ? NearKernel == FarKernel : null;
        public bool Agrees => CharpolySymmetric == Rule && (RouteAPairs is null || RouteAPairs == Rule);
    }

    /// <summary>The Dirac-algebra reading at spin I = TwiceSpin/2 through the rational similarity:
    /// whether the five components anticommute pairwise, whether each squares to a multiple of 1,
    /// and if so the five constants.</summary>
    public sealed record SpinDiracReading(int TwiceSpin, bool Anticommute, bool SquaresScalar,
        IReadOnlyList<GaussianRational>? Squares)
    {
        public bool IsDirac => Anticommute && SquaresScalar;
    }

    /// <summary>The fixed points (J₀, J₁, J₂, q, ω): the page's counting point first (J₁ = J₂ = 0,
    /// J₀ = 1, ω = √6·ω_Q, i.e. q = 3, ω = 6), then breaks with J₀ > 0 and J₁·J₂ = 0 away from
    /// ω_Q = ω = 0, J₀ = 0 with both frequencies on, the boundary pairings at ω_Q = ω = 0, every
    /// density positive, and four negative-density points read by route B alone.</summary>
    public static readonly IReadOnlyList<(long J0, long J1, long J2, long Q, long Omega)> FixedPoints = new[]
    {
        (1L, 0L, 0L, 3L, 6L),     // the counting point: breaks, 4 against 2
        (1L, 0L, 2L, 1L, 0L),     // J₁ = 0, splitting on: breaks
        (2L, 1L, 0L, 0L, 1L),     // J₂ = 0, offset on: breaks
        (0L, 1L, 2L, 2L, 3L),     // J₀ = 0, both frequencies on: pairs
        (1L, 2L, 3L, 1L, 1L),     // every density positive: breaks
        (1L, 0L, 2L, 0L, 0L),     // J₁ = 0 at ω_Q = ω = 0: pairs
        (2L, 1L, 0L, 0L, 0L),     // J₂ = 0 at ω_Q = ω = 0: pairs
        (3L, 0L, 0L, 0L, 0L),     // J₀ alone at ω_Q = ω = 0: pairs, 8 against 8
        (0L, 0L, 0L, 1L, 2L),     // no relaxation: L = −i[H, ·]
        (1L, 1L, 1L, 0L, 0L),     // extreme narrowing: breaks
        (0L, -1L, 2L, 1L, 2L),    // negative density, J₀ = 0: pairs (route B only)
        (1L, -1L, 2L, 0L, 0L),    // negative density, J₁·J₂ ≠ 0: breaks (route B only)
        (-1L, 0L, 2L, 0L, 0L),    // negative J₀, J₁ = 0 at ω_Q = ω = 0: pairs (route B only)
        (-1L, 0L, 2L, 1L, 0L),    // negative J₀, J₁ = 0, splitting on: breaks (route B only)
    };

    /// <summary>The control point for the two broken maps (H ≠ 0 and J₁, J₂ ≠ 0, so each control
    /// has something to break).</summary>
    public static readonly (long J0, long J1, long J2, long Q, long Omega) ControlPoint = (1, 2, 3, 1, 1);

    public long J0 { get; }
    public long J1 { get; }
    public long J2 { get; }

    /// <summary>q = ω_Q·√6/2, the coefficient of diag(1, −1, −1, 1) in H.</summary>
    public long Q { get; }
    public long Omega { get; }

    // gate 1: the Dirac algebra
    public bool SqrtThreeFactorisationHolds { get; }
    public bool AdjointRelationsHold { get; }
    public bool SplittingParameterHolds { get; }
    public bool DiracAlgebraHolds { get; }
    public IReadOnlyList<GaussianRational> SquareConstants { get; }
    public bool PauliStringsMatch { get; }
    public IReadOnlyList<string> PauliLabels { get; } = new[] { "ZZ", "ZX", "ZY", "XI", "YI" };
    public IReadOnlyList<GaussianRational> PauliFactors { get; }
    public IReadOnlyList<SpinDiracReading> SpinReadings { get; }

    // gate 2: the one map
    public bool InverseHolds { get; }
    public bool MirrorIdentityHolds { get; }
    public bool TransposeDroppedBreaks { get; }
    public bool UnitInPlaceOfGammaBreaks { get; }

    // gate 4: the two sectors at the user point
    public bool ProjectorCommutesWithL { get; }
    public bool ProjectorCommutesWithS { get; }
    public bool J0PartIsMinusThreeJ0P { get; }
    public bool DarkSectorPairs { get; }
    public bool BrightSectorPairs { get; }
    public bool SectorMirrorsHold => ProjectorCommutesWithL && ProjectorCommutesWithS && J0PartIsMinusThreeJ0P && DarkSectorPairs && BrightSectorPairs;

    // gate 5: the handshake with NMR
    public bool LongitudinalSplitHolds { get; }
    public IReadOnlyList<GaussianRational> LongitudinalWeights { get; }
    public bool TransverseSplitHolds { get; }
    public IReadOnlyList<GaussianRational> TransverseWeights { get; }
    public bool SplittingDetunesSatellites { get; }
    public bool SwappedAssignmentBreaks { get; }
    public bool SplittingFeedsPartner { get; }
    public IReadOnlyList<GaussianRational> TripletIntensities { get; }
    public bool DarkLinesBlindToJ0AndSplitting { get; }
    public bool TripleQuantumOnDarkCentre { get; }
    public bool PopulationsPairAboutDarkCentre { get; }
    public bool ExtremeNarrowingOneRate { get; }
    public bool SatellitePartnerOnBrightMirror { get; }
    /// <summary>The eigenoperator equations and their controls; the weights and the triplet are
    /// readings of the spin matrices (no generator input) and stand beside it.</summary>
    public bool HandshakeHolds => LongitudinalSplitHolds && TransverseSplitHolds && SplittingDetunesSatellites && SwappedAssignmentBreaks
        && SplittingFeedsPartner
        && DarkLinesBlindToJ0AndSplitting && TripleQuantumOnDarkCentre && PopulationsPairAboutDarkCentre
        && ExtremeNarrowingOneRate && SatellitePartnerOnBrightMirror;

    // gate 3: the global-centre theorem
    public IReadOnlyList<CentrePoint> Points { get; }
    public CentrePoint UserPoint { get; }
    public bool TheoremHolds => Points.All(p => p.Agrees) && UserPoint.Agrees;

    private readonly Mat _d;                         // D = 3I_z² − 15/4
    private readonly Mat _g0;                        // Γ̂₀ = D/3 = diag(1, −1, −1, 1)
    private readonly Mat _iz;
    private readonly Dictionary<int, Mat> _r = new(); // T₂,ₘ = √3·R₂,ₘ, m = ±1, ±2
    private readonly Mat _s, _sInv, _p;

    public SodiumTwoMirrorsWitness(long j0 = 1, long j1 = 2, long j2 = 3, long q = 1, long omega = 1)
    {
        foreach (var (v, name) in new[] { (j0, nameof(j0)), (j1, nameof(j1)), (j2, nameof(j2)), (q, nameof(q)), (omega, nameof(omega)) })
            if (v > MaxMagnitude || v < -MaxMagnitude)
                throw new ArgumentOutOfRangeException(name, v, $"|input| ≤ {MaxMagnitude}");
        (J0, J1, J2, Q, Omega) = (j0, j1, j2, q, omega);

        // ---- spin 3/2 from the spin formula, entries in ℚ(√3)
        var x = new BigRational(15, 4);
        _iz = Zero(Levels);
        var ipA = Zero(Levels);
        var ipB = Zero(Levels);
        for (int k = 0; k < Levels; k++)
        {
            var m = new BigRational(3 - 2 * k, 2);
            _iz[k, k] = m;
            if (k == 0) continue;
            var (a, b) = SqrtInQ3(x - m * (m + 1));       // I_+|m_k⟩ = √(X − m(m + 1))·|m_{k−1}⟩
            ipA[k - 1, k] = a;
            ipB[k - 1, k] = b;
        }
        var ip = new Root3(ipA, ipB);
        var im = new Root3(Transpose(ipA), Transpose(ipB)); // I_− = I_+ᵀ (real)
        var izR = new Root3(_iz, Zero(Levels));
        var t = new Dictionary<int, Root3>
        {
            [1] = (izR * ip + ip * izR).Scale(new BigRational(-1, 2)),
            [-1] = (izR * im + im * izR).Scale(new BigRational(1, 2)),
            [2] = (ip * ip).Scale(new BigRational(1, 2)),
            [-2] = (im * im).Scale(new BigRational(1, 2)),
        };
        SqrtThreeFactorisationHolds = t.Values.All(tm => IsZero(tm.A)) && t.Values.All(tm => !IsZero(tm.B));
        foreach (var (mm, tm) in t) _r[mm] = tm.B;
        AdjointRelationsHold = Equal(Adjoint(_r[1]), Scale(_r[-1], -1)) && Equal(Adjoint(_r[2]), _r[-2]);

        _d = Sub(Scale(Mul(_iz, _iz), 3), Scale(Identity(Levels), x));
        _g0 = Diagonal(1, -1, -1, 1);
        SplittingParameterHolds = Equal(_d, Scale(_g0, 3));

        // the five components as (scale², rational matrix): Γ₀ = D/√6, Γ₁..Γ₄ = √3·(…)
        var gammas = new (BigRational Scale2, Mat G)[]
        {
            (new BigRational(1, 6), _d),
            (3, Add(_r[1], Adjoint(_r[1]))),
            (3, Scale(Sub(_r[1], Adjoint(_r[1])), GaussianRational.I)),
            (3, Add(_r[2], Adjoint(_r[2]))),
            (3, Scale(Sub(_r[2], Adjoint(_r[2])), GaussianRational.I)),
        };
        bool hermitian = gammas.All(g => Equal(Adjoint(g.G), g.G));
        var (anti, squares) = ReadClifford(gammas);
        DiracAlgebraHolds = hermitian && anti && squares is not null;
        SquareConstants = squares ?? Array.Empty<GaussianRational>();

        var pauli = new[] { Kron(PZ, PZ), Kron(PZ, PX), Kron(PZ, PY), Kron(PX, PI), Kron(PY, PI) };
        var factors = new List<GaussianRational>();
        bool match = true;
        for (int a = 0; a < 5; a++)
        {
            var f = ProportionalFactor(gammas[a].G, pauli[a]);
            match &= f is not null;
            factors.Add(f ?? GaussianRational.Zero);
        }
        PauliStringsMatch = match;
        PauliFactors = factors;

        SpinReadings = new[] { 2, 3, 4, 5 }.Select(ReadSpinBySimilarity).ToArray();

        // ---- the one map and its inverse, as superoperators
        _s = Super(X => Mul(_g0, Transpose(X)));
        _sInv = Super(Y => Mul(Transpose(Y), _g0));
        InverseHolds = Equal(Mul(_s, _sInv), Identity(Dim)) && Equal(Mul(_sInv, _s), Identity(Dim));
        bool controlInverses = Equal(Mul(Super(X => Mul(_g0, X)), Super(Y => Mul(_g0, Y))), Identity(Dim))
                               && Equal(Mul(Super(Transpose), Super(Transpose)), Identity(Dim));
        InverseHolds &= controlInverses;
        _p = Zero(Dim);
        for (int i = 0; i < Levels; i++)
            for (int j = 0; j < Levels; j++)
                if (IsOuter(i) != IsOuter(j)) _p[Levels * i + j, Levels * i + j] = GaussianRational.One;
        ProjectorCommutesWithS = Equal(Mul(_p, _s), Mul(_s, _p));

        // ---- gate 2 controls at the fixed control point
        var (cl, cj0) = Generator(ControlPoint.J0, ControlPoint.J1, ControlPoint.J2, ControlPoint.Q, ControlPoint.Omega);
        var cRhs = MirrorRight(cl, cj0, ControlPoint.J1, ControlPoint.J2);
        var noT = Super(X => Mul(_g0, X));
        var noTInv = Super(Y => Mul(_g0, Y));
        TransposeDroppedBreaks = !Equal(Mul(Mul(noT, cl), noTInv), cRhs);
        var unit = Super(Transpose);
        UnitInPlaceOfGammaBreaks = !Equal(Mul(Mul(unit, cl), unit), cRhs);

        // ---- gate 3 at every fixed point and at the user point
        Points = FixedPoints.Select(pt => ReadPoint(pt.J0, pt.J1, pt.J2, pt.Q, pt.Omega)).ToArray();
        UserPoint = ReadPoint(j0, j1, j2, q, omega);

        // ---- gates 2 and 4 at the user point
        var (l, lj0) = Generator(j0, j1, j2, q, omega);
        MirrorIdentityHolds = UserPoint.MirrorIdentityHolds;
        ProjectorCommutesWithL = Equal(Mul(_p, l), Mul(l, _p));
        J0PartIsMinusThreeJ0P = Equal(lj0, Scale(_p, -3 * j0));
        DarkSectorPairs = SectorPairs(l, dark: true, 2 * Spin32QuadrupolarMirrorsClaim.DarkCentre(j1, j2));
        BrightSectorPairs = SectorPairs(l, dark: false, 2 * Spin32QuadrupolarMirrorsClaim.BrightCentre(j0, j1, j2));

        // ---- gate 5: the handshake with NMR, the textbook relaxation read as strings and counted
        long c1 = Spin32QuadrupolarMirrorsClaim.DarkCentre(j1, j2), c2 = Spin32QuadrupolarMirrorsClaim.BrightCentre(j0, j1, j2);
        var zi = Diagonal(1, 1, -1, -1);                   // ZI on |00⟩, |01⟩, |10⟩, |11⟩
        var iz2 = Diagonal(1, -1, 1, -1);                  // IZ
        var central = ipA;                                 // the CT part of I_+: 2·|1/2⟩⟨−1/2|
        var satellites = ipB;                              // I_+ = central + √3·satellites
        var partner = Zero(Levels);                        // the satellites' antisymmetric partner
        partner[0, 1] = GaussianRational.One;
        partner[2, 3] = -GaussianRational.One;
        var tq = Zero(Levels);
        tq[0, 3] = GaussianRational.One;
        GaussianRational Rate(long re, long freq) => new(-re, -freq);   // −re − i·freq
        var l0 = Generator(j0, j1, j2, 0, omega).L;                       // the user's rates, no splitting

        // T₁: I_z = ZI + ½·IZ, the parts at 6J₂ and 6J₁, the weights their norms squared
        LongitudinalSplitHolds = Equal(_iz, Add(zi, Scale(iz2, new BigRational(1, 2))))
                                 && Eigen(l, zi, Rate(6 * j2, 0)) && Eigen(l, iz2, Rate(6 * j1, 0));
        LongitudinalWeights = new[] { Weight(zi, _iz), Weight(Scale(iz2, new BigRational(1, 2)), _iz) };
        // T₂ at ω_Q = 0: I_+ = CT + √3·satellites, the CT part at c₁, the satellites at 3(J₀ + J₁)
        TransverseSplitHolds = Eigen(l, central, Rate(c1, omega)) && Eigen(l0, satellites, Rate(3 * (j0 + j1), omega));
        var sNorm = 3 * Trace(Mul(Adjoint(satellites), satellites));
        var cNorm = Trace(Mul(Adjoint(central), central));
        TransverseWeights = new[] { sNorm / (sNorm + cNorm), cNorm / (sNorm + cNorm) };
        // controls: a splitting detunes the two satellite cells (ω ± 2q); at the control point, where J₁ ≠ J₂ and
        // J₀ ≠ J₂, the swapped rate assignments fail
        SplittingDetunesSatellites = !Eigen(Generator(j0, j1, j2, 1, omega).L, satellites, Rate(3 * (j0 + j1), omega));
        var (cj0c, cj1c, cj2c, _, cwc) = ControlPoint;
        var cl0 = Generator(cj0c, cj1c, cj2c, 0, cwc).L;
        SwappedAssignmentBreaks = !Eigen(cl0, zi, Rate(6 * cj1c, 0)) && !Eigen(cl0, iz2, Rate(6 * cj2c, 0))
                                  && !Eigen(cl0, satellites, Rate(Spin32QuadrupolarMirrorsClaim.DarkCentre(cj1c, cj2c), cwc))
                                  && !Eigen(cl0, central, Rate(3 * (cj0c + cj1c), cwc));
        // 3:4:3, the squared entries of I_+
        TripletIntensities = new[] { 3 * satellites[0, 1] * satellites[0, 1], central[1, 2] * central[1, 2], 3 * satellites[2, 3] * satellites[2, 3] };
        // the central transition and the triple quantum: dark to J₀ and to the splitting alike, on −c₁
        var lj0Unit = Generator(1, 0, 0, 0, 0).LJ0;
        var splittingUnit = Generator(0, 0, 0, 1, 0).L;
        DarkLinesBlindToJ0AndSplitting = IsZero(Apply(lj0Unit, central)) && IsZero(Apply(splittingUnit, central))
                                         && IsZero(Apply(lj0Unit, tq)) && IsZero(Apply(splittingUnit, tq))
                                         && !IsZero(Apply(lj0Unit, satellites)) && !IsZero(Apply(splittingUnit, satellites));
        // the splitting reaches the partner and nothing else: −i[Γ̂₀, satellites] = −2i·partner
        SplittingFeedsPartner = Equal(Apply(splittingUnit, satellites), Scale(partner, new GaussianRational(0, -2)));
        TripleQuantumOnDarkCentre = Eigen(l, tq, Rate(c1, 3 * omega));
        // the populations pair about −c₁: quadrupolar order Γ̂₀ at −2c₁ against equilibrium at 0 (L is unital), and
        // the two T₁ modes 6J₁ + 6J₂ = 2c₁ with each other (their eigen equations above)
        PopulationsPairAboutDarkCentre = Eigen(l, _g0, Rate(2 * c1, 0)) && IsZero(Apply(l, Identity(Levels)));
        // extreme narrowing: one density for all three, one rate for the four observable parts
        var ln = Generator(1, 1, 1, 0, omega).L;
        ExtremeNarrowingOneRate = Eigen(ln, zi, Rate(6, 0)) && Eigen(ln, iz2, Rate(6, 0))
                                  && Eigen(ln, central, Rate(6, omega)) && Eigen(ln, satellites, Rate(6, omega));
        // at ω_Q = 0: the satellites' mirror image −S(satellites)† (S keeps the bright sector, † brings the order back) is
        // the rank-2 coherence T₂,₁ (up to a factor), an eigenoperator at 2c₂ − 3(J₀ + J₁), orthogonal to I_+
        SatellitePartnerOnBrightMirror = Equal(Scale(Adjoint(Apply(_s, satellites)), -1), partner)
                                         && ProportionalFactor(partner, _r[1]) is not null
                                         && Trace(Mul(Adjoint(partner), satellites)).IsZero && Trace(Mul(Adjoint(partner), central)).IsZero
                                         && Eigen(l0, partner, Rate(2 * c2 - 3 * (j0 + j1), omega));
    }

    // ===== the handshake =====

    /// <summary>The superoperator applied to a 4 × 4 operator, vec index 4i + j.</summary>
    private static Mat Apply(Mat super, Mat op)
    {
        var r = Zero(Levels);
        for (int a = 0; a < Levels; a++)
            for (int b = 0; b < Levels; b++)
                for (int i = 0; i < Levels; i++)
                    for (int j = 0; j < Levels; j++)
                        if (!op[i, j].IsZero) r[a, b] += super[Levels * a + b, Levels * i + j] * op[i, j];
        return r;
    }

    private static bool Eigen(Mat super, Mat op, GaussianRational lambda) => Equal(Apply(super, op), Scale(op, lambda));

    private static GaussianRational Trace(Mat a)
    {
        var t = GaussianRational.Zero;
        for (int i = 0; i < a.GetLength(0); i++) t += a[i, i];
        return t;
    }

    /// <summary>Tr(part†part)/Tr(whole†whole), the weight of an orthogonal part in the decay of the whole.</summary>
    private static GaussianRational Weight(Mat part, Mat whole) =>
        Trace(Mul(Adjoint(part), part)) / Trace(Mul(Adjoint(whole), whole));

    // ===== the generator =====

    /// <summary>L and its J₀ part L_{J₀}, as 16 × 16 superoperators on vec index 4i + j of |i⟩⟨j|.</summary>
    private (Mat L, Mat LJ0) Generator(long j0, long j1, long j2, long q, long omega)
    {
        var h = Add(Scale(_g0, q), Scale(_iz, omega));
        var d2 = Mul(_d, _d);
        Mat Dissipate(Mat c, Mat cDagC, Mat rho) =>
            Sub(Mul(Mul(c, rho), Adjoint(c)), Scale(Add(Mul(cDagC, rho), Mul(rho, cDagC)), new BigRational(1, 2)));
        Mat J0Part(Mat rho) => Scale(Dissipate(_d, d2, rho), new BigRational(j0, 6));
        Mat Full(Mat rho)
        {
            var acc = Scale(Sub(Mul(h, rho), Mul(rho, h)), -GaussianRational.I);
            acc = Add(acc, J0Part(rho));
            foreach (var (m, r) in _r)
            {
                long rate = Math.Abs(m) == 1 ? j1 : j2;
                if (rate == 0) continue;
                acc = Add(acc, Scale(Dissipate(r, Mul(Adjoint(r), r), rho), 3 * rate)); // T = √3·R
            }
            return acc;
        }
        return (Super(Full), Super(J0Part));
    }

    /// <summary>−L + 2·L_{J₀} − 6(J₁ + J₂)·1, the right side of the mirror identity.</summary>
    private static Mat MirrorRight(Mat l, Mat lj0, long j1, long j2) =>
        Sub(Add(Scale(l, -1), Scale(lj0, 2)), Scale(Identity(Dim), 6 * (j1 + j2)));

    private CentrePoint ReadPoint(long j0, long j1, long j2, long q, long omega)
    {
        var (l, lj0) = Generator(j0, j1, j2, q, omega);
        long twoS = Spin32QuadrupolarMirrorsClaim.TwiceGlobalCentre(j0, j1, j2);
        bool rule = Spin32QuadrupolarMirrorsClaim.HasGlobalCentre(j0, j1, j2, q, omega);

        int? near = null, far = null;
        if (j0 >= 0 && j1 >= 0 && j2 >= 0 && twoS != 0)
        {
            near = Dim - Rank(l);
            far = Dim - Rank(Add(l, Scale(Identity(Dim), twoS)));
        }
        bool charpoly = CharpolySymmetric(l, Enumerable.Range(0, Dim).ToArray(), twoS);

        bool mirror = Equal(Mul(Mul(_s, l), _sInv), MirrorRight(l, lj0, j1, j2));
        bool sectors = Equal(Mul(_p, l), Mul(l, _p))
                       && Equal(lj0, Scale(_p, -3 * j0))
                       && SectorPairs(l, dark: true, 2 * Spin32QuadrupolarMirrorsClaim.DarkCentre(j1, j2))
                       && SectorPairs(l, dark: false, 2 * Spin32QuadrupolarMirrorsClaim.BrightCentre(j0, j1, j2));
        return new CentrePoint(j0, j1, j2, q, omega, twoS, rule, near, far, charpoly, mirror, sectors);
    }

    private static bool IsOuter(int level) => level == 0 || level == Levels - 1;

    private static bool SectorPairs(Mat l, bool dark, long twoC)
    {
        var idx = new List<int>();
        for (int i = 0; i < Levels; i++)
            for (int j = 0; j < Levels; j++)
                if ((IsOuter(i) == IsOuter(j)) == dark) idx.Add(Levels * i + j);
        return CharpolySymmetric(l, idx.ToArray(), twoC);
    }

    /// <summary>Route B: χ of the principal submatrix on <paramref name="idx"/> satisfies
    /// χ(−x − twoC) = χ(x) identically. The matrix is scaled by the lcm c of its denominators to
    /// Gaussian integers (Berkowitz, division free); the eigenvalues scale by c, so the shift does
    /// too. The degree is even here, so no sign enters.</summary>
    private static bool CharpolySymmetric(Mat l, int[] idx, long twoC)
    {
        int n = idx.Length;
        BigInteger c = BigInteger.One;
        foreach (int a in idx)
            foreach (int b in idx)
                c = BigRational.Lcm(BigRational.Lcm(c, l[a, b].Re.Denominator), l[a, b].Im.Denominator);
        var m = new GaussianInteger[n, n];
        for (int a = 0; a < n; a++)
            for (int b = 0; b < n; b++)
            {
                var re = l[idx[a], idx[b]].Re * c;
                var im = l[idx[a], idx[b]].Im * c;
                m[a, b] = new GaussianInteger(re.Numerator, im.Numerator);  // integral by the choice of c
            }
        var chi = GaussianPolynomial.Trim(GaussianMatrixCharpoly.Characteristic(m));
        var mirrored = GaussianPolynomial.ComposeLinear(chi, -GaussianInteger.One, new GaussianInteger(-c * twoC, BigInteger.Zero));
        return chi.Length == mirrored.Length && chi.Zip(mirrored).All(z => z.First.Equals(z.Second));
    }

    // ===== the Dirac algebra =====

    private static (bool Anticommute, IReadOnlyList<GaussianRational>? Squares) ReadClifford((BigRational Scale2, Mat G)[] gammas)
    {
        bool anti = true;
        for (int a = 0; a < gammas.Length; a++)
            for (int b = a + 1; b < gammas.Length; b++)
                anti &= IsZero(Add(Mul(gammas[a].G, gammas[b].G), Mul(gammas[b].G, gammas[a].G)));
        var squares = new List<GaussianRational>();
        foreach (var (s2, g) in gammas)
        {
            var sq = Mul(g, g);
            if (!IsScalar(sq, out var c)) return (anti, null);
            squares.Add((GaussianRational)s2 * c);
        }
        return (anti, squares);
    }

    /// <summary>The control route at spin I = twiceSpin/2: the diagonal similarity that makes I_±
    /// rational (I_+ ↦ X − m(m + 1) on the superdiagonal, I_− ↦ 1 below it, X = I(I + 1)) keeps
    /// anticommutators and squares, and builds the components from T₂,±ₘ alone (Γ₁ = T₂,₁ − T₂,₋₁,
    /// Γ₂ = i(T₂,₁ + T₂,₋₁), Γ₃ = T₂,₂ + T₂,₋₂, Γ₄ = i(T₂,₂ − T₂,₋₂), using T₂,₁† = −T₂,₋₁ and
    /// T₂,₂† = T₂,₋₂, which the √3 route checks at 3/2).</summary>
    private static SpinDiracReading ReadSpinBySimilarity(int twiceSpin)
    {
        int d = twiceSpin + 1;
        var x = new BigRational(twiceSpin * (twiceSpin + 2), 4);
        var iz = Zero(d);
        var ip = Zero(d);
        var im = Zero(d);
        for (int k = 0; k < d; k++)
        {
            var m = new BigRational(twiceSpin - 2 * k, 2);
            iz[k, k] = m;
            if (k == 0) continue;
            ip[k - 1, k] = x - m * (m + 1);
            im[k, k - 1] = GaussianRational.One;
        }
        var half = new BigRational(1, 2);
        var t1 = Scale(Add(Mul(iz, ip), Mul(ip, iz)), -half);
        var tm1 = Scale(Add(Mul(iz, im), Mul(im, iz)), half);
        var t2 = Scale(Mul(ip, ip), half);
        var tm2 = Scale(Mul(im, im), half);
        var dq = Sub(Scale(Mul(iz, iz), 3), Scale(Identity(d), x));
        var gammas = new (BigRational, Mat)[]
        {
            (new BigRational(1, 6), dq),
            (1, Sub(t1, tm1)),
            (1, Scale(Add(t1, tm1), GaussianRational.I)),
            (1, Add(t2, tm2)),
            (1, Scale(Sub(t2, tm2), GaussianRational.I)),
        };
        var (anti, squares) = ReadClifford(gammas);
        return new SpinDiracReading(twiceSpin, anti, squares is not null, squares);
    }

    private static GaussianRational? ProportionalFactor(Mat g, Mat p)
    {
        GaussianRational? f = null;
        for (int i = 0; i < Levels && f is null; i++)
            for (int j = 0; j < Levels && f is null; j++)
                if (!p[i, j].IsZero) f = g[i, j] / p[i, j];
        if (f is null || f.Value.IsZero) return null;
        return Equal(g, Scale(p, f.Value)) ? f : null;
    }

    private static readonly Mat PI = Diagonal(1, 1);
    private static readonly Mat PZ = Diagonal(1, -1);
    private static readonly Mat PX = new GaussianRational[,] { { 0, 1 }, { 1, 0 } };
    private static readonly Mat PY = new GaussianRational[,] { { 0, -GaussianRational.I }, { GaussianRational.I, 0 } };

    // ===== ℚ(√3) for the spin matrices =====

    /// <summary>A + √3·B with A, B rational matrices.</summary>
    private readonly record struct Root3(Mat A, Mat B)
    {
        public static Root3 operator *(Root3 u, Root3 v) =>
            new(Add(Mul(u.A, v.A), SodiumTwoMirrorsWitness.Scale(Mul(u.B, v.B), 3)), Add(Mul(u.A, v.B), Mul(u.B, v.A)));
        public static Root3 operator +(Root3 u, Root3 v) => new(Add(u.A, v.A), Add(u.B, v.B));
        public Root3 Scale(BigRational s) => new(SodiumTwoMirrorsWitness.Scale(A, s), SodiumTwoMirrorsWitness.Scale(B, s));
    }

    /// <summary>√r in ℚ(√3) for a rational r ≥ 0: (√r, 0) when r is a rational square, (0, √(r/3))
    /// when r/3 is; anything else is outside the field and throws.</summary>
    private static (GaussianRational A, GaussianRational B) SqrtInQ3(BigRational r)
    {
        if (RationalSqrt(r) is { } a) return (a, GaussianRational.Zero);
        if (RationalSqrt(r / 3) is { } b) return (GaussianRational.Zero, b);
        throw new InvalidOperationException($"√{r} is not in ℚ(√3)");
    }

    private static BigRational? RationalSqrt(BigRational r)
    {
        if (r.Sign < 0) return null;
        var n = ISqrt(r.Numerator);
        var d = ISqrt(r.Denominator);
        return n * n == r.Numerator && d * d == r.Denominator ? new BigRational(n, d) : null;
    }

    private static BigInteger ISqrt(BigInteger v)
    {
        if (v.IsZero) return v;
        var x = v;
        var y = (x + 1) / 2;
        while (y < x) { x = y; y = (x + v / x) / 2; }
        return x;
    }

    // ===== exact matrices over ℚ(i) =====

    private static Mat Zero(int n)
    {
        var m = new GaussianRational[n, n];
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++) m[i, j] = GaussianRational.Zero;
        return m;
    }

    private static Mat Identity(int n)
    {
        var m = Zero(n);
        for (int i = 0; i < n; i++) m[i, i] = GaussianRational.One;
        return m;
    }

    private static Mat Diagonal(params int[] entries)
    {
        var m = Zero(entries.Length);
        for (int i = 0; i < entries.Length; i++) m[i, i] = entries[i];
        return m;
    }

    private static Mat Mul(Mat a, Mat b)
    {
        int n = a.GetLength(0), k = a.GetLength(1), p = b.GetLength(1);
        var r = new GaussianRational[n, p];
        for (int i = 0; i < n; i++)
            for (int j = 0; j < p; j++)
            {
                var s = GaussianRational.Zero;
                for (int t = 0; t < k; t++)
                    if (!a[i, t].IsZero && !b[t, j].IsZero) s += a[i, t] * b[t, j];
                r[i, j] = s;
            }
        return r;
    }

    private static Mat Add(Mat a, Mat b) => Map2(a, b, (u, v) => u + v);
    private static Mat Sub(Mat a, Mat b) => Map2(a, b, (u, v) => u - v);
    private static Mat Scale(Mat a, GaussianRational s) => Map2(a, a, (u, _) => s * u);
    private static Mat Scale(Mat a, BigRational s) => Scale(a, (GaussianRational)s);

    private static Mat Map2(Mat a, Mat b, Func<GaussianRational, GaussianRational, GaussianRational> f)
    {
        int n = a.GetLength(0), p = a.GetLength(1);
        var r = new GaussianRational[n, p];
        for (int i = 0; i < n; i++)
            for (int j = 0; j < p; j++) r[i, j] = f(a[i, j], b[i, j]);
        return r;
    }

    private static Mat Transpose(Mat a)
    {
        int n = a.GetLength(0), p = a.GetLength(1);
        var r = new GaussianRational[p, n];
        for (int i = 0; i < n; i++)
            for (int j = 0; j < p; j++) r[j, i] = a[i, j];
        return r;
    }

    private static Mat Adjoint(Mat a)
    {
        var r = Transpose(a);
        for (int i = 0; i < r.GetLength(0); i++)
            for (int j = 0; j < r.GetLength(1); j++) r[i, j] = r[i, j].Conjugate;
        return r;
    }

    private static Mat Kron(Mat a, Mat b)
    {
        int n = a.GetLength(0), m = b.GetLength(0);
        var r = new GaussianRational[n * m, n * m];
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++)
                for (int k = 0; k < m; k++)
                    for (int l = 0; l < m; l++) r[i * m + k, j * m + l] = a[i, j] * b[k, l];
        return r;
    }

    private static bool Equal(Mat a, Mat b)
    {
        if (a.GetLength(0) != b.GetLength(0) || a.GetLength(1) != b.GetLength(1)) return false;
        for (int i = 0; i < a.GetLength(0); i++)
            for (int j = 0; j < a.GetLength(1); j++)
                if (a[i, j] != b[i, j]) return false;
        return true;
    }

    private static bool IsZero(Mat a)
    {
        foreach (var z in a) if (!z.IsZero) return false;
        return true;
    }

    private static bool IsScalar(Mat a, out GaussianRational c)
    {
        c = a[0, 0];
        return Equal(a, Scale(Identity(a.GetLength(0)), c));
    }

    /// <summary>The superoperator of a linear map on 4 × 4 matrices, column 4i + j holding the image
    /// of |i⟩⟨j| in the same vec order.</summary>
    private static Mat Super(Func<Mat, Mat> f)
    {
        var s = Zero(Dim);
        for (int i = 0; i < Levels; i++)
            for (int j = 0; j < Levels; j++)
            {
                var e = Zero(Levels);
                e[i, j] = GaussianRational.One;
                var y = f(e);
                for (int a = 0; a < Levels; a++)
                    for (int b = 0; b < Levels; b++) s[Levels * a + b, Levels * i + j] = y[a, b];
            }
        return s;
    }

    /// <summary>Exact rank over ℚ(i) by Gaussian elimination (route A).</summary>
    private static int Rank(Mat a)
    {
        var m = (Mat)a.Clone();
        int rows = m.GetLength(0), cols = m.GetLength(1), rank = 0;
        for (int col = 0; col < cols && rank < rows; col++)
        {
            int piv = -1;
            for (int r = rank; r < rows; r++) if (!m[r, col].IsZero) { piv = r; break; }
            if (piv < 0) continue;
            for (int c = 0; c < cols; c++) (m[rank, c], m[piv, c]) = (m[piv, c], m[rank, c]);
            for (int r = rank + 1; r < rows; r++)
            {
                if (m[r, col].IsZero) continue;
                var f = m[r, col] / m[rank, col];
                for (int c = col; c < cols; c++) m[r, c] -= f * m[rank, c];
            }
            rank++;
        }
        return rank;
    }

    // ===== IInspectable =====

    public string DisplayName => $"SodiumTwoMirrorsWitness (J₀={J0}, J₁={J1}, J₂={J2}, q={Q}, ω={Omega})";

    private bool DiracGate => SqrtThreeFactorisationHolds && AdjointRelationsHold && SplittingParameterHolds && DiracAlgebraHolds && PauliStringsMatch;
    private bool OnlyThreeHalves => SpinReadings.All(r => r.IsDirac == (r.TwiceSpin == 3));
    private bool MirrorEverywhere => InverseHolds && Points.All(p => p.MirrorIdentityHolds) && MirrorIdentityHolds;
    private bool ControlsBreak => TransposeDroppedBreaks && UnitInPlaceOfGammaBreaks;
    private bool SectorsEverywhere => Points.All(p => p.SectorMirrorsHold) && SectorMirrorsHold;

    private static string Verdict(bool ok, string yes, string no) => ok ? yes : no;

    public string Summary =>
        $"spin 3/2 built from the spin matrices, exact over ℚ(i): the five quadrupole components " +
        $"{Verdict(DiracGate, "are Dirac matrices", "FAIL to be Dirac matrices")} (squares {string.Join(", ", SquareConstants)}; ZZ, ZX, ZY, XI, YI), " +
        $"{Verdict(OnlyThreeHalves, "at spin 3/2 only among 1, 3/2, 2, 5/2", "with the spin control FAILING")}; " +
        $"S·L·S⁻¹ = −L + 2·L_{{J₀}} − 6(J₁ + J₂) {Verdict(MirrorEverywhere, "holds", "FAILS")} at the {Points.Count} fixed points and the user point, " +
        $"its two controls {Verdict(ControlsBreak, "break it", "FAIL to break it")}; the dark and bright sectors " +
        $"{Verdict(SectorsEverywhere, "pair about −c₁ and −c₂", "FAIL to pair")}; the global-centre rule " +
        $"{Verdict(TheoremHolds, "holds", "FAILS")} at {Points.Count(p => p.Rule)} pairing and " +
        $"{Points.Count(p => !p.Rule)} breaking fixed points and at the user point, route A (F158's kernels) and route B (χ) " +
        $"{Verdict(TheoremHolds, "agreeing with it and with each other", "DISAGREEING")}; at the user point the spectrum " +
        $"{(UserPoint.CharpolySymmetric ? "has" : "has no")} global centre −s, 2s = {UserPoint.TwoS}; the handshake with NMR " +
        $"{Verdict(HandshakeHolds, "holds (the textbook parts of I_z and I_+ are eigenoperators at the textbook rates)", "FAILS")}.";

    private static string Show(CentrePoint p) =>
        $"(J₀, J₁, J₂, q, ω) = ({p.J0}, {p.J1}, {p.J2}, {p.Q}, {p.Omega}), 2s = {p.TwoS}: rule {(p.Rule ? "pairs" : "breaks")}, " +
        (p.RouteAApplies ? $"route A dim ker L = {p.NearKernel} against dim ker(L + 2s) = {p.FarKernel}, "
            : p.TwoS == 0 ? "route A empty (2s = 0, every density zero), " : "route A outside F158 (a negative density), ") +
        $"route B χ(−x − 2s) {(p.CharpolySymmetric ? "=" : "≠")} χ(x){(p.Agrees ? "" : " (DISAGREEMENT)")}";

    public IEnumerable<IInspectable> Children
    {
        get
        {
            yield return new InspectableNode(
                displayName: "the gate: five Dirac matrices at spin 3/2",
                summary: $"from I_± in ℚ(√3): every T₂,ₘ with m ≠ 0 {Verdict(SqrtThreeFactorisationHolds, "is √3 times a rational matrix", "is NOT √3 times a rational matrix")}, " +
                         $"T₂,₁† = −T₂,₋₁ and T₂,₂† = T₂,₋₂ {Verdict(AdjointRelationsHold, "hold", "FAIL")}, D = 3I_z² − 15/4 " +
                         $"{Verdict(SplittingParameterHolds, "is 3·diag(1, −1, −1, 1)", "is NOT 3·diag(1, −1, −1, 1)")} (so ω_Q·T₂₀ = q·diag(1, −1, −1, 1), q = ω_Q·√6/2); " +
                         $"Γ₀..Γ₄ {Verdict(DiracAlgebraHolds, "are Hermitian, anticommute pairwise and square to", "FAIL the Clifford relations; squares")} " +
                         $"{string.Join(", ", SquareConstants)} times 1, and their rational parts (Γ₀ = G₀/√6, Γₐ = √3·Gₐ) {Verdict(PauliStringsMatch, "are", "are NOT")} " +
                         $"{string.Join(", ", PauliLabels.Zip(PauliFactors, (s, f) => $"({f})·{s}"))}.",
                provenance: NodeProvenance.Live);

            yield return new InspectableNode(
                displayName: "the spin control",
                summary: "through the rational similarity I_+ ↦ X − m(m + 1), I_− ↦ 1: " +
                         string.Join("; ", SpinReadings.Select(r =>
                             $"I = {(r.TwiceSpin % 2 == 0 ? $"{r.TwiceSpin / 2}" : $"{r.TwiceSpin}/2")}: {(r.Anticommute ? "anticommute" : "do NOT anticommute")}, " +
                             $"{(r.SquaresScalar ? $"square to {string.Join(", ", r.Squares!)}" : "squares NOT multiples of 1")}")) +
                         $"; {Verdict(OnlyThreeHalves, "Dirac at 3/2 alone, so the gate can fail", "the control is COMPROMISED")}" +
                         $"{Verdict(SpinReadings.Single(r => r.TwiceSpin == 3).Squares?.SequenceEqual(SquareConstants) == true, ", and at 3/2 the two routes read the same squares", ", the two routes read DIFFERENT squares at 3/2")}.",
                provenance: NodeProvenance.Live);

            yield return new InspectableNode(
                displayName: "one map, two mirrors",
                summary: $"S(X) = Γ̂₀·Xᵀ, S⁻¹(Y) = Yᵀ·Γ̂₀ {Verdict(InverseHolds, "(inverse checked)", "(inverse FAILS)")}: S·L·S⁻¹ = −L + 2·L_{{J₀}} − 6(J₁ + J₂)·1 " +
                         $"{Verdict(MirrorEverywhere, "holds exactly", "FAILS")} at the user point and at all {Points.Count} fixed points{(Points.Any(p => p.J0 < 0 || p.J1 < 0 || p.J2 < 0) ? ", negative densities included" : "")}; " +
                         $"controls at (J₀, J₁, J₂, q, ω) = {ControlPoint}: Γ̂₀·X without the transpose {Verdict(TransposeDroppedBreaks, "breaks it", "does NOT break it")}, " +
                         $"Xᵀ with 1 in place of Γ̂₀ {Verdict(UnitInPlaceOfGammaBreaks, "breaks it", "does NOT break it")}.",
                provenance: NodeProvenance.Live);

            yield return new InspectableNode(
                displayName: "the two sectors",
                summary: $"P, the projector on the bright cells |m⟩⟨m′| (|m| ≠ |m′|), {Verdict(ProjectorCommutesWithL, "commutes", "does NOT commute")} with L and " +
                         $"{Verdict(ProjectorCommutesWithS, "with S", "NOT with S")}; L_{{J₀}} {Verdict(J0PartIsMinusThreeJ0P, "=", "≠")} −3J₀·P; the dark part " +
                         $"{Verdict(DarkSectorPairs, "pairs", "does NOT pair")} about −c₁ = −{Spin32QuadrupolarMirrorsClaim.DarkCentre(J1, J2)} and the bright part " +
                         $"{Verdict(BrightSectorPairs, "pairs", "does NOT pair")} about −c₂ = −{Spin32QuadrupolarMirrorsClaim.BrightCentre(J0, J1, J2)} (χ of each 8 × 8 block); " +
                         $"{Verdict(Points.All(p => p.SectorMirrorsHold), "the same at every fixed point", "NOT at every fixed point")}.",
                provenance: NodeProvenance.Live);

            yield return new InspectableNode(
                displayName: "the handshake with NMR",
                summary: $"the textbook relaxation of a spin 3/2, read as strings: I_z = ZI + ½·IZ with ZI at 6J₂ and IZ at 6J₁ " +
                         $"{Verdict(LongitudinalSplitHolds, "(eigenoperators)", "(NOT eigenoperators)")}, squared norms {string.Join(", ", LongitudinalWeights)}; " +
                         $"I_+ = CT + √3·satellites with the CT part at c₁ = {Spin32QuadrupolarMirrorsClaim.DarkCentre(J1, J2)} and, at ω_Q = 0, the satellites at " +
                         $"3(J₀ + J₁) = {3 * (J0 + J1)} {Verdict(TransverseSplitHolds, "(eigenoperators)", "(NOT eigenoperators)")}, squared norms {string.Join(", ", TransverseWeights)}" +
                         $"{(J1 == J2 || J0 == J2 ? " (at this point two of the rates coincide; the control point below separates them)" : "")}; " +
                         $"controls: a splitting {Verdict(SplittingDetunesSatellites, "detunes the satellites", "does NOT detune the satellites")}, " +
                         $"{Verdict(SplittingFeedsPartner, "feeding exactly the partner below (−i[Γ̂₀, satellites] = −2i·partner)", "NOT feeding the partner alone")}, and at {ControlPoint} the swapped " +
                         $"assignments {Verdict(SwappedAssignmentBreaks, "fail", "do NOT fail")}; the squared entries of I_+ {string.Join(" : ", TripletIntensities)}; " +
                         $"the central transition and the triple quantum {Verdict(DarkLinesBlindToJ0AndSplitting, "blind to J₀ and to the splitting alike, the satellites not", "NOT blind as stated")}; " +
                         $"the triple quantum {Verdict(TripleQuantumOnDarkCentre, "on −c₁, the mean of the two T₁ rates", "NOT on −c₁")}; quadrupolar order Γ̂₀ " +
                         $"{Verdict(PopulationsPairAboutDarkCentre, "at −2c₁, the partner of equilibrium about −c₁, as the two T₁ modes are partners", "NOT at −2c₁")}; extreme narrowing " +
                         $"{Verdict(ExtremeNarrowingOneRate, "gives one rate", "does NOT give one rate")}; R₂′ − R₂″ = 3(J₀ − J₂) = (c₂ − c₁) − 3J₂ (arithmetic); " +
                         $"at ω_Q = 0 the satellites' mirror image −S(satellites)† {Verdict(SatellitePartnerOnBrightMirror, $"is the rank-2 coherence T₂,₁, orthogonal to I_+, and sits at 2c₂ − 3(J₀ + J₁) = {3 * (J0 + J1 + 2 * J2)}, the satellites' partner about −c₂", "is NOT T₂,₁ at its rate as stated")}.",
                provenance: NodeProvenance.Live);

            yield return new InspectableNode(
                displayName: "the global-centre theorem, two routes",
                summary: $"one centre −s exactly when J₀ = 0, or ω_Q = ω = 0 and J₁·J₂ = 0: {Verdict(TheoremHolds, "both routes meet the rule at every point", "a route or the rule DISAGREES")}.",
                children: Points.Append(UserPoint).Select((p, i) => (IInspectable)new InspectableNode(
                    displayName: i == Points.Count ? "the user point" : $"point {i + 1}",
                    summary: Show(p),
                    provenance: NodeProvenance.Live)).ToArray(),
                provenance: NodeProvenance.Live);
        }
    }

    public InspectablePayload Payload => InspectablePayload.Empty;
}
