using System.Numerics;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Numerics;
using RCPsiSquared.Core.Symmetry;
using Mat = RCPsiSquared.Core.Numerics.GaussianRational[,];

namespace RCPsiSquared.Diagnostics.Foundation;

/// <summary>The live lab for F167 (<c>experiments/THE_SODIUM_NUCLEUS_TWO_MIRRORS.md</c>): one spin
/// 3/2, levels m = 3/2, 1/2, −1/2, −3/2, under quadrupolar relaxation secular in the Larmor frequency
/// (the splitting kept whole) at infinite temperature, dissipator Σ_m J_|m|·(T₂ₘρT₂ₘ† − ½{T₂ₘ†T₂ₘ, ρ}), Hamiltonian ω_Q·T₂₀ + ω·I_z.
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
/// splitting and swapped rate assignments as the controls; and the bright pairs' exceptional point:
/// the satellite and double-quantum pairs as closed 2 × 2 blocks whose two rates sum to 2c₂ at every
/// splitting, with the discriminants 4(9J₂² − 4q²) and 4(9J₁² − 4q²), each read ON its exceptional
/// point (2q = 3J₂, 2q = 3J₁, the density non-zero) as a 2 × 2 Jordan block of the full generator
/// (geometric multiplicity 1, algebraic 2, no eigenvalue there one step off in q) beside the scalar
/// control at J₂ = q = 0, the coalesced vector self-orthogonal, and the satellite part of I_+ in the rates a
/// measurement reads: its initial rate R₂′ = 3(J₀ + J₁) at the five splittings read, below the exceptional point two simple
/// rates with the weights 4/3 and −1/3 (at <see cref="BelowEpPoint"/>), above it two mode frequencies
/// pulled in to D_ω² = D² − R₁(4/5)², and the merging where the nominal separation D = 4q equals the
/// 4/5 longitudinal rate; and the dynamic frequency shift (gate 7, at the fixed <see cref="ShiftPoint"/>, K₁ and K₂
/// integers): the Lamb-shift Hamiltonian from the spin matrices against its closed form, the mirror identity with it,
/// the bright pairs moving as wholes, the Jordan block at the shifted eigenvalue, and the shifted global-centre rule
/// on both routes, the kernel dimensions and the characteristic polynomial, at eleven points; and a spin-lock (gate 8): the
/// split of I_x into its central-transition and satellite parts read off the spin matrices, the mirror identity and the
/// sectors kept by the central-transition part and broken by the satellite part at one point, and the global centre under
/// a lock on both routes at ten points (integer mixtures of the two parts; the spectral break is the verifier's); and the
/// irreducible tensors in strings (gate 9, at <see cref="TensorPoint"/>): the page's sixteen rows read for rank and order
/// off the spin matrices, the ten that are eigenoperators and the six turned between two string classes with their
/// weights, the triple-quantum path's two equal and opposite overlaps, the one map on the tensors, and one string class
/// split by the order. These are statements on eigenvalues, not on the lines of a spectrum; which come
/// first for an observable reading a pair is the page's section, gated by its verifier.</para>
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

    /// <summary>The fixed points of gate 6, (J₀, J₁, J₂, q, ω), integers chosen so that every root read is a
    /// Gaussian integer: the satellite pair ON its exceptional point (2q = 3J₂), the double-quantum pair ON
    /// its own (2q = 3J₁), the diabolic control (J₂ = 0 at q = 0, where the satellite block is a scalar), one
    /// point below the satellite EP (√(9J₂² − 4q²) = 9) and one above it (√(4q² − 9J₂²) = 16). J₀ > 0 keeps
    /// the central transition off the bright centre and ω ≠ 0 keeps the orders ±p apart, so no other
    /// eigenvalue meets the one read.</summary>
    public static readonly (long J0, long J1, long J2, long Q, long Omega) SatelliteEpPoint = (1, 4, 2, 3, 1);
    public static readonly (long J0, long J1, long J2, long Q, long Omega) DoubleQuantumEpPoint = (1, 2, 4, 3, 1);
    public static readonly (long J0, long J1, long J2, long Q, long Omega) DiabolicPoint = (1, 4, 0, 0, 1);
    public static readonly (long J0, long J1, long J2, long Q, long Omega) BelowEpPoint = (1, 3, 5, 6, 1);
    public static readonly (long J0, long J1, long J2, long Q, long Omega) AboveEpPoint = (1, 3, 4, 10, 1);

    /// <summary>One bright pair read ON its exceptional point, λ* = −c₂ − ip·ω (order p = 1 the satellites, 2
    /// the double quantum): the block's part N = M − λ*·1, non-zero with N² = 0; the full generator's ranks
    /// of (L − λ*)^k, k = 1, 2, 3, so the geometric multiplicity 16 − rank₁ and the algebraic 16 − rank₃;
    /// and the rank of L − λ* one step off in q on either side, where λ* is no eigenvalue.</summary>
    public sealed record EpReading(long J0, long J1, long J2, long Q, long Omega, int Order, GaussianRational Lambda,
        bool PartNonZero, bool PartSquaresToZero, int Rank1, int Rank2, int Rank3, int RankStepDown, int RankStepUp)
    {
        public int Geometric => Dim - Rank1;
        public int Algebraic => Dim - Rank3;
        public bool IsJordanPair => PartNonZero && PartSquaresToZero && Geometric == 1 && Algebraic == 2
                                    && RankStepDown == Dim && RankStepUp == Dim;
    }

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

    // gate 6: the bright pairs' exceptional point
    /// <summary>Both bright pairs, the satellites |3/2⟩⟨1/2|, |−1/2⟩⟨−3/2| and the double quantum
    /// |3/2⟩⟨−1/2|, |1/2⟩⟨−3/2|, span L-invariant sets of two cells at the user point and at the two EP
    /// points; the control set, a satellite with the central transition in place of its partner, does not
    /// (at the satellite EP point, where J₂ ≠ 0).</summary>
    public bool BrightPairsClosed { get; }
    /// <summary>The sum rule: at the user point each pair's block has the trace −2c₂ − 2ip·ω, so its two rates
    /// sum to 2c₂ there (at every splitting in the verifier's symbolic blocks); the control set's trace differs (in its real part where J₀ ≠ 0, and in its imaginary
    /// part by 2q).</summary>
    public bool BrightPairTracesHold { get; }
    /// <summary>At the user point the discriminants (tr M)² − 4·det M are 4(9J₂² − 4q²) (satellites) and
    /// 4(9J₁² − 4q²) (double quantum); at the two EP points, where J₁ ≠ J₂, the swapped assignment fails.</summary>
    public bool BrightPairDiscriminantsHold { get; }
    public EpReading SatelliteEp { get; }
    public EpReading DoubleQuantumEp { get; }
    /// <summary>The same-object control: at J₂ = 0 and q = 0 the satellite block is the scalar −c₂ − iω, and in
    /// the full generator that eigenvalue is semisimple (L − λ and (L − λ)² both of rank 14).</summary>
    public bool DiabolicControlSemisimple { get; }
    /// <summary>At the satellite EP the kernel vector of L − λ*, read off N as (N₀₁, −N₀₀), is non-zero, lies on
    /// the two satellite cells and is self-orthogonal (vᵀv = 0); with two components that gives equal weight on
    /// the satellites' sum, the observable part of I_+, and on their mirror image, the partner. The control: below
    /// the exceptional point (<see cref="BelowEpPoint"/>) the slow eigenvector of the same pair is not
    /// self-orthogonal.</summary>
    public bool EpVectorSelfOrthogonal { get; }
    /// <summary>A reading: the generator is complex symmetric in the cell basis at the satellite EP point (H
    /// diagonal, and the jumps T₂,±ₘ in pairs of equal rate, so that the dissipator is the Γₐ's, each symmetric or
    /// antisymmetric), the setting in which a coalesced vector is
    /// self-orthogonal.</summary>
    public bool GeneratorComplexSymmetric { get; }
    /// <summary>The satellite part s of I_+ has ⟨s|L|s⟩/⟨s|s⟩ = −3(J₀ + J₁) − iω, the initial decay rate R₂′,
    /// at the user point, at ω_Q = 0 and at the three points around the satellite EP. The control: one satellite
    /// cell alone does not, its quotient carrying the full −c₂ and the splitting (at the satellite EP point).</summary>
    public bool InitialSlopeSplittingFree { get; }
    /// <summary>Ω = √(9J₂² − 4q²) at <see cref="BelowEpPoint"/> and Ω′ = √(4q² − 9J₂²) at
    /// <see cref="AboveEpPoint"/>, from these closed forms at the declared points (9 and 16), the values the rank
    /// readings below test against the generator and the node prints.</summary>
    public long BelowEpOmega { get; }
    public long AboveEpOmega { get; }
    /// <summary>Below the satellite EP (<see cref="BelowEpPoint"/>): the weights of the satellite part on its
    /// slow and its fast mode, at the rates c₂ ∓ Ω.</summary>
    public IReadOnlyList<GaussianRational> BelowEpWeights { get; }
    /// <summary>c₂ ∓ Ω are simple eigenvalues of L (L − λ and (L − λ)² both of rank 15) and the slow weight is
    /// 4/3; the fast weight, −1/3, is one minus it (the two spectral projectors sum to the identity), shown beside
    /// the gate rather than in it.</summary>
    public bool BelowEpTwoRates { get; }
    /// <summary>Above the satellite EP (<see cref="AboveEpPoint"/>): the two modes oscillate at ±Ω′ about the order
    /// frame and not at the nominal ±2q (neither nominal position is an eigenvalue), and their separation
    /// obeys D_ω² + R₁(4/5)² = D², with D = 4q and R₁(4/5) read off L as the rate of ZI; the 1/5 mode's rate in its
    /// place fails. Mode frequencies, not line positions (the lines are the page's verifier's).</summary>
    public bool AboveEpFrequenciesPulled { get; }
    /// <summary>The merging in the rates a measurement reads: the satellite pair is a Jordan block at the satellite EP
    /// point (the reading <see cref="SatelliteEp"/>, which carries the location), and there the nominal separation 4q is
    /// the rate of the 4/5 longitudinal mode ZI read off L, not the 1/5 mode's IZ; the double-quantum pair is one at
    /// the double-quantum EP point, where 4q is IZ's rate and not ZI's.</summary>
    public bool EpAtTheLongitudinalRates { get; }
    /// <summary>Readings at the user point: where each bright pair stands against its exceptional point, from
    /// the sign of its discriminant and from its coupling, both read off the generator's block.</summary>
    public string SatelliteRegime { get; }
    public string DoubleQuantumRegime { get; }
    public bool ExceptionalPointHolds => BrightPairsClosed && BrightPairTracesHold && BrightPairDiscriminantsHold
        && SatelliteEp.IsJordanPair && DoubleQuantumEp.IsJordanPair && DiabolicControlSemisimple && EpVectorSelfOrthogonal
        && InitialSlopeSplittingFree && BelowEpTwoRates && AboveEpFrequenciesPulled && EpAtTheLongitudinalRates;

    // gate 7: the dynamic frequency shift
    /// <summary>The point (J₀, J₁, J₂, q, ω, K₁, K₂) at which the shift is read, the shift added as the Lamb-shift
    /// Hamiltonian K₁(T₂,₁†T₂,₁ − T₂,₋₁†T₂,₋₁) + K₂(T₂,₂†T₂,₂ − T₂,₋₂†T₂,₋₂), diagonal: −3K₂·ZI − 3K₁·IZ in strings.
    /// The satellite exceptional point is read with the same K₁, K₂. At (1, 3), 2K₁ ≠ K₂, so the part of the shift that
    /// no frame removes is present (at 2K₁ = K₂ the shift is only a turn of the frame).</summary>
    public static readonly (long J0, long J1, long J2, long Q, long Omega, long K1, long K2) ShiftPoint = (1, 2, 3, 1, 1, 1, 3);
    /// <summary>The mirror identity S·L·S⁻¹ = −L + 2·L_{J₀} − 6(J₁ + J₂) holds with the shift (it is real and
    /// diagonal, so Γ₀ commutes with it and the transpose reverses it).</summary>
    public bool ShiftKeepsMirror { get; }
    /// <summary>Each bright pair stays closed and moves as a whole, both cells by 6i·K₁ (satellites) or 6i·K₂ (double
    /// quanta), so its discriminant is unchanged; the central transition moves by 6i(K₂ − K₁).</summary>
    public bool ShiftMovesPairsWhole { get; }
    /// <summary>At the satellite exceptional point with the shift, the eigenvalue λ* + 6i·K₁ has the ranks of a
    /// 2 × 2 Jordan block in L (15, 14, 14), and one step off in q it is no eigenvalue.</summary>
    public bool ShiftKeepsJordan { get; }
    /// <summary>The Lamb-shift Hamiltonian built from the spin matrices is the diagonal
    /// (−3(K₁+K₂), 3(K₁−K₂), −3(K₁−K₂), 3(K₁+K₂)) on m = 3/2 … −3/2, that is −3K₂·ZI − 3K₁·IZ, read at
    /// (K₁, K₂) = (1, 0), (0, 1) and the shift point, which pins the linear map.</summary>
    public bool ShiftClosedForm { get; }
    public sealed record ShiftCentreCase(long J0, long J1, long J2, long Q, long Omega, long K1, long K2);
    /// <summary>One shifted point read on both routes: route A, dim ker L against dim ker(L + 2s) (every case has
    /// non-negative densities and 2s ≠ 0, inside F158's class), and route B, χ(−x − 2s) = χ(x).</summary>
    public sealed record ShiftCentreReading(ShiftCentreCase Case, bool Rule, int NearKernel, int FarKernel, bool CharpolySymmetric)
    {
        public bool Agrees => Rule == CharpolySymmetric && Rule == (NearKernel == FarKernel);
    }
    /// <summary>Points for the global-centre rule with the shift: J₁ = 0 at ω = 6K₁ and J₂ = 0 at ω = 3K₂ (pairing), also
    /// with the other K non-zero (the satellites resonate at 6K₁ whatever K₂ is); the same at the unshifted resonance
    /// ω = 0, at an ω that a rule mixing the two K would accept, and with a splitting (breaking); J₀ = 0 with every
    /// term on (pairing); all densities on (breaking).</summary>
    public static readonly IReadOnlyList<ShiftCentreCase> ShiftCentreCases = new[]
    {
        new ShiftCentreCase(1, 0, 3, 0, 6, 1, 0), new ShiftCentreCase(1, 0, 3, 0, 0, 1, 0),
        new ShiftCentreCase(1, 2, 0, 0, 3, 0, 1), new ShiftCentreCase(1, 2, 0, 0, 0, 0, 1),
        new ShiftCentreCase(1, 0, 3, 0, 6, 1, 1), new ShiftCentreCase(1, 0, 3, 0, 12, 1, 1),
        new ShiftCentreCase(1, 2, 0, 0, 3, 1, 1), new ShiftCentreCase(1, 2, 0, 0, 6, 1, 1),
        new ShiftCentreCase(1, 0, 3, 1, 6, 1, 1),
        new ShiftCentreCase(0, 2, 3, 1, 1, 1, 2), new ShiftCentreCase(1, 2, 3, 0, 0, 1, 1),
    };
    public IReadOnlyList<ShiftCentreReading> ShiftCentreReadings { get; }
    public bool DynamicShiftHolds => ShiftClosedForm && ShiftKeepsMirror && ShiftMovesPairsWhole && ShiftKeepsJordan
        && ShiftCentreReadings.All(r => r.Agrees) && ShiftCentreReadings.Any(r => r.Rule) && ShiftCentreReadings.Any(r => !r.Rule);

    // gate 8: a spin-lock, ω₁·I_x in the rotating frame, split by the spin matrices into its two parts
    /// <summary>I_x = (I_+ + I_−)/2 splits along I_+'s rational and √3 parts: the rational part is the central
    /// transition, (XX + YY)/2, real and commuting with Γ₀; the √3 part is the satellites, (√3/2)·IX, anticommuting.</summary>
    public bool LockSplitHolds { get; }
    /// <summary>The lock strength of gate 8.</summary>
    public const long LockNu = 2;
    /// <summary>At <see cref="ShiftPoint"/> with LockNu times the central-transition part added: the mirror identity
    /// holds and the bright projector still commutes with L.</summary>
    public bool SelectiveLockKeepsMirror { get; }
    /// <summary>With the satellite direction added instead, the mirror identity fails and the projector no longer
    /// commutes. Both residuals are linear in the coefficient, so the √3 the exact field
    /// carries changes neither; the spectral break itself is the verifier's (x times an irreducible factor of degree 15,
    /// no eigenvalue paired at either centre).</summary>
    public bool SatelliteLockBreaksMirror { get; }
    public sealed record LockCentreCase(long A, long B, long J0, long J1, long J2, long Q, long Omega, long K1, long K2);
    public sealed record LockCentreReading(LockCentreCase Case, bool Rule, int NearKernel, int FarKernel, bool CharpolySymmetric)
    {
        public bool Agrees => Rule == CharpolySymmetric && Rule == (NearKernel == FarKernel);
    }
    /// <summary>The global centre under a lock a·(satellite direction) + b·(central transition), a and b non-zero
    /// integers in place of the field's √3 and 1: only J₁ on, q = 0 and ω = 3K₂ pairs (also with K₁ off the satellites'
    /// resonance, and at ω = 0 with no shift); the satellites' resonance ω = 6K₁ with K₂ = 0, ω off, a splitting, J₂ or J₀
    /// switched on break, and so do the two branches the lock removes, J₂ alone at the satellites' resonance (the shift's
    /// J₁ = 0 branch) and J₀ alone at q = ω = 0 (the theorem's J₁·J₂ = 0 branch).</summary>
    public static readonly IReadOnlyList<LockCentreCase> LockCentreCases = new[]
    {
        new LockCentreCase(2, 1, 0, 2, 0, 0, 3, 0, 1), new LockCentreCase(5, 3, 0, 2, 0, 0, 3, 5, 1),
        new LockCentreCase(2, 1, 0, 2, 0, 0, 0, 0, 0), new LockCentreCase(2, 1, 0, 2, 0, 0, 6, 1, 0),
        new LockCentreCase(2, 1, 0, 2, 0, 0, 1, 0, 0), new LockCentreCase(2, 1, 0, 2, 0, 1, 3, 0, 1),
        new LockCentreCase(2, 1, 0, 2, 3, 0, 3, 0, 1), new LockCentreCase(2, 1, 1, 2, 0, 0, 3, 0, 1),
        new LockCentreCase(2, 1, 0, 0, 3, 0, 6, 1, 0), new LockCentreCase(2, 1, 1, 0, 0, 0, 0, 0, 0),
    };
    public IReadOnlyList<LockCentreReading> LockCentreReadings { get; }
    public bool SpinLockHolds => LockSplitHolds && SelectiveLockKeepsMirror && SatelliteLockBreaksMirror
        && LockCentreReadings.All(r => r.Agrees) && LockCentreReadings.Any(r => r.Rule) && LockCentreReadings.Any(r => !r.Rule);

    // gate 9: the irreducible tensors in strings, at the fixed TensorPoint
    /// <summary>(J₀, J₁, J₂, ω) of gate 9, with J₀ ≠ J₂ and J₁ ≠ J₂ so that the turned pairs' two string classes decay
    /// at different rates.</summary>
    public static readonly (long J0, long J1, long J2, long Omega) TensorPoint = (1, 2, 3, 1);
    /// <summary>Each of the page's sixteen string combinations (rank l, order p) satisfies Σₐ[Iₐ, [Iₐ, M]] = l(l + 1)·M
    /// and [I_z, M] = p·M with the spin matrices in ℚ(√3); the rank-3 row at order 1 read as rank 1 fails.</summary>
    public bool TensorRowsHold { get; }
    /// <summary>At ω_Q = 0 the ten rows of rank 0, rank 2 and rank 3 at |p| = 2, 3 are eigenoperators at the string
    /// classes' rates; each of the six rank-1 and rank-3 rows at |p| ≤ 1 splits into its parts on two string classes,
    /// read off the row, each an eigenoperator, at two different rates, with the squared-norm weights 4/5 : 1/5 (p = 0)
    /// and 3/5 : 2/5 (p = ±1) and their reverse.</summary>
    public bool TensorBasesTurned { get; }
    /// <summary>A consistency reading of the weights, not an independent route: the rank-3 row at order 1 meets the
    /// satellite part of the rank-1 row, read off the row, with a squared normalized overlap of exactly 6/25, and the two
    /// rows are orthogonal, so the central-transition part's overlap is its negative; with the parts' rates of
    /// <see cref="TensorBasesTurned"/> the rank-3 content of e^{Lt}·T₁,₁ at ω = 0 (otherwise times e^{−iωt}) is
    /// ±(√6/5)·(e^{−3(J₀ + J₁)t} − e^{−c₁t}), the sign the phases' (fixed by the verifier). Both numbers follow from the
    /// 3/5 : 2/5 weights once each order-1 part is one direction.</summary>
    public bool TripleQuantumPathHolds { get; }
    /// <summary>S(X) = Γ₀·Xᵀ: T₀,₀ ↔ T₂,₀, T₂,±₂ onto T₃,∓₂, T₃,±₃ onto T₃,∓₃ (onto: proportional, the signs in the
    /// ladder's phases are the verifier's), and the reflection [[4/5, 3/5], [3/5, −4/5]] on the order-0 rows of rank 1
    /// and 3 (equal norms).</summary>
    public bool MirrorOnTensorsHolds { get; }
    /// <summary>XX + YY lives on the central-transition cells (p = ±1), XX − YY on the triple-quantum cells (p = ±3),
    /// and at ω = 0 both are eigenoperators at c₁: one string class, two orders; and the strings inside the cells of
    /// order |p| fill them only at p = 0 and |p| = 2, so only there is an order filter a projection onto strings.</summary>
    public bool OrderSplitsAStringClass { get; }
    public bool TensorDictionaryHolds => TensorRowsHold && TensorBasesTurned && TripleQuantumPathHolds
        && MirrorOnTensorsHolds && OrderSplitsAStringClass;

    // gate 3: the global-centre theorem
    public IReadOnlyList<CentrePoint> Points { get; }
    public CentrePoint UserPoint { get; }
    public bool TheoremHolds => Points.All(p => p.Agrees) && UserPoint.Agrees;

    private readonly Mat _d;                         // D = 3I_z² − 15/4
    private readonly Mat _g0;                        // Γ̂₀ = D/3 = diag(1, −1, −1, 1)
    private readonly Mat _iz;
    private readonly Mat _ixCt, _ixSat;               // I_x = _ixCt + √3·_ixSat: the central transition and the satellites
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
        _ixCt = Scale(Add(ipA, Transpose(ipA)), new BigRational(1, 2));
        _ixSat = Scale(Add(ipB, Transpose(ipB)), new BigRational(1, 2));
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

        // ---- gate 6: the bright pairs' exceptional point
        var (sj0, sj1, sj2, sq, sw) = SatelliteEpPoint;
        var (dj0, dj1, dj2, dq, dw) = DoubleQuantumEpPoint;
        var lSat = Generator(sj0, sj1, sj2, sq, sw).L;
        var lDq = Generator(dj0, dj1, dj2, dq, dw).L;
        // both pairs are closed sets of two cells; a satellite with the central transition in place of its partner is not
        BrightPairsClosed = Closed(l, SatelliteIdx) && Closed(l, DoubleQuantumIdx)
                            && Closed(lSat, SatelliteIdx) && Closed(lDq, DoubleQuantumIdx)
                            && !Closed(lSat, SatelliteWithCentralIdx);
        // the sum rule: the trace of each pair's block, −2c₂ − 2ip·ω at every splitting; the control set's is not
        long c2Sat = Spin32QuadrupolarMirrorsClaim.BrightCentre(sj0, sj1, sj2);
        BrightPairTracesHold = Trace(Block(l, SatelliteIdx)) == new GaussianRational(-2 * c2, -2 * omega)
                               && Trace(Block(l, DoubleQuantumIdx)) == new GaussianRational(-2 * c2, -4 * omega)
                               && Trace(Block(lSat, SatelliteWithCentralIdx)) != new GaussianRational(-2 * c2Sat, -2 * sw);
        // the discriminants: the satellites carry J₂, the double quantum J₁; the swapped assignment fails where J₁ ≠ J₂
        GaussianRational DiscFormula(long jp, long qq) => new(4 * (9 * jp * jp - 4 * qq * qq), 0);
        BrightPairDiscriminantsHold = Disc2(Block(l, SatelliteIdx)) == DiscFormula(j2, q)
                                      && Disc2(Block(l, DoubleQuantumIdx)) == DiscFormula(j1, q)
                                      && Disc2(Block(lSat, SatelliteIdx)) != DiscFormula(sj1, sq)
                                      && Disc2(Block(lDq, DoubleQuantumIdx)) != DiscFormula(dj2, dq);
        SatelliteRegime = Regime(Block(l, SatelliteIdx));
        DoubleQuantumRegime = Regime(Block(l, DoubleQuantumIdx));
        // ON the exceptional points: a 2 × 2 Jordan block of the whole generator
        SatelliteEp = ReadEp(lSat, SatelliteEpPoint, 1);
        DoubleQuantumEp = ReadEp(lDq, DoubleQuantumEpPoint, 2);
        // the same-object control: no coupling and no splitting leave the satellite block a scalar, semisimple
        var (aj0, aj1, aj2, aq, aw) = DiabolicPoint;
        var lDia = Generator(aj0, aj1, aj2, aq, aw).L;
        var lambdaDia = new GaussianRational(-Spin32QuadrupolarMirrorsClaim.BrightCentre(aj0, aj1, aj2), -aw);
        var shiftedDia = Shift(lDia, lambdaDia);
        DiabolicControlSemisimple = IsScalar(Block(lDia, SatelliteIdx), out var scalarDia) && scalarDia == lambdaDia
                                    && Rank(shiftedDia) == Dim - 2 && Rank(Mul(shiftedDia, shiftedDia)) == Dim - 2;
        // the coalesced vector, read off the rank-one part N = [[a, b], [c, d]] as (b, −a): non-zero, an eigenvector of
        // the whole generator, and self-orthogonal (for two components that already puts as much on the satellites'
        // sum (1, 1) as on their mirror image, the partner (1, −1))
        var partSat = Sub(Block(lSat, SatelliteIdx), Scale(Identity(2), SatelliteEp.Lambda));
        var (va, vb) = (partSat[0, 1], -partSat[0, 0]);
        var coalesced = Zero(Levels);
        coalesced[0, 1] = va;
        coalesced[2, 3] = vb;
        bool coalescedSelfOrthogonal = !(va.IsZero && vb.IsZero) && Eigen(lSat, coalesced, SatelliteEp.Lambda)
                                       && (va * va + vb * vb).IsZero;
        GeneratorComplexSymmetric = Equal(Transpose(lSat), lSat);
        // the satellite part of I_+ starts at R₂′ = 3(J₀ + J₁) whatever the splitting: its Rayleigh quotient
        var (bj0, bj1, bj2, bq, bw) = BelowEpPoint;
        var (uj0, uj1, uj2, uq, uw) = AboveEpPoint;
        var lBelow = Generator(bj0, bj1, bj2, bq, bw).L;
        var lAbove = Generator(uj0, uj1, uj2, uq, uw).L;
        var satNorm = Trace(Mul(Adjoint(satellites), satellites));
        GaussianRational Rayleigh(Mat lx) => Trace(Mul(Adjoint(satellites), Apply(lx, satellites))) / satNorm;
        // the control: one satellite cell alone carries the full −c₂ and the splitting in its quotient
        var oneCell = Zero(Levels);
        oneCell[0, 1] = GaussianRational.One;
        var oneCellQuotient = Trace(Mul(Adjoint(oneCell), Apply(lSat, oneCell)));
        InitialSlopeSplittingFree = Rayleigh(l) == new GaussianRational(-3 * (j0 + j1), -omega)
                                    && Rayleigh(l0) == new GaussianRational(-3 * (j0 + j1), -omega)
                                    && Rayleigh(lSat) == new GaussianRational(-3 * (sj0 + sj1), -sw)
                                    && Rayleigh(lBelow) == new GaussianRational(-3 * (bj0 + bj1), -bw)
                                    && Rayleigh(lAbove) == new GaussianRational(-3 * (uj0 + uj1), -uw)
                                    && oneCellQuotient != new GaussianRational(-3 * (sj0 + sj1), -sw);
        // below the satellite EP: two real rates c₂ ∓ Ω, and the weights of the satellite part on them from the
        // block's spectral projectors (M − λ∓)/(λ± − λ∓), read on s = (1, 1)
        long c2Below = Spin32QuadrupolarMirrorsClaim.BrightCentre(bj0, bj1, bj2);
        long omBelow = ExactSqrt(9 * bj2 * bj2 - 4 * bq * bq);
        BelowEpOmega = omBelow;
        var lamSlow = new GaussianRational(-c2Below + omBelow, -bw);
        var lamFast = new GaussianRational(-c2Below - omBelow, -bw);
        var mBelow = Block(lBelow, SatelliteIdx);
        GaussianRational OnSum(Mat m) => m[0, 0] + m[0, 1] + m[1, 0] + m[1, 1];
        var wSlow = OnSum(Sub(mBelow, Scale(Identity(2), lamFast))) / (2 * (lamSlow - lamFast));
        var wFast = OnSum(Sub(mBelow, Scale(Identity(2), lamSlow))) / (2 * (lamFast - lamSlow));
        BelowEpWeights = new[] { wSlow, wFast };
        bool Simple(Mat lx, GaussianRational lam)
        {
            var sh = Shift(lx, lam);
            return Rank(sh) == Dim - 1 && Rank(Mul(sh, sh)) == Dim - 1;
        }
        BelowEpTwoRates = Simple(lBelow, lamSlow) && Simple(lBelow, lamFast) && wSlow == new BigRational(4, 3);
        // the self-orthogonality control: below the point the slow eigenvector of the same pair, a column of
        // M − λ_fast (whose columns span the slow eigenspace), is not self-orthogonal
        var (sa, sb) = (mBelow[0, 0] - lamFast, mBelow[1, 0]);
        EpVectorSelfOrthogonal = coalescedSelfOrthogonal && !(sa.IsZero && sb.IsZero) && !(sa * sa + sb * sb).IsZero;
        // above it: the two modes at ±Ω′ about the order frame, at neither nominal position ±2q; their separation 2Ω′
        // against D = 4q and the 4/5 longitudinal rate, read off L as the rate of ZI (the 1/5 rate, IZ's, in its place
        // fails). Mode frequencies, not the maxima of a spectrum.
        long c2Above = Spin32QuadrupolarMirrorsClaim.BrightCentre(uj0, uj1, uj2);
        long omAbove = ExactSqrt(4 * uq * uq - 9 * uj2 * uj2);
        AboveEpOmega = omAbove;
        BigRational dOmega = 2 * omAbove, dNom = 4 * uq;
        AboveEpFrequenciesPulled = Rank(Shift(lAbove, new GaussianRational(-c2Above, -uw + omAbove))) == Dim - 1
                                   && Rank(Shift(lAbove, new GaussianRational(-c2Above, -uw - omAbove))) == Dim - 1
                                   && Rank(Shift(lAbove, new GaussianRational(-c2Above, -uw + 2 * uq))) == Dim
                                   && Rank(Shift(lAbove, new GaussianRational(-c2Above, -uw - 2 * uq))) == Dim
                                   && RealRateOf(lAbove, zi) is { } r45 && dOmega * dOmega + r45 * r45 == dNom * dNom
                                   && RealRateOf(lAbove, iz2) is { } r15 && dOmega * dOmega + r15 * r15 != dNom * dNom;
        // the merging in the rates a measurement reads: where each pair is a Jordan block, the nominal separation 4q meets the
        // 4/5 rate for the satellites and the 1/5 rate for the double quantum, each read off L
        EpAtTheLongitudinalRates = SatelliteEp.IsJordanPair && DoubleQuantumEp.IsJordanPair
                                   && RealRateOf(lSat, zi) is { } sat45 && sat45 == 4 * sq
                                   && RealRateOf(lSat, iz2) is { } sat15 && sat15 != 4 * sq
                                   && RealRateOf(lDq, iz2) is { } dq15 && dq15 == 4 * dq
                                   && RealRateOf(lDq, zi) is { } dq45 && dq45 != 4 * dq;

        // gate 7: the dynamic frequency shift, the Lamb-shift Hamiltonian of the shifts K(m·ω₀) = Im Γ, odd in m
        var (hj0, hj1, hj2, hq, hw, hk1, hk2) = ShiftPoint;
        var (lK, lKj0) = Generator(hj0, hj1, hj2, hq, hw, hk1, hk2);
        var lNoK = Generator(hj0, hj1, hj2, hq, hw).L;
        ShiftKeepsMirror = Equal(Mul(Mul(_s, lK), _sInv), MirrorRight(lK, lKj0, hj1, hj2));
        bool MovesWhole(int[] idx, long by) => Closed(lK, idx)
            && Equal(Sub(Block(lK, idx), Block(lNoK, idx)), Scale(Identity(2), new GaussianRational(0, by)));
        int ct = Levels * 1 + 2;
        ShiftMovesPairsWhole = MovesWhole(SatelliteIdx, 6 * hk1) && MovesWhole(DoubleQuantumIdx, 6 * hk2)
                               && lK[ct, ct] - lNoK[ct, ct] == new GaussianRational(0, 6 * (hk2 - hk1));
        // the satellite exceptional point with the same shift: the ranks of the Jordan block at λ* + 6i·K₁
        var (ej0, ej1, ej2, eq, ew) = SatelliteEpPoint;
        var shK = Shift(Generator(ej0, ej1, ej2, eq, ew, hk1, hk2).L, SatelliteEp.Lambda + new GaussianRational(0, 6 * hk1));
        var shK2 = Mul(shK, shK);
        var lamK = SatelliteEp.Lambda + new GaussianRational(0, 6 * hk1);
        ShiftKeepsJordan = Rank(shK) == Dim - 1 && Rank(shK2) == Dim - 2 && Rank(Mul(shK2, shK)) == Dim - 2
                           && Rank(Shift(Generator(ej0, ej1, ej2, eq - 1, ew, hk1, hk2).L, lamK)) == Dim
                           && Rank(Shift(Generator(ej0, ej1, ej2, eq + 1, ew, hk1, hk2).L, lamK)) == Dim;
        // the shift read off the spin matrices against its closed form diag(−3(K₁+K₂), 3(K₁−K₂), −3(K₁−K₂), 3(K₁+K₂))
        // read at (1, 0), (0, 1) and the shift point: the term is linear in (K₁, K₂), so the two unit points pin it
        bool ClosedAt(long a1, long a2)
        {
            var read = LambShift(a1, a2);
            long[] closed = { -3 * (a1 + a2), 3 * (a1 - a2), -3 * (a1 - a2), 3 * (a1 + a2) };
            return Enumerable.Range(0, Levels).All(i => Enumerable.Range(0, Levels).All(j =>
                read[i, j] == (i == j ? new GaussianRational(closed[i], 0) : GaussianRational.Zero)));
        }
        ShiftClosedForm = ClosedAt(1, 0) && ClosedAt(0, 1) && ClosedAt(hk1, hk2);
        // the global centre with the shift: both routes against the rule, where it pairs and where it breaks
        ShiftCentreReadings = ShiftCentreCases.Select(c =>
        {
            var lc = Generator(c.J0, c.J1, c.J2, c.Q, c.Omega, c.K1, c.K2).L;
            long two = Spin32QuadrupolarMirrorsClaim.TwiceGlobalCentre(c.J0, c.J1, c.J2);
            return new ShiftCentreReading(c, Spin32QuadrupolarMirrorsClaim.HasGlobalCentreWithShift(c.J0, c.J1, c.J2, c.Q, c.Omega, c.K1, c.K2),
                Dim - Rank(lc), Dim - Rank(Add(lc, Scale(Identity(Dim), two))),
                CharpolySymmetric(lc, Enumerable.Range(0, Dim).ToArray(), two));
        }).ToArray();

        // gate 8: the spin-lock
        var half = new GaussianRational(new BigRational(1, 2), BigRational.Zero);
        bool Pattern(Mat m, (int, int)[] cells, GaussianRational v) => Enumerable.Range(0, Levels).All(i =>
            Enumerable.Range(0, Levels).All(j => m[i, j] == (cells.Contains((i, j)) ? v : GaussianRational.Zero)));
        LockSplitHolds = Pattern(_ixCt, new[] { (1, 2), (2, 1) }, GaussianRational.One)
                         && Pattern(_ixSat, new[] { (0, 1), (1, 0), (2, 3), (3, 2) }, half)
                         && Equal(Mul(_ixCt, _g0), Mul(_g0, _ixCt)) && Equal(Mul(_ixSat, _g0), Scale(Mul(_g0, _ixSat), -1));
        var (lockL, lockJ0) = Generator(hj0, hj1, hj2, hq, hw, hk1, hk2, Scale(_ixCt, LockNu));
        SelectiveLockKeepsMirror = Equal(Mul(Mul(_s, lockL), _sInv), MirrorRight(lockL, lockJ0, hj1, hj2))
                                   && Equal(Mul(_p, lockL), Mul(lockL, _p));
        var (satL, satJ0) = Generator(hj0, hj1, hj2, hq, hw, hk1, hk2, _ixSat);
        SatelliteLockBreaksMirror = !Equal(Mul(Mul(_s, satL), _sInv), MirrorRight(satL, satJ0, hj1, hj2))
                                    && !Equal(Mul(_p, satL), Mul(satL, _p));
        LockCentreReadings = LockCentreCases.Select(c =>
        {
            var lc = Generator(c.J0, c.J1, c.J2, c.Q, c.Omega, c.K1, c.K2, Add(Scale(_ixSat, c.A), Scale(_ixCt, c.B))).L;
            long two = Spin32QuadrupolarMirrorsClaim.TwiceGlobalCentre(c.J0, c.J1, c.J2);
            return new LockCentreReading(c,
                Spin32QuadrupolarMirrorsClaim.HasGlobalCentreUnderLock(c.J0, c.J1, c.J2, c.Q, c.Omega, c.K1, c.K2),
                Dim - Rank(lc), Dim - Rank(Add(lc, Scale(Identity(Dim), two))),
                CharpolySymmetric(lc, Enumerable.Range(0, Dim).ToArray(), two));
        }).ToArray();

        (TensorRowsHold, TensorBasesTurned, TripleQuantumPathHolds, MirrorOnTensorsHolds, OrderSplitsAStringClass) =
            ReadTensors(ip, im);
    }

    // ===== the irreducible tensors in strings =====

    private (bool Rows, bool Turned, bool Path, bool Mirror, bool Order) ReadTensors(Root3 ip, Root3 im)
    {
        Mat S2(Mat a, Mat b) => Kron(a, b);
        Mat Lin(params (GaussianRational C, Mat M)[] terms) =>
            terms.Aggregate(Zero(Levels), (acc, t) => Add(acc, Scale(t.M, t.C)));
        GaussianRational G(long v) => new(v, 0);
        Root3 R(Mat a) => new(a, Zero(Levels));
        bool Zero3(Root3 u) => IsZero(u.A) && IsZero(u.B);
        Root3 Comm(Root3 u, Root3 v) => u * v + (v * u).Scale(-1);

        var (ii, zi, iz, zz) = (S2(PI, PI), S2(PZ, PI), S2(PI, PZ), S2(PZ, PZ));
        Mat Sat(int sg) => Lin((G(1), S2(PI, PX)), (new GaussianRational(0, sg), S2(PI, PY)));
        Mat Ct(int sg) => Lin((G(1), S2(PX, PX)), (G(1), S2(PY, PY)), (new GaussianRational(0, sg), S2(PY, PX)), (new GaussianRational(0, -sg), S2(PX, PY)));
        var rows = new Dictionary<(int L, int P), Root3>
        {
            [(0, 0)] = R(ii), [(1, 0)] = R(Lin((G(2), zi), (G(1), iz))),
            [(2, 0)] = R(zz), [(3, 0)] = R(Lin((G(2), iz), (G(-1), zi))),
        };
        foreach (int sg in new[] { 1, -1 })
        {
            rows[(1, sg)] = new Root3(Ct(sg), Sat(sg));                                // (XX + YY …) + √3·(IX ± i·IY)
            rows[(3, sg)] = new Root3(Scale(Sat(sg), 2), Scale(Ct(sg), -1));            // 2·(IX ± i·IY) − √3·(XX + YY …)
            rows[(2, sg)] = R(Lin((G(1), S2(PZ, PX)), (new GaussianRational(0, sg), S2(PZ, PY))));
            rows[(2, 2 * sg)] = R(Lin((G(1), S2(PX, PI)), (new GaussianRational(0, sg), S2(PY, PI))));
            rows[(3, 2 * sg)] = R(Lin((G(1), S2(PX, PZ)), (new GaussianRational(0, sg), S2(PY, PZ))));
            rows[(3, 3 * sg)] = R(Lin((G(1), S2(PX, PX)), (G(-1), S2(PY, PY)), (new GaussianRational(0, sg), S2(PX, PY)), (new GaussianRational(0, sg), S2(PY, PX))));
        }

        var ix = (ip + im).Scale(new BigRational(1, 2));
        var iyTimesI = (ip + im.Scale(-1)).Scale(new BigRational(1, 2));          // i·I_y = (I₊ − I₋)/2
        var izR = R(_iz);
        bool Carries(Root3 m, int l, int p)
        {
            // Σ[Iₐ, [Iₐ, M]] with [I_y, [I_y, M]] = −[i·I_y, [i·I_y, M]]
            var cas = Comm(ix, Comm(ix, m)) + Comm(iyTimesI, Comm(iyTimesI, m)).Scale(-1) + Comm(izR, Comm(izR, m));
            return Zero3(cas + m.Scale(-l * (l + 1))) && Zero3(Comm(izR, m) + m.Scale(-p));
        }
        bool rowsHold = rows.Count == 16 && rows.All(kv => Carries(kv.Value, kv.Key.L, kv.Key.P))
                        && !Carries(rows[(3, 1)], 1, 1);

        var (j0, j1, j2, omega) = TensorPoint;
        var l = Generator(j0, j1, j2, 0, omega).L;
        long c1 = Spin32QuadrupolarMirrorsClaim.DarkCentre(j1, j2);
        GaussianRational Rate(long re, long freq) => new(-re, -freq);
        var eigenRates = new Dictionary<(int, int), long>
        {
            [(0, 0)] = 0, [(2, 0)] = 6 * (j1 + j2), [(2, 1)] = 3 * (j0 + j1 + 2 * j2), [(2, -1)] = 3 * (j0 + j1 + 2 * j2),
            [(2, 2)] = 3 * (j0 + 2 * j1 + j2), [(2, -2)] = 3 * (j0 + 2 * j1 + j2),
            [(3, 2)] = 3 * (j0 + j2), [(3, -2)] = 3 * (j0 + j2), [(3, 3)] = c1, [(3, -3)] = c1,
        };
        bool ten = eigenRates.All(kv => Eigen(l, rows[kv.Key].A, Rate(kv.Value, kv.Key.Item2 * omega)));
        GaussianRational Norm(Mat a) => Trace(Mul(Adjoint(a), a));
        // a row's part on a set of strings, read off the row itself: Σ Tr(s·M)/4 · s on both of its ℚ(√3) components
        Mat OnStrings(Mat m, Mat[] strings) => strings.Aggregate(Zero(Levels),
            (acc, st) => Add(acc, Scale(st, Trace(Mul(st, m)) / G(4))));
        Root3 Part(Root3 r, Mat[] strings) => new(OnStrings(r.A, strings), OnStrings(r.B, strings));
        // ⟨u, v⟩ = rational + √3·root, with ⟨A + √3B, C + √3D⟩ = ⟨A, C⟩ + 3⟨B, D⟩ + √3(⟨A, D⟩ + ⟨B, C⟩)
        (GaussianRational Rational, GaussianRational Root) Inner(Root3 u, Root3 v) =>
            (Trace(Mul(Adjoint(u.A), v.A)) + G(3) * Trace(Mul(Adjoint(u.B), v.B)),
             Trace(Mul(Adjoint(u.A), v.B)) + Trace(Mul(Adjoint(u.B), v.A)));
        bool EigenPart(Root3 part, GaussianRational lambda) => !Zero3(part)
            && (IsZero(part.A) || Eigen(l, part.A, lambda)) && (IsZero(part.B) || Eigen(l, part.B, lambda));
        var (zis, izs) = (new[] { zi }, new[] { iz });
        var sats = new[] { S2(PI, PX), S2(PI, PY) };
        var cts = new[] { S2(PX, PX), S2(PX, PY), S2(PY, PX), S2(PY, PY) };
        // each turned row: its part on class a (weight wa) and on class b, both eigenoperators at two different rates,
        // the two parts summing to the row
        bool TurnedRow((int L, int P) key, Mat[] a, long rateA, Mat[] b, long rateB, BigRational wa)
        {
            var row = rows[key];
            var (pa, pb) = (Part(row, a), Part(row, b));
            var (na, nb, nr) = (Inner(pa, pa), Inner(pb, pb), Inner(row, row));
            return Zero3(row + (pa + pb).Scale(-1)) && rateA != rateB
                   && EigenPart(pa, Rate(rateA, key.P * omega)) && EigenPart(pb, Rate(rateB, key.P * omega))
                   && na.Root.IsZero && nb.Root.IsZero && nr.Root.IsZero
                   && na.Rational / nr.Rational == new GaussianRational(wa, BigRational.Zero)
                   && nb.Rational / nr.Rational == new GaussianRational(1 - wa, BigRational.Zero);
        }
        long satRate = 3 * (j0 + j1);
        bool sixTurned = TurnedRow((1, 0), zis, 6 * j2, izs, 6 * j1, new BigRational(4, 5))
                         && TurnedRow((3, 0), zis, 6 * j2, izs, 6 * j1, new BigRational(1, 5))
                         && new[] { 1, -1 }.All(sg => TurnedRow((1, sg), sats, satRate, cts, c1, new BigRational(3, 5))
                                                      && TurnedRow((3, sg), sats, satRate, cts, c1, new BigRational(2, 5)));
        bool turned = ten && sixTurned;

        // the rank-3 row meets the satellite part of the rank-1 row at order 1 in a pure √3 overlap whose square, over
        // the two rows' squared norms, is 6/25; the CT part's overlap is its negative (the rows are orthogonal), so with
        // the parts' two rates (TensorBasesTurned) the rank-3 content of e^{Lt}·T₁,₁ is ±(√6/5)·(e^{−3(J₀+J₁)t} − e^{−c₁t})
        var row11 = rows[(1, 1)];
        var row31 = rows[(3, 1)];
        var oSat = Inner(row31, Part(row11, sats));
        var (n11, n31) = (Inner(row11, row11), Inner(row31, row31));
        bool path = oSat.Rational.IsZero && n11.Root.IsZero && n31.Root.IsZero
                    && G(3) * oSat.Root * oSat.Root / (n11.Rational * n31.Rational) == G(6) / G(25)
                    && Inner(row31, row11).Rational.IsZero && Inner(row31, row11).Root.IsZero;

        Mat Mir(Mat x) => Mul(_g0, Transpose(x));
        bool Onto(Mat x, Mat y) => ProportionalFactor(Mir(x), y) is not null;
        var (t10, t30) = (rows[(1, 0)].A, rows[(3, 0)].A);
        bool mirror = Equal(Mir(ii), zz) && Equal(Mir(zz), ii)
                      && new[] { 1, -1 }.All(sg => Onto(rows[(2, 2 * sg)].A, rows[(3, -2 * sg)].A)
                                                  && Onto(rows[(3, 2 * sg)].A, rows[(2, -2 * sg)].A)
                                                  && Onto(rows[(3, 3 * sg)].A, rows[(3, -3 * sg)].A))
                      && Norm(t10) == Norm(t30)
                      && Equal(Mir(t10), Add(Scale(t10, new BigRational(4, 5)), Scale(t30, new BigRational(3, 5))))
                      && Equal(Mir(t30), Add(Scale(t10, new BigRational(3, 5)), Scale(t30, new BigRational(-4, 5))));

        var plus = Add(S2(PX, PX), S2(PY, PY));
        var minus = Sub(S2(PX, PX), S2(PY, PY));
        bool OnCells(Mat m, params (int, int)[] cells) => Enumerable.Range(0, Levels).All(a =>
            Enumerable.Range(0, Levels).All(b => m[a, b].IsZero != cells.Contains((a, b))));
        var l0 = Generator(j0, j1, j2, 0, 0).L;
        // an order filter projects onto a span of strings exactly where the strings inside its cells fill them:
        // at p = 0 and |p| = 2, not at |p| = 1 (four strings, six cells) and |p| = 3 (none, two cells)
        var paulis = new[] { PI, PX, PY, PZ };
        var all = paulis.SelectMany(a => paulis.Select(b => S2(a, b))).ToArray();
        int OrderOf(int a, int b) => Math.Abs(b - a);                     // |m_a − m_b| with m_k = 3/2 − k
        (int Cells, int Strings) Fill(int absP) => (
            Enumerable.Range(0, Levels).Sum(a => Enumerable.Range(0, Levels).Count(b => OrderOf(a, b) == absP)),
            all.Count(m => Enumerable.Range(0, Levels).All(a => Enumerable.Range(0, Levels)
                .All(b => m[a, b].IsZero || OrderOf(a, b) == absP))));
        bool order = OnCells(plus, (1, 2), (2, 1)) && OnCells(minus, (0, 3), (3, 0))
                     && Eigen(l0, plus, Rate(c1, 0)) && Eigen(l0, minus, Rate(c1, 0))
                     && Fill(0) == (4, 4) && Fill(2) == (4, 4) && Fill(1) == (6, 4) && Fill(3) == (2, 0);
        return (rowsHold, turned, path, mirror, order);
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

    // ===== the bright pairs' exceptional point =====

    private static readonly int[] SatelliteIdx = { 1, 11 };             // |3/2⟩⟨1/2|, |−1/2⟩⟨−3/2| (vec index 4i + j)
    private static readonly int[] DoubleQuantumIdx = { 2, 7 };          // |3/2⟩⟨−1/2|, |1/2⟩⟨−3/2|
    private static readonly int[] SatelliteWithCentralIdx = { 1, 6 };   // the control: |3/2⟩⟨1/2| with |1/2⟩⟨−1/2|

    /// <summary>Whether the cells span an L-invariant subspace: no column of a cell has an entry outside the cells.</summary>
    private static bool Closed(Mat l, int[] idx)
    {
        foreach (int c in idx)
            for (int r = 0; r < Dim; r++)
                if (Array.IndexOf(idx, r) < 0 && !l[r, c].IsZero) return false;
        return true;
    }

    private static Mat Block(Mat l, int[] idx)
    {
        var b = Zero(idx.Length);
        for (int a = 0; a < idx.Length; a++)
            for (int c = 0; c < idx.Length; c++) b[a, c] = l[idx[a], idx[c]];
        return b;
    }

    /// <summary>(tr M)² − 4·det M of a 2 × 2 block, zero exactly where its two eigenvalues meet.</summary>
    private static GaussianRational Disc2(Mat m)
    {
        var tr = m[0, 0] + m[1, 1];
        return tr * tr - 4 * (m[0, 0] * m[1, 1] - m[0, 1] * m[1, 0]);
    }

    private static Mat Shift(Mat l, GaussianRational lambda) => Sub(l, Scale(Identity(l.GetLength(0)), lambda));

    /// <summary>The Jordan reading at any point, for the tests that feed it points where the pair is NOT at its
    /// exceptional point and watch it fail.</summary>
    internal EpReading ReadEpAt((long J0, long J1, long J2, long Q, long Omega) pt, int order) =>
        ReadEp(Generator(pt.J0, pt.J1, pt.J2, pt.Q, pt.Omega).L, pt, order);

    private EpReading ReadEp(Mat l, (long J0, long J1, long J2, long Q, long Omega) pt, int order)
    {
        var lambda = new GaussianRational(-Spin32QuadrupolarMirrorsClaim.BrightCentre(pt.J0, pt.J1, pt.J2), -order * pt.Omega);
        var part = Shift(Block(l, order == 1 ? SatelliteIdx : DoubleQuantumIdx), lambda);
        var a = Shift(l, lambda);
        var a2 = Mul(a, a);
        int RankOffBy(long step) => Rank(Shift(Generator(pt.J0, pt.J1, pt.J2, pt.Q + step, pt.Omega).L, lambda));
        return new EpReading(pt.J0, pt.J1, pt.J2, pt.Q, pt.Omega, order, lambda, !IsZero(part), IsZero(Mul(part, part)),
            Rank(a), Rank(a2), Rank(Mul(a2, a)), RankOffBy(-1), RankOffBy(1));
    }

    /// <summary>r with L(op) = −r·op for a real r, read off the first non-zero entry of op and checked on every
    /// entry; null when op is no eigenoperator at a real rate.</summary>
    private static BigRational? RealRateOf(Mat super, Mat op)
    {
        var image = Apply(super, op);
        for (int i = 0; i < Levels; i++)
            for (int j = 0; j < Levels; j++)
                if (!op[i, j].IsZero)
                {
                    var lambda = image[i, j] / op[i, j];
                    return lambda.Im.IsZero && Equal(image, Scale(op, lambda)) ? -lambda.Re : null;
                }
        return null;
    }

    /// <summary>√v for a perfect square v ≥ 0; the gate-6 points are chosen so, anything else throws.</summary>
    private static long ExactSqrt(long v)
    {
        if (v < 0) throw new InvalidOperationException($"{v} < 0 has no real square root");
        var r = (long)ISqrt(v);
        return r * r == v ? r : throw new InvalidOperationException($"{v} is not a perfect square");
    }

    /// <summary>Where a bright pair stands, from its discriminant read off the generator (real for real inputs).</summary>
    private static string Regime(Mat block)
    {
        var disc = Disc2(block);
        bool coupled = !block[0, 1].IsZero;
        return !disc.Im.IsZero ? "with a complex discriminant (NOT expected for real inputs)"
            : disc.Re.Sign > 0 ? "below its exceptional point (two real rates about the bright centre)"
            : disc.Re.Sign < 0 ? (coupled ? "past its exceptional point (two mode frequencies on one width c₂)" : "uncoupled (two detuned cells, no exceptional point at this density zero)")
            : coupled ? "on its exceptional point (a 2 × 2 Jordan block)" : "a scalar (no coupling and no splitting: diabolic)";
    }

    /// <summary>The exact generator at a point, for a caller that reads it with another instrument (the tests read
    /// the exceptional point once more with <c>EpCharacter</c>, numerically).</summary>
    public GaussianRational[,] GeneratorMatrix(long j0, long j1, long j2, long q, long omega) => Generator(j0, j1, j2, q, omega).L;

    // ===== the generator =====

    /// <summary>L and its J₀ part L_{J₀}, as 16 × 16 superoperators on vec index 4i + j of |i⟩⟨j|.</summary>
    private (Mat L, Mat LJ0) Generator(long j0, long j1, long j2, long q, long omega, long k1 = 0, long k2 = 0, Mat? extra = null)
    {
        var h = Add(Scale(_g0, q), Scale(_iz, omega));
        if (k1 != 0 || k2 != 0) h = Add(h, LambShift(k1, k2));
        if (extra is not null) h = Add(h, extra);
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

    /// <summary>The Lamb-shift Hamiltonian K₁(T₂,₁†T₂,₁ − T₂,₋₁†T₂,₋₁) + K₂(T₂,₂†T₂,₂ − T₂,₋₂†T₂,₋₂), built from the
    /// spin matrices (T = √3·R, so T†T = 3·R†R).</summary>
    private Mat LambShift(long k1, long k2)
    {
        Mat TdT(int m) => Scale(Mul(Adjoint(_r[m]), _r[m]), 3);
        return Add(Scale(Sub(TdT(1), TdT(-1)), k1), Scale(Sub(TdT(2), TdT(-2)), k2));
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
        $"{Verdict(HandshakeHolds, "holds (the textbook parts of I_z and I_+ are eigenoperators at the textbook rates)", "FAILS")}; " +
        $"the bright pairs' exceptional point {Verdict(ExceptionalPointHolds, "holds (the satellites and the double quantum each a 2 × 2 Jordan block on Re λ = −c₂, at 2q = 3J₂ and at 2q = 3J₁ respectively, the density non-zero; the scalar control semisimple)", "FAILS")}, " +
        $"and at the user point the satellite pair is {SatelliteRegime}, the double-quantum pair {DoubleQuantumRegime}; " +
        $"with the dynamic frequency shift {Verdict(DynamicShiftHolds, "both mirrors and the exceptional point hold and the global centre follows the shifted rule", "FAILS")}; " +
        $"under a spin-lock {Verdict(SpinLockHolds, "the central-transition part keeps both mirrors, the satellite part breaks the identity and mixes the sectors, and the global centre follows the lock's rule", "FAILS")}; " +
        $"the irreducible tensors in strings {Verdict(TensorDictionaryHolds, "hold (sixteen rows of rank and order; ten eigenoperators, six turned between two string classes)", "FAIL")}.";

    private static string Show(EpReading e, string name, string condition) =>
        $"the {name} at (J₀, J₁, J₂, q, ω) = ({e.J0}, {e.J1}, {e.J2}, {e.Q}, {e.Omega}), {condition}: λ* = {e.Lambda}, the block's part N " +
        $"{(e.PartNonZero ? "non-zero" : "ZERO")} with N² {(e.PartSquaresToZero ? "= 0" : "≠ 0")}, geometric multiplicity {e.Geometric} and " +
        $"algebraic {e.Algebraic} in L (ranks {e.Rank1}, {e.Rank2}, {e.Rank3}), and one step off in q no eigenvalue there (ranks " +
        $"{e.RankStepDown}, {e.RankStepUp}): {(e.IsJordanPair ? "a 2 × 2 Jordan block" : "NOT a Jordan pair")}";

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
                displayName: "the bright pairs' exceptional point",
                summary: $"the satellites (|3/2⟩⟨1/2|, |−1/2⟩⟨−3/2|) and the double quantum (|3/2⟩⟨−1/2|, |1/2⟩⟨−3/2|) " +
                         $"{Verdict(BrightPairsClosed, "are closed 2 × 2 blocks (a satellite with the central transition is not)", "are NOT closed as stated")}, " +
                         $"their traces {Verdict(BrightPairTracesHold, "are −2c₂ − 2ip·ω at the user point, the two rates of a pair summing to 2c₂ there (at every splitting in the verifier's symbolic blocks)", "are NOT −2c₂ − 2ip·ω")}, " +
                         $"their discriminants {Verdict(BrightPairDiscriminantsHold, "are 4(9J₂² − 4q²) and 4(9J₁² − 4q²) (the swapped assignment fails)", "are NOT as stated")}; " +
                         $"at the user point the satellite pair is {SatelliteRegime}, the double-quantum pair {DoubleQuantumRegime}; " +
                         $"ON the exceptional point: {Show(SatelliteEp, "satellites", $"2q = 3J₂ {(2 * SatelliteEp.Q == 3 * SatelliteEp.J2 ? "there" : "NOT there")}")}; " +
                         $"{Show(DoubleQuantumEp, "double quantum", $"2q = 3J₁ {(2 * DoubleQuantumEp.Q == 3 * DoubleQuantumEp.J1 ? "there" : "NOT there")}")}; " +
                         $"the control at {DiabolicPoint}, J₂ = q = 0, {Verdict(DiabolicControlSemisimple, "is a scalar block and semisimple in L (ranks 14, 14)", "is NOT semisimple as stated")}; " +
                         $"the coalesced vector {Verdict(EpVectorSelfOrthogonal, "is self-orthogonal, vᵀv = 0 (with two components: as much on the satellites' sum as on their mirror image), the pair's slow eigenvector below the point not", "is NOT as stated")}" +
                         $"{(GeneratorComplexSymmetric ? " (L is complex symmetric in the cell basis)" : " (L is NOT complex symmetric in the cell basis)")}; " +
                         $"the satellite part of I_+ {Verdict(InitialSlopeSplittingFree, "starts at R₂′ = 3(J₀ + J₁) at each of the five splittings read (every splitting in the verifier's symbolic series), one satellite cell alone not", "does NOT start at R₂′")}, " +
                         $"below the exceptional point at {BelowEpPoint} {Verdict(BelowEpTwoRates, $"decays at the two rates c₂ ∓ {BelowEpOmega} with the weights {string.Join(" and ", BelowEpWeights)}", "does NOT decay as stated")}, " +
                         $"above it at {AboveEpPoint} {Verdict(AboveEpFrequenciesPulled, $"has its two modes at ±{AboveEpOmega} about the order frame and at neither ±2q = ±{2 * AboveEpPoint.Q}, with D_ω² + R₁(4/5)² = D² for D = 4q and R₁(4/5) the rate of ZI read off L, the 1/5 rate failing in its place (mode frequencies, not line positions)", "is NOT as stated")}; " +
                         $"the merging {Verdict(EpAtTheLongitudinalRates, "sits where the nominal separation 4q equals the 4/5 longitudinal rate (satellites) or the 1/5 one (double quantum), both read off L", "is NOT at the longitudinal rates")}.",
                provenance: NodeProvenance.Live);

            yield return new InspectableNode(
                displayName: "the dynamic frequency shift",
                summary: $"the Lamb-shift Hamiltonian of the shifts K(m·ω₀) = Im Γ, odd in m, −3K₂·ZI − 3K₁·IZ, added at {ShiftPoint}: " +
                         $"the mirror identity {Verdict(ShiftKeepsMirror, "holds", "FAILS")}; the bright pairs " +
                         $"{Verdict(ShiftMovesPairsWhole, $"move as wholes, the satellites by 6i·K₁ = {6 * ShiftPoint.K1}i and the double quantum by 6i·K₂ = {6 * ShiftPoint.K2}i, the central transition by 6i(K₂ − K₁) = {6 * (ShiftPoint.K2 - ShiftPoint.K1)}i, so the satellites sit 6(2K₁ − K₂) = {6 * (2 * ShiftPoint.K1 - ShiftPoint.K2)} off it", "do NOT move as wholes")}; " +
                         $"the satellite exceptional point {Verdict(ShiftKeepsJordan, "stays a 2 × 2 Jordan block, at λ* + 6i·K₁ (ranks 15, 14, 14)", "is NOT a Jordan block there")}; " +
                         $"the global centre with the shift (J₀ = 0, or q = 0 and J₁ = 0 with ω = 6K₁, or q = 0 and J₂ = 0 with ω = 3K₂) " +
                         $"{Verdict(ShiftCentreReadings.All(r => r.Agrees), $"agrees with dim ker L = dim ker(L + 2s) and with χ(−x − 2s) = χ(x) at all {ShiftCentreReadings.Count} points ({ShiftCentreReadings.Count(r => r.Rule)} pairing, {ShiftCentreReadings.Count(r => !r.Rule)} breaking)", "DISAGREES with a route")}.",
                provenance: NodeProvenance.Live);

            yield return new InspectableNode(
                displayName: "a spin-lock",
                summary: $"ω₁·I_x in the rotating frame, split by the spin matrices: {Verdict(LockSplitHolds, "the central transition (XX + YY)/2 commutes with Γ₀, the satellites (√3/2)·IX anticommute", "the split is NOT as stated")}; " +
                         $"at {ShiftPoint} the selective lock ω₁ = {LockNu} {Verdict(SelectiveLockKeepsMirror, "keeps the mirror identity and the sectors", "BREAKS the mirror identity")}, the satellite part " +
                         $"{Verdict(SatelliteLockBreaksMirror, "breaks the identity and mixes the sectors", "does NOT break it")}; the global centre under a lock (only J₁, q = 0, ω = 3K₂; read with integer mixtures a:b of the two parts in place of √3:1, whose lit element is ZZ − (a/b)·XI) " +
                         $"{Verdict(LockCentreReadings.All(r => r.Agrees), $"agrees with both routes at all {LockCentreReadings.Count} points ({LockCentreReadings.Count(r => r.Rule)} pairing, {LockCentreReadings.Count(r => !r.Rule)} breaking)", "DISAGREES with a route")}.",
                provenance: NodeProvenance.Live);

            yield return new InspectableNode(
                displayName: "the irreducible tensors in strings",
                summary: $"the page's sixteen string combinations {Verdict(TensorRowsHold, "carry rank l and order p (Casimir l(l + 1), [I_z, ·] = p, from the spin matrices in ℚ(√3)), the rank-3 row at order 1 failing as rank 1", "do NOT carry their rank and order")}; " +
                         $"at (J₀, J₁, J₂, ω) = {TensorPoint}, ω_Q = 0, {Verdict(TensorBasesTurned, "rank 0, rank 2 and rank 3 at |p| = 2, 3 are eigenoperators, and the six rank-1 and rank-3 rows at |p| ≤ 1 are turned between two string classes at two rates, 4/5 : 1/5 (ZI, IZ) and 3/5 : 2/5 (satellites, central transition) and the reverse", "the two bases are NOT as stated")}; " +
                         $"the rank-3 content of e^(Lt)·T₁,₁ (off resonance times e^(−iωt)) {Verdict(TripleQuantumPathHolds, "is ±(√6/5)(e^(−3(J₀ + J₁)t) − e^(−c₁t)), the squared overlap 6/25 a consistency reading of the weights", "is NOT that")}; " +
                         $"S(X) = Γ₀·Xᵀ {Verdict(MirrorOnTensorsHolds, "carries T₀,₀ ↔ T₂,₀, T₂,±₂ onto T₃,∓₂ and back, T₃,±₃ onto T₃,∓₃ (up to factors; the signs are the verifier's), and reflects the order-0 rows of rank 1 and 3 by [[4/5, 3/5], [3/5, −4/5]]", "does NOT act as stated")}; " +
                         $"XX + YY and XX − YY {Verdict(OrderSplitsAStringClass, "sit at orders ±1 and ±3 and decay in one string class at c₁", "do NOT split as stated")}.",
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
