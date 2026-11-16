using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Numerics;
using RCPsiSquared.Core.Quadrupolar;
using K = RCPsiSquared.Core.Numerics.Cyclotomic24;
using KM = RCPsiSquared.Core.Numerics.Cyclotomic24Matrix;
using GM = RCPsiSquared.Core.Numerics.GaussianRationalMatrix;
using Mat = RCPsiSquared.Core.Numerics.Cyclotomic24[,];

namespace RCPsiSquared.Diagnostics.Foundation;

/// <summary>The experiment in strings, live: the section of that name in
/// <c>experiments/THE_SODIUM_NUCLEUS_TWO_MIRRORS.md</c> (F167), read for one spin 3/2 in Souza et al.'s labelling
/// |3/2⟩, |1/2⟩, |−1/2⟩, |−3/2⟩ = |00⟩, |01⟩, |10⟩, |11⟩ (the first bit high), exactly over ℚ(ζ₂₄)
/// (<see cref="Cyclotomic24"/>: i, √2, √3 and every phase that is a multiple of 15°), the spin operators built from the
/// spin formula first (the turn about z from I_z) and the generator taken from the shared builder in Core. A child of
/// <see cref="SodiumTwoMirrorsWitness"/>; the page's verifier (<c>simulations/sodium_two_mirrors.py experiment</c>)
/// holds the twelve statements these gates read, in sympy (φ a symbol, time functions, the Theodolite's echo).
///
/// <para>How "every φ" is read without a symbol: each identity in a pulse's phase φ is a Laurent polynomial in
/// e^{iφ} (the turn about z alone in e^{iφ/2}, read in ψ = φ/2) whose exponents span at most 13 values, and the
/// witness checks it at 24 points of the circle (φ at the multiples of 15°, for the turn about z ψ there and so φ at the
/// multiples of 30°), more points than the span, so it holds at every φ.
/// The turn's global phase e^{−3iφ/2} is dropped where only a conjugation by it enters. How "every splitting and
/// offset" is read: the shared builder is linear in (J₀, J₁, J₂, q, ω) by construction (H = q·ZZ + ω·I_z, each jump's
/// dissipator scaled by its J), guarded at two generic points, so a linear statement read at the five unit points holds
/// at every point. Those statements involve only the diagonal H and the jumps, so they do not test the builder's vec
/// order or the sign of H; the parent witness and Core's tests do. Time-dependent forms are read as amplitudes on the rates, not as
/// functions of t. Ranks over the field are read through the regular representation (<see cref="Cyclotomic24Matrix.Rank"/>).</para>
///
/// <para>Twelve gates, as the verifier's twelve claims: the labelling, the coherence order, unitality and the thermal
/// source, the inversion recovery, the turn about z, the phase cycle, the 180° pulse, the echo, the 180° pulse in the
/// page's bits, the 90° pulse, the receiver, the pools. One of the verifier's readings stays there: the Theodolite's
/// letter layer reading the echo (its Python twin, <c>cube_letters.py</c>; Diagnostics reaches neither MirrorWorld nor the
/// Python cube tools), which this witness replaces by superoperator identities on the strings. The pools' heights and arrow are named as
/// <c>LetterCube</c>'s reading of a one-site jump on the page and run by neither. The verifier's phases π/5 and π/7 lie
/// in no field here; the witness reads the 24 phases instead.</para>
///
/// <para>Typed home <see cref="Core.Symmetry.Spin32QuadrupolarMirrorsClaim"/>.</para></summary>
public sealed class SodiumExperimentInStringsWitness : IInspectable
{
    private const int Levels = 4;
    private const int Phases = K.Order;                      // the 24 phases k·15°

    private static readonly string[] Labels =
        "IXYZ".SelectMany(a => "IXYZ".Select(b => $"{a}{b}")).ToArray();

    // ===== gate 1: the labelling =====
    /// <summary>The page's bits a (outer/inner) and b (the sign of m) are Souza et al.'s through a CNOT and a swap,
    /// b = bit 1, a = bit 1 XOR bit 2, a written first: I_z = Z_b + Z_a·Z_b/2 and ZZ (∝ T₂₀) goes to Z_a.</summary>
    public bool LabellingReproducesIz { get; }
    /// <summary>Every string goes to one string with a sign ±1 (as under any relabelling of the four levels), the seven
    /// images the page lists among them.</summary>
    public bool LabellingSignedImages { get; }
    /// <summary>How many of the 15 strings besides II land on another k_Z (8).</summary>
    public int StringsMovedInKz { get; }
    /// <summary>The control: with a and b assigned the other way round I_z is not reproduced.</summary>
    public bool WrongAssignmentFails { get; }
    public bool LabellingHolds => LabellingReproducesIz && LabellingSignedImages && StringsMovedInKz == 8 && WrongAssignmentFails;

    // ===== gate 2: the coherence order =====
    /// <summary>On all sixteen cells |m⟩⟨m′|: p = 2v₁ + v₂ (p read off I_z), every string of the cell at the letter cube's
    /// k_Z = |v₁| + |v₂|, and the bits' popcount difference (bra minus ket) v₁ + v₂ (an identity of z = 1 − 2·bit), with
    /// v_i = (z_i − z_i′)/2.</summary>
    public bool OrderIsWeightedBitOrders { get; }
    /// <summary>The control: the weights the other way round, p = v₁ + 2v₂, fail.</summary>
    public bool SwappedWeightsFail { get; }
    /// <summary>In the page's bits the bits' orders do not fix p: the v carrying two orders are exactly (±1, 0) and (0, ±1).</summary>
    public bool PageBitsDoNotFixOrder { get; }
    public bool CoherenceOrderHolds => OrderIsWeightedBitOrders && SwappedWeightsFail && PageBitsDoNotFixOrder;

    // ===== gate 3: unitality and the thermal source =====
    /// <summary>The generator is linear in (J₀, J₁, J₂, q, ω): at two generic points it is the sum of the unit points' (a
    /// guard on the builder, which is linear by construction).</summary>
    public bool GeneratorLinear { get; }
    /// <summary>L(II) = 0 and Tr L(S) = 0 for every string, at the five unit points, so at every point.</summary>
    public bool UnitalAndTracePreserving { get; }
    /// <summary>The control: amplitude damping on bit 1, read through the same superoperator route, keeps every string's
    /// trace and moves the identity.</summary>
    public bool AmplitudeDampingMovesIdentity { get; }
    /// <summary>−L(I_z) = 6J₂·ZI + 3J₁·IZ at every point, so a relaxation toward (1 + ε·I_z)/4 has the source
    /// (ε/4)(6J₂·ZI + 3J₁·IZ), F120's form.</summary>
    public bool ThermalSourceHolds { get; }
    public bool UnitalityHolds => GeneratorLinear && UnitalAndTracePreserving && AmplitudeDampingMovesIdentity && ThermalSourceHolds;

    // ===== gate 4: the inversion recovery =====
    /// <summary>ZI and IZ are eigenoperators at −6J₂ and −6J₁ at every point.</summary>
    public bool LongitudinalEigenoperators { get; }
    /// <summary>The rank-3 population T₃₀ ∝ 2·IZ − ZI in units of the equilibrium magnetization's rank-1 amplitude, as
    /// amplitudes on (1, e^{−6J₂t}, e^{−6J₁t}): from −I_z toward +I_z, and for the deviation −I_z (a saturation
    /// recovery, or the same inverted start at the model's infinite temperature).</summary>
    public IReadOnlyList<BigRational> InversionAmplitudes { get; }
    public IReadOnlyList<BigRational> SaturationAmplitudes { get; }
    public bool InversionRecoveryHolds { get; }

    // ===== gate 5: the turn about z =====
    /// <summary>e^{−iφI_z} = e^{−iφZ} ⊗ e^{−iφZ/2} at every φ, the cell |m⟩⟨m′| turned by e^{−i(2v₁ + v₂)φ}, and at 90°
    /// the Weyl clock matrix Z ⊗ diag(1, i) times e^{−3iπ/4}; the control (weights swapped) fails.</summary>
    public bool TurnFactorsInSouzaBits { get; }
    /// <summary>In the page's bits the turn is e^{−iφZ_b}·(cos(φ/2) − i·sin(φ/2)·Z_a·Z_b) at every φ.</summary>
    public bool TurnFormInPageBits { get; }
    /// <summary>Its operator Schmidt rank in the page's bits at φ = 30°·k, k = 0 … 23.</summary>
    public IReadOnlyList<int> TurnRanksInPageBits { get; }
    public bool TurnHolds => TurnFactorsInSouzaBits && TurnFormInPageBits
        && TurnRanksInPageBits.Select((r, k) => r == (k % 6 == 0 ? 1 : 2)).All(b => b);

    // ===== gate 6: the phase cycle =====
    /// <summary>n uniform steps 2πk/n weighted by e^{ip₀φ_k} keep the pathways with p′ − p ≡ p₀ (mod n) and remove the
    /// rest, n = 1 … 8, read in ℤ[ζ_n] (<see cref="CyclotomicRing"/>).</summary>
    public bool CycleKeepsResidue { get; }
    /// <summary>A pulse's phase turned by φ multiplies each pathway p → p′ by e^{−i(p′ − p)φ} (the 90° pulse, every φ, an
    /// identity of the conjugation by a diagonal turn; the control with the opposite sign fails), and the turned pulse is
    /// the pulse built about its own axis by Sylvester's formula, at all 24 phases.</summary>
    public bool PathwayPhase { get; }
    /// <summary>Run on the spin: eight steps of the 90° pulse on T₃₀ keep exactly the part at p₀ = 3, 1 or 2, an equality
    /// any pulse obeys by the selection rule; what the spin adds is which parts are there, none at order 2 and some at 3
    /// and 1. The control, two steps, keeps more than the order 1.</summary>
    public bool EightStepsOnTheSpin { get; }
    public bool PhaseCycleHolds => CycleKeepsResidue && PathwayPhase && EightStepsOnTheSpin;

    // ===== gate 7: the 180° pulse =====
    /// <summary>The ideal 180° pulse at phase φ is i·(X·e^{2iφZ}) ⊗ (X·e^{iφZ}) at every φ, and at all 24 phases the pulse
    /// built about its own axis; anti-diagonal with entries of modulus 1, so |m⟩⟨m′| goes to |−m⟩⟨−m′| at every φ.</summary>
    public bool Pulse180IsAProduct { get; }
    /// <summary>It keeps every string on its k_Z at every φ.</summary>
    public bool Pulse180KeepsKz { get; }
    /// <summary>At 0° and 90° every string goes to ± itself; the control, 45°, moves strings onto other strings.</summary>
    public bool Pulse180SignedAtQuarterTurns { get; }
    /// <summary>The 90° pulse does not keep k_Z: IZ goes to XY.</summary>
    public bool Pulse90MovesKz { get; }
    public bool Pulse180Holds => Pulse180IsAProduct && Pulse180KeepsKz && Pulse180SignedAtQuarterTurns && Pulse90MovesKz;

    // ===== gate 8: the echo =====
    /// <summary>The 180° pulse at 0°, as the sandwich ρ ↦ RρR†, reverses the offset ZI + IZ/2 and keeps the splitting ZZ
    /// (both read from Core's <see cref="Spin32Strings.Hamiltonian"/>).</summary>
    public bool EchoReversesOffsetKeepsSplitting { get; }
    /// <summary>None of the sixteen plain sandwiches ρ ↦ AρB with A, B products of one-site invertible diagonal or
    /// anti-diagonal factors reverses offset and splitting together (on a diagonal H such a factor acts by its kind alone, so I and X
    /// represent every kind).</summary>
    public bool NoPlainSandwichReversesBoth { get; }
    /// <summary>ρ ↦ ρᵀ reverses the offset and the splitting, each and together (as it reverses every real symmetric H,
    /// (−i[H, ρᵀ])ᵀ = i[Hᵀ, ρ]), and is no pulse: its Choi matrix, built from it, has the eigenvalue −1 (a pulse's is
    /// positive).</summary>
    public bool TransposeReversesBothAndIsNoPulse { get; }
    public bool EchoHolds => EchoReversesOffsetKeepsSplitting && NoPlainSandwichReversesBoth && TransposeReversesBothAndIsNoPulse;

    // ===== gate 9: the 180° pulse in the page's bits =====
    /// <summary>Operator Schmidt ranks in the page's bits at 0°, 90°, 30°, 45°; the gate also reads the pulse there as
    /// i·IX at 0° and −i·ZY at 90°.</summary>
    public IReadOnlyList<int> Pulse180RanksInPageBits { get; }
    public bool Pulse180PageBitsHolds { get; }

    // ===== gate 10: the 90° pulse =====
    /// <summary>I_z → −I_y, ZI → −(√3/2)·IY − YX/2, and T₃₀ onto the orders +3 and −3 with the squared fraction 5/16 at each;
    /// with them the turn without its phase factoring in Souza et al.'s bits at the 24 phases, which carries the pulse's
    /// operator Schmidt rank to every phase.</summary>
    public bool Pulse90Images { get; }
    /// <summary>Operator Schmidt rank in Souza et al.'s bits, at 0° and 15° (every phase by the factoring turn).</summary>
    public IReadOnlyList<int> Pulse90RanksInSouzaBits { get; }
    /// <summary>In the page's bits at 0°, 30°, 45°, 90°.</summary>
    public IReadOnlyList<int> Pulse90RanksInPageBits { get; }
    public bool Pulse90Holds => Pulse90Images && Pulse90RanksInSouzaBits.SequenceEqual(new[] { 4, 4 })
        && Pulse90RanksInPageBits.SequenceEqual(new[] { 2, 4, 4, 4 });

    // ===== gate 11: the receiver =====
    /// <summary>I₊ = √3·(I ⊗ |0⟩⟨1|) + 2·(|0⟩⟨1| ⊗ |1⟩⟨0|), its cells at v = (0, 1) (entry √3) and (1, −1) (entry 2),
    /// squared norms 6 : 4.</summary>
    public bool ReceiverCells { get; }
    /// <summary>Tr(I₊·ρ) sees I₋ and the rank-1 part a 90° pulse makes of I_z; T₂,₋₁ and T₃,₋₁ (rank and order read
    /// here), and all a 90° pulse carries from T₃₀, are invisible.</summary>
    public bool ReceiverReadsRankOne { get; }
    /// <summary>At zero splitting and offset the satellite and central parts of T₃,₋₁ are eigenoperators at
    /// −3(J₀ + J₁) and −3(J₁ + J₂) at every density, with the amplitudes ±√6/5 on T₁,₋₁ (both unit norm), cancelling
    /// where the rates meet, J₀ = J₂.</summary>
    public bool RelaxationCarriesRankThree { get; }
    /// <summary>The splitting ZZ carries T₃,₋₁ into ZX, ZY at first order (invisible) and into view at second; the dynamic
    /// shift −3K₂·ZI − 3K₁·IZ (Core's <see cref="Spin32Strings.Hamiltonian"/>) into view at first order,
    /// 48√3·(2K₁ − K₂).</summary>
    public bool SplittingAndShiftInView { get; }
    /// <summary>Operator Schmidt rank of I₊: in Souza et al.'s bits, in the page's bits.</summary>
    public IReadOnlyList<int> ReceiverRanks { get; }
    public bool ReceiverHolds => ReceiverCells && ReceiverReadsRankOne && RelaxationCarriesRankThree && SplittingAndShiftInView
        && ReceiverRanks.SequenceEqual(new[] { 2, 3 });

    // ===== gate 12: the pools =====
    /// <summary>The label's Pauli-transfer matrix under |1⟩⟨0| at (1 − f₀)·k and |0⟩⟨1| at f₀·k: X and Y at −k/2, Z at
    /// −k, the arrow from I to Z k·(2f₀ − 1), at every f₀ and k (linear in the two rates, read at unit rates).</summary>
    public bool PoolTransferMatrix { get; }
    /// <summary>The fractions are the steady state (read at five f₀, a quadratic in f₀), so ⟨Z_label⟩ = 2f₀ − 1; the rates
    /// swapped do not keep them.</summary>
    public bool PoolSteadyState { get; }
    public bool PoolsHold => PoolTransferMatrix && PoolSteadyState;

    /// <summary>The field's own control: [[1, √2], [√2, 2]] and [[√3, 3], [1, √3]] have rank 1 over ℚ(ζ₂₄), [[1, √2],
    /// [√2, 3]] rank 2.</summary>
    public bool FieldRankControl { get; }

    public bool AllHold => LabellingHolds && CoherenceOrderHolds && UnitalityHolds && InversionRecoveryHolds && TurnHolds
        && PhaseCycleHolds && Pulse180Holds && EchoHolds && Pulse180PageBitsHolds && Pulse90Holds && ReceiverHolds && PoolsHold
        && FieldRankControl;

    private readonly Dictionary<string, Mat> _s = new();
    private readonly Mat _iz, _ip, _im, _ix, _iy, _pab, _pwrong;

    public SodiumExperimentInStringsWitness()
    {
        foreach (var l in Labels) _s[l] = KM.FromGaussian(Spin32Strings.Matrix(l));

        // ---- spin 3/2 from the spin formula, I₊|m_k⟩ = a(m_k)·|m_{k−1}⟩, a² = I(I + 1) − m(m + 1)
        var space = SpinSpace.Of(3);
        _iz = KM.Zero(Levels);
        _ip = KM.Zero(Levels);
        for (int k = 0; k < Levels; k++)
        {
            _iz[k, k] = K.FromRational(space.M(k));
            if (k > 0) _ip[k - 1, k] = SqrtOf(space.RaiseSquared(k));
        }
        _im = KM.Transpose(_ip);
        _ix = KM.Scale(KM.Add(_ip, _im), new BigRational(1, 2));
        _iy = KM.Scale(KM.Subtract(_ip, _im), K.I * new BigRational(-1, 2));         // (I₊ − I₋)/(2i)

        _pab = Relabel((b1, b2) => b1 ^ b2, (b1, _) => b1);
        _pwrong = Relabel((b1, _) => b1, (b1, b2) => b1 ^ b2);

        FieldRankControl = KM.Rank(Mat2(K.One, K.Sqrt2, K.Sqrt2, R(2))) == 1
                           && KM.Rank(Mat2(K.Sqrt3, R(3), K.One, K.Sqrt3)) == 1
                           && KM.Rank(Mat2(K.One, K.Sqrt2, K.Sqrt2, R(3))) == 2;

        (LabellingReproducesIz, LabellingSignedImages, StringsMovedInKz, WrongAssignmentFails) = ReadLabelling();
        (OrderIsWeightedBitOrders, SwappedWeightsFail, PageBitsDoNotFixOrder) = ReadOrder();
        (GeneratorLinear, UnitalAndTracePreserving, AmplitudeDampingMovesIdentity, ThermalSourceHolds) = ReadUnitality();
        (LongitudinalEigenoperators, InversionAmplitudes, SaturationAmplitudes, InversionRecoveryHolds) = ReadInversion();
        (TurnFactorsInSouzaBits, TurnFormInPageBits, TurnRanksInPageBits) = ReadTurn();
        (CycleKeepsResidue, PathwayPhase, EightStepsOnTheSpin) = ReadPhaseCycle();
        (Pulse180IsAProduct, Pulse180KeepsKz, Pulse180SignedAtQuarterTurns, Pulse90MovesKz) = ReadPulse180();
        (EchoReversesOffsetKeepsSplitting, NoPlainSandwichReversesBoth, TransposeReversesBothAndIsNoPulse) = ReadEcho();
        (Pulse180RanksInPageBits, Pulse180PageBitsHolds) = ReadPulse180InPageBits();
        (Pulse90Images, Pulse90RanksInSouzaBits, Pulse90RanksInPageBits) = ReadPulse90();
        (ReceiverCells, ReceiverReadsRankOne, RelaxationCarriesRankThree, SplittingAndShiftInView, ReceiverRanks) = ReadReceiver();
        (PoolTransferMatrix, PoolSteadyState) = ReadPools();
    }

    // ===== gate 1 =====

    private (bool, bool, int, bool) ReadLabelling()
    {
        var za = _s["ZI"];                                          // in the page's order a is the high bit
        var zb = _s["IZ"];
        var izTarget = KM.Add(zb, KM.Scale(_s["ZZ"], new BigRational(1, 2)));    // Z_b + Z_a·Z_b/2
        bool iz = KM.AreEqual(Conj(_pab, _iz), izTarget) && KM.AreEqual(Conj(_pab, _s["ZZ"]), za);
        var images = Labels.ToDictionary(l => l, l => StringsOf(Conj(_pab, _s[l])));
        bool signed = images.Values.All(v => v.Count == 1 && (v.Values.Single() == K.One || v.Values.Single() == -K.One));
        var listed = new (string From, string To, int Sign)[]
        {
            ("ZI", "IZ", 1), ("IZ", "ZZ", 1), ("ZZ", "ZI", 1), ("IX", "XI", 1), ("XX", "IX", 1), ("XY", "ZY", 1), ("YY", "ZX", -1),
        };
        signed &= listed.All(t => images[t.From].Count == 1 && images[t.From].TryGetValue(t.To, out var c) && c == R(t.Sign));
        int moved = Labels.Count(l => signed && Kz(images[l].Keys.Single()) != Kz(l));
        bool wrong = !KM.AreEqual(Conj(_pwrong, _iz), izTarget);
        return (iz, signed, moved, wrong);
    }

    // ===== gate 2 =====

    private (bool, bool, bool) ReadOrder()
    {
        static int Z(int bit) => 1 - 2 * bit;
        static (int B1, int B2) Souza(int k) => (k >> 1, k & 1);
        (int A, int B) Page(int k)                                   // read off the relabelling of gate 1
        {
            int row = Enumerable.Range(0, Levels).Single(r => _pab[r, k] == K.One);
            return (row >> 1, row & 1);
        }
        (int, int) V((int, int) x, (int, int) y) => ((Z(x.Item1) - Z(y.Item1)) / 2, (Z(x.Item2) - Z(y.Item2)) / 2);
        bool weighted = true, swapped = true;
        for (int a = 0; a < Levels; a++)
            for (int b = 0; b < Levels; b++)
            {
                int p = Order(a, b);
                var (v1, v2) = V(Souza(a), Souza(b));
                int popDiff = (Souza(b).B1 + Souza(b).B2) - (Souza(a).B1 + Souza(a).B2);
                var cell = KM.Zero(Levels);
                cell[a, b] = K.One;
                weighted &= p == 2 * v1 + v2 && popDiff == v1 + v2
                            && StringsOf(cell).Keys.All(u => Kz(u) == Math.Abs(v1) + Math.Abs(v2));
                swapped &= p == v1 + 2 * v2;
            }
        var orders = new Dictionary<(int, int), HashSet<int>>();
        for (int a = 0; a < Levels; a++)
            for (int b = 0; b < Levels; b++)
            {
                var v = V(Page(a), Page(b));
                if (!orders.TryGetValue(v, out var set)) orders[v] = set = new HashSet<int>();
                set.Add(Order(a, b));
            }
        var clash = orders.Where(kv => kv.Value.Count > 1).Select(kv => kv.Key).ToHashSet();
        bool page = clash.SetEquals(new[] { (1, 0), (-1, 0), (0, 1), (0, -1) });
        return (weighted, !swapped, page);
    }

    // ===== gates 3 and 4: the generator =====

    private static readonly (long J0, long J1, long J2, long Q, long Omega)[] UnitPoints =
        { (1, 0, 0, 0, 0), (0, 1, 0, 0, 0), (0, 0, 1, 0, 0), (0, 0, 0, 1, 0), (0, 0, 0, 0, 1) };

    /// <summary>The generic points for the linearity guard.</summary>
    public static readonly IReadOnlyList<(long J0, long J1, long J2, long Q, long Omega)> GenericPoints = new[] { (2L, 3L, 5L, 7L, 11L), (-3L, 4L, -1L, 2L, 9L) };

    private static Mat Generator((long J0, long J1, long J2, long Q, long Omega) p) =>
        KM.FromGaussian(new QuadrupolarGenerator(SpinSpace.Of(3), new QuadrupolarRates(p.J0, p.J1, p.J2),
            SecularHamiltonian.Of(p.Q, p.Omega)).Physical());

    private (bool, bool, bool, bool) ReadUnitality()
    {
        var units = UnitPoints.Select(Generator).ToArray();
        bool linear = GenericPoints.All(pt =>
        {
            long[] g = { pt.J0, pt.J1, pt.J2, pt.Q, pt.Omega };
            var sum = KM.Zero(16);
            for (int e = 0; e < 5; e++) sum = KM.Add(sum, KM.Scale(units[e], new BigRational(g[e])));
            return KM.AreEqual(Generator(pt), sum);
        });
        bool unital = units.All(l => KM.IsZero(Apply(l, _s["II"])) && Labels.All(s => KM.Trace(Apply(l, _s[s])).IsZero));
        // amplitude damping toward bit value 1 on bit 1 as a superoperator, Sandwich(J, J†) − ½(J†J ⊗ 1 + 1 ⊗ J†J),
        // J = |1⟩⟨0| ⊗ I, read through the same Apply as the generator
        var jump = GM.Zero(4);
        jump[2, 0] = GaussianRational.One;
        jump[3, 1] = GaussianRational.One;
        var jdj = GM.Multiply(GM.Adjoint(jump), jump);
        var one = GM.Identity(4);
        var damp = KM.FromGaussian(GM.Subtract(GM.Sandwich(jump, GM.Adjoint(jump)),
            GM.Scale(GM.Add(GM.Sandwich(jdj, one), GM.Sandwich(one, jdj)), new BigRational(1, 2))));
        bool damping = Labels.All(s => KM.Trace(Apply(damp, _s[s])).IsZero) && !KM.IsZero(Apply(damp, _s["II"]));
        // −L_e(I_z) for the five unit points: (0, 3·IZ, 6·ZI, 0, 0)
        var expected = new[] { KM.Zero(4), KM.Scale(_s["IZ"], new BigRational(3)), KM.Scale(_s["ZI"], new BigRational(6)), KM.Zero(4), KM.Zero(4) };
        bool source = Enumerable.Range(0, 5).All(e => KM.AreEqual(KM.Scale(Apply(units[e], _iz), new BigRational(-1)), expected[e]));
        return (linear, unital, damping, source);
    }

    private (bool, IReadOnlyList<BigRational>, IReadOnlyList<BigRational>, bool) ReadInversion()
    {
        var units = UnitPoints.Select(Generator).ToArray();
        long[] lamZi = { 0, 0, -6, 0, 0 }, lamIz = { 0, -6, 0, 0, 0 };
        bool eigen = Enumerable.Range(0, 5).All(e =>
            KM.AreEqual(Apply(units[e], _s["ZI"]), KM.Scale(_s["ZI"], new BigRational(lamZi[e])))
            && KM.AreEqual(Apply(units[e], _s["IZ"]), KM.Scale(_s["IZ"], new BigRational(lamIz[e]))));
        var t10 = KM.Add(KM.Scale(_s["ZI"], new BigRational(2)), _s["IZ"]);
        var t30 = KM.Subtract(KM.Scale(_s["IZ"], new BigRational(2)), _s["ZI"]);
        bool equalNorms = Ip(t10, t10) == Ip(t30, t30);
        // the rank-1 amplitude of the start −I_z, |⟨T₁₀, −I_z⟩|/|T₁₀|; with |T₁₀| = |T₃₀| the norms cancel
        var rank1 = Rational(Ip(t10, KM.Scale(_iz, new BigRational(-1))));
        BigRational Unit(Mat m) => Rational(Ip(t30, m)) / (rank1.Sign < 0 ? -rank1 : rank1);
        Mat Along(string s, Mat m) => KM.Scale(_s[s], Ip(_s[s], m) * new BigRational(1, 4));
        // the deviation from equilibrium on span(ZI, IZ), which it lies in entirely
        IReadOnlyList<BigRational> Amplitudes(Mat equilibrium, Mat deviation) => new[]
        {
            Unit(equilibrium), Unit(Along("ZI", deviation)), Unit(Along("IZ", deviation)),
        };
        var minus2Iz = KM.Scale(_iz, new BigRational(-2));
        bool spans = KM.AreEqual(KM.Add(Along("ZI", minus2Iz), Along("IZ", minus2Iz)), minus2Iz);
        var inversion = Amplitudes(_iz, minus2Iz);                                     // −I_z relaxing toward +I_z
        var saturation = Amplitudes(KM.Zero(4), KM.Scale(_iz, new BigRational(-1)));   // the deviation −I_z
        bool holds = eigen && equalNorms && spans && Carries(t30, 3, 0) && Carries(t10, 1, 0)
                     && inversion.SequenceEqual(new[] { BigRational.Zero, new BigRational(4, 5), new BigRational(-4, 5) })
                     && saturation.SequenceEqual(new[] { BigRational.Zero, new BigRational(2, 5), new BigRational(-2, 5) });
        return (eigen, inversion, saturation, holds);
    }

    // ===== gate 5 =====

    private (bool, bool, IReadOnlyList<int>) ReadTurn()
    {
        bool factors = true, page = true;
        var ranks = new int[Phases];
        static K W(int k, int power) => K.RootPower(-k * power);                  // w^power, w = e^{−iψ}, ψ = 15°·k
        for (int k = 0; k < Phases; k++)
        {
            var turn = TurnHalfAngle(k);                                           // e^{−imφ} = w^{2m}, φ = 2ψ
            factors &= KM.AreEqual(turn, KM.Kron(Diag(W(k, 2), W(k, -2)), Diag(W(k, 1), W(k, -1))));
            for (int a = 0; a < Levels; a++)
                for (int b = 0; b < Levels; b++)
                {
                    var (v1, v2) = SouzaV(a, b);
                    factors &= turn[a, a] * turn[b, b].Conjugate == W(k, 2 * (2 * v1 + v2));
                }
            var cos = K.Cos(k);
            var sin = K.Sin(k);
            var form = KM.Multiply(KM.Kron(KM.Identity(2), Diag(W(k, 2), W(k, -2))),
                KM.Subtract(KM.Scale(KM.Identity(4), cos), KM.Scale(_s["ZZ"], K.I * sin)));
            var inPage = Conj(_pab, turn);
            page &= KM.AreEqual(inPage, form);
            ranks[k] = SchmidtRank(inPage);
        }
        // at φ = 90° (ψ = 45°, k = 3): the clock matrix Z ⊗ diag(1, i) times e^{−3iπ/4}
        var clock = KM.Scale(KM.Kron(Diag(K.One, -K.One), Diag(K.One, K.I)), K.RootPower(-9));
        factors &= KM.AreEqual(TurnHalfAngle(3), clock);
        // the control: the cell phase with the weights swapped fails at φ = 30° (half-angle 15°)
        var t1 = TurnHalfAngle(1);
        bool control = Enumerable.Range(0, Levels).Any(a => Enumerable.Range(0, Levels).Any(b =>
        {
            var (v1, v2) = SouzaV(a, b);
            return t1[a, a] * t1[b, b].Conjugate != W(1, 2 * (v1 + 2 * v2));
        }));
        return (factors && control, page, ranks);
    }

    // ===== gate 6 =====

    /// <summary>e^{−iθA} for θ = tTheta·π/12 (tTheta even) and A with the four eigenvalues ±1/2, ±3/2, by Sylvester's
    /// formula on the spectral projectors; null where Π(A − λ) ≠ 0, where the formula would not apply.</summary>
    private static Mat? Expo(Mat a, int tTheta)
    {
        if (tTheta % 2 != 0) throw new ArgumentException("θ must be an even multiple of 15°.", nameof(tTheta));
        int[] twice = { 3, 1, -1, -3 };
        var one = KM.Identity(Levels);
        Mat Shifted(int l) => KM.Subtract(a, KM.Scale(one, new BigRational(l, 2)));
        if (!KM.IsZero(twice.Select(Shifted).Aggregate(KM.Multiply))) return null;
        var r = KM.Zero(Levels);
        foreach (int l in twice)
        {
            var proj = one;
            foreach (int lp in twice.Where(x => x != l))
                proj = KM.Scale(KM.Multiply(proj, Shifted(lp)), new BigRational(2, l - lp));
            r = KM.Add(r, KM.Scale(proj, K.RootPower(-tTheta * l / 2)));
        }
        return r;
    }

    /// <summary>cos φ·I_x + sin φ·I_y at φ = 15°·kPhi.</summary>
    private Mat Axis(int kPhi) => KM.Add(KM.Scale(_ix, K.Cos(kPhi)), KM.Scale(_iy, K.Sin(kPhi)));

    /// <summary>The turn about z by φ = 15°·kPhi without its global phase, diag(1, e^{iφ}, e^{2iφ}, e^{3iφ}); a conjugation
    /// by it is a conjugation by e^{−iφI_z}.</summary>
    private Mat TurnUpToPhase(int kPhi) => Diag(Enumerable.Range(0, Levels).Select(l => K.RootPower(kPhi * Order(0, l))).ToArray());

    /// <summary>e^{−iφI_z} at φ = 2ψ, ψ = 15°·k: the entries w^{2m} with w = e^{−iψ} = ζ^{−k}, 2m read off I_z.</summary>
    private Mat TurnHalfAngle(int k) => Diag(Enumerable.Range(0, Levels).Select(l => K.RootPower(-k * TwiceM(l))).ToArray());

    /// <summary>2m at level l, read off I_z (an integer at spin 3/2).</summary>
    private int TwiceM(int level) => Rational(_iz[level, level] * new BigRational(2)) is { IsInteger: true } t
        ? (int)t.Numerator : throw new InvalidOperationException("2m is not an integer.");

    /// <summary>m_a − m_b, the coherence order of the cell |a⟩⟨b|, read off I_z.</summary>
    private int Order(int a, int b) => (TwiceM(a) - TwiceM(b)) is var d && d % 2 == 0
        ? d / 2 : throw new InvalidOperationException("m_a − m_b is not an integer.");

    private Mat Phased(Mat pulse, int kPhi) => Conj(TurnUpToPhase(kPhi), pulse, unitary: true);

    private Mat? _r90, _r180;
    private Mat R90 => _r90 ??= Expo(_ix, 6) ?? throw new InvalidOperationException("I_x has not the minimal polynomial of spin 3/2.");
    private Mat R180 => _r180 ??= Expo(_ix, 12) ?? throw new InvalidOperationException("I_x has not the minimal polynomial of spin 3/2.");

    private Mat Part(Mat m, int p0)
    {
        var r = KM.Zero(Levels);
        for (int a = 0; a < Levels; a++)
            for (int b = 0; b < Levels; b++)
                if (Order(a, b) == p0) r[a, b] = m[a, b];
        return r;
    }

    private Mat T30 => KM.Subtract(KM.Scale(_s["IZ"], new BigRational(2)), _s["ZI"]);

    private (bool, bool, bool) ReadPhaseCycle()
    {
        bool keep = true;
        for (int n = 1; n <= 8; n++)
            for (int p = -6; p <= 6; p++)
                for (int p0 = -3; p0 <= 3; p0++)
                {
                    var acc = CyclotomicRing.Zero(n);
                    for (int k = 0; k < n; k++) CyclotomicRing.AddRootPower(acc, n, (long)(p0 - p) * k, 1);
                    var expected = CyclotomicRing.Zero(n);
                    if (((p - p0) % n + n) % n == 0) CyclotomicRing.AddRootPower(expected, n, 0, n);
                    keep &= CyclotomicRing.AreEqual(acc, expected);
                }

        var r = R90;
        bool path = true;
        for (int kPhi = 0; kPhi < Phases; kPhi++)
        {
            var rph = Phased(r, kPhi);
            for (int a = 0; a < Levels; a++)
                for (int c = 0; c < Levels; c++)
                    for (int b = 0; b < Levels; b++)
                        for (int e = 0; e < Levels; e++)
                            path &= rph[a, c] * rph[b, e].Conjugate
                                    == r[a, c] * r[b, e].Conjugate * K.RootPower(-(Order(a, b) - Order(c, e)) * kPhi);
        }
        var rph1 = Phased(r, 1);
        bool wrongSign = Enumerable.Range(0, Levels).Any(a => Enumerable.Range(0, Levels).Any(c =>
            Enumerable.Range(0, Levels).Any(b => Enumerable.Range(0, Levels).Any(e =>
                rph1[a, c] * rph1[b, e].Conjugate != r[a, c] * r[b, e].Conjugate * K.RootPower(Order(a, b) - Order(c, e))))));
        bool indep = Enumerable.Range(0, Phases).All(k => Expo(Axis(k), 6) is { } direct && KM.AreEqual(direct, Phased(r, k)));

        var t30 = T30;
        var after = Conj(r, t30, unitary: true);
        Mat Cycle(int steps, int p0)
        {
            if (Phases % steps != 0) throw new ArgumentException($"{steps} steps do not divide the {Phases} phases.");
            var acc = KM.Zero(Levels);
            for (int k = 0; k < steps; k++)
            {
                int kPhi = Phases / steps * k;                     // φ_k = 2πk/steps
                acc = KM.Add(acc, KM.Scale(Conj(Phased(r, kPhi), t30, unitary: true), K.RootPower(p0 * kPhi)));
            }
            return KM.Scale(acc, new BigRational(1, steps));
        }
        bool physical = new[] { 3, 1, 2 }.All(p0 => KM.AreEqual(Cycle(8, p0), Part(after, p0)))
                        && KM.IsZero(Part(after, 2)) && !KM.IsZero(Part(after, 3)) && !KM.IsZero(Part(after, 1))
                        && !KM.AreEqual(Cycle(2, 1), Part(after, 1));
        return (keep, path && wrongSign && indep, physical);
    }

    // ===== gate 7 =====

    private static Mat Product180(int kPhi)
    {
        var x = Mat2(K.Zero, K.One, K.One, K.Zero);
        var left = KM.Multiply(x, Diag(K.RootPower(2 * kPhi), K.RootPower(-2 * kPhi)));
        var right = KM.Multiply(x, Diag(K.RootPower(kPhi), K.RootPower(-kPhi)));
        return KM.Scale(KM.Kron(left, right), K.I);
    }

    private (bool, bool, bool, bool) ReadPulse180()
    {
        var r = R180;
        bool product = true, kz = true;
        for (int kPhi = 0; kPhi < Phases; kPhi++)
        {
            var rph = Phased(r, kPhi);
            product &= KM.AreEqual(rph, Product180(kPhi));
            for (int a = 0; a < Levels; a++)
                for (int b = 0; b < Levels; b++)
                    product &= b == Levels - 1 - a ? rph[a, b] * rph[a, b].Conjugate == K.One : rph[a, b].IsZero;
            foreach (var l in Labels)
            {
                var img = StringsOf(Conj(rph, _s[l], unitary: true));
                kz &= img.Count > 0 && img.Keys.All(u => Kz(u) == Kz(l));
            }
        }
        product &= Enumerable.Range(0, Phases).All(k => Expo(Axis(k), 12) is { } direct && KM.AreEqual(direct, Product180(k)));
        bool ToItself(int kPhi) => Labels.All(l =>
        {
            var img = StringsOf(Conj(Phased(r, kPhi), _s[l], unitary: true));
            return img.Count == 1 && img.TryGetValue(l, out var c) && (c == K.One || c == -K.One);
        });
        bool signed = ToItself(0) && ToItself(6) && !ToItself(3);
        var iz90 = StringsOf(Conj(R90, _s["IZ"], unitary: true));
        bool moves = iz90.Count == 1 && iz90.TryGetValue("XY", out var cx) && cx == K.One;
        return (product, kz, signed, moves);
    }

    // ===== gate 8 =====

    private (bool, bool, bool) ReadEcho()
    {
        GaussianRational[,] S(string l) => Spin32Strings.Matrix(l);
        var hOff = FromStrings(SecularHamiltonian.Of(0, 1));                         // ω·I_z at ω = 1: ZI + IZ/2
        var hSplit = FromStrings(SecularHamiltonian.Of(1, 0));                       // q·ZZ at q = 1
        var hBoth = GM.Add(hOff, hSplit);
        GaussianRational[,] L(GaussianRational[,] h) => QuadrupolarGenerator.CommutatorPart(h);
        bool Reverses(GaussianRational[,] s, GaussianRational[,] sInv, GaussianRational[,] h) =>
            GM.AreEqual(GM.Multiply(GM.Multiply(s, L(h)), sInv), GM.Scale(L(h), -1));
        bool Keeps(GaussianRational[,] s, GaussianRational[,] sInv, GaussianRational[,] h) =>
            GM.AreEqual(GM.Multiply(GM.Multiply(s, L(h)), sInv), L(h));
        var pulse = ToGaussian(Phased(R180, 0));                                    // the 180° pulse at 0°, i·XX by gate 7
        var echo = GM.Sandwich(pulse, GM.Adjoint(pulse));
        var echoInv = GM.Sandwich(GM.Adjoint(pulse), pulse);
        bool offsetSplit = GM.AreEqual(hOff, GM.Add(S("ZI"), GM.Scale(S("IZ"), new BigRational(1, 2)))) && GM.AreEqual(hSplit, S("ZZ"))
                           && Reverses(echo, echoInv, hOff) && Keeps(echo, echoInv, hSplit);

        string[] letters = { "I", "X" };                                        // one representative per kind of site factor
        bool none = true;
        foreach (var a1 in letters) foreach (var a2 in letters) foreach (var b1 in letters) foreach (var b2 in letters)
        {
            var a = S(a1 + a2);
            var b = S(b1 + b2);
            none &= !Reverses(GM.Sandwich(a, b), GM.Sandwich(a, b), hBoth);   // I and X are their own inverses
        }

        // the transpose as a superoperator, and its Choi matrix Σ|i⟩⟨j| ⊗ T(|i⟩⟨j|) = Σ|i⟩⟨j| ⊗ |j⟩⟨i|, the swap
        var t = GM.Zero(16);
        for (int x = 0; x < Levels; x++)
            for (int y = 0; y < Levels; y++) t[Levels * y + x, Levels * x + y] = GaussianRational.One;
        bool transposed = Reverses(t, t, hBoth) && Reverses(t, t, hOff) && Reverses(t, t, hSplit);
        // the Choi matrix Σ|i⟩⟨j| ⊗ T(|i⟩⟨j|) read off the superoperator t: entry [4i + a, 4j + b] = T(|i⟩⟨j|)[a, b]
        var choi = GM.Zero(16);
        for (int i = 0; i < Levels; i++)
            for (int j = 0; j < Levels; j++)
                for (int a = 0; a < Levels; a++)
                    for (int b = 0; b < Levels; b++) choi[Levels * i + a, Levels * j + b] = t[Levels * a + b, Levels * i + j];
        var v = GM.Zero(16, 1);
        v[Levels * 0 + 1, 0] = GaussianRational.One;
        v[Levels * 1 + 0, 0] = -GaussianRational.One;
        bool noPulse = GM.AreEqual(GM.Multiply(choi, v), GM.Scale(v, -1));
        return (offsetSplit, none, transposed && noPulse);
    }

    // ===== gate 9 =====

    private (IReadOnlyList<int>, bool) ReadPulse180InPageBits()
    {
        var ranks = new[] { 0, 6, 2, 3 }.Select(kPhi => SchmidtRank(Conj(_pab, Phased(R180, kPhi)))).ToArray();
        // the pulse itself in the page's bits: at 0° i·IX, a flip of b alone; at 90° −i·ZY
        var at0 = StringsOf(Conj(_pab, Phased(R180, 0)));
        var at90 = StringsOf(Conj(_pab, Phased(R180, 6)));
        bool holds = ranks.SequenceEqual(new[] { 1, 1, 2, 2 })
                     && at0.Count == 1 && at0.TryGetValue("IX", out var c1) && c1 == K.I
                     && at90.Count == 1 && at90.TryGetValue("ZY", out var c2) && c2 == -K.I;
        return (ranks, holds);
    }

    // ===== gate 10 =====

    private (bool, IReadOnlyList<int>, IReadOnlyList<int>) ReadPulse90()
    {
        var r = R90;
        bool iz = KM.AreEqual(Conj(r, _iz, unitary: true), KM.Scale(_iy, new BigRational(-1)));
        var zi = StringsOf(Conj(r, _s["ZI"], unitary: true));
        bool ziImage = zi.Count == 2 && zi.TryGetValue("IY", out var cy) && cy == K.Sqrt3 * new BigRational(-1, 2)
                       && zi.TryGetValue("YX", out var cyx) && cyx == K.FromRational(new BigRational(-1, 2));
        var t30 = T30;
        var b30 = Conj(r, t30, unitary: true);
        BigRational Fraction(int p0)
        {
            var s = K.Zero;
            for (int a = 0; a < Levels; a++)
                for (int b = 0; b < Levels; b++)
                    if (Order(a, b) == p0) s += b30[a, b] * b30[a, b].Conjugate;
            return Rational(s) / Rational(Ip(t30, t30));
        }
        bool orders = Fraction(3) == new BigRational(5, 16) && Fraction(-3) == new BigRational(5, 16);
        // the turn without its phase is a product in Souza et al.'s bits at every sampled phase, so a pulse's operator
        // Schmidt rank there does not depend on its phase
        bool turnProduct = Enumerable.Range(0, Phases).All(k =>
            KM.AreEqual(TurnUpToPhase(k), KM.Kron(Diag(K.One, K.RootPower(2 * k)), Diag(K.One, K.RootPower(k)))));
        var souza = new[] { 0, 1 }.Select(k => SchmidtRank(Phased(r, k))).ToArray();
        var page = new[] { 0, 2, 3, 6 }.Select(k => SchmidtRank(Conj(_pab, Phased(r, k)))).ToArray();
        return (iz && ziImage && orders && turnProduct, souza, page);
    }

    // ===== gate 11 =====

    private (bool, bool, bool, bool, IReadOnlyList<int>) ReadReceiver()
    {
        var raise = Mat2(K.Zero, K.One, K.Zero, K.Zero);           // |0⟩⟨1|
        var lower = Mat2(K.Zero, K.Zero, K.One, K.Zero);           // |1⟩⟨0|
        bool form = KM.AreEqual(_ip, KM.Add(KM.Scale(KM.Kron(KM.Identity(2), raise), K.Sqrt3),
            KM.Scale(KM.Kron(raise, lower), R(2))));
        var byV = new Dictionary<(int, int), List<K>>();
        for (int a = 0; a < Levels; a++)
            for (int b = 0; b < Levels; b++)
                if (!_ip[a, b].IsZero)
                {
                    var v = SouzaV(a, b);
                    if (!byV.TryGetValue(v, out var list)) byV[v] = list = new List<K>();
                    list.Add(_ip[a, b]);
                }
        BigRational NormAt((int, int) v) => Rational(byV[v].Aggregate(K.Zero, (s, z) => s + z * z.Conjugate));
        bool cells = form && byV.Count == 2 && byV.ContainsKey((0, 1)) && byV.ContainsKey((1, -1))
                     && byV[(0, 1)].All(z => z == K.Sqrt3) && byV[(1, -1)].All(z => z == R(2))
                     && NormAt((0, 1)) == 6 && NormAt((1, -1)) == 4;

        Mat Lin(params (K C, string S)[] terms) => terms.Aggregate(KM.Zero(Levels), (acc, t) => KM.Add(acc, KM.Scale(_s[t.S], t.C)));
        var minusI = -K.I;
        var satDir = Lin((K.One, "IX"), (minusI, "IY"));                                        // IX − i·IY
        var t3m1 = KM.Subtract(KM.Scale(satDir, R(2)),
            KM.Scale(Lin((K.One, "XX"), (K.One, "YY"), (minusI, "YX"), (K.I, "XY")), K.Sqrt3));  // ∝ T₃,₋₁
        var t2m1 = Lin((K.One, "ZX"), (minusI, "ZY"));                                           // ∝ T₂,₋₁
        Mat Comm(Mat x, Mat y) => KM.Subtract(KM.Multiply(x, y), KM.Multiply(y, x));
        K Read(Mat m) => KM.Trace(KM.Multiply(_ip, m));
        var r90 = R90;
        bool rankOne = Carries(t3m1, 3, -1) && Carries(t2m1, 2, -1) && Carries(_im, 1, -1) && !Carries(t3m1, 1, -1)
                       && Read(t3m1).IsZero && !Read(_im).IsZero && Read(t2m1).IsZero
                       && !Read(Conj(r90, _iz, unitary: true)).IsZero && Read(Conj(r90, T30, unitary: true)).IsZero;

        // relaxation at zero splitting and offset: the satellite part (along IX − i·IY) and the central part
        var satPart = KM.Scale(satDir, Ip(satDir, t3m1) * Inverse(Ip(satDir, satDir)));
        var ctPart = KM.Subtract(t3m1, satPart);
        var densities = new (long, long, long, long, long)[] { (1, 0, 0, 0, 0), (0, 1, 0, 0, 0), (0, 0, 1, 0, 0) };
        long[] lamSat = { -3, -3, 0 }, lamCt = { 0, -3, -3 };
        bool eigen = Enumerable.Range(0, 3).All(e =>
        {
            var l = Generator(densities[e]);
            return KM.AreEqual(Apply(l, satPart), KM.Scale(satPart, new BigRational(lamSat[e])))
                   && KM.AreEqual(Apply(l, ctPart), KM.Scale(ctPart, new BigRational(lamCt[e])));
        });
        var norm = K.Sqrt2 * new BigRational(20);                                               // |I₋|·|T₃,₋₁| = √10·√80
        bool normRight = Ip(_im, _im) == R(10) && Ip(t3m1, t3m1) == R(80) && norm * norm == R(800);
        // the amplitudes Tr(I₊·part)/norm against ±√6/5, compared with the norm multiplied through
        var sqrt6Over5 = K.Sqrt2 * K.Sqrt3 * new BigRational(1, 5);
        bool relax = eigen && normRight && Read(satPart) == sqrt6Over5 * norm && Read(ctPart) == -sqrt6Over5 * norm;

        // the splitting at first and second order, the dynamic shift's two coefficients
        var once = Comm(_s["ZZ"], t3m1);
        var onceStrings = StringsOf(once);
        var hK1 = KM.FromGaussian(FromStrings(SecularHamiltonian.None.WithLambShift(1, 0)));    // the K₁ term, −3K₁·IZ
        var hK2 = KM.FromGaussian(FromStrings(SecularHamiltonian.None.WithLambShift(0, 1)));    // the K₂ term, −3K₂·ZI
        var shiftK1 = Read(Comm(hK1, t3m1));
        var shiftK2 = Read(Comm(hK2, t3m1));
        bool view = Read(once).IsZero && onceStrings.Count > 0 && onceStrings.Keys.All(s => s is "ZX" or "ZY")
                    && !Read(Comm(_s["ZZ"], once)).IsZero
                    && KM.AreEqual(hK1, KM.Scale(_s["IZ"], R(-3))) && KM.AreEqual(hK2, KM.Scale(_s["ZI"], R(-3)))
                    && shiftK1 == K.Sqrt3 * new BigRational(96) && shiftK2 == K.Sqrt3 * new BigRational(-48);
        var ranks = new[] { SchmidtRank(_ip), SchmidtRank(Conj(_pab, _ip)) };
        return (cells, rankOne, relax, view, ranks);
    }

    // ===== gate 12 =====

    private static (bool, bool) ReadPools()
    {
        GaussianRational[,] P(char c) => c switch
        {
            'I' => new GaussianRational[,] { { 1, 0 }, { 0, 1 } },
            'X' => new GaussianRational[,] { { 0, 1 }, { 1, 0 } },
            'Y' => new GaussianRational[,] { { 0, -GaussianRational.I }, { GaussianRational.I, 0 } },
            _ => new GaussianRational[,] { { 1, 0 }, { 0, -1 } },
        };
        var lower = new GaussianRational[,] { { 0, 0 }, { 1, 0 } };                // |1⟩⟨0|: out of pool 0
        var raise = new GaussianRational[,] { { 0, 1 }, { 0, 0 } };                // |0⟩⟨1|: back
        GaussianRational[,] D(GaussianRational[,] j, BigRational rate, GaussianRational[,] x)
        {
            var jdj = GM.Multiply(GM.Adjoint(j), j);
            var body = GM.Subtract(GM.Multiply(GM.Multiply(j, x), GM.Adjoint(j)),
                GM.Scale(GM.Add(GM.Multiply(jdj, x), GM.Multiply(x, jdj)), new BigRational(1, 2)));
            return GM.Scale(body, rate);
        }
        GaussianRational Tr(GaussianRational[,] a) => a[0, 0] + a[1, 1];
        GaussianRational[,] Ptm(Func<GaussianRational[,], GaussianRational[,]> map)
        {
            var m = GM.Zero(4);
            for (int r = 0; r < 4; r++)
                for (int c = 0; c < 4; c++)
                    m[r, c] = Tr(GM.Multiply(P("IXYZ"[r]), map(P("IXYZ"[c])))) / 2;
            return m;
        }
        GaussianRational[,] Expected(BigRational heightScale, BigRational arrow)
        {
            var m = GM.Zero(4);
            m[1, 1] = -heightScale / 2;
            m[2, 2] = -heightScale / 2;
            m[3, 3] = -heightScale;
            m[3, 0] = arrow;
            return m;
        }
        bool ptm = GM.AreEqual(Ptm(x => D(lower, 1, x)), Expected(1, -1)) && GM.AreEqual(Ptm(x => D(raise, 1, x)), Expected(1, 1));
        // with the rates (1 − f₀)k and f₀k the sum is Expected(k, k(2f₀ − 1)): read at two points of (f₀, k)
        bool atPoints = new[] { (new BigRational(1, 3), new BigRational(2)), (new BigRational(3, 4), new BigRational(5)) }.All(pt =>
            GM.AreEqual(Ptm(x => GM.Add(D(lower, (1 - pt.Item1) * pt.Item2, x), D(raise, pt.Item1 * pt.Item2, x))),
                Expected(pt.Item2, pt.Item2 * (2 * pt.Item1 - 1))));
        var fs = new[] { BigRational.Zero, new BigRational(1, 4), new BigRational(1, 2), new BigRational(2, 3), BigRational.One };
        GaussianRational[,] Fractions(BigRational f0) => new GaussianRational[,] { { f0, 0 }, { 0, 1 - f0 } };
        bool steady = fs.All(f0 =>
            GM.AreEqual(GM.Add(D(lower, 1 - f0, Fractions(f0)), D(raise, f0, Fractions(f0))), GM.Zero(2)));
        var quarter = new BigRational(1, 4);
        bool swappedFails = !GM.AreEqual(GM.Add(D(lower, quarter, Fractions(quarter)), D(raise, 1 - quarter, Fractions(quarter))), GM.Zero(2));
        return (ptm && atPoints, steady && swappedFails);
    }

    // ===== helpers =====

    private static K R(long n) => K.FromRational(new BigRational(n));

    /// <summary>Σₐ[Iₐ, [Iₐ, M]] = l(l + 1)·M and [I_z, M] = p·M: M carries rank l and order p.</summary>
    private bool Carries(Mat m, int l, int p)
    {
        Mat Comm(Mat x, Mat y) => KM.Subtract(KM.Multiply(x, y), KM.Multiply(y, x));
        var cas = KM.Add(KM.Add(Comm(_ix, Comm(_ix, m)), Comm(_iy, Comm(_iy, m))), Comm(_iz, Comm(_iz, m)));
        return KM.AreEqual(cas, KM.Scale(m, new BigRational(l * (l + 1)))) && KM.AreEqual(Comm(_iz, m), KM.Scale(m, new BigRational(p)));
    }

    /// <summary>A Hamiltonian as Core writes it in Souza et al.'s strings.</summary>
    private static GaussianRational[,] FromStrings(SecularHamiltonian h) => Spin32Strings.Hamiltonian(h)
        .Aggregate(GM.Zero(4), (acc, t) => GM.Add(acc, GM.Scale(Spin32Strings.Matrix(t.Label), t.Coefficient)));

    /// <summary>A matrix whose entries lie in ℚ(i), as Gaussian rationals; anything else throws.</summary>
    private static GaussianRational[,] ToGaussian(Mat m)
    {
        var r = GM.Zero(m.GetLength(0), m.GetLength(1));
        for (int i = 0; i < m.GetLength(0); i++)
            for (int j = 0; j < m.GetLength(1); j++)
            {
                var z = m[i, j];
                if (Enumerable.Range(1, K.Degree - 1).Any(e => e != 6 && !z[e].IsZero))
                    throw new InvalidOperationException($"{z} is not in ℚ(i).");
                r[i, j] = new GaussianRational(z[0], z[6]);                         // i = ζ⁶
            }
        return r;
    }

    /// <summary>The value of a number the gate requires to be rational; anything else throws rather than being truncated.</summary>
    private static BigRational Rational(K z) =>
        z.TryRational(out var q) ? q : throw new InvalidOperationException($"{z} is not rational.");

    private static K Inverse(K z)
    {
        if (!z.TryRational(out var q) || q.IsZero) throw new InvalidOperationException($"{z} is not a non-zero rational.");
        return K.FromRational(BigRational.One / q);
    }

    /// <summary>√r for a rational r ≥ 0 whose root lies in ℚ(√2, √3): r, r/2, r/3 or r/6 a rational square.</summary>
    private static K SqrtOf(BigRational r)
    {
        foreach (var (div, root) in new[] { (1, K.One), (2, K.Sqrt2), (3, K.Sqrt3), (6, K.Sqrt2 * K.Sqrt3) })
            if (RationalSqrt(r / div) is { } q) return root * q;
        throw new InvalidOperationException($"√{r} is not in ℚ(√2, √3).");
    }

    private static BigRational? RationalSqrt(BigRational r)
    {
        if (r.Sign < 0) return null;
        var n = ISqrt(r.Numerator);
        var d = ISqrt(r.Denominator);
        return n * n == r.Numerator && d * d == r.Denominator ? new BigRational(n, d) : null;
    }

    private static System.Numerics.BigInteger ISqrt(System.Numerics.BigInteger v)
    {
        if (v.IsZero) return v;
        var x = v;
        var y = (x + 1) / 2;
        while (y < x) { x = y; y = (x + v / x) / 2; }
        return x;
    }

    private static Mat Diag(params K[] entries)
    {
        var m = KM.Zero(entries.Length);
        for (int i = 0; i < entries.Length; i++) m[i, i] = entries[i];
        return m;
    }

    private static Mat Mat2(K a, K b, K c, K d)
    {
        var m = KM.Zero(2);
        (m[0, 0], m[0, 1], m[1, 0], m[1, 1]) = (a, b, c, d);
        return m;
    }

    /// <summary>U·M·U⁻¹ for a permutation (U⁻¹ = Uᵀ) or, with <paramref name="unitary"/>, a unitary (U⁻¹ = U†).</summary>
    private static Mat Conj(Mat u, Mat m, bool unitary = false) =>
        KM.Multiply(u, m, unitary ? KM.Adjoint(u) : KM.Transpose(u));

    /// <summary>The permutation from Souza et al.'s index 2·bit₁ + bit₂ to the page's 2a + b.</summary>
    private static Mat Relabel(Func<int, int, int> a, Func<int, int, int> b)
    {
        var p = KM.Zero(Levels);
        for (int k = 0; k < Levels; k++)
        {
            int b1 = k >> 1, b2 = k & 1;
            p[2 * a(b1, b2) + b(b1, b2), k] = K.One;
        }
        return p;
    }

    private static (int V1, int V2) SouzaV(int a, int b)
    {
        static int Z(int bit) => 1 - 2 * bit;
        return ((Z(a >> 1) - Z(b >> 1)) / 2, (Z(a & 1) - Z(b & 1)) / 2);
    }

    private static int Kz(string label) => label.Count(c => c is 'X' or 'Y');

    /// <summary>⟨A, B⟩ = Tr(A†B).</summary>
    private static K Ip(Mat a, Mat b) => KM.Trace(KM.Multiply(KM.Adjoint(a), b));

    /// <summary>The string coefficients Tr(S·M)/4 that are not zero.</summary>
    private Dictionary<string, K> StringsOf(Mat m)
    {
        var r = new Dictionary<string, K>();
        foreach (var l in Labels)
        {
            var s = _s[l];
            var t = K.Zero;
            for (int i = 0; i < Levels; i++)
                for (int j = 0; j < Levels; j++)
                    if (!s[i, j].IsZero && !m[j, i].IsZero) t += s[i, j] * m[j, i];
            if (!t.IsZero) r[l] = t * new BigRational(1, 4);
        }
        return r;
    }

    /// <summary>The superoperator applied to a 4 × 4 operator, vec index 4i + j.</summary>
    private static Mat Apply(Mat super, Mat op)
    {
        var r = KM.Zero(Levels);
        for (int a = 0; a < Levels; a++)
            for (int b = 0; b < Levels; b++)
                for (int i = 0; i < Levels; i++)
                    for (int j = 0; j < Levels; j++)
                        if (!op[i, j].IsZero && !super[Levels * a + b, Levels * i + j].IsZero)
                            r[a, b] += super[Levels * a + b, Levels * i + j] * op[i, j];
        return r;
    }

    /// <summary>The operator Schmidt rank of a two-qubit operator across its two bits: the rank of the realigned
    /// matrix R[(i₁, j₁), (i₂, j₂)] = U[(i₁, i₂), (j₁, j₂)] over ℚ(ζ₂₄).</summary>
    private static int SchmidtRank(Mat u)
    {
        var r = KM.Zero(Levels);
        for (int row = 0; row < Levels; row++)
            for (int col = 0; col < Levels; col++)
                r[row, col] = u[2 * (row / 2) + col / 2, 2 * (row % 2) + col % 2];
        return KM.Rank(r);
    }

    // ===== IInspectable =====

    public string DisplayName => "the experiment in strings (spin 3/2, Souza et al.'s bits, exact over ℚ(ζ₂₄))";

    private static string Verdict(bool ok, string yes, string no) => ok ? yes : no;

    public string Summary =>
        $"twelve gates on preparation, pulses, phase cycle, receiver and pools in strings and the field's own control, {(AllHold ? "all holding" : "NOT all holding")}: " +
        $"the labelling {Verdict(LabellingHolds, "holds", "FAILS")}, the coherence order {Verdict(CoherenceOrderHolds, "holds", "FAILS")}, " +
        $"unitality and the thermal source {Verdict(UnitalityHolds, "hold", "FAIL")}, the inversion recovery {Verdict(InversionRecoveryHolds, "holds", "FAILS")}, " +
        $"the turn about z {Verdict(TurnHolds, "holds", "FAILS")}, the phase cycle {Verdict(PhaseCycleHolds, "holds", "FAILS")}, " +
        $"the 180° pulse {Verdict(Pulse180Holds, "holds", "FAILS")}, the echo {Verdict(EchoHolds, "holds", "FAILS")}, " +
        $"the 180° pulse in the page's bits {Verdict(Pulse180PageBitsHolds, "holds", "FAILS")}, the 90° pulse {Verdict(Pulse90Holds, "holds", "FAILS")}, " +
        $"the receiver {Verdict(ReceiverHolds, "holds", "FAILS")}, the pools {Verdict(PoolsHold, "hold", "FAIL")}; " +
        $"every φ of a pulse read at the 24 phases k·15° (the turn about z at k·30°, its half-angle at k·15°), every point at the five unit points.";

    public IEnumerable<IInspectable> Children
    {
        get
        {
            yield return Node("the labelling",
                $"the page's bits through a CNOT and a swap (b = bit 1, a = bit 1 XOR bit 2): I_z = Z_b + Z_a·Z_b/2 and ZZ → Z_a " +
                $"{Verdict(LabellingReproducesIz, "hold", "FAIL")}; every string to ± one string, the page's seven images among them, " +
                $"{Verdict(LabellingSignedImages, "holds", "FAILS")}; {StringsMovedInKz} of the 15 strings besides II land on another k_Z; " +
                $"a and b the other way round {Verdict(WrongAssignmentFails, "do not reproduce I_z", "DO reproduce I_z (FAILS)")}.");
            yield return Node("the coherence order",
                $"p = 2v₁ + v₂, k_Z = |v₁| + |v₂| and the popcount difference v₁ + v₂ on all sixteen cells " +
                $"{Verdict(OrderIsWeightedBitOrders, "hold", "FAIL")}; p = v₁ + 2v₂ {Verdict(SwappedWeightsFail, "fails", "does NOT fail")}; " +
                $"in the page's bits v = (±1, 0), (0, ±1) {Verdict(PageBitsDoNotFixOrder, "each carry two orders", "do NOT each carry two orders")}.");
            yield return Node("unitality and the thermal source",
                $"the generator {Verdict(GeneratorLinear, "is linear in (J₀, J₁, J₂, q, ω)", "is NOT linear")} (guarded at {string.Join(" and ", GenericPoints)}); " +
                $"L(II) = 0 and Tr L(S) = 0 {Verdict(UnitalAndTracePreserving, "at every point", "FAIL")}; amplitude damping on bit 1 " +
                $"{Verdict(AmplitudeDampingMovesIdentity, "keeps the traces and moves the identity", "does NOT behave as the control")}; " +
                $"−L(I_z) = 6J₂·ZI + 3J₁·IZ {Verdict(ThermalSourceHolds, "at every point, so the thermal source is (ε/4)(6J₂·ZI + 3J₁·IZ)", "FAILS")}.");
            yield return Node("the inversion recovery",
                $"ZI at −6J₂ and IZ at −6J₁ {Verdict(LongitudinalEigenoperators, "at every point", "FAIL")}; " +
                $"the rank-3 population in units of the equilibrium's rank-1 amplitude, on (1, e^(−6J₂t), e^(−6J₁t)): " +
                $"inversion {string.Join(", ", InversionAmplitudes)}, saturation {string.Join(", ", SaturationAmplitudes)} " +
                $"{Verdict(InversionRecoveryHolds, "(4/5 against 2/5, zero when J₁ = J₂)", "(NOT as stated)")}.");
            yield return Node("the turn about z",
                $"e^(−iφI_z) = e^(−iφZ) ⊗ e^(−iφZ/2), the cells turned by e^(−i(2v₁ + v₂)φ), and Z ⊗ diag(1, i)·e^(−3iπ/4) at 90° " +
                $"{Verdict(TurnFactorsInSouzaBits, "hold", "FAIL")}; in the page's bits e^(−iφZ_b)(cos(φ/2) − i·sin(φ/2)·Z_aZ_b) " +
                $"{Verdict(TurnFormInPageBits, "holds", "FAILS")}, operator Schmidt ranks at φ = 30°·k: {string.Join(", ", TurnRanksInPageBits)} " +
                $"{Verdict(TurnHolds, "(a product exactly at the multiples of 180° among them)", "(NOT as stated)")}.");
            yield return Node("the phase cycle",
                $"n steps keep p′ − p ≡ p₀ (mod n), n = 1 … 8, {Verdict(CycleKeepsResidue, "exactly", "FAIL")}; a pulse's phase on a pathway " +
                $"{Verdict(PathwayPhase, "is e^(−i(p′ − p)φ) at every φ, the turned pulse the pulse about its own axis at the 24 phases", "is NOT as stated")}; " +
                $"eight steps of the 90° pulse on T₃₀ {Verdict(EightStepsOnTheSpin, "keep exactly the order 3, 1 or 2, of which the spin has none at 2; two steps keep more", "do NOT keep one order")}.");
            yield return Node("the 180° pulse",
                $"i·(X·e^(2iφZ)) ⊗ (X·e^(iφZ)) {Verdict(Pulse180IsAProduct, "at every φ, |m⟩⟨m′| → |−m⟩⟨−m′|", "FAILS")}; every string on its k_Z " +
                $"{Verdict(Pulse180KeepsKz, "at every φ", "NOT kept")}; ± itself at 0° and 90° {Verdict(Pulse180SignedAtQuarterTurns, "(not at 45°)", "FAILS")}; " +
                $"the 90° pulse {Verdict(Pulse90MovesKz, "sends IZ to XY", "does NOT send IZ to XY")}.");
            yield return Node("the echo",
                $"the 180° pulse at 0° (ρ ↦ RρR†, R = i·XX) reverses ZI + IZ/2 and keeps ZZ{Verdict(EchoReversesOffsetKeepsSplitting, "", " (FAILS)")}; no plain one-site monomial " +
                $"sandwich reverses both {Verdict(NoPlainSandwichReversesBoth, "(all sixteen kinds read)", "(one DOES: FAILS)")}; ρ ↦ ρᵀ " +
                $"{Verdict(TransposeReversesBothAndIsNoPulse, "reverses each and both and is no pulse (its Choi matrix, read off it, has −1)", "FAILS")}; " +
                $"the Theodolite's reading of the same is the verifier's.");
            yield return Node("the 180° pulse in the page's bits",
                $"operator Schmidt ranks at 0°, 90°, 30°, 45°: {string.Join(", ", Pulse180RanksInPageBits)}, the pulse i·IX at 0° (a flip of b alone) and −i·ZY at 90°" +
                $"{Verdict(Pulse180PageBitsHolds, "", " (NOT as stated)")}.");
            yield return Node("the 90° pulse",
                $"I_z → −I_y, ZI → −(√3/2)·IY − YX/2, T₃₀ onto ±3 with 5/16 at each, and the turn's factoring{Verdict(Pulse90Images, "", " (FAILS)")}; operator Schmidt " +
                $"rank in Souza et al.'s bits {string.Join(", ", Pulse90RanksInSouzaBits)} at 0° and 15° (every phase by the factoring turn), " +
                $"in the page's bits {string.Join(", ", Pulse90RanksInPageBits)} at 0°, 30°, 45°, 90° {Verdict(Pulse90Holds, "(as stated)", "(NOT as stated)")}.");
            yield return Node("the receiver",
                $"I₊'s cells at v = (0, 1) (√3) and (1, −1) (2), squared norms 6 : 4{Verdict(ReceiverCells, "", " (FAIL)")}; Tr(I₊·ρ) " +
                $"{Verdict(ReceiverReadsRankOne, "reads only rank 1 at order −1 (T₂,₋₁, T₃,₋₁ and a 90° pulse's image of T₃₀ invisible)", "does NOT read as stated")}; " +
                $"relaxation at zero splitting and offset {Verdict(RelaxationCarriesRankThree, "carries a unit T₃,₋₁ into T₁,₋₁ with ±√6/5 on e^(−3(J₀ + J₁)t) and e^(−3(J₁ + J₂)t)", "does NOT carry as stated")}; " +
                $"the splitting and the dynamic shift {Verdict(SplittingAndShiftInView, "bring it into view at second and at first order, Tr(I₊·[H_LS, T]) = 48√3·(2K₁ − K₂) for T = 2(IX − i·IY) − √3(XX + YY − i(YX − XY)) ∝ T₃,₋₁", "do NOT as stated")}; " +
                $"operator Schmidt rank {ReceiverRanks[0]} in Souza et al.'s bits, {ReceiverRanks[1]} in the page's {Verdict(ReceiverRanks.SequenceEqual(new[] { 2, 3 }), "(as stated)", "(NOT as stated)")}.");
            yield return Node("the pools",
                $"|1⟩⟨0| at (1 − f₀)k and |0⟩⟨1| at f₀k on a label bit: X, Y at −k/2, Z at −k, the arrow k(2f₀ − 1) " +
                $"{Verdict(PoolTransferMatrix, "at every f₀ and k", "FAIL")}; the fractions {Verdict(PoolSteadyState, "are the steady state, the swapped rates not", "are NOT the steady state")}; " +
                $"the heights and the arrow are what LetterCube reads for a one-site jump (named, not run).");
            yield return Node("the field",
                $"ranks over ℚ(ζ₂₄) through the regular representation: the controls {Verdict(FieldRankControl, "read 1, 1, 2", "FAIL")}.");
        }
    }

    private static IInspectable Node(string name, string summary) =>
        new InspectableNode(displayName: name, summary: summary, provenance: NodeProvenance.Live);

    public InspectablePayload Payload => InspectablePayload.Empty;
}
