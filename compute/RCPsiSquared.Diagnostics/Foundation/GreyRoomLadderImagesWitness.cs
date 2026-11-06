using System.Globalization;
using System.Numerics;
using RCPsiSquared.Core.BlockSpectrum;
using RCPsiSquared.Core.CoherenceBlocks;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Pauli;
using ComplexMatrix = MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex>;

namespace RCPsiSquared.Diagnostics.Foundation;

/// <summary>The live lab for F171, the grey rooms inside the coloured blocks (claim
/// <see cref="GreyRoomLadderImagesClaim"/>, proof <c>docs/proofs/PROOF_GREY_GLASS_LADDER_IMAGES.md</c>, gates
/// <c>simulations/grey_glass_gate.py</c> and <c>simulations/grey_glass_ring_rooms.py</c>, the page
/// <c>docs/THE_GREY_GLASS.md</c>): under Z dephasing at one rate on every site the dissipator is one number
/// −2γk on each Hamming plane of a joint-popcount block, and the η ladder Φ and the spin ladder S⁻ carry the
/// boundary block (k, 0) of the XY chain, one-to-one, into every block (p, q) whose window holds k, a grey room
/// of dimension C(N, k) at −2γk, kept by the Hamiltonian part. This witness builds those images on every block
/// of its window and reads, at inspect time, where they lie, whether they are kept, and how large the room
/// around them is.
///
/// <para><b>What is built independently, and what is shared.</b> The chain or ring is a
/// <see cref="PauliHamiltonian"/> of XX + YY bonds, the wrap bond one more Pauli term, materialised as the full
/// 2^N matrix, and every block (p, q) is cut out of the Liouvillian by
/// <see cref="PerBlockLiouvillianBuilder.BuildBlockZ"/>, the committed road the F140 and F170 witnesses take.
/// The ladders are committed objects consumed as they stand: Φ is <see cref="SpectatorIntertwiner.BuildW"/>, the
/// spectator W of the typed parent F125, and S⁻ is the transpose of <see cref="SpectatorIntertwiner.BuildSPlus"/>,
/// the sideways spin ladder whose typed claim is adjacent to F171's and not a parent (its tier is Candidate for
/// a multiplet half F171 does not use). So the ladders, the Hamiltonian and the blocks share nothing with the
/// Python gate, which builds its ladders, sectors and generators in its own code, or with the sober twin
/// <c>MirrorWorld.GreyRoom</c>, which builds the Hamiltonian by hand in integers. The ROOM is not a second
/// route: it is the same incremental Krylov elimination the twin runs, written apart in each project,
/// <see cref="FreeFermionFloorWitness.KrylovRowRankModP"/> with the height a parameter, here at the free-fermion
/// floor's pair of primes, one of which differs from the twin's; for the rooms this road checks the construction
/// of H and the blocks, not the rank. The twin counts and does not build the ladders; this witness builds
/// them, and the keeping and the losses are its own reads.</para>
///
/// <para><b>What is exact and in which direction.</b> Every bond coupling and γ sit on the dyadic grid
/// 1/<see cref="Grid"/>, so Grid·L is an exact Gaussian-integer matrix, its diagonal the dissipator and its
/// off-diagonal −i times the integer commutator A = Grid·ad_H, read off the builder and checked symmetric. Four
/// readings per image, each with its own direction of certainty. (1) The image lies on the plane k: every
/// nonzero entry at Hamming distance k, exact. (2) The image is one-to-one: its rank modulo a prime equals
/// C(N, k); a rank modulo p can only drop, so the rational rank is at least that, and it is at most the number
/// of columns: exact. (3) Kept: the intertwining identity A_target·M = M·A_source holds entry for entry in
/// integers, M the image of the whole block (k, 0) by S⁻^s Φ^r, so the image of a kept block is kept: exact.
/// On the ring the source's wrap bond may carry the opposite sign from the target's, since each ladder step
/// turns the flux the particle number sets (PROOF_RING_FROZEN_BAND Lemma 1, typed in F170); both signs are
/// tried and the one that closes the identity is recorded, <see cref="RoomRead.SourceWrapSign"/>, and read
/// against the law of <see cref="WrapSignLaw"/>. Not kept: the rank of the image beside its motion exceeds
/// C(N, k) modulo a prime, which can only undercount, so the rational rank exceeds it too and the image is not
/// invariant: exact. (4) The room: dim(block) minus the rank of the Krylov rows the cells off the plane
/// generate under A, modulo two primes, an UPPER bound on the rational room. Where a kept, one-to-one image of
/// dimension C(N, k) meets a room bound of C(N, k) the rational room is exactly C(N, k):
/// <see cref="RoomRead.Certified"/>, exact at the coupling read.</para>
///
/// <para><b>The couplings read.</b> One bond profile, uniform by default and any dyadic profile on request
/// (<see cref="Bonds"/>, the ring's last entry its wrap bond); every read is invariant under a common scale of
/// the couplings, and γ enters only the dissipator's diagonal, so a loss or a room is read at the profile's shape.
/// The chain's statements hold at any real couplings and are read here at the profile given; the even ring's
/// losses at odd k, which the proof's §6 computes at N = 4, 6, 8, hold there at every nonzero wrap bond, since
/// the deciding map holds no coupling; the odd ring's losses are measured in the
/// proof at the uniform ring and at couplings of both signs, and a loss read here at another profile is a read
/// at that profile. The block (1, 1) is certified against the proof's C(N, 2) + 1 on the uniform open chain
/// alone; at another profile its certification, where it happens, is exact at that input and no theorem.</para>
///
/// <para><b>Two-sided, and theorems apart from reads.</b> The grey-room lemma's premise is read off the committed
/// dissipator rather than assumed: every diagonal entry of every block is −2γ times the cell's Hamming distance
/// (<see cref="LightIsOneNumberPerPlane"/>). The two ladders commute as integer matrices on every block where
/// both can step (<see cref="LaddersCommute"/>). What F171 proves at the profile read stands in
/// <see cref="Failures"/> when it fails and <see cref="Holds"/> says whether it did: every image on its plane and
/// one-to-one, every proved-kept image kept with its room at least C(N, k), the even ring's losses at odd k off
/// the corners (a theorem at N = 4 and 6 at every nonzero wrap bond, the deciding map holding no coupling), the
/// wrap-sign law, the uniform open chain's (1, 1), and on a palindromic chain the Heisenberg (1, 1) at least
/// ⌊N/2⌋ (<see cref="HeisenbergFloorHolds"/>: F140's frozen modes at exactly −4γ lie, under uniform light, on
/// the plane of distance 2 in the kernel of the Hamiltonian part, so they are in the room; where the mod-p room
/// meets that floor it is certified, <see cref="HeisenbergBlock11Certified"/>). What the page and the proof
/// MEASURE stands in <see cref="OpenReads"/> when it comes out otherwise and <see cref="ReadsAsThePage"/> says
/// whether it did; a read that differs is an open row at this profile, not a finding against F171: the odd
/// ring's losses off the corner routes (measured in the proof at the uniform ring and at couplings of both
/// signs), the controls that must fail, read as losses by rank and not as residuals since a residual can be a
/// block-constant shift with the image still kept: the spin ladder with its stagger dropped, the V of
/// PROOF_CODIM1_BY_ADDITIVITY §6, on every proved-kept image with a spin step that does not fill the corner
/// (0, k) (s ∉ {0, k}; where the spin steps fill it the unstaggered map fills the same block)
/// (<see cref="StaggerDroppedImagesLost"/> of <see cref="StaggerDroppedImagesRead"/>), and on the Heisenberg
/// chain at the same profile, where it has no zero bond, the ZZ term losing both ladders' images from N = 3 on
/// (<see cref="HeisenbergEtaImageLost"/>, <see cref="HeisenbergSpinImageLost"/>: on a connected chain the
/// one-particle ZZ diagonal differs at the ends; at N = 2 the ZZ energy is a function of the particle number
/// alone and both images survive, read the other way round; with zero bonds making equal dimers the diagonal is
/// constant and the spin image survives, so nothing is asked there), and on the uniform chain the Heisenberg
/// (1, 1) reading exactly ⌊N/2⌋, the page's read at N = 3 to 8 (<see cref="HeisenbergBlock11Room"/>). On the
/// ring the reduction of each popcount sector to a fixed-flux fermion ring is docked through
/// <see cref="RingFrozenBandWitness.SectorIdentityHolds"/>, which runs on F170's witness's own uniform ring and
/// not on these blocks (<see cref="FluxRuleDock"/>); the flux on THESE blocks is read by the wrap-sign law.</para>
///
/// <para><b>What stays a read.</b> The room's excess over the ladder image, where the image is kept and the
/// bound is larger (the chain's middle height k = N/2 strictly inside a window at couplings the reflection
/// keeps, 8 against 6 at N = 4 and 22 against 20 at N = 6, the reflection's products; the rings' extra rooms),
/// is an upper bound with the proved lower bound beside it and nothing certified between. The mirror products
/// themselves, the half-chain reading, the antipodal pairs and the sum-of-full-pairs criterion are the proof's
/// §5 and §8 and are not built here.</para>
///
/// <para><b>Running it.</b> <c>inspect --root greyroom [--N n] [--ring] [--J j | --bonds a,b,...] [--gamma g]</c>;
/// the sweep at N = 6, chain or ring, takes a few seconds. The tests,
/// <c>dotnet test compute/RCPsiSquared.Diagnostics.Tests --filter FullyQualifiedName~GreyRoomLadderImagesClaimTests</c>,
/// run in about half a minute, the single block (4, 2) at N = 8 included. When <see cref="Holds"/> is false,
/// <see cref="Failures"/> names every theorem that failed; when <see cref="ReadsAsThePage"/> is false,
/// <see cref="OpenReads"/> names every read that came out otherwise.</para></summary>
public sealed class GreyRoomLadderImagesWitness : IInspectable
{
    private static readonly CultureInfo Inv = CultureInfo.InvariantCulture;

    /// <summary>A chain needs two sites.</summary>
    public const int MinN = 2;

    /// <summary>A ring needs three: at N = 2 the wrap bond would double the one bond.</summary>
    public const int MinRingN = 3;

    /// <summary>Largest N whose every block is swept at construction; the middle block has 400 cells.</summary>
    public const int MaxSweepN = 6;

    /// <summary>Largest N a single block is read at (<see cref="ReadRoom"/>); the middle block at N = 8 has
    /// 4900 cells and the Krylov closure is cubic in it.</summary>
    public const int MaxN = 8;

    /// <summary>The dyadic grid every coupling and γ are snapped onto.</summary>
    public const long Grid = 1024;

    /// <summary>|J_b| ≤ 2 on every bond: the identities multiply grid-unit couplings by the images' integer
    /// entries, and this keeps every product far below 2^63.</summary>
    public const double MaxCoupling = 2.0;

    /// <summary>The free-fermion floor witness's pair of primes, below 2^30 so a product of two residues fits a
    /// long; one of them differs from the sober twin's pair. The commutator is a real integer matrix, so no
    /// square root of −1 is needed.</summary>
    private static readonly long[] Primes = { 1_000_000_007L, 998_244_353L };

    /// <summary>One image, the ladder image of (k, 0) in the block (p, q), with the room around it.
    /// <paramref name="SpinSteps"/> and <paramref name="EtaSteps"/> are s and r. <paramref name="RuleProvedKept"/>
    /// is F171's closed-form rule (<see cref="GreyRoomLadderImagesClaim.ImageProvedKept"/>);
    /// <paramref name="KeptExactly"/> the identity read here from the start (k, 0), and true without an identity read
/// where the image fills its block; <paramref name="SourceWrapSign"/>
    /// the sign of the source ring's wrap bond relative to the target's that closes it (+1 on the chain; +1 where
    /// either sign closes it, the starts (0, 0) and (N, 0), whose kets cannot hop, and a ring whose wrap bond is
    /// zero; 0 when no identity was read, the image filling its block or no sign closing it);
    /// <paramref name="LostByRank"/> the exact loss. <paramref name="RoomUpperBound"/> is the mod-p room,
    /// <paramref name="Certified"/> the meeting of the two bounds.</summary>
    public readonly record struct RoomRead(
        int P, int Q, int K, int BlockDimension, int PlaneDimension, int LadderDimension,
        int SpinSteps, int EtaSteps, bool FillsBlock, bool RuleProvedKept,
        bool ImageOnPlane, int ImageRankModP, bool KeptExactly, int SourceWrapSign,
        bool LostByRank, int RoomUpperBound, bool Certified);

    public int N { get; }
    public bool Ring { get; }
    public double Gamma { get; }

    /// <summary>The bond couplings read, snapped to the grid: N − 1 on the chain, N on the ring with the wrap
    /// bond last.</summary>
    public IReadOnlyList<double> Bonds { get; }

    /// <summary>Every bond equal, the wrap bond included.</summary>
    public bool Uniform { get; }

    /// <summary>The ring's wrap bond is nonzero; on the chain false (there is none). With it zero the ring is
    /// the open chain and F171 keeps every image.</summary>
    public bool WrapBondNonzero { get; }

    /// <summary>On the chain: the profile is a palindrome, so the chain's reflection keeps the couplings and
    /// the reflection's products at the middle height exist. Not asked on the ring.</summary>
    public bool ChainReflectionKeepsCouplings { get; }

    /// <summary>No bond is zero: the chain or ring is connected.</summary>
    public bool Connected { get; }

    /// <summary>Every (p, q, k) of every window at N, p and q both ways round (the transposes are read, not
    /// inferred).</summary>
    public IReadOnlyList<RoomRead> Rooms { get; }

    /// <summary>The grey-room lemma's premise off the builder: on every cell of every block the dissipator is
    /// −2γ times the Hamming distance, exactly.</summary>
    public bool LightIsOneNumberPerPlane { get; }

    /// <summary>Φ·S⁻ = S⁻·Φ as integer matrices on every block where both can step.</summary>
    public bool LaddersCommute { get; }

    /// <summary>On the ring: every popcount sector of the spin ring equals the fixed-flux fermion ring with
    /// sign w·(−1)^(p−1), both wraps, entry for entry (PROOF_RING_FROZEN_BAND Lemma 1, F170), docked through
    /// <see cref="RingFrozenBandWitness.SectorIdentityHolds"/> on that witness's own uniform ring; a dock of the
    /// lemma, not a read of these blocks. True on the chain, where there is nothing to reduce.</summary>
    public bool FluxRuleDock { get; }

    /// <summary>On the ring with a nonzero wrap bond: every kept image whose start can hop (0 &lt; k &lt; N) closed
    /// its identity at the sign <see cref="WrapSignLaw"/> gives. True on the chain.</summary>
    public bool WrapSignLawHolds { get; }

    /// <summary>The proved-kept images rebuilt with the spin ladder's stagger dropped (PROOF_CODIM1_BY_ADDITIVITY
    /// §6's V), over the routes with a spin step that does not fill the corner (0, k), s ∉ {0, k} (where the spin
    /// steps fill it the unstaggered map fills the same block and the image is kept anyway); and how many of them
    /// are lost by rank. A read: every one was lost at every profile tried, and the stagger is what makes the spin
    /// ladder commute with the Hamiltonian part.</summary>
    public int StaggerDroppedImagesRead { get; }
    public int StaggerDroppedImagesLost { get; }

    /// <summary>The mod-p room of the open Heisenberg chain's block (1, 1) at the distance 2, same profile: an
    /// upper bound. On a palindromic chain it is at least ⌊N/2⌋ (<see cref="HeisenbergFloorHolds"/>), and on the
    /// uniform chain the page reads exactly ⌊N/2⌋ at N = 3 to 8 (a mirror-symmetric non-uniform chain of six reads
    /// 3 too, a chain of four whose profile the reflection does not keep reads 1 against 2).</summary>
    public int HeisenbergBlock11Room { get; }

    /// <summary>On a palindromic chain: the Heisenberg (1, 1) room is at least ⌊N/2⌋, F140's frozen modes at
    /// exactly −4γ, which under uniform light lie on the plane of distance 2 in the kernel of the Hamiltonian part
    /// (a mode at Re λ = −4γ has all its weight at distance 2 by the Absorption cell cost, and then L_H v = 0). A
    /// theorem at every reversal-invariant coupling; true where nothing is asked (the ring, a profile the
    /// reflection does not keep).</summary>
    public bool HeisenbergFloorHolds { get; }

    /// <summary>On a palindromic chain: the mod-p room of the Heisenberg (1, 1) meets F140's floor ⌊N/2⌋, so the
    /// rational room is exactly ⌊N/2⌋ at this profile.</summary>
    public bool HeisenbergBlock11Certified { get; }

    /// <summary>On the Heisenberg chain at the same profile, Φ of (1, 0) into (2, 1) is lost by rank: true from
    /// N = 3 on, the ZZ term being quartic in the fermions; false at N = 2, where the ZZ energy is a function
    /// of the particle number alone and only shifts the identity by a constant.</summary>
    public bool HeisenbergEtaImageLost { get; }

    /// <summary>On the Heisenberg chain at the same profile, S⁻ of (2, 0) into (1, 1) is lost by rank: true
    /// from N = 3 on, the ZZ diagonal breaking ΣhΣ = −h; false at N = 2 for the same reason as Φ's.</summary>
    public bool HeisenbergSpinImageLost { get; }

    /// <summary>Every theorem of F171 at this profile that failed, named; empty when <see cref="Holds"/>.</summary>
    public IReadOnlyList<string> Failures { get; }

    /// <summary>Every read the page or the proof measured that came out otherwise here, named; empty when
    /// <see cref="ReadsAsThePage"/>. An entry is an open row at this profile, not a finding against F171.</summary>
    public IReadOnlyList<string> OpenReads { get; }

    public GreyRoomLadderImagesWitness(int n = 5, bool ring = false, double j = 0.75, double gamma = 0.25,
        IReadOnlyList<double>? bondCouplings = null)
    {
        int min = ring ? MinRingN : MinN;
        if (n < min || n > MaxSweepN)
            throw new ArgumentOutOfRangeException(nameof(n), n,
                $"the sweep reads {min} ≤ N ≤ {MaxSweepN}; a single block is read to N = {MaxN} by ReadRoom");
        var profile = Profile(n, ring, j, bondCouplings);
        long gNum = SnapGamma(gamma);

        N = n;
        Ring = ring;
        Gamma = gNum / (double)Grid;
        Bonds = profile.Select(b => b / (double)Grid).ToArray();
        Uniform = profile.All(b => b == profile[0]);
        WrapBondNonzero = ring && profile[n - 1] != 0;
        ChainReflectionKeepsCouplings = !ring && Enumerable.Range(0, n - 1).All(b => profile[b] == profile[n - 2 - b]);
        Connected = profile.All(b => b != 0);
        int wrap = ring ? 1 : 0;

        var cache = new BlockCache(n, profile, gNum);
        var rooms = new List<RoomRead>();
        bool light = true;
        for (int p = 0; p <= n; p++)
            for (int q = 0; q <= n; q++)
                foreach (int k in GreyRoomLadderImagesClaim.Window(n, p, q))
                {
                    rooms.Add(ReadRoom(cache, wrap, p, q, k, out bool lightHere));
                    light &= lightHere;
                }
        Rooms = rooms;
        LightIsOneNumberPerPlane = light;
        LaddersCommute = LaddersCommuteAt(n);
        FluxRuleDock = !ring || (RingFrozenBandWitness.SectorIdentityHolds(n, Math.Abs(profile[0]) == 0 ? Grid : Math.Abs(profile[0]), +1) &&
                                 RingFrozenBandWitness.SectorIdentityHolds(n, Math.Abs(profile[0]) == 0 ? Grid : Math.Abs(profile[0]), -1));
        WrapSignLawHolds = !WrapBondNonzero || rooms.All(r =>
            r.FillsBlock || !r.KeptExactly || r.K == 0 || r.K == n || r.SourceWrapSign == WrapSignLaw(n, r.SpinSteps, r.EtaSteps));

        int read = 0, lost = 0;
        foreach (var r in rooms)
        {
            if (!r.RuleProvedKept || r.SpinSteps == 0 || r.SpinSteps == r.K) continue;
            var target = cache.Get(wrap, r.P, r.Q, heisenberg: false);
            var unstaggered = LadderImage(n, (r.K, 0), r.SpinSteps, r.EtaSteps, dropStagger: true).Columns;
            read++;
            if (IsLostByRank(target, unstaggered, r.LadderDimension)) lost++;
        }
        StaggerDroppedImagesRead = read;
        StaggerDroppedImagesLost = lost;

        if (!ring)
        {
            var block11 = cache.Get(0, 1, 1, heisenberg: true);
            HeisenbergBlock11Room = RoomUpperBound(block11, 2);
            HeisenbergEtaImageLost = IsLostByRank(cache.Get(0, 2, 1, heisenberg: true),
                LadderImage(n, (1, 0), 0, 1, dropStagger: false).Columns, Choose(n, 1));
            HeisenbergSpinImageLost = IsLostByRank(block11,
                LadderImage(n, (2, 0), 1, 0, dropStagger: false).Columns, Choose(n, 2));
        }
        HeisenbergFloorHolds = ring || !ChainReflectionKeepsCouplings || HeisenbergBlock11Room >= n / 2;
        HeisenbergBlock11Certified = !ring && ChainReflectionKeepsCouplings && HeisenbergBlock11Room == n / 2;

        var failures = new List<string>();
        var reads = new List<string>();
        if (!LightIsOneNumberPerPlane) failures.Add("the dissipator is not −2γ·k on every cell of every block");
        if (!LaddersCommute) failures.Add("Φ and S⁻ do not commute as integer matrices");
        if (!FluxRuleDock) failures.Add("the sector reduction to fixed-flux rings fails on F170's witness's ring");
        if (!WrapSignLawHolds) failures.Add("a kept image closed its identity at a wrap sign the law does not give");
        foreach (var r in rooms)
        {
            string at = $"({r.P}, {r.Q}) k={r.K}";
            if (!r.ImageOnPlane) failures.Add($"{at}: the image leaves its plane");
            if (r.ImageRankModP != r.LadderDimension) failures.Add($"{at}: the image is not one-to-one (rank {r.ImageRankModP})");
            if (r.RuleProvedKept)
            {
                if (!r.FillsBlock && !r.KeptExactly) failures.Add($"{at}: proved kept, but no identity closes");
                if (r.RoomUpperBound < r.LadderDimension) failures.Add($"{at}: the room bound {r.RoomUpperBound} is below the kept image's C(N, k)");
            }
            else if (!r.LostByRank)
            {
                // the even ring's losses are a theorem at every nonzero wrap bond at the N swept; the odd ring's are measured
                if (n % 2 == 0) failures.Add($"{at}: not proved kept, and not shown lost by rank");
                else reads.Add($"{at}: not proved kept, and not shown lost by rank at this profile (the proof measured the loss at the uniform ring)");
            }
            if (!ring && Uniform && r.P == 1 && r.Q == 1 && !r.Certified)
                failures.Add($"{at}: the uniform open chain's (1, 1) is not certified at C(N, k)");
        }
        if (!HeisenbergFloorHolds)
            failures.Add($"the palindromic Heisenberg (1, 1) reads {HeisenbergBlock11Room} below F140's floor ⌊N/2⌋ = {n / 2}");
        if (StaggerDroppedImagesLost != StaggerDroppedImagesRead)
            reads.Add($"the stagger dropped keeps {StaggerDroppedImagesRead - StaggerDroppedImagesLost} of {StaggerDroppedImagesRead} images");
        if (!ring)
        {
            if (Uniform && HeisenbergBlock11Room != n / 2) reads.Add($"the uniform Heisenberg (1, 1) reads {HeisenbergBlock11Room} against ⌊N/2⌋ = {n / 2}");
            if (Connected)
            {
                bool expectLost = n >= 3;
                if (HeisenbergEtaImageLost != expectLost) reads.Add($"the connected Heisenberg chain {(expectLost ? "keeps" : "loses")} Φ's image at N = {n}");
                if (HeisenbergSpinImageLost != expectLost) reads.Add($"the connected Heisenberg chain {(expectLost ? "keeps" : "loses")} S⁻'s image at N = {n}");
            }
        }
        Failures = failures;
        OpenReads = reads;
    }

    /// <summary>Every theorem of F171 at this N and profile holds: every image on its plane and one-to-one, kept
    /// exactly with its room at least C(N, k) where proved kept, the even ring's losses at odd k, the uniform open
    /// chain's (1, 1) certified at both heights, the lemma's premise, the commuting ladders, the flux dock, the
    /// wrap-sign law, the palindromic Heisenberg floor. A false is a finding against the page at this N and
    /// profile, or an error; <see cref="Failures"/> says which read.</summary>
    public bool Holds => Failures.Count == 0;

    /// <summary>Every read the page or the proof measured came out as there: the odd ring's losses, the controls
    /// lost by rank, the uniform Heisenberg (1, 1) at ⌊N/2⌋. A false is an open row at this profile;
    /// <see cref="OpenReads"/> says which.</summary>
    public bool ReadsAsThePage => OpenReads.Count == 0;

    /// <summary>The number of (p, q, k) whose rational room is certified exactly C(N, k).</summary>
    public int CertifiedCount => Rooms.Count(r => r.Certified);

    /// <summary>The sign of the source ring's wrap bond, relative to the target's, at which the identity from
    /// the start (k, 0) closes on a ring with a nonzero wrap bond, for a start that can hop (0 &lt; k &lt; N):
    /// (−1)^(s+r), each ladder step turning the flux the particle number sets, on the even ring and on the odd
    /// ring's single-ladder routes; on the odd ring's mixed routes, which fill a corner on the way in one order of
    /// the two ladders, −(−1)^(s+r),
    /// since through the corner the comparison crosses from the start's ket to the target's bra, and on the odd
    /// ring Σh_fΣ = −h_{−f} adds one turn.</summary>
    public static int WrapSignLaw(int n, int spinSteps, int etaSteps)
    {
        int sign = (spinSteps + etaSteps) % 2 == 0 ? 1 : -1;
        return n % 2 == 1 && spinSteps > 0 && etaSteps > 0 ? -sign : sign;
    }

    // ================= the reads =================

    /// <summary>One block, one height, at any N up to <see cref="MaxN"/>: the ladder image of (k, 0) in
    /// (p, q) and the room around it, at the uniform coupling <paramref name="j"/> or the profile
    /// <paramref name="bondCouplings"/>.</summary>
    public static RoomRead ReadRoom(int n, bool ring, int p, int q, int k, double j = 0.75, double gamma = 0.25,
        IReadOnlyList<double>? bondCouplings = null)
    {
        int min = ring ? MinRingN : MinN;
        if (n < min || n > MaxN)
            throw new ArgumentOutOfRangeException(nameof(n), n, $"{min} ≤ N ≤ {MaxN}");
        if (!GreyRoomLadderImagesClaim.Window(n, p, q).Contains(k))
            throw new ArgumentOutOfRangeException(nameof(k), k, $"the block ({p}, {q}) of {n} sites holds no cell at Hamming distance {k}");
        var cache = new BlockCache(n, Profile(n, ring, j, bondCouplings), SnapGamma(gamma));
        return ReadRoom(cache, ring ? 1 : 0, p, q, k, out _);
    }

    private static RoomRead ReadRoom(BlockCache cache, int wrap, int p, int q, int k, out bool lightHere)
    {
        int n = cache.N;
        bool ring = wrap != 0;
        var target = cache.Get(wrap, p, q, heisenberg: false);
        lightHere = target.LightIsOneNumber;
        int dim = target.Dim;
        int plane = target.Hd.Count(h => h == k);
        int ladder = Choose(n, k);
        var (s, r) = GreyRoomLadderImagesClaim.LadderSteps(n, p, q, k);
        bool fills = GreyRoomLadderImagesClaim.IsCorner(n, p, q);
        bool rule = GreyRoomLadderImagesClaim.ImageProvedKept(n, p, q, k, ring, wrapBondNonzero: ring && cache.Bonds[n - 1] != 0);

        var image = LadderImage(n, (k, 0), s, r, dropStagger: false);
        if (image.TargetDim != dim)
            throw new InvalidOperationException($"the ladder image of ({k}, 0) by S⁻^{s}Φ^{r} does not land in ({p}, {q})");
        bool onPlane = true;
        foreach (var col in image.Columns)
            for (int c = 0; c < dim && onPlane; c++)
                if (col[c] != 0 && target.Hd[c] != k) onPlane = false;
        int rank = int.MinValue;
        foreach (long prime in Primes) rank = Math.Max(rank, RankModP(image.Columns, prime));

        // the exact identity from the start (k, 0), the source's wrap bond at either sign on the ring
        bool kept = fills;
        int wrapSign = 0;
        if (!fills)
            foreach (int w in ring ? new[] { 1, -1 } : new[] { 0 })
            {
                var source = cache.Get(w, k, 0, heisenberg: false);
                if (Residual(target.A, image.Columns, source.A) == 0)
                {
                    kept = true;
                    wrapSign = ring ? w : 1;
                    break;
                }
            }

        bool lost = !kept && IsLostByRank(target, image.Columns, ladder);
        int room = RoomUpperBound(target, k);
        bool certified = kept && onPlane && rank == ladder && room == ladder;
        return new RoomRead(p, q, k, dim, plane, ladder, s, r, fills, rule, onPlane, rank, kept, wrapSign, lost, room, certified);
    }

    /// <summary>The exact loss: the image beside its motion has rank above <paramref name="ladder"/> modulo a
    /// prime, which can only undercount, so over the rationals too, and the image is not invariant.</summary>
    private static bool IsLostByRank(BlockData target, long[][] columns, int ladder)
    {
        var moved = new long[2 * columns.Length][];
        for (int i = 0; i < columns.Length; i++)
        {
            moved[i] = columns[i];
            moved[columns.Length + i] = Apply(target.A, columns[i], target.Dim);
        }
        foreach (long prime in Primes)
            if (RankModP(moved, prime) > ladder) return true;
        return false;
    }

    private static int RoomUpperBound(BlockData block, int k)
    {
        int room = int.MaxValue;
        foreach (long prime in Primes)
            room = Math.Min(room, block.Dim - FreeFermionFloorWitness.KrylovRowRankModP(block.Hd, k, block.A, block.Dim, prime));
        return room;
    }

    /// <summary>Φ·S⁻ against S⁻·Φ, integer matrices from (p, q) to (p, q + 2), on every block where both can
    /// step.</summary>
    public static bool LaddersCommuteAt(int n)
    {
        for (int p = 1; p <= n - 1; p++)
            for (int q = 0; q <= n - 2; q++)
            {
                int dim = Choose(n, p) * Choose(n, q);
                var spinFirst = LadderImage(n, (p, q), 1, 1, dropStagger: false, spinFirst: true);
                var etaFirst = LadderImage(n, (p, q), 1, 1, dropStagger: false, spinFirst: false);
                for (int c = 0; c < dim; c++)
                    if (!spinFirst.Columns[c].SequenceEqual(etaFirst.Columns[c])) return false;
            }
        return true;
    }

    // ================= the atoms =================

    private static long Snap(double v) => (long)Math.Round(v * Grid);

    private static long SnapGamma(double gamma)
    {
        long gNum = Snap(gamma);
        if (gNum <= 0 || gamma > 1.0)
            throw new ArgumentOutOfRangeException(nameof(gamma), gamma,
                "γ must lie in (0, 1] and not snap to 0: the planes are the dissipator's eigenspaces at every γ > 0, " +
                "and γ = 0 has no light to be one number per plane");
        return gNum;
    }

    /// <summary>The bond profile on the grid: uniform <paramref name="j"/> on every bond, or the given
    /// couplings, N − 1 on the chain and N on the ring with the wrap bond last, each within ±<see cref="MaxCoupling"/>,
    /// not all zero.</summary>
    private static long[] Profile(int n, bool ring, double j, IReadOnlyList<double>? bondCouplings)
    {
        int bonds = ring ? n : n - 1;
        var profile = new long[bonds];
        if (bondCouplings is null)
        {
            if (!(Math.Abs(j) <= MaxCoupling) || Snap(j) == 0)
                throw new ArgumentOutOfRangeException(nameof(j), j, $"J must lie within ±{MaxCoupling} and not snap to 0");
            Array.Fill(profile, Snap(j));
            return profile;
        }
        if (bondCouplings.Count != bonds)
            throw new ArgumentException($"a {(ring ? "ring" : "chain")} on {n} sites carries {bonds} bonds; got {bondCouplings.Count}", nameof(bondCouplings));
        for (int b = 0; b < bonds; b++)
        {
            if (!(Math.Abs(bondCouplings[b]) <= MaxCoupling))
                throw new ArgumentOutOfRangeException(nameof(bondCouplings), bondCouplings[b], $"every coupling must lie within ±{MaxCoupling}");
            profile[b] = Snap(bondCouplings[b]);
        }
        if (profile.All(b => b == 0))
            throw new ArgumentException("every coupling snaps to 0: nothing moves", nameof(bondCouplings));
        return profile;
    }

    private static int Choose(int n, int k) => GreyRoomLadderImagesClaim.LadderDimension(n, k);

    /// <summary>A sparse integer matrix by column: Columns[c] lists the rows and values of column c.</summary>
    private sealed record SparseMap(List<(int Row, long Val)>[] Columns, int Rows);

    private sealed class BlockData
    {
        public required int Dim;
        public required int[] Hd;
        public required List<(int Row, long Val)>[] A;
        public required bool LightIsOneNumber;
    }

    /// <summary>The blocks of one N at one profile and γ, built once each, keyed by (wrap sign, p, q, book): the
    /// wrap sign multiplies the ring's last bond, 0 on the chain.</summary>
    private sealed class BlockCache
    {
        public int N { get; }
        public long[] Bonds { get; }
        private readonly long gNum;
        private readonly Dictionary<(int, int, int, bool), BlockData> blocks = new();
        private readonly Dictionary<(int, bool), ComplexMatrix> hamiltonians = new();

        public BlockCache(int n, long[] bonds, long gNum) { N = n; Bonds = bonds; this.gNum = gNum; }

        public BlockData Get(int wrapSign, int p, int q, bool heisenberg)
        {
            if (blocks.TryGetValue((wrapSign, p, q, heisenberg), out var b)) return b;
            if (!hamiltonians.TryGetValue((wrapSign, heisenberg), out var h))
            {
                h = Hamiltonian(N, Bonds, wrapSign, heisenberg);
                hamiltonians[(wrapSign, heisenberg)] = h;
            }
            b = Build(N, h, gNum, p, q);
            blocks[(wrapSign, p, q, heisenberg)] = b;
            return b;
        }
    }

    /// <summary>The chain or ring as Pauli terms: (J_b/2)(XX + YY) on every bond, the ring's wrap bond (N − 1, 0)
    /// at wrapSign·(J_w/2), and on the Heisenberg book (J_b/2)·ZZ beside them (the isotropic point, hop J_b).</summary>
    private static ComplexMatrix Hamiltonian(int n, long[] bonds, int wrapSign, bool heisenberg)
    {
        var terms = new List<PauliTerm>(3 * bonds.Length);
        for (int b = 0; b < bonds.Length; b++)
        {
            int s1 = b, s2 = (b + 1) % n;
            Complex c = (b == n - 1 ? wrapSign : 1) * bonds[b] / (double)Grid / 2.0;
            if (c == Complex.Zero) continue;
            terms.Add(PauliTerm.TwoSite(n, s1, PauliLetter.X, s2, PauliLetter.X, c));
            terms.Add(PauliTerm.TwoSite(n, s1, PauliLetter.Y, s2, PauliLetter.Y, c));
            if (heisenberg) terms.Add(PauliTerm.TwoSite(n, s1, PauliLetter.Z, s2, PauliLetter.Z, c));
        }
        return new PauliHamiltonian(n, terms).ToMatrix();
    }

    /// <summary>The block (p, q) cut from the Liouvillian, cells ket-outer and bra-inner over
    /// <see cref="BlockBasis.PopcountStates"/> ascending (the ladders' basis), its Hamming distances, and the
    /// integer commutator A = Grid·ad_H read off the off-diagonal (−i times an integer, checked symmetric) and,
    /// on the Heisenberg book, the diagonal's imaginary part; the real diagonal read against −2γ·k.</summary>
    private static BlockData Build(int n, ComplexMatrix h, long gNum, int p, int q)
    {
        var kets = BlockBasis.PopcountStates(n, p);
        var bras = BlockBasis.PopcountStates(n, q);
        int d = 1 << n;
        int dim = kets.Count * bras.Count;
        var flat = new int[dim];
        var hd = new int[dim];
        int idx = 0;
        foreach (long ket in kets)
            foreach (long bra in bras)
            {
                flat[idx] = (int)(ket * d + bra);
                hd[idx++] = BitOperations.PopCount((ulong)(ket ^ bra));
            }
        var gamma = Enumerable.Repeat(gNum / (double)Grid, n).ToArray();
        var block = PerBlockLiouvillianBuilder.BuildBlockZ(h, gamma, flat);

        var a = new List<(int Row, long Val)>[dim];
        for (int c = 0; c < dim; c++) a[c] = new List<(int, long)>();
        bool light = true;
        for (int col = 0; col < dim; col++)
            for (int row = 0; row < dim; row++)
            {
                Complex z = block[row, col];
                if (z == Complex.Zero) continue;
                long re = Exact(z.Real), im = Exact(z.Imaginary);
                if (row == col)
                {
                    if (re != -2 * gNum * hd[col]) light = false;
                    if (im != 0) a[col].Add((col, -im));
                }
                else
                {
                    if (re != 0) throw new InvalidOperationException($"entry ({row}, {col}) has a real off-diagonal part: H is not real");
                    a[col].Add((row, -im));
                }
            }
        for (int c = 0; c < dim; c++)
            foreach (var (row, val) in a[c])
                if (!a[row].Contains((c, val)))
                    throw new InvalidOperationException($"the commutator is not symmetric at ({row}, {c})");
        return new BlockData { Dim = dim, Hd = hd, A = a, LightIsOneNumber = light };
    }

    private static long Exact(double v)
    {
        double s = v * Grid;
        double r = Math.Round(s);
        if (s != r)
            throw new InvalidOperationException($"block entry {v.ToString("R", Inv)} is not a multiple of 1/{Grid}: the dyadic grid was broken");
        return (long)r;
    }

    private static SparseMap Sparse(ComplexMatrix m, bool transpose)
    {
        int rows = transpose ? m.ColumnCount : m.RowCount;
        int cols = transpose ? m.RowCount : m.ColumnCount;
        var columns = new List<(int Row, long Val)>[cols];
        for (int c = 0; c < cols; c++) columns[c] = new List<(int, long)>();
        for (int i = 0; i < m.RowCount; i++)
            for (int k = 0; k < m.ColumnCount; k++)
            {
                Complex z = m[i, k];
                if (z == Complex.Zero) continue;
                if (z.Imaginary != 0 || z.Real != Math.Round(z.Real))
                    throw new InvalidOperationException($"ladder entry {z} is not an integer");
                long v = (long)z.Real;
                if (transpose) columns[i].Add((k, v)); else columns[k].Add((i, v));
            }
        return new SparseMap(columns, rows);
    }

    /// <summary>S⁻ from (p, q) to (p − 1, q + 1): the transpose of the committed S⁺ from (p − 1, q + 1), the
    /// Hilbert-Schmidt adjoint in the orthonormal cell basis since the entries are real. With
    /// <paramref name="dropStagger"/> the same map without the (−1)^l, the transpose of PROOF_CODIM1_BY_ADDITIVITY
    /// §6's V, the control that must fail.</summary>
    private static SparseMap SpinDown(int n, int p, int q, bool dropStagger)
    {
        if (!dropStagger) return Sparse(SpectatorIntertwiner.BuildSPlus(n, p - 1, q + 1), transpose: true);
        var inP = BlockBasis.PopcountStates(n, p);
        var inQ = BlockBasis.PopcountStates(n, q);
        var outP = BlockBasis.PopcountStates(n, p - 1);
        var outQ = BlockBasis.PopcountStates(n, q + 1);
        var outPIdx = new Dictionary<long, int>();
        for (int i = 0; i < outP.Count; i++) outPIdx[outP[i]] = i;
        var outQIdx = new Dictionary<long, int>();
        for (int i = 0; i < outQ.Count; i++) outQIdx[outQ[i]] = i;
        var columns = new List<(int Row, long Val)>[inP.Count * inQ.Count];
        for (int i = 0; i < inP.Count; i++)
            for (int j = 0; j < inQ.Count; j++)
            {
                long a = inP[i], b = inQ[j];
                var col = columns[i * inQ.Count + j] = new List<(int, long)>();
                for (int l = 0; l < n; l++)
                {
                    // the site l leaves the ket and joins the bra; the JW string sign counts the sites before l,
                    // so it is the same on the states before and after the move
                    if (!SpectatorIntertwiner.Occupied(n, a, l) || SpectatorIntertwiner.Occupied(n, b, l)) continue;
                    long a2 = a & ~SpectatorIntertwiner.SiteMask(n, l), b2 = b | SpectatorIntertwiner.SiteMask(n, l);
                    int row = outPIdx[a2] * outQ.Count + outQIdx[b2];
                    col.Add((row, SpectatorIntertwiner.JwSign(n, a, l) * SpectatorIntertwiner.JwSign(n, b, l)));
                }
            }
        return new SparseMap(columns, outP.Count * outQ.Count);
    }

    private static SparseMap EtaUp(int n, int p, int q) => Sparse(SpectatorIntertwiner.BuildW(n, p, q), transpose: false);

    private static long[] Apply(List<(int Row, long Val)>[] byColumn, long[] v, int rows)
    {
        var outp = new long[rows];
        for (int c = 0; c < v.Length; c++)
        {
            if (v[c] == 0) continue;
            foreach (var (row, val) in byColumn[c]) outp[row] = checked(outp[row] + val * v[c]);
        }
        return outp;
    }

    /// <summary>The ladder image of the whole block <paramref name="source"/> after <paramref name="spin"/>
    /// steps of S⁻ and <paramref name="eta"/> steps of Φ, as integer columns over the target's cells, one per
    /// cell of the source. The spin steps first unless <paramref name="spinFirst"/> is false (the two commute,
    /// <see cref="LaddersCommuteAt"/>).</summary>
    private static (long[][] Columns, int TargetDim) LadderImage(int n, (int P, int Q) source, int spin, int eta,
        bool dropStagger, bool spinFirst = true)
    {
        int dim = Choose(n, source.P) * Choose(n, source.Q);
        var columns = new long[dim][];
        for (int c = 0; c < dim; c++) { columns[c] = new long[dim]; columns[c][c] = 1; }
        int p = source.P, q = source.Q;
        int targetDim = dim;
        void Step(SparseMap map)
        {
            for (int c = 0; c < dim; c++) columns[c] = Apply(map.Columns, columns[c], map.Rows);
            targetDim = map.Rows;
        }
        void SpinSteps() { for (int i = 0; i < spin; i++) { Step(SpinDown(n, p, q, dropStagger)); p--; q++; } }
        void EtaSteps() { for (int i = 0; i < eta; i++) { Step(EtaUp(n, p, q)); p++; q++; } }
        if (spinFirst) { SpinSteps(); EtaSteps(); } else { EtaSteps(); SpinSteps(); }
        return (columns, targetDim);
    }

    /// <summary>The largest entry of A_target·M − M·A_source, M the image's columns indexed by the source's
    /// cells; exact in integers.</summary>
    private static long Residual(List<(int Row, long Val)>[] aTarget, long[][] m, List<(int Row, long Val)>[] aSource)
    {
        int rows = m[0].Length;
        long worst = 0;
        for (int j = 0; j < m.Length; j++)
        {
            var left = Apply(aTarget, m[j], rows);
            var right = new long[rows];
            foreach (var (i, val) in aSource[j])
                for (int x = 0; x < rows; x++) right[x] = checked(right[x] + val * m[i][x]);
            for (int x = 0; x < rows; x++) worst = Math.Max(worst, Math.Abs(left[x] - right[x]));
        }
        return worst;
    }

    private static long Mod(long a, long p) { long m = a % p; return m < 0 ? m + p : m; }

    private static long ModPow(long b, long e, long p)
    {
        long acc = 1;
        b = Mod(b, p);
        while (e > 0)
        {
            if ((e & 1) == 1) acc = acc * b % p;
            b = b * b % p;
            e >>= 1;
        }
        return acc;
    }

    /// <summary>The rank over GF(p) of the vectors (rows), by forward elimination on copies.</summary>
    private static int RankModP(long[][] vectors, long p)
    {
        var a = vectors.Select(v => v.Select(x => Mod(x, p)).ToArray()).ToArray();
        int cols = a.Length == 0 ? 0 : a[0].Length;
        int rank = 0;
        for (int col = 0; col < cols && rank < a.Length; col++)
        {
            int piv = -1;
            for (int row = rank; row < a.Length; row++) if (a[row][col] != 0) { piv = row; break; }
            if (piv < 0) continue;
            (a[rank], a[piv]) = (a[piv], a[rank]);
            var pr = a[rank];
            long inv = ModPow(pr[col], p - 2, p);
            for (int row = rank + 1; row < a.Length; row++)
            {
                var ar = a[row];
                if (ar[col] == 0) continue;
                long f = ar[col] * inv % p;
                for (int k = col; k < cols; k++) if (pr[k] != 0) ar[k] = Mod(ar[k] - f * pr[k] % p, p);
            }
            rank++;
        }
        return rank;
    }

    // ================= the reading =================

    private static string Fmt(double v) => v.ToString("0.######", Inv);

    private string ProfileText => Uniform ? $"J={Fmt(Bonds[0])} on every bond" : $"J=[{string.Join(", ", Bonds.Select(Fmt))}]";

    public string DisplayName =>
        $"GreyRoomLadderImagesWitness (F171: the ladder images of the boundary blocks as grey rooms, live; N={N} {(Ring ? "ring" : "chain")})";

    public string Summary =>
        $"N={N} {(Ring ? "ring" : "open chain")}, {ProfileText}, γ={Fmt(Gamma)}: {Rooms.Count} (block, height) entries; " +
        $"every image on its plane: {Rooms.All(r => r.ImageOnPlane)}, one-to-one: {Rooms.All(r => r.ImageRankModP == r.LadderDimension)}; " +
        $"kept exactly {Rooms.Count(r => r.KeptExactly)} (proved kept {Rooms.Count(r => r.RuleProvedKept)}), " +
        $"lost by rank {Rooms.Count(r => r.LostByRank)} (not proved kept {Rooms.Count(r => !r.RuleProvedKept)}); " +
        $"rooms certified {CertifiedCount}; the light one number per plane: {LightIsOneNumberPerPlane}; " +
        $"the ladders commute: {LaddersCommute}; " +
        (Ring ? $"the flux dock: {FluxRuleDock}, the wrap-sign law: {WrapSignLawHolds}; "
              : $"(1, 1) at the distance 2: room ≤ {Rooms.Single(r => r.P == 1 && r.Q == 1 && r.K == 2).RoomUpperBound} against C(N, 2) = {Choose(N, 2)}; ") +
        (Holds ? "every theorem as F171 says" : $"{Failures.Count} theorem(s) not as F171 says: a finding at this N and profile, or an error") +
        (ReadsAsThePage ? "; every read as the page." : $"; {OpenReads.Count} read(s) not as the page: open rows at this profile.");

    public IEnumerable<IInspectable> Children
    {
        get
        {
            yield return new InspectableNode(
                $"the lemma's premise: the light is one number per plane = {LightIsOneNumberPerPlane}",
                summary: "every diagonal entry of every block cut by the committed builder is −2γ times the cell's " +
                         "Hamming distance, read exactly; the planes are the dissipator's eigenspaces, so a subspace " +
                         "the Hamiltonian part keeps inside one is grey (PROOF_GREY_GLASS_LADDER_IMAGES §2)",
                provenance: NodeProvenance.Live);
            yield return new InspectableNode(
                $"the two ladders commute = {LaddersCommute} (a theorem); the stagger dropped loses {StaggerDroppedImagesLost} of {StaggerDroppedImagesRead} images (a read: all, at every profile tried)",
                summary: "Φ·S⁻ = S⁻·Φ as integer matrices on every block where both can step (the two SU(2)s of the " +
                         "doubled picture); the spin ladder without its (−1)^l, PROOF_CODIM1_BY_ADDITIVITY §6's V, " +
                         "rebuilt on every proved-kept image with a spin step that does not fill the corner (0, k), " +
                         "is read lost by rank (a residual alone could be a constant shift with the image still kept)",
                provenance: NodeProvenance.Live);
            if (Ring)
            {
                yield return new InspectableNode(
                    $"the flux dock (PROOF_RING_FROZEN_BAND Lemma 1, F170): sectors equal the fixed-flux rings = {FluxRuleDock}",
                    summary: "every popcount sector of the spin ring with wrap w equals the fermion ring with wrap sign " +
                             "w·(−1)^(p−1), both wraps, entry for entry, read on F170's witness's own uniform ring (a dock " +
                             "of the lemma, not a read of these blocks); each ladder step turns that flux",
                    provenance: NodeProvenance.Live);
                yield return new InspectableNode(
                    $"the wrap-sign law on these blocks = {WrapSignLawHolds}",
                    summary: "every kept image whose start can hop closed A·M = M·A with the source ring's wrap bond at " +
                             "(−1)^(s+r) times the target's, each step turning the flux; on the odd ring's mixed routes, " +
                             "which pass through a corner, at −(−1)^(s+r), the comparison crossing from the start's ket " +
                             "to the target's bra where Σh_fΣ = −h_{−f} adds one turn",
                    provenance: NodeProvenance.Live);
            }
            for (int p = 0; p <= N; p++)
                for (int q = 0; q <= N; q++)
                {
                    var reads = Rooms.Where(r => r.P == p && r.Q == q).ToList();
                    yield return new InspectableNode(
                        $"block ({p}, {q}), {reads[0].BlockDimension} cells{(reads[0].FillsBlock ? ", a corner" : "")}: " +
                        string.Join("; ", reads.Select(Describe)),
                        provenance: NodeProvenance.Live);
                }
            if (!Ring)
                yield return new InspectableNode(
                    $"the Heisenberg control: (1, 1) at the distance 2 reads ≤ {HeisenbergBlock11Room} against ⌊N/2⌋ = {N / 2}" +
                    (ChainReflectionKeepsCouplings ? $", F140's floor on this palindromic chain (holds: {HeisenbergFloorHolds}; certified: {HeisenbergBlock11Certified})" : ", not a floor at a profile the reflection does not keep") +
                    $"; Φ's image lost {HeisenbergEtaImageLost}, S⁻'s lost {HeisenbergSpinImageLost} (a read: lost from N = 3 on where no bond is zero)",
                    summary: "on the open Heisenberg chain at the same profile the ZZ term loses both ladders' images by rank " +
                             "(Φ of (1, 0) into (2, 1), S⁻ of (2, 0) into (1, 1)) from N = 3 on where the chain has no zero bond, " +
                             "the one-particle ZZ diagonal differing at the ends; at N = 2 the ZZ energy is a function of the " +
                             "particle number and both survive, and with zero bonds making equal dimers the spin image survives; " +
                             "the block (1, 1) holds at the distance 2 an upper bound, at least F140's ⌊N/2⌋ on a palindromic " +
                             "chain (a theorem: a frozen mode at exactly −4γ lies on the plane of distance 2 in the kernel of " +
                             "the Hamiltonian part), exactly ⌊N/2⌋ on the uniform chain, the page's read at N = 3 to 8",
                    provenance: NodeProvenance.Live);
            yield return new InspectableNode(
                Holds ? "where the theorems fail: nowhere" : $"where the theorems fail: {Failures.Count} read(s)",
                summary: Holds ? "every theorem of F171 holds at this N and profile" : string.Join("; ", Failures),
                provenance: NodeProvenance.Live);
            yield return new InspectableNode(
                ReadsAsThePage ? "where the reads differ from the page: nowhere" : $"where the reads differ from the page: {OpenReads.Count} read(s), open rows at this profile",
                summary: ReadsAsThePage ? "the odd ring's losses, the controls and the uniform Heisenberg (1, 1) read as the page" : string.Join("; ", OpenReads),
                provenance: NodeProvenance.Live);
            yield return new InspectableNode("what stays a read, and what stays outside",
                summary: "a room above its ladder image is an upper bound with the proved C(N, k) below it and nothing " +
                         "certified between (the chain's middle height k = N/2 strictly inside a window at couplings the " +
                         "reflection keeps, 8 against 6 at N = 4 and 22 against 20 at N = 6, the reflection's products; " +
                         "the rings' rooms); the mirror products, the half-chain reading, the antipodal pairs and the " +
                         "sum-of-full-pairs criterion (the proof's §5 and §8) are not built here; a loss is exact at the " +
                         "N and profile read and claimed at no other");
        }
    }

    private string Describe(RoomRead r)
    {
        string kept;
        if (r.FillsBlock) kept = "fills the block (kept)";
        else if (r.KeptExactly)
        {
            string sign = !Ring ? "" : !WrapBondNonzero || r.K == 0 || r.K == N ? " (either wrap sign on the start)"
                : r.SourceWrapSign == -1 ? " with J_w → −J_w on the start" : " at the same J_w";
            kept = $"kept: A·M = M·A exact from ({r.K}, 0){sign}";
        }
        else if (r.LostByRank) kept = "LOST: the image beside its motion has rank > C(N, k) mod p, so A does not keep it (exact)";
        else kept = "neither read (an open row at this profile)";
        string rule = r.RuleProvedKept == (r.KeptExactly || r.FillsBlock) ? "" : " [AGAINST THE RULE]";
        string room = r.Certified ? $"room = {r.RoomUpperBound} of the plane's {r.PlaneDimension} (certified)"
            : r.KeptExactly || r.FillsBlock ? $"room in [{r.LadderDimension}, {r.RoomUpperBound}] of the plane's {r.PlaneDimension}"
            : $"room ≤ {r.RoomUpperBound} of the plane's {r.PlaneDimension}, no ladder bound";
        return $"k={r.K}: image C(N, {r.K}) = {r.LadderDimension}, rank {r.ImageRankModP}, on the plane {r.ImageOnPlane}, {kept}{rule}; {room}";
    }

    public InspectablePayload Payload => InspectablePayload.Empty;
}
