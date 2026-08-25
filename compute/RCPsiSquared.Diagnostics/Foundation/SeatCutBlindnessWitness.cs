using System.Globalization;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Symmetry;

namespace RCPsiSquared.Diagnostics.Foundation;

/// <summary>The live lab for F157, the blind seat (claim <see cref="SeatCutBlindnessClaim"/>,
/// proof <c>docs/proofs/PROOF_BLIND_SEAT_SPAN_AND_NODE_LEMMA.md</c>, gates
/// <c>simulations/blind_seat_span_proof.py</c> and <c>simulations/seat_cut_blindness.py</c>): put
/// the Z-dephasing on ONE seat of a chain and count, at inspect time, what it cannot touch.
///
/// <code>
///     blind(j)  =  N − rank[ e_j, H e_j, H² e_j, ..., H^N e_j ]
///     dim ker L_SE(j)  =  1 + dim commutant(H restricted to the Krylov complement)
/// </code>
///
/// <para><b>What is recomputed and what is written down.</b> Everything numeric here is
/// recomputed: the integer single-excitation matrix is built from the bond profile, the Krylov
/// rank is taken by elimination over GF(p) at two primes, and the kernel dimension is a second,
/// independent elimination on the masked commutator. What is written down is the CLOSED FORM the
/// count is checked against, <see cref="SeatCutBlindnessClaim.BlindHeisenberg"/> and
/// <see cref="SeatCutBlindnessClaim.BlindXy"/>, which applies only on a uniform profile and is
/// reported as "not applicable" otherwise. Two independent computations meeting, in the house
/// pattern of <see cref="QuditPartialPalindromeWitness"/>.</para>
///
/// <para><b>Ranks, never a spectrum, and the one-sidedness is disclosed.</b> A rank over GF(p)
/// can only DROP relative to the rank over ℚ, so a blind count read this way can only ever be too
/// LARGE and a kernel dimension likewise. Two primes are taken and the LARGER rank is kept for the
/// count, the SMALLER nullity for the span, which is the tightening direction in both cases. Note
/// what that leaves standing: the span and the count are each a best-available upper bound, they
/// may come from DIFFERENT primes, and both can in principle be too large together, which is why
/// the proof file's Theorem A gate takes both sides over ℚ instead and this witness is a live
/// reading rather than the theorem's evidence. There is no tolerance anywhere in this file because no quantity
/// compared is a float: a float SVD rank of the same kernel returns 21 instead of 6 at N = 11,
/// J = 10⁻⁵ while reporting a singular-value gap of 5.95·10³.</para>
///
/// <para><b>Two guards, for two different costs.</b> The COUNT is an N × N elimination and walks
/// past the 2^N wall, so it is allowed to <see cref="MaxN"/> = 200. The SPAN eliminates on a
/// system with about N² columns, so it is computed only up to <see cref="MaxSpanN"/> = 12 and
/// reported as not computed above that. Neither guard is physics.</para>
///
/// <para><b>Falsifiers are printed beside the counts, not implied.</b> Every reading here sits
/// beside a zero the same routine produces: both END seats of a zero-free chain are blind 0 by
/// Lemma J1, and on the Heisenberg book a chain whose length shares no odd divisor with any
/// 2j+1 has no blind seat at all. A count without its paired zero would not distinguish a working
/// rank routine from one that returns N.</para>
///
/// <para>Children: the per-seat table (live), the closed-form comparison (live), the span identity
/// (live, only at N ≤ <see cref="MaxSpanN"/>), the falsifiers (live), the parity forcing (live), and
/// the scope fences. Root: <c>inspect --root blind</c>.</para></summary>
public sealed class SeatCutBlindnessWitness : IInspectable
{
    private static readonly CultureInfo Inv = CultureInfo.InvariantCulture;

    /// <summary>Two primes for the rank. The larger rank is kept, a rank mod p being able only to
    /// drop. Both are below 2^31 so that a product of two reduced residues stays inside
    /// <see cref="long"/> with room to spare.
    ///
    /// <para>A rank over GF(p) is blind to a factor of p in every entry, so a coupling divisible by
    /// a ranking prime would report the whole space blind if only that prime were used. TWO primes
    /// is the defence, and it is exact rather than probabilistic here: divisibility by BOTH at once
    /// needs |J| ≈ 2147483647·999999937, which <see cref="MaxCoupling"/> refuses outright.</para></summary>
    private static readonly long[] Primes = { 2147483647L, 999999937L };

    /// <summary>The largest chain the live COUNT will run on. The count is an N × N elimination,
    /// so this is generous on purpose: the law is checked past the 2^N wall.</summary>
    public const int MaxN = 200;

    /// <summary>The largest |J| a bond may carry. <see cref="BuildH"/> doubles a coupling for the hop
    /// and sums up to N−1 of them on the ZZ diagonal, so this keeps the matrix well inside int64 with
    /// room for the largest guarded N. Scale, not physics: the commutant does not move under it.</summary>
    public const long MaxCoupling = 1L << 40;

    /// <summary>The largest chain the live SPAN will run on. The span eliminates on a system with
    /// about N² columns, so it is capped far lower than the count.</summary>
    public const int MaxSpanN = 12;

    public int N { get; }
    public SeatCutBook Book { get; }

    /// <summary>The integer bond profile, N−1 entries. Zero entries are permitted and are where
    /// the two-halves gcd phrasing stops being valid; the count and the span theorem do not care.</summary>
    public IReadOnlyList<long> Bonds { get; }

    /// <summary>Whether every bond carries the same nonzero coupling, which is the hypothesis of
    /// the two closed forms. False for an irregular or a zero-bearing profile.</summary>
    public bool IsUniform { get; }

    private readonly long[,] _h;
    private readonly int[] _blind;

    public SeatCutBlindnessWitness(int n = 7, SeatCutBook book = SeatCutBook.Heisenberg,
                                   IReadOnlyList<long>? bonds = null)
    {
        if (n < 2)
            throw new ArgumentOutOfRangeException(nameof(n), $"a chain needs at least 2 sites; got N = {n}.");
        if (n > MaxN)
            throw new ArgumentOutOfRangeException(nameof(n),
                $"N = {n} exceeds the live-count guard {MaxN}; the count is an N x N elimination and the guard " +
                "is cost, not physics.");
        if (bonds is not null && bonds.Count != n - 1)
            throw new ArgumentException($"a chain on {n} sites carries {n - 1} bonds; got {bonds.Count}.", nameof(bonds));
        if (bonds is not null)
            foreach (long b in bonds)
                if (Math.Abs(b) > MaxCoupling)
                    throw new ArgumentOutOfRangeException(nameof(bonds),
                        $"|J| = {Math.Abs(b)} exceeds {MaxCoupling}: the ZZ diagonal sums up to N-1 couplings and " +
                        "the hop doubles one, so a larger magnitude could wrap int64 silently. The commutant is " +
                        "scale-invariant, so divide the profile through instead.");

        N = n;
        Book = book;
        Bonds = bonds is null ? Enumerable.Repeat(1L, n - 1).ToArray() : bonds.ToArray();
        IsUniform = Bonds.Count > 0 && Bonds[0] != 0 && Bonds.All(b => b == Bonds[0]);

        _h = BuildH();
        _blind = new int[N];
        for (int seat = 0; seat < N; seat++) _blind[seat] = ComputeBlind(seat);
    }

    /// <summary>The site-indexed integer single-excitation matrix, in the repo's integer
    /// convention (<c>simulations/seat_cut_blindness.py</c>'s <c>se_hamiltonian_int</c>): a hop of
    /// 2J on each bond, and on the Heisenberg book a diagonal paying −J at a bond's own two ends
    /// and +J elsewhere, bond by bond. On the XY book the diagonal is zero, which is what makes the
    /// parity forcing possible.</summary>
    private long[,] BuildH()
    {
        var h = new long[N, N];
        for (int b = 0; b < N - 1; b++)
        {
            h[b, b + 1] += 2 * Bonds[b];
            h[b + 1, b] += 2 * Bonds[b];
        }
        if (Book == SeatCutBook.Heisenberg)
            for (int s = 0; s < N; s++)
                for (int b = 0; b < N - 1; b++)
                    h[s, s] += (s == b || s == b + 1) ? -Bonds[b] : Bonds[b];
        return h;
    }

    /// <summary>blind(seat), recomputed: N minus the rank of the seat's Krylov matrix, the larger
    /// rank over the two primes.</summary>
    public int Blind(int seat)
    {
        if (seat < 0 || seat >= N)
            throw new ArgumentOutOfRangeException(nameof(seat), $"seat must lie in 0..{N - 1}; got {seat}.");
        return _blind[seat];
    }

    private int ComputeBlind(int seat)
    {
        int rank = 0;
        foreach (long p in Primes) rank = Math.Max(rank, KrylovRankModP(seat, p));
        return N - rank;
    }

    private int KrylovRankModP(int seat, long p)
    {
        var rows = new List<long[]>();
        var vec = new long[N];
        vec[seat] = 1;
        for (int k = 0; k <= N; k++)
        {
            rows.Add((long[])vec.Clone());
            var next = new long[N];
            for (int a = 0; a < N; a++)
            {
                long s = 0;
                for (int b = 0; b < N; b++)
                    s = (s + Mod(_h[a, b], p) * vec[b]) % p;
                next[a] = s;
            }
            vec = next;
        }
        return RankModP(rows, N, p);
    }

    /// <summary>dim ker L_SE(seat), recomputed by a second, independent elimination: the kernel is
    /// { rho : [H, rho] = 0 and rho clear on the seat's off-diagonal row and column }, so neither
    /// γ's magnitude nor a scale of H enters. Returns null above <see cref="MaxSpanN"/>.</summary>
    public int? Span(int seat)
    {
        if (seat < 0 || seat >= N)
            throw new ArgumentOutOfRangeException(nameof(seat), $"seat must lie in 0..{N - 1}; got {seat}.");
        if (N > MaxSpanN) return null;
        int nullity = int.MaxValue;
        foreach (long p in Primes) nullity = Math.Min(nullity, SpanNullityModP(seat, p));
        return nullity;
    }

    private int SpanNullityModP(int seat, long p)
    {
        var allowed = new List<(int I, int J)>();
        var index = new Dictionary<(int, int), int>();
        for (int i = 0; i < N; i++)
            for (int j = 0; j < N; j++)
                if ((i == seat) == (j == seat))
                {
                    index[(i, j)] = allowed.Count;
                    allowed.Add((i, j));
                }

        var rows = new List<long[]>();
        for (int a = 0; a < N; a++)
            for (int b = 0; b < N; b++)
            {
                var row = new long[allowed.Count];
                bool any = false;
                foreach (var (i, j) in allowed)
                {
                    long v = 0;
                    if (j == b) v += _h[a, i];
                    if (i == a) v -= _h[j, b];
                    if (v != 0) { row[index[(i, j)]] = Mod(v, p); any = true; }
                }
                if (any) rows.Add(row);
            }
        return allowed.Count - RankModP(rows, allowed.Count, p);
    }

    /// <summary>The closed form for this seat, or null when the profile is not uniform and the
    /// closed form therefore does not apply.</summary>
    public int? ClosedForm(int seat) =>
        IsUniform ? SeatCutBlindnessClaim.Blind(N, seat, Book) : null;

    /// <summary>How many seats the live count and the closed form agree on. Meaningless off a
    /// uniform profile, where it returns null.</summary>
    public int? ClosedFormAgreements()
    {
        if (!IsUniform) return null;
        int ok = 0;
        for (int seat = 0; seat < N; seat++)
            if (_blind[seat] == SeatCutBlindnessClaim.Blind(N, seat, Book)) ok++;
        return ok;
    }

    /// <summary>The seats this profile leaves blind, live.</summary>
    public IReadOnlyList<int> BlindSeats() =>
        Enumerable.Range(0, N).Where(s => _blind[s] > 0).ToArray();

    private static long Mod(long x, long p)
    {
        long r = x % p;
        return r < 0 ? r + p : r;
    }

    private static long ModPow(long b, long e, long p)
    {
        long r = 1;
        b = Mod(b, p);
        while (e > 0)
        {
            if ((e & 1) == 1) r = r * b % p;
            b = b * b % p;
            e >>= 1;
        }
        return r;
    }

    private static long ModInverse(long a, long p) => ModPow(a, p - 2, p);

    /// <summary>Rank over GF(p) by Gauss-Jordan. Entries are reduced into [0, p) on entry and kept
    /// there, so every product of two of them stays below p² and inside <see cref="long"/>.</summary>
    private static int RankModP(List<long[]> rows, int cols, long p)
    {
        var work = rows.Select(r => r.Select(x => Mod(x, p)).ToArray()).ToList();
        int rank = 0;
        for (int c = 0; c < cols && rank < work.Count; c++)
        {
            int piv = -1;
            for (int r = rank; r < work.Count; r++)
                if (work[r][c] != 0) { piv = r; break; }
            if (piv < 0) continue;
            (work[rank], work[piv]) = (work[piv], work[rank]);
            long inv = ModInverse(work[rank][c], p);
            for (int k = c; k < cols; k++) work[rank][k] = work[rank][k] * inv % p;
            for (int r = 0; r < work.Count; r++)
            {
                if (r == rank) continue;
                long f = work[r][c];
                if (f == 0) continue;
                for (int k = c; k < cols; k++)
                    work[r][k] = ((work[r][k] - f * work[rank][k] % p) % p + p) % p;
            }
            rank++;
        }
        return rank;
    }

    private string BookLabel => Book == SeatCutBook.Heisenberg ? "Heisenberg (ZZ on)" : "XY (ZZ off)";

    private string ProfileLabel =>
        IsUniform ? $"uniform J = {Bonds[0].ToString(Inv)}"
                  : "[" + string.Join(", ", Bonds.Select(b => b.ToString(Inv))) + "]";

    public string DisplayName =>
        $"SeatCutBlindnessWitness (F157 live lab, N={N}, {BookLabel}, {ProfileLabel})";

    public string Summary
    {
        get
        {
            var blindSeats = BlindSeats();
            string seats = blindSeats.Count == 0
                ? "no seat is blind"
                : "blind seats " + string.Join(", ", blindSeats.Select(s => $"{s}:{_blind[s]}"));
            int? agree = ClosedFormAgreements();
            string law = agree is null
                ? "closed form not applicable (profile is not uniform)"
                : $"closed form agrees at {agree}/{N} seats";
            return $"N={N} {BookLabel}, {ProfileLabel}: {seats}; {law}";
        }
    }

    public IEnumerable<IInspectable> Children
    {
        get
        {
            // 1. The per-seat table: the live rank beside the written closed form.
            var rows = new List<IInspectable>();
            for (int seat = 0; seat < N && seat < 24; seat++)
            {
                int? law = ClosedForm(seat);
                int? span = Span(seat);
                string lawPart = law is null ? "closed form n/a" :
                    (law == _blind[seat] ? $"closed form {law} MATCH" : $"closed form {law} MISMATCH");
                string spanPart = span is null
                    ? $"span not computed (N > {MaxSpanN})"
                    : (span == 1 + _blind[seat]
                        ? $"dim ker = {span} = 1 + blind"
                        : $"dim ker = {span}, which is NOT 1 + blind = {1 + _blind[seat]}; H is degenerate on the Krylov complement");
                rows.Add(new InspectableNode($"seat {seat}",
                    summary: $"blind = {_blind[seat]} (Krylov rank at two primes); {lawPart}; {spanPart}",
                    provenance: NodeProvenance.Live));
            }
            yield return new InspectableNode("per seat: the live count, the closed form, the span",
                summary: N > 24
                    ? $"the first 24 of {N} seats; the count itself runs to N = {MaxN}"
                    : $"all {N} seats, recomputed at inspect time",
                children: rows,
                provenance: NodeProvenance.Live);

            // 2. The two independent computations meeting (or not).
            int? agreements = ClosedFormAgreements();
            yield return new InspectableNode("the live count vs the closed form",
                summary: agreements is null
                    ? $"the profile {ProfileLabel} is not uniform, so the closed forms " +
                      "(gcd(2j+1, N) - 1)/2 and gcd(j+1, N+1) - 1 do not apply; the Krylov count still does, " +
                      "and the criterion deg gcd(chi(H), chi(H struck at j)) still gives it"
                    : (agreements == N
                        ? $"MATCH at all {N} seats: two independent computations meet, an exact GF(p) Krylov rank " +
                          $"and the closed form of the {BookLabel} book"
                        : $"MISMATCH: the closed form agrees at only {agreements} of {N} seats, which on a uniform " +
                          "zero-free chain would be a finding about the construction and not a tolerance"),
                provenance: NodeProvenance.Live);

            // 3. The span identity, and the criterion behind it.
            if (N <= MaxSpanN)
            {
                int held = 0, broke = 0;
                for (int seat = 0; seat < N; seat++)
                {
                    int? span = Span(seat);
                    if (span is null) continue;
                    if (span == 1 + _blind[seat]) held++; else broke++;
                }
                bool zeroFree = Bonds.All(b => b != 0);
                yield return new InspectableNode("the span identity (Theorem A, Corollary C)",
                    summary: $"dim ker L_SE(j) = 1 + blind(j) holds at {held} of {N} seats and breaks at {broke}. " +
                             (zeroFree
                                ? "This profile has no zero bond, so the Jacobi node lemma makes H simple and " +
                                  "Corollary B says it must hold everywhere; a break here would be a finding."
                                : "This profile HAS a zero bond, so the chain is cut and its pieces may repeat " +
                                  "each other; where that puts a repeated eigenvalue inside the Krylov complement the " +
                                  "identity breaks, and the kernel is 1 + dim commutant(H|_W) instead. The zero " +
                                  "bond is not itself the discriminator: [1,1,0,1,1] breaks it and [1,0,1] holds."),
                    provenance: NodeProvenance.Live);
            }

            // 4. The falsifiers: zeros the same routine has to produce.
            yield return new InspectableNode("the falsifiers (zeros from the same routine)",
                summary: $"end seats, live: blind(0) = {_blind[0]}, blind({N - 1}) = {_blind[N - 1]}. " +
                         (Bonds.All(b => b != 0)
                            ? "On a zero-free chain both must be 0 by Lemma J1, no eigenvector of an unreduced " +
                              "Jacobi matrix vanishing at an end."
                            : "This profile has a zero bond, so an end may be cut off and Lemma J1 does not apply.") +
                         $" Seats blind on this profile: {(BlindSeats().Count == 0 ? (Book == SeatCutBook.Heisenberg ? "none, which is the generic case on the Heisenberg book" : "none, which on the XY book means N is even, since an odd chain has its odd seats forced") : string.Join(", ", BlindSeats()))}.",
                provenance: NodeProvenance.Live);

            // 5. The parity forcing, which is a statement about the OTHER book.
            var forced = Enumerable.Range(0, N)
                .Where(s => SeatCutBlindnessClaim.ParityForcedXy(N, s, Book)).ToArray();
            yield return new InspectableNode("the parity-forced third kind (XY only)",
                summary: Book == SeatCutBook.Xy
                    ? (N % 2 == 1
                        ? $"N is odd, so every odd seat is blind at EVERY zero-free profile: {string.Join(", ", forced)}. " +
                          $"Live check on this profile: {string.Join(", ", forced.Select(s => $"{s}:{_blind[s]}"))}."
                        : "N is even, so parity forces nothing here; blindness on this book is then the ordinary " +
                          "coincidence the gcd law counts.")
                    : "not applicable on the Heisenberg book: its diagonal is nonzero, so the odd-size singularity " +
                      "argument does not run, and blindness off the uniform chain is a coincidence a mirror can " +
                      "force and disorder can only meet.",
                provenance: NodeProvenance.Live);

            // 6. Scope, written down rather than recomputed.
            yield return new InspectableNode("scope and fences",
                summary: "ONE seat and the single-excitation sector: the span theorem's cyclic-vector step needs " +
                         "the projector to come from the same vector that generates the Krylov space, and above " +
                         "popcount 1 the Wedderburn blocks are not forced. At a seat whose own ray is H-invariant " +
                         "(an isolated seat, or its bonds detuned to zero) the count stops equalling the blind " +
                         "SUBSPACE by one, though the Krylov-gcd identity and the span theorem hold anyway. The " +
                         "two-halves gcd phrasing needs a chain with no zero bond. Proof: " +
                         "docs/proofs/PROOF_BLIND_SEAT_SPAN_AND_NODE_LEMMA.md.");
        }
    }

    public InspectablePayload Payload => InspectablePayload.Empty;
}
