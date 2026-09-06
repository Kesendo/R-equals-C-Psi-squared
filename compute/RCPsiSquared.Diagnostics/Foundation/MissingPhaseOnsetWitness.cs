using System;
using System.Collections.Generic;
using System.Linq;
using System.Numerics;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Numerics;

namespace RCPsiSquared.Diagnostics.Foundation;

/// <summary>The lab for §8 of <c>experiments/THE_MOTION_AND_THE_MISSING_PHASE.md</c>: the three readings
/// of a broken mirror do not share an onset, and each one's first nonzero Taylor order and leading
/// coefficient is an exact Gaussian rational. §9 of that page states the recurrence and names this
/// object as the home of its exact-arithmetic checks at N = 5, 7, 9, 11, 13; they run here, at inspect
/// time, rather than being re-derived by hand.
/// <para>What is DERIVED (on the page, not here): the closed forms of §8, the (−2i)^m hop-parity
/// argument of §8.2, the rank-two output structure of §8.3, and the Dyson integral of §8.4. What is
/// RECOMPUTED LIVE (here): the §9 coefficient recurrence in exact ℚ(i), the first nonzero order of each
/// reading, each leading coefficient, and the comparison of every one of them against the page's closed
/// form. The comparison is exact equality in ℚ(i); a nonzero residual is a finding about the derivation
/// or about this code, never a tolerance to widen.</para>
/// <para>The method has a home: <c>docs/proofs/PROOF_F94_BORN_DOMINANT_FOUR_THIRDS.md</c> and
/// <c>docs/proofs/PROOF_F96_BORN_SUBDOMINANT_SLOPES.md</c> take exact short-time Lindblad coefficients
/// and apply the partial trace to the COEFFICIENT before the reading; F96 states the recipe generally.
/// The block this runs on is the (1,1) corner with a site-resolved rate, the object of
/// <c>docs/proofs/PROOF_R90_FROZEN_DIVISOR.md</c> §1, here with all the light on the centre seat and one
/// end bond detuned. The B block's rule is GammaFold's veil at Σ_l γ_l = γ (MirrorWorld
/// <c>GammaFold</c>); B is not a density matrix, and its negative-rate solution is algebra, not a
/// dephasing experiment.</para>
/// <para>Cost: the recurrence runs to order 2m+2 on N×N matrices over ℚ(i), which is cheap; what grows
/// is the height of the rationals, so <see cref="MaxSites"/> caps N. The finite-time companion at N = 7
/// is the committed producer <c>simulations/missing_phase_long_time.py</c>, a different route (spectral,
/// floating) to different quantities; nothing here reads its output.</para></summary>
public sealed class MissingPhaseOnsetWitness : IInspectable
{
    /// <summary>The cost guard, a policy cap rather than a measured ceiling: the page's checks stop at
    /// N = 13 (see <see cref="PageSites"/>), where the deepest reading runs to order 13. What grows with N
    /// is the height of the rationals rather than the dimension, so a caller wanting more should raise
    /// this deliberately and measure, not assume the gap between 13 and 21 has been walked.</summary>
    public const int MaxSites = 21;

    public BigRational Epsilon { get; }
    public BigRational Gamma { get; }

    /// <summary>The five odd chain lengths §9 names, at the couplings §9 names.</summary>
    public static readonly int[] PageSites = { 5, 7, 9, 11, 13 };

    public MissingPhaseOnsetWitness() : this(new BigRational(1, 3), new BigRational(2, 5)) { }

    public MissingPhaseOnsetWitness(BigRational epsilon, BigRational gamma)
    {
        if (gamma.Sign <= 0)
            throw new ArgumentOutOfRangeException(nameof(gamma), gamma.ToString(),
                "§8's two distance forms carry |ε|γ and the page sets γ > 0; at γ < 0 that oracle changes sign " +
                "while the computed distance, a norm, does not, so the comparison would be between two " +
                "different claims. The leakage form is linear in γ and flips WITH the computed value, so it " +
                "is not what this guard is about, and at γ = 0 there is nothing to read at all. " +
                "The recurrence itself accepts any γ through Coefficients, which is how the veil at −γ is run.");
        Epsilon = epsilon;
        Gamma = gamma;
    }

    // ---- the system of §8.1 -------------------------------------------------------------------

    /// <summary>The N × N hopping matrix h of §8.1: zero diagonal, h[0,1] = h[1,0] = 2(1+ε), every other
    /// neighbour pair 2. The page's convention puts the matrix element at 2J, and only the left end bond
    /// carries the defect.</summary>
    public static GaussianRational[,] Hopping(int n, BigRational epsilon)
    {
        Guard(n);
        var h = new GaussianRational[n, n];
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++)
                h[i, j] = GaussianRational.Zero;
        var end = new BigRational(2) * (BigRational.One + epsilon);
        h[0, 1] = end; h[1, 0] = end;
        for (int j = 1; j < n - 1; j++)
        {
            h[j, j + 1] = new BigRational(2);
            h[j + 1, j] = new BigRational(2);
        }
        return h;
    }

    /// <summary>A(0) = B(0) = |d₀⟩⟨d₀| with |d₀⟩ = (|0⟩ − |N−1⟩)/√2, the reflection-odd end state of §1.
    /// The √2 squares away, so the block is rational.</summary>
    public static GaussianRational[,] InitialBlock(int n)
    {
        Guard(n);
        var a = Zeros(n);
        var half = new BigRational(1, 2);
        a[0, 0] = half;
        a[n - 1, n - 1] = half;
        a[0, n - 1] = -half;
        a[n - 1, 0] = -half;
        return a;
    }

    /// <summary>The §9 recurrence, run in exact ℚ(i) to <paramref name="maxOrder"/>:
    /// A_{k+1} = [−i(hA_k − A_k h) + γ(zA_k z − A_k)]/(k+1) and the same with −γ(zB_k z + B_k) for B,
    /// plus the matched γ = 0 reference. Index k of each returned list is the coefficient of t^k.</summary>
    public static (GaussianRational[][,] A, GaussianRational[][,] B, GaussianRational[][,] Reference)
        Coefficients(int n, BigRational epsilon, BigRational gamma, int maxOrder)
        => Coefficients(n, epsilon, gamma, maxOrder, (n - 1) / 2);

    /// <summary>The same recurrence with the light on an arbitrary seat rather than on the centre. The
    /// page's system is the centre case; the other seats exist so that a reading can be shown to depend
    /// on where the light sits, which is what makes the centre gates say something.</summary>
    public static (GaussianRational[][,] A, GaussianRational[][,] B, GaussianRational[][,] Reference)
        Coefficients(int n, BigRational epsilon, BigRational gamma, int maxOrder, int seat)
    {
        Guard(n);
        if (maxOrder < 0) throw new ArgumentOutOfRangeException(nameof(maxOrder));
        if (seat < 0 || seat >= n) throw new ArgumentOutOfRangeException(nameof(seat));
        var h = Hopping(n, epsilon);
        int c = seat;

        var a = new GaussianRational[maxOrder + 1][,];
        var b = new GaussianRational[maxOrder + 1][,];
        var reference = new GaussianRational[maxOrder + 1][,];
        a[0] = InitialBlock(n);
        b[0] = InitialBlock(n);
        reference[0] = InitialBlock(n);

        for (int k = 0; k < maxOrder; k++)
        {
            a[k + 1] = Step(h, a[k], c, gamma, plus: true, k);
            b[k + 1] = Step(h, b[k], c, gamma, plus: false, k);
            reference[k + 1] = Step(h, reference[k], c, BigRational.Zero, plus: true, k);
        }
        return (a, b, reference);
    }

    /// <summary>One recurrence step. <paramref name="plus"/> selects the A rule γ(zXz − X) over the B
    /// rule −γ(zXz + X); on a cell (i,j) the two read −2γ·[exactly one index is the centre] and
    /// −2γ·[the two indices AGREE at the centre], the cell (c,c) belonging to the second and not to the
    /// first. The B rule is the Lattice's turned rule read site by site.</summary>
    private static GaussianRational[,] Step(GaussianRational[,] h, GaussianRational[,] x, int c,
                                            BigRational gamma, bool plus, int k)
    {
        int n = h.GetLength(0);
        var next = Zeros(n);
        var minusI = new GaussianRational(BigRational.Zero, -BigRational.One);
        var divisor = new GaussianRational(new BigRational(k + 1), BigRational.Zero);

        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++)
            {
                var commutator = GaussianRational.Zero;
                for (int s = 0; s < n; s++)
                {
                    if (!h[i, s].IsZero) commutator += h[i, s] * x[s, j];
                    if (!h[s, j].IsZero) commutator -= x[i, s] * h[s, j];
                }
                // z = I − 2|c⟩⟨c|, so (z X z)[i,j] = z_i z_j X[i,j] with z_i = −1 exactly at i = c.
                int zz = ((i == c ? -1 : 1) * (j == c ? -1 : 1));
                var factor = plus ? new BigRational(zz - 1) : new BigRational(-(zz + 1));
                next[i, j] = (minusI * commutator + new GaussianRational(gamma * factor, BigRational.Zero) * x[i, j])
                             / divisor;
            }
        return next;
    }

    // ---- the three readings of §7, order by order ----------------------------------------------

    /// <summary>What one (N, ε, γ) returns: the first nonzero Taylor order of each reading and its exact
    /// leading coefficient. <c>DOutSquared</c> is the SQUARE of d_out's coefficient, because d_out
    /// carries a √2 and the square is what stays inside ℚ.</summary>
    public sealed record Reading(int Sites, int M, int DOutOrder, BigRational? DOutCoefficientSquared,
                                 int DTwoOrder, BigRational DTwoCoefficient,
                                 int LeakageOrder, BigRational LeakageCoefficient,
                                 GaussianRational CentreEndA, GaussianRational CentreEndB);

    /// <summary>Runs the recurrence and reads the three onsets. Every number returned is exact.</summary>
    public static Reading Read(int n, BigRational epsilon, BigRational gamma)
        => ReadWithLightOnSeat(n, epsilon, gamma, (n - 1) / 2);

    /// <summary>The three readings with the light on <paramref name="seat"/>. The ENCODING does not move
    /// with the light: the decoder still traces out the chain centre, because the code of §1 is built
    /// from the chain's reflection and not from where the dephasing happens to sit.</summary>
    public static Reading ReadWithLightOnSeat(int n, BigRational epsilon, BigRational gamma, int seat)
    {
        Guard(n);
        if (seat < 0 || seat >= n) throw new ArgumentOutOfRangeException(nameof(seat));
        int m = (n - 1) / 2;
        int c = m;
        // The deepest reading the page names is δp_out at order 2m+1; one spare order lets a wrong
        // onset show up as a LATER first-nonzero rather than as a silent truncation.
        int maxOrder = 2 * m + 2;
        var (a, b, reference) = Coefficients(n, epsilon, gamma, maxOrder, seat);

        var deltaA = new GaussianRational[maxOrder + 1][,];
        var deltaB = new GaussianRational[maxOrder + 1][,];
        for (int k = 0; k <= maxOrder; k++)
        {
            deltaA[k] = Subtract(a[k], reference[k]);
            deltaB[k] = Subtract(b[k], reference[k]);
        }

        int dTwoOrder = -1;
        var dTwo = BigRational.Zero;
        for (int k = 0; k <= maxOrder && dTwoOrder < 0; k++)
        {
            var value = TwoSiteDistance(deltaA[k], n);
            if (!value.IsZero) { dTwoOrder = k; dTwo = value; }
        }

        int dOutOrder = -1;
        BigRational? dOutSquared = null;
        for (int k = 0; k <= maxOrder && dOutOrder < 0; k++)
            if (!DecodedDifferenceVanishes(deltaA[k], deltaB[k], n, c))
            {
                dOutOrder = k;
                // The order is exact whatever the shape of the difference; the VALUE exists only in the
                // rank-two case of §8.3, and a trace norm that is not representable is reported as absent
                // rather than skipped, which would move the onset to a later order that is not the onset.
                dOutSquared = DecodedIsRankTwo(deltaA[k], c, n)
                            ? DecodedDistanceSquared(deltaA[k], deltaB[k], n, c)
                            : (BigRational?)null;
            }

        int leakageOrder = -1;
        var leakage = BigRational.Zero;
        for (int k = 0; k <= maxOrder && leakageOrder < 0; k++)
        {
            var value = LeakageDifference(deltaA[k], n);
            if (!value.IsZero) { leakageOrder = k; leakage = value; }
        }

        // maxOrder is 2m+2 and m is at least 2, so order m+1 is always inside the computed range.
        var centreEndA = deltaA[m + 1][c, 0];
        var centreEndB = deltaB[m + 1][c, 0];

        return new Reading(n, m, dOutOrder, dOutSquared, dTwoOrder, dTwo,
                           leakageOrder, leakage, centreEndA, centreEndB);
    }

    /// <summary>d_2's coefficient at one Taylor order: max over site pairs of ½‖Tr_outside(δρ)‖₁, which
    /// §8.2 reduces exactly to |s|/2 + max(|s|/2, |v|) with s = δA_aa + δA_bb and v = Re δA_ab.
    /// PRECONDITION, satisfied by every physical coefficient beyond order 0 and to be given deliberately
    /// to a hand-built input: Tr δA = 0. The reduced pair's |00⟩⟨00| and |11⟩⟨11| entries are −s/2 only
    /// because the sites outside the pair carry −s between them. The B
    /// blocks drop out of a two-site trace by F70, and the imaginary part of A_ab cancels between the
    /// two copies; both are why this is a rational and not a norm needing an eigensolver.</summary>
    public static BigRational TwoSiteDistance(GaussianRational[,] deltaA, int n)
    {
        var best = BigRational.Zero;
        for (int p = 0; p < n; p++)
            for (int q = p + 1; q < n; q++)
            {
                var s = deltaA[p, p].Re + deltaA[q, q].Re;
                var v = deltaA[p, q].Re;
                var halfS = Abs(s) / new BigRational(2);
                var value = halfS + Max(halfS, Abs(v));
                if (Greater(value, best)) best = value;
            }
        return best;
    }

    /// <summary>d_out's coefficient SQUARED at one Taylor order, or null when the decoded difference at
    /// that order is not the rank-two object §8.3 derives. The decoded output of §8.3 puts δA on the
    /// non-centre block and on |f⟩⟨f|, and δB on the coupling between |f⟩ and the non-centre states; at
    /// the leading order only the latter survives, as |f⟩⟨w| + |w⟩⟨f| with w the two end amplitudes, so
    /// the trace norm is 2‖w‖ and the reading is ‖w‖. Returning null rather than a number is the point:
    /// it says this order is NOT the rank-two case, and the caller must not read a coefficient off it.
    /// </summary>
    public static BigRational DecodedDistanceSquared(GaussianRational[,] deltaA, GaussianRational[,] deltaB,
                                                     int n, int c)
    {
        if (!DecodedIsRankTwo(deltaA, c, n))
            throw new InvalidOperationException(
                "The decoded difference at this order is not the rank-two object of §8.3: its A part does " +
                "not vanish, so its trace norm needs the eigenvalues of a matrix this method does not " +
                "diagonalise. Ask DecodedIsRankTwo first.");

        // The B part: the coupling |f⟩⟨j| for j ≠ c. Its squared norm is Σ_j |δB_{c,j}|², and the trace
        // norm of |f⟩⟨w| + |w⟩⟨f| is 2‖w‖, so the reading ½‖·‖₁ is ‖w‖ itself.
        var normSquared = BigRational.Zero;
        for (int j = 0; j < n; j++)
            if (j != c) normSquared += deltaB[c, j].NormSquared;
        return normSquared;
    }

    /// <summary>Whether the decoded difference at one order is the rank-two object §8.3 derives: its A
    /// part, the non-centre block together with the |f⟩⟨f| entry, vanishes, leaving only the coupling
    /// between |f⟩ and the non-centre states. On the page's own system this holds at the leading order;
    /// off the centre seat it does not, which is the case that must not be silently passed over.
    /// <para>This is a SUFFICIENT condition for the closed form below, not a characterisation of when an
    /// exact value exists. With the non-centre block zero and δA_cc = a ≠ 0 the difference is
    /// a|f⟩⟨f| + |f⟩⟨w| + |w⟩⟨f|, whose ½‖·‖₁ is ½√(a² + 4‖w‖²) and whose square is rational too; that
    /// case is excluded here rather than handled, because it does not arise on the page's system and an
    /// untested branch is worth less than an absent one.</para></summary>
    public static bool DecodedIsRankTwo(GaussianRational[,] deltaA, int c, int n)
    {
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++)
                if (i != c && j != c && !deltaA[i, j].IsZero)
                    return false;
        return deltaA[c, c].IsZero;
    }

    /// <summary>Whether the decoded difference at one order is exactly zero. This is decidable at every
    /// order and for every seat, which is why the ONSET is exact even where the trace norm is not
    /// representable.</summary>
    public static bool DecodedDifferenceVanishes(GaussianRational[,] deltaA, GaussianRational[,] deltaB,
                                                 int n, int c)
    {
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++)
                if (i != c && j != c && !deltaA[i, j].IsZero)
                    return false;
        if (!deltaA[c, c].IsZero) return false;
        for (int j = 0; j < n; j++)
            if (j != c && !deltaB[c, j].IsZero) return false;
        return true;
    }

    /// <summary>δp_out's coefficient at one Taylor order. §8.4 has p_out = Tr(QA) with Q = (I+R)/2 and
    /// R the chain reflection, so the coefficient is Tr(Q δA) and is read off the coefficient matrix
    /// directly, with no propagation and no integral.</summary>
    public static BigRational LeakageDifference(GaussianRational[,] deltaA, int n)
    {
        Guard(n);
        // Tr(Q δA) = ½·Tr(δA) + ½·Σ_j δA[j, N−1−j].
        var total = GaussianRational.Zero;
        for (int j = 0; j < n; j++) total += deltaA[j, j] + deltaA[n - 1 - j, j];
        if (!total.Im.IsZero)
            throw new InvalidOperationException(
                $"Tr(Q·δA) came out with an imaginary part {total.Im}. Both factors are Hermitian, so this " +
                "is a finding about the recurrence, not a part to discard.");
        return total.Re * new BigRational(1, 2);
    }

    // ---- the reduction itself, against the physical generator ------------------------------------

    /// <summary>The cost guard for the full-space route: 2^N × 2^N exact matrices, so N = 7 is already
    /// 128 × 128 with a commutator at every order. The reduction is a claim about the generator and not
    /// about N, so one N where it is checked against the physical Lindbladian settles it.</summary>
    public const int MaxSitesFullSpace = 7;

    /// <summary>The Taylor coefficients of the FULL 2^N × 2^N density matrix under the page's own
    /// Lindbladian, L(ρ) = −i[H,ρ] + γ(Z_c ρ Z_c − ρ) with H = Σ_l J_l (X_l X_{l+1} + Y_l Y_{l+1}) at Δ = 0,
    /// seeded with the physical state ½·W[[A₀,B₀],[B₀,A₀]]W†. This is the route §8.1's reduction is a
    /// claim ABOUT: the two-block form is worth nothing unless it reproduces this, and the reduced
    /// generator agreeing with a doubled generator built from the same h and z would not show it.</summary>
    public static GaussianRational[][,] FullLindbladCoefficients(int n, BigRational epsilon, BigRational gamma,
                                                                 int maxOrder)
    {
        Guard(n);
        if (maxOrder < 0) throw new ArgumentOutOfRangeException(nameof(maxOrder));
        if (n > MaxSitesFullSpace)
            throw new ArgumentOutOfRangeException(nameof(n), n,
                $"the full-space route stops at N = {MaxSitesFullSpace}; 2^N × 2^N exact matrices grow past use.");
        int dim = 1 << n;
        int c = (n - 1) / 2;

        // H as its nonzero list. The XX+YY term on a bond exchanges an excited and an empty site at
        // amplitude 2J_l, and nothing else moves at Δ = 0.
        var hops = new List<(int Row, int Col, BigRational Value)>();
        for (int l = 0; l < n - 1; l++)
        {
            var coupling = l == 0 ? new BigRational(2) * (BigRational.One + epsilon) : new BigRational(2);
            if (coupling.IsZero) continue;
            for (int state = 0; state < dim; state++)
            {
                bool low = (state & (1 << l)) != 0, high = (state & (1 << (l + 1))) != 0;
                if (low == high) continue;
                int moved = state ^ (1 << l) ^ (1 << (l + 1));
                hops.Add((moved, state, coupling));
            }
        }

        var rho = new GaussianRational[maxOrder + 1][,];
        rho[0] = PhysicalSeed(n);
        var minusI = new GaussianRational(BigRational.Zero, -BigRational.One);

        for (int k = 0; k < maxOrder; k++)
        {
            var next = new GaussianRational[dim, dim];
            for (int i = 0; i < dim; i++)
                for (int j = 0; j < dim; j++)
                {
                    int parity = (((i >> c) & 1) ^ ((j >> c) & 1)) == 0 ? 1 : -1;
                    next[i, j] = new GaussianRational(gamma * new BigRational(parity - 1), BigRational.Zero)
                                 * rho[k][i, j];
                }
            foreach (var (row, col, value) in hops)
            {
                var weight = new GaussianRational(value, BigRational.Zero);
                for (int j = 0; j < dim; j++) next[row, j] += minusI * weight * rho[k][col, j];   // −i·Hρ
                for (int i = 0; i < dim; i++) next[i, col] -= minusI * weight * rho[k][i, row];   // +i·ρH
            }
            var divisor = new GaussianRational(new BigRational(k + 1), BigRational.Zero);
            for (int i = 0; i < dim; i++)
                for (int j = 0; j < dim; j++)
                    next[i, j] /= divisor;
            rho[k + 1] = next;
        }
        return rho;
    }

    /// <summary>The physical initial state, built WITHOUT <see cref="Embed"/> and without
    /// <see cref="InitialBlock"/>, so that the k = 0 comparison of the reduction gate is a check rather
    /// than a tautology. It is ½|Φ⟩⟨Φ| with |Φ⟩ = |d₀⟩ + F|d₀⟩, the §2 preparation |+⟩⟨+|_q ⊗ |d₀⟩⟨d₀|
    /// read in the computational basis: four amplitudes ±1/√2 at e₀, e_{N−1} and their global flips, so
    /// every entry of the state is ±¼ and the √2 never appears.</summary>
    public static GaussianRational[,] PhysicalSeed(int n)
    {
        Guard(n);
        int dim = 1 << n, all = dim - 1;
        var rho = new GaussianRational[dim, dim];
        for (int i = 0; i < dim; i++)
            for (int j = 0; j < dim; j++)
                rho[i, j] = GaussianRational.Zero;

        var support = new (int Index, int Sign)[]
        {
            (1, +1), (1 << (n - 1), -1), (all ^ 1, +1), (all ^ (1 << (n - 1)), -1),
        };
        var quarter = new BigRational(1, 4);
        foreach (var (a, sa) in support)
            foreach (var (b, sb) in support)
                rho[a, b] = new GaussianRational(quarter * new BigRational(sa * sb), BigRational.Zero);
        return rho;
    }

    /// <summary>The same embedding applied to an arbitrary block pair, so a reduced coefficient can be
    /// lifted back into the full space and compared there.</summary>
    public static GaussianRational[,] Embed(GaussianRational[,] a, GaussianRational[,] b, int n)
    {
        Guard(n);
        int dim = 1 << n, all = dim - 1;
        var rho = new GaussianRational[dim, dim];
        for (int i = 0; i < dim; i++)
            for (int j = 0; j < dim; j++)
                rho[i, j] = GaussianRational.Zero;
        var half = new GaussianRational(new BigRational(1, 2), BigRational.Zero);
        for (int j = 0; j < n; j++)
            for (int k = 0; k < n; k++)
            {
                int ket = 1 << j, bra = 1 << k;
                rho[ket, bra] += half * a[j, k];
                rho[all ^ ket, all ^ bra] += half * a[j, k];
                rho[ket, all ^ bra] += half * b[j, k];
                rho[all ^ ket, bra] += half * b[j, k];
            }
        return rho;
    }

    /// <summary>How many cells of how many orders the two-block reduction and the full Lindbladian
    /// disagree on, at one N. Zero is the claim of §8.1. The full-space arm is seeded by
    /// <see cref="PhysicalSeed"/>, which is built from the state rather than from the blocks, so the
    /// k = 0 comparison already tests the embedding instead of restating it.</summary>
    public static int ReductionMismatchCount(int n, BigRational epsilon, BigRational gamma, int maxOrder)
    {
        var full = FullLindbladCoefficients(n, epsilon, gamma, maxOrder);
        var (a, b, _) = Coefficients(n, epsilon, gamma, maxOrder);
        int dim = 1 << n, mismatches = 0;
        for (int k = 0; k <= maxOrder; k++)
        {
            var lifted = Embed(a[k], b[k], n);
            for (int i = 0; i < dim; i++)
                for (int j = 0; j < dim; j++)
                    if (lifted[i, j] != full[k][i, j]) mismatches++;
        }
        return mismatches;
    }

    // ---- the closed forms of §8, for the comparison ---------------------------------------------

    /// <summary>§8's closed form for d_out's leading coefficient, squared: [√2·m·2^m/(m+1)!·|ε|γ]².</summary>
    public static BigRational PageDOutSquared(int m, BigRational epsilon, BigRational gamma)
    {
        var coefficient = new BigRational(m) * Pow2(m) / Factorial(m + 1) * Abs(epsilon) * gamma;
        return new BigRational(2) * coefficient * coefficient;
    }

    /// <summary>§8's closed form for d_2's leading coefficient: 2^m/(m+1)!·|ε|γ at even m, and
    /// 2^(m+1)/(m+1)!·|ε|γ·max(1,|1+ε|) at odd m, where the order is one higher.</summary>
    public static BigRational PageDTwo(int m, BigRational epsilon, BigRational gamma)
    {
        if (m % 2 == 0)
            return Pow2(m) / Factorial(m + 1) * Abs(epsilon) * gamma;
        var stretch = Max(BigRational.One, Abs(BigRational.One + epsilon));
        return Pow2(m + 1) / Factorial(m + 1) * Abs(epsilon) * gamma * stretch;
    }

    /// <summary>§8's closed form for δp_out's leading coefficient: (−1)^m·2^(2m+1)/(2m+1)!·ε²γ.</summary>
    public static BigRational PageLeakage(int m, BigRational epsilon, BigRational gamma)
    {
        var magnitude = Pow2(2 * m + 1) / Factorial(2 * m + 1) * epsilon * epsilon * gamma;
        return m % 2 == 0 ? magnitude : -magnitude;
    }

    /// <summary>§8.2's closed form for the first dephasing-sensitive centre/end entry of δA:
    /// −γε(−2i)^m/(m+1)!, at order m+1.</summary>
    public static GaussianRational PageCentreEndA(int m, BigRational epsilon, BigRational gamma)
    {
        var scalar = new GaussianRational(-gamma * epsilon / Factorial(m + 1), BigRational.Zero);
        return scalar * PowMinusTwoI(m);
    }

    /// <summary>§8.3's closed form for the same entry of δB: −mγε(−2i)^m/(m+1)!. The extra factor m is
    /// the veil's, B_γ = e^(−2γt)·A_{−γ}, differentiated m+1 times at the origin.</summary>
    public static GaussianRational PageCentreEndB(int m, BigRational epsilon, BigRational gamma)
        => new GaussianRational(new BigRational(m), BigRational.Zero) * PageCentreEndA(m, epsilon, gamma);

    // ---- helpers ---------------------------------------------------------------------------------

    private static void Guard(int n)
    {
        if (n < 5 || n % 2 == 0)
            throw new ArgumentOutOfRangeException(nameof(n), n, "N must be odd and at least 5 (N = 2m+1, m ≥ 2).");
        if (n > MaxSites)
            throw new ArgumentOutOfRangeException(nameof(n), n, $"N exceeds the cost guard MaxSites = {MaxSites}.");
    }

    private static GaussianRational[,] Zeros(int n)
    {
        var x = new GaussianRational[n, n];
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++)
                x[i, j] = GaussianRational.Zero;
        return x;
    }

    private static GaussianRational[,] Subtract(GaussianRational[,] x, GaussianRational[,] y)
    {
        int n = x.GetLength(0);
        var d = new GaussianRational[n, n];
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++)
                d[i, j] = x[i, j] - y[i, j];
        return d;
    }

    private static GaussianRational PowMinusTwoI(int m)
    {
        var value = GaussianRational.One;
        var minusTwoI = new GaussianRational(BigRational.Zero, new BigRational(-2));
        for (int k = 0; k < m; k++) value *= minusTwoI;
        return value;
    }

    private static BigRational Pow2(int k)
    {
        var value = BigRational.One;
        for (int i = 0; i < k; i++) value *= new BigRational(2);
        return value;
    }

    private static BigRational Factorial(int k)
    {
        BigInteger f = BigInteger.One;
        for (int i = 2; i <= k; i++) f *= i;
        return new BigRational(f);
    }

    // ℚ has an order but BigRational does not expose one; these three are local rather than a widening
    // of that type, so the promotion of ℚ(i) into Core stays the only arithmetic change this makes.
    private static BigRational Abs(BigRational q) => q.Sign < 0 ? -q : q;
    private static bool Greater(BigRational a, BigRational b) => (a - b).Sign > 0;
    private static BigRational Max(BigRational a, BigRational b) => Greater(a, b) ? a : b;

    // ---- the live tree -----------------------------------------------------------------------------

    public string DisplayName => "Missing-phase onsets, recomputed in exact ℚ(i)";

    /// <summary>One reading per chain length, computed once per instance: the tree asks for them twice,
    /// in the summary and in the children, and the recurrence is the same run both times.</summary>
    private IReadOnlyList<Reading> Readings =>
        _readings ??= PageSites.Select(n => Read(n, Epsilon, Gamma)).ToList();

    private IReadOnlyList<Reading>? _readings;

    public string Summary
    {
        get
        {
            var rows = Readings;
            var agree = rows.Count(r => Agrees(r, Epsilon, Gamma));
            var lag = rows.All(r => r.DTwoOrder - r.DOutOrder == (r.M % 2 == 0 ? 0 : 1))
                    ? "d₂ lags d_out by exactly one order at odd m and starts with it at even m"
                    : "the d₂/d_out order gap is NOT the alternating one";
            var alternates = rows.All(r => r.LeakageCoefficient.Sign == (r.M % 2 == 0 ? 1 : -1))
                    ? "the leakage difference alternates in sign with m"
                    : "the leakage difference does NOT alternate in sign with m";
            return $"The three readings of a broken mirror at ε = {Epsilon}, γ = {Gamma}: " +
                   $"{agree} of {rows.Count} chain lengths reproduce every closed form of §8 exactly. " +
                   $"Onsets (N: d_out, d₂, leakage) = " +
                   string.Join("; ", rows.Select(r => $"{r.Sites}: t^{r.DOutOrder}, t^{r.DTwoOrder}, t^{r.LeakageOrder}"))
                   + $". {lag}, and {alternates}.";
        }
    }

    public IEnumerable<IInspectable> Children
    {
        get
        {
            foreach (var r in Readings)
            {
                int n = r.Sites;
                var sign = r.LeakageCoefficient.Sign < 0 ? "−" : "+";
                var dOut = r.DOutOrder < 0
                    ? "no onset within the orders searched"
                    : r.DOutCoefficientSquared.HasValue
                        ? $"opens at t^{r.DOutOrder}, coefficient² = {r.DOutCoefficientSquared.Value}"
                        : $"opens at t^{r.DOutOrder}, coefficient not representable here (the difference is " +
                          "not the rank-two case of §8.3)";
                yield return new InspectableNode(
                    $"N = {n} (m = {r.M})",
                    $"d_out {dOut}, " +
                    $"d₂ at t^{r.DTwoOrder} (coefficient = {r.DTwoCoefficient}), " +
                    $"leakage at {sign}t^{r.LeakageOrder} (coefficient = {r.LeakageCoefficient}); " +
                    $"centre/end entries δA = {r.CentreEndA}, δB = {r.CentreEndB}; " +
                    $"page closed forms {(Agrees(r, Epsilon, Gamma) ? "met exactly" : "NOT met")}.",
                    provenance: NodeProvenance.Live);
            }

            int mismatches = ReductionMismatchCount(ReductionSites, Epsilon, Gamma, ReductionOrders);
            yield return new InspectableNode(
                $"the §8.1 reduction against the full Lindbladian (N = {ReductionSites})",
                mismatches == 0
                    ? $"½·W[[A,B],[B,A]]W† equals the {1 << ReductionSites}×{1 << ReductionSites} state under " +
                      $"L(ρ) = −i[H,ρ] + γ(Z_c ρ Z_c − ρ) in every cell of every order through t^{ReductionOrders}. " +
                      "The full arm is seeded from the state, not from the blocks, so order 0 is a check of " +
                      "the embedding and not a restatement of it."
                    : $"{mismatches} cells disagree with the physical Lindbladian; the reduction of §8.1 " +
                      "does not hold as implemented.",
                provenance: NodeProvenance.Live);
        }
    }

    /// <summary>Where the reduction is checked against the physical generator. One N settles it: §8.1 is
    /// a claim about the generator, not about the chain length.</summary>
    public const int ReductionSites = 5;
    private const int ReductionOrders = 6;

    /// <summary>Whether one reading reproduces every closed form of §8 exactly. Exact equality: these
    /// are two routes to the same rational, so a difference is a finding, not a tolerance.</summary>
    public static bool Agrees(Reading r, BigRational epsilon, BigRational gamma)
        => r.DOutOrder == r.M + 1
        && r.DOutCoefficientSquared.HasValue
        && r.DOutCoefficientSquared.Value == PageDOutSquared(r.M, epsilon, gamma)
        && r.DTwoOrder == (r.M % 2 == 0 ? r.M + 1 : r.M + 2)
        && r.DTwoCoefficient == PageDTwo(r.M, epsilon, gamma)
        && r.LeakageOrder == 2 * r.M + 1
        && r.LeakageCoefficient == PageLeakage(r.M, epsilon, gamma)
        && r.CentreEndA == PageCentreEndA(r.M, epsilon, gamma)
        && r.CentreEndB == PageCentreEndB(r.M, epsilon, gamma);

    public InspectablePayload Payload => InspectablePayload.Empty;
}
