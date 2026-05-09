using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>F4 closed form (Tier 1, Clebsch-Gordan decomposition; ANALYTICAL_FORMULAS line 150):
///
/// <code>
///   Stat(N) = Σ_J  m(J, N) · (2J + 1)²
///
///   m(J, N) = C(N, N/2 − J) − C(N, N/2 − J − 1)
///
///   where J ranges over allowed total spins for N spin-1/2 particles:
///     N even: J ∈ {0, 1, 2, ..., N/2}
///     N odd:  J ∈ {1/2, 3/2, ..., N/2}
/// </code>
///
/// <para>F4 counts the number of stationary modes (kernel of the purely-Hamiltonian
/// Liouvillian L = −i[H, ·] when Σγ = 0) for an N-spin Heisenberg chain. Each
/// total-spin-J multiplet of multiplicity m(J, N) contributes (2J+1)² operators
/// that commute with the Heisenberg H within their J-subspace.</para>
///
/// <para><b>Verified table:</b></para>
/// <list type="bullet">
///   <item>N=2: Stat = 10 (J=0 contributes 1, J=1 contributes 9)</item>
///   <item>N=3: Stat = 24 (J=1/2 contributes 8, J=3/2 contributes 16)</item>
///   <item>N=4: Stat = 54 (J=0 contributes 2, J=1 contributes 27, J=2 contributes 25)</item>
///   <item>N=5: Stat = 120 (J=1/2: 20, J=3/2: 64, J=5/2: 36)</item>
/// </list>
///
/// <para>Pi2-Foundation anchor:</para>
/// <list type="bullet">
///   <item><b>IrrepDimensionCoefficient = 2 = a_0</b>: in 2J+1. Live from
///         <see cref="Pi2DyadicLadderClaim.Term"/>(0). Same anchor as F1's
///         TwoFactor, F50's DegeneracyFactor, F66's UpperPoleCoefficient,
///         F43's XorRateCoefficient. The (2J+1)² in F4 is structurally
///         a_0²·J² + 2·a_0·J + 1 = "twice anchor", appearing squared.</item>
/// </list>
///
/// <para>The Clebsch-Gordan multiplicity m(J, N) is a combinatorial difference
/// of binomial coefficients (no Pi2 anchor); the structural identity Σ
/// m(J)·(2J+1) = 2^N (sum of irrep dimensions = total Hilbert dimension)
/// makes the closed form rigorous via Schur-Weyl duality.</para>
///
/// <para>Tier1Derived: F4 is Tier 1 via Clebsch-Gordan SU(2)-irrep decomposition
/// (Schur-Weyl); valid for Heisenberg Hamiltonian, Σγ = 0, all N. Pi2-Foundation
/// anchoring is composition through Pi2DyadicLadder.</para>
///
/// <para>Anchors: <c>docs/ANALYTICAL_FORMULAS.md</c> F4 (line 150) +
/// <c>experiments/CAVITY_MODES_FORMULA.md</c> +
/// <c>compute/RCPsiSquared.Core/Symmetry/Pi2DyadicLadderClaim.cs</c>.</para></summary>
public sealed class F4StationaryModeCountPi2Inheritance : Claim
{
    private readonly Pi2DyadicLadderClaim _ladder;

    /// <summary>The "2" coefficient in 2J+1 (irrep dimension formula). Live from
    /// Pi2DyadicLadder a_0. Same anchor as F1/F50/F66/F43.</summary>
    public double IrrepDimensionCoefficient => _ladder.Term(0);

    /// <summary>Clebsch-Gordan multiplicity of total spin J in N spin-1/2 particles:
    /// <c>m(J, N) = C(N, N/2 − J) − C(N, N/2 − J − 1)</c>. The 2J = twoJ
    /// parameter is integer (twoJ takes values N, N−2, ..., 0 or 1 depending on
    /// parity of N).</summary>
    public int SpinMultiplicity(int N, int twoJ)
    {
        if (N < 1) throw new ArgumentOutOfRangeException(nameof(N), N, "F4 requires N ≥ 1.");
        if (twoJ < 0) throw new ArgumentOutOfRangeException(nameof(twoJ), twoJ, "twoJ must be ≥ 0.");
        if (twoJ > N) throw new ArgumentOutOfRangeException(nameof(twoJ), twoJ, $"twoJ must be ≤ N; got twoJ={twoJ}, N={N}.");
        if ((N - twoJ) % 2 != 0)
            throw new ArgumentException($"N − twoJ must be even (N={N}, twoJ={twoJ}); twoJ has wrong parity for N spin-1/2.");

        int kHigh = (N - twoJ) / 2;
        int kLow = kHigh - 1;
        return Binomial(N, kHigh) - Binomial(N, kLow);
    }

    /// <summary>Irrep dimension <c>(2J+1)</c> from twoJ = 2J. The "2" is
    /// IrrepDimensionCoefficient = a_0; the +1 is from the discreteness of integer
    /// dimensions.</summary>
    public int IrrepDimension(int twoJ)
    {
        if (twoJ < 0) throw new ArgumentOutOfRangeException(nameof(twoJ), twoJ, "twoJ must be ≥ 0.");
        return twoJ + 1;
    }

    /// <summary>F4's stationary mode count <c>Stat(N) = Σ_J m(J, N) · (2J + 1)²</c>.
    /// Sum over allowed twoJ values: N, N−2, N−4, ..., 0 (even N) or 1 (odd N).</summary>
    public int StationaryModeCount(int N)
    {
        if (N < 1) throw new ArgumentOutOfRangeException(nameof(N), N, "F4 requires N ≥ 1.");
        int total = 0;
        // twoJ ranges from N downward by 2 to 0 (even N) or 1 (odd N).
        for (int twoJ = N; twoJ >= 0; twoJ -= 2)
        {
            int m = SpinMultiplicity(N, twoJ);
            int dim = IrrepDimension(twoJ);
            total += m * dim * dim;
        }
        return total;
    }

    /// <summary>Drift check: Σ_J m(J, N) · (2J+1) = 2^N (Schur-Weyl identity, the
    /// total Hilbert dimension equals the sum of irrep dimensions weighted by
    /// multiplicities).</summary>
    public bool SchurWeylDimensionIdentityHolds(int N)
    {
        if (N < 1) throw new ArgumentOutOfRangeException(nameof(N), N, "F4 requires N ≥ 1.");
        int total = 0;
        for (int twoJ = N; twoJ >= 0; twoJ -= 2)
        {
            int m = SpinMultiplicity(N, twoJ);
            int dim = IrrepDimension(twoJ);
            total += m * dim;
        }
        return total == (1 << N);
    }

    /// <summary>List of (twoJ, multiplicity, irrepDimension, contribution) for each
    /// allowed J, useful for inspecting the F4 sum decomposition.</summary>
    public IReadOnlyList<(int TwoJ, int Multiplicity, int IrrepDim, int Contribution)> SpinDecomposition(int N)
    {
        if (N < 1) throw new ArgumentOutOfRangeException(nameof(N), N, "F4 requires N ≥ 1.");
        var result = new List<(int, int, int, int)>();
        for (int twoJ = N; twoJ >= 0; twoJ -= 2)
        {
            int m = SpinMultiplicity(N, twoJ);
            int dim = IrrepDimension(twoJ);
            result.Add((twoJ, m, dim, m * dim * dim));
        }
        return result;
    }

    private static int Binomial(int N, int k)
    {
        if (k < 0 || k > N) return 0;
        if (k == 0 || k == N) return 1;
        if (k > N - k) k = N - k;
        long c = 1;
        for (int i = 0; i < k; i++)
        {
            c = c * (N - i) / (i + 1);
        }
        return (int)c;
    }

    public F4StationaryModeCountPi2Inheritance(Pi2DyadicLadderClaim ladder)
        : base("F4 stationary mode count Stat(N) = Σ_J m(J)·(2J+1)² (Clebsch-Gordan SU(2)-irrep decomposition); 2 = a_0 in 2J+1; verified N=2..5: 10, 24, 54, 120",
               Tier.Tier1Derived,
               "docs/ANALYTICAL_FORMULAS.md F4 + " +
               "experiments/CAVITY_MODES_FORMULA.md + " +
               "compute/RCPsiSquared.Core/Symmetry/Pi2DyadicLadderClaim.cs")
    {
        _ladder = ladder ?? throw new ArgumentNullException(nameof(ladder));
    }

    public override string DisplayName =>
        "F4 stationary mode count as Pi2-Foundation a_0 (in 2J+1) inheritance";

    public override string Summary =>
        $"Stat(N) = Σ_J m(J)·(2J+1)²; Schur-Weyl Clebsch-Gordan; 2 = a_0 (= {IrrepDimensionCoefficient}); Heisenberg H, Σγ=0; N=2:10, N=3:24, N=4:54, N=5:120 ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("F4 closed form",
                summary: "Stat(N) = Σ_J m(J,N)·(2J+1)² where m(J,N) = C(N, N/2−J) − C(N, N/2−J−1); kernel dimension of L_H at Σγ=0 for Heisenberg chain");
            yield return InspectableNode.RealScalar("IrrepDimensionCoefficient (= a_0 = 2)", IrrepDimensionCoefficient);
            yield return new InspectableNode("Schur-Weyl identity",
                summary: "Σ_J m(J,N)·(2J+1) = 2^N (drift check via SchurWeylDimensionIdentityHolds); the total Hilbert dimension equals the sum of irrep dimensions weighted by multiplicities");
            yield return new InspectableNode("verified table",
                summary: $"Stat(2) = {StationaryModeCount(2)}, Stat(3) = {StationaryModeCount(3)}, Stat(4) = {StationaryModeCount(4)}, Stat(5) = {StationaryModeCount(5)}, Stat(6) = {StationaryModeCount(6)}");
            yield return new InspectableNode("N=4 spin decomposition",
                summary: string.Join("; ", SpinDecomposition(4).Select(t => $"twoJ={t.TwoJ}: m={t.Multiplicity}, (2J+1)²={t.IrrepDim*t.IrrepDim}, c={t.Contribution}")));
        }
    }
}
