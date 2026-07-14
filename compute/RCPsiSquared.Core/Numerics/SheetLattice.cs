namespace RCPsiSquared.Core.Numerics;

/// <summary>The F127 sheet lattice (PROOF_F127_RESIDUE_COLLAPSE §2), exact integer combinatorics.
/// The cross form 𝔉 partial-fractions into 72 atoms indexed (i, j, ξ, υ, e); each atom carries two
/// simple poles (the pole sign s = ±1 of cos(a_i ± b_j)) and each pole meets a sheet on two sides
/// (the sheet sign τ = ±1), giving 72·2·2 = 288 pole events. Each event's pole condition is a linear
/// form L·(a₁,a₂,a₃,b₁,b₂,b₃) ≡ 0 whose six coefficients are all ±1; canonicalizing by the sign of
/// the first coefficient collapses the 288 events onto exactly 32 sheets, nine events each, one per
/// (i,j) block. This is the exact integer twin of the committed
/// <c>simulations/f127_sheet_lattice.py</c> (<c>enum_events</c>/<c>canon</c>): no floats, no
/// randomness, so the live witness recomputes the lattice's combinatorics from scratch at inspect
/// time. The residue that every sheet's nine events share is the §3 core function T
/// (<see cref="CrossFormCertificate.EvaluateCoreT"/>).</summary>
public static class SheetLattice
{
    public const int AtomCount = 72;
    public const int EventCount = 288;
    public const int SheetCount = 32;
    public const int EventsPerSheet = 9;

    /// <summary>One pole event: its atom coordinates (i, j, ξ, υ, e), its pole sign s and sheet
    /// sign τ, its raw linear form L and the sign-canonicalized sheet key <c>Canon</c>.</summary>
    public readonly record struct Event(
        int I, int J, string Xi, string Up, int E, int S, int Tau,
        IReadOnlyList<int> L, IReadOnlyList<int> Canon);

    /// <summary>The live-verifiable combinatorics of the lattice.</summary>
    public readonly record struct Report(
        int Atoms, int Events, bool AllCoefficientsPmOne,
        int DistinctSheets, int MinEventsPerSheet, int MaxEventsPerSheet,
        bool EachSheetOneEventPerBlock);

    private static (int Lo, int Hi) Complement(int i)
    {
        int lo = -1, hi = -1;
        for (int t = 0; t < 3; t++)
            if (t != i) { if (lo < 0) lo = t; else hi = t; }
        return (lo, hi);   // lo < hi
    }

    /// <summary>Sign so the first coefficient is +1 (0 canonicalizes to itself; only reached if a
    /// coefficient is off the ±1 lattice, which <see cref="Report.AllCoefficientsPmOne"/> flags).</summary>
    private static int[] Canonicalize(int[] l)
    {
        int sign = l[0] > 0 ? 1 : l[0] < 0 ? -1 : 1;
        var c = new int[6];
        for (int k = 0; k < 6; k++) c[k] = sign * l[k];
        return c;
    }

    /// <summary>The 288 pole events, built exactly. ξ (psum: +1,+1; pdif: −1 on the smaller, +1 on
    /// the larger complement index) on the two non-i a-slots; e·υ on the two non-j b-slots; −τ on
    /// slot i; −τ·s on slot 3+j.</summary>
    public static List<Event> Enumerate()
    {
        var events = new List<Event>(EventCount);
        string[] kinds = { "psum", "pdif" };
        for (int i = 0; i < 3; i++)
        {
            var (ja, la) = Complement(i);
            for (int j = 0; j < 3; j++)
            {
                var (jb, lb) = Complement(j);
                foreach (var xi in kinds)
                    foreach (var up in kinds)
                        foreach (int e in new[] { 1, -1 })
                        {
                            var vec = new int[6];
                            if (xi == "psum") { vec[ja] += 1; vec[la] += 1; }
                            else { vec[ja] += -1; vec[la] += 1; }
                            if (up == "psum") { vec[3 + jb] += e; vec[3 + lb] += e; }
                            else { vec[3 + jb] += -e; vec[3 + lb] += e; }
                            foreach (int s in new[] { 1, -1 })
                                foreach (int tau in new[] { 1, -1 })
                                {
                                    var l = (int[])vec.Clone();
                                    l[i] -= tau;
                                    l[3 + j] -= tau * s;
                                    events.Add(new Event(i, j, xi, up, e, s, tau, l, Canonicalize(l)));
                                }
                        }
            }
        }
        return events;
    }

    private static int CanonKey(IReadOnlyList<int> canon)
    {
        int key = 0;
        for (int k = 0; k < 6; k++) key |= (canon[k] > 0 ? 1 : 0) << k;
        return key;
    }

    public static Report Analyze()
    {
        var events = Enumerate();
        int atoms = events.Select(ev => (ev.I, ev.J, ev.Xi, ev.Up, ev.E)).Distinct().Count();
        bool allPm = events.All(ev => ev.L.All(c => c == 1 || c == -1));
        var sheets = events.GroupBy(ev => CanonKey(ev.Canon)).ToList();
        int min = sheets.Min(g => g.Count());
        int max = sheets.Max(g => g.Count());
        bool oneEach = sheets.All(g => g.Select(ev => (ev.I, ev.J)).Distinct().Count() == g.Count());
        return new Report(atoms, events.Count, allPm, sheets.Count, min, max, oneEach);
    }
}
