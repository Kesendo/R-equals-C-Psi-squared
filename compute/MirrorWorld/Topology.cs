namespace MirrorWorld;

// The geometry (the doc's "a little geometry: who can reach whom"). The handshake rides these bonds. Every
// one is excitation-conserving, so they all share the same joint-popcount block structure (F63): the geometry
// shapes how novelty is born and spreads, not which cells are forbidden. Standalone, re-derived from the atom.
public static class Topology
{
    // a line: 0-1-2-...-(N-1).
    public static (int a, int b)[] Chain(int n)
    {
        var bonds = new List<(int, int)>();
        for (int i = 0; i + 1 < n; i++) bonds.Add((i, i + 1));
        return bonds.ToArray();
    }

    // the line closed into a loop (for N>2; a 2-ring is just the one bond).
    public static (int a, int b)[] Ring(int n)
    {
        var bonds = new List<(int, int)>(Chain(n));
        if (n > 2) bonds.Add((n - 1, 0));
        return bonds.ToArray();
    }

    // a hub (site 0) joined to every spoke.
    public static (int a, int b)[] Star(int n)
    {
        var bonds = new List<(int, int)>();
        for (int i = 1; i < n; i++) bonds.Add((0, i));
        return bonds.ToArray();
    }

    // all-to-all.
    public static (int a, int b)[] Complete(int n)
    {
        var bonds = new List<(int, int)>();
        for (int a = 0; a < n; a++)
            for (int b = a + 1; b < n; b++) bonds.Add((a, b));
        return bonds.ToArray();
    }

    public static (int a, int b)[] Named(string name, int n) => name switch
    {
        "ring" => Ring(n),
        "star" => Star(n),
        "complete" => Complete(n),
        _ => Chain(n),
    };
}
