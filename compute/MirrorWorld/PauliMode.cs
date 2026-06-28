namespace MirrorWorld;

// A Pauli string: the symmetry-adapted basis of a disagreement-count sector, the superposition.
// Each Pauli string is a fixed superposition of the bare pairs (Pair) at the same disagreement
// count, so the two bases have equal size; the Pauli strings are the ones that also diagonalize the
// mirror group. Its XY-weight k = n_X + n_Y IS the disagreement count (the only thing the watching
// reads), rate Re lambda = -2*gamma*k. The diagonal mirrors {(-1)^n_Y, (-1)^n_Z} quarter a sector
// into four equal Klein cells by (n_Y, n_Z) parity.
public sealed class PauliMode : GameObject
{
    public char[] Letters { get; }   // I / X / Y / Z per site
    public double Gamma { get; }

    public PauliMode(World world, char[] letters, double gamma) : base(world)
    {
        Letters = letters;
        Gamma = gamma;
    }

    public int Nx => Letters.Count(c => c == 'X');
    public int Ny => Letters.Count(c => c == 'Y');
    public int Nz => Letters.Count(c => c == 'Z');

    public int K => Nx + Ny;                          // XY-weight = disagreement count
    public double Rate => -2.0 * Gamma * K;
    public (int Y, int Z) Klein => (Ny % 2, Nz % 2);  // the quartering

    public override IReadOnlyList<string> Own => new[] { "letters", "k", "rate", "klein" };

    // All 4^N Pauli strings (site 0 = least-significant base-4 digit): the symmetry-adapted basis.
    // Shared by the sim and the tests.
    public static IEnumerable<PauliMode> Enumerate(World world, int n, double gamma)
    {
        char[] alphabet = { 'I', 'X', 'Y', 'Z' };
        int total = 1 << (2 * n);
        for (int idx = 0; idx < total; idx++)
        {
            var l = new char[n];
            int x = idx;
            for (int s = 0; s < n; s++) { l[s] = alphabet[x & 3]; x >>= 2; }
            yield return new PauliMode(world, l, gamma);
        }
    }
}
