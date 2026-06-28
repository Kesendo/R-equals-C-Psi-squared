using System.Globalization;
using MirrorWorld;

CultureInfo.CurrentCulture = CultureInfo.InvariantCulture;   // dots, not commas

const double gamma = 0.5;
var world = new World();
Console.WriteLine($"empty world (Z-dephasing, no Hamiltonian), gamma={gamma}");
Console.WriteLine($"World.Own (left, the frame): {string.Join(",", world.Own)}");
Console.WriteLine();

foreach (int N in new[] { 2, 3, 4 })
{
    var modes = EnumeratePauli(N).Select(l => new PauliMode(world, l, gamma)).ToList();
    Console.WriteLine($"==== N={N}  ({modes.Count} Pauli-string modes = {1 << N}*sum C(N,k); the superposition basis) ====");

    for (int k = 0; k <= N; k++)
    {
        int count = modes.Count(m => m.K == k);
        string tag = k == 0 ? "  [populations / kernel]"
                   : k == N ? "  [drain]"
                   : k == N - k ? "  [SELF-MIRROR center]" : "";
        Console.WriteLine($"  k={k}: {count,3} modes   rate {-2.0 * gamma * k,5:0.0}   mirror k={N - k}{tag}");
    }

    Console.WriteLine(N % 2 == 0
        ? $"  fold k<->N-k about k=N/2={N / 2}: even N, the self-mirror sector EXISTS (k={N / 2})"
        : $"  fold k<->N-k about k=N/2={(N / 2.0):0.0}: odd N, axis BETWEEN rungs, NO self-mirror sector");

    int kc = (N % 2 == 0) ? N / 2 : (N - 1) / 2;   // central (even) or inner (odd)
    var cells = modes.Where(m => m.K == kc)
                     .GroupBy(m => m.Klein)
                     .OrderBy(g => g.Key.Y).ThenBy(g => g.Key.Z)
                     .Select(g => $"(nY%2={g.Key.Y},nZ%2={g.Key.Z}):{g.Count()}");
    Console.WriteLine($"  superposition / Klein quartering of sector k={kc}: {string.Join("  ", cells)}   [truly = (0,0)]");

    // Grading B (T1): the joint-popcount blocks a number-conserving Hamiltonian respects.
    Console.WriteLine($"  Grading B: {(N + 1) * (N + 1)} joint-popcount blocks (p,q), size C(N,p)*C(N,q):");
    long total = 0;
    for (int p = 0; p <= N; p++)
    {
        var sizes = new List<string>();
        for (int q = 0; q <= N; q++)
        {
            var b = new Block(world, N, p, q);
            total += b.Size;
            sizes.Add($"{b.Size,3}");
        }
        Console.WriteLine($"    p={p}: {string.Join(" ", sizes)}");
    }
    Console.WriteLine($"    sum = {total} = 4^{N}; the diagonal (p,p) carries k=0, the {N + 1} populations (the kernel)");
    Console.WriteLine();
}

// The special cases as T1 closed forms , adopted, not re-derived (all already Tier-1 in the repo).
Console.WriteLine("==== special cases (T1 closed forms, adopted) ====");
for (int n = 2; n <= 6; n++)
    Console.WriteLine($"  N={n}: V-effect 1+cos(pi/N) = {1.0 + Math.Cos(Math.PI / n):0.0000}   ceiling g2(K_N)=4/N = {4.0 / n:0.0000}");
Console.WriteLine($"  N=3: minimal EP, Q_EP = 2/g_eff (the whole rate ladder is HD=1 and HD=3)");
Console.WriteLine($"  N=4: self-fold diabolic q_EP = sqrt((-1+sqrt13)/6) = {Math.Sqrt((-1.0 + Math.Sqrt(13.0)) / 6.0):0.0000}  (real q, N=4 only); (2,2) ceiling K_4 = 2-2/sqrt3 = {2.0 - 2.0 / Math.Sqrt(3.0):0.0000}");
Console.WriteLine($"  N=5: diabolics leave the real axis into complex-conjugate pairs (only loud defective EPs remain real)");

// all 4^N Pauli strings, site 0 the least-significant base-4 digit
static IEnumerable<char[]> EnumeratePauli(int N)
{
    char[] alphabet = { 'I', 'X', 'Y', 'Z' };
    int total = 1 << (2 * N);
    for (int idx = 0; idx < total; idx++)
    {
        var l = new char[N];
        int x = idx;
        for (int s = 0; s < N; s++) { l[s] = alphabet[x & 3]; x >>= 2; }
        yield return l;
    }
}
