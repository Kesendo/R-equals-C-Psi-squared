namespace MirrorWorld;

// A joint-popcount block (p, q) = (popcount of bra, popcount of ket): the second grading. A
// number-conserving Hamiltonian (XY/Heisenberg/XXZ) conserves the two popcounts separately
// (U(1)xU(1)), so the Liouvillian is block-diagonal in (p, q): (N+1)^2 blocks, block (p,q) of size
// C(N,p)*C(N,q), summing to 4^N. Inside a block the disagreement count k runs over
// {|p-q|, |p-q|+2, ...} (fixed parity p+q); the Hamiltonian mixes those rungs (the superposition
// leaving the integer grid), the dephasing keeps them apart. T1: JointPopcountSectors,
// LiouvillianBlockSpectrum (per-block eig = full spectrum, bit-exact N=3,4,5).
public sealed class Block : GameObject
{
    public int N { get; }
    public int P { get; }   // popcount of the bra
    public int Q { get; }   // popcount of the ket

    public Block(World world, int n, int p, int q) : base(world)
    {
        N = n;
        P = p;
        Q = q;
    }

    public long Size => Binomial(N, P) * Binomial(N, Q);
    public bool Diagonal => P == Q;            // (p,p) carry even k incl k=0, the populations
    public int MinK => Math.Abs(P - Q);        // the lowest disagreement rung living in this block

    public override IReadOnlyList<string> Own => new[] { "p", "q", "size", "minK" };

    public static long Binomial(int n, int k)
    {
        if (k < 0 || k > n) return 0;
        long r = 1;
        for (int i = 1; i <= k; i++) r = r * (n - k + i) / i;
        return r;
    }
}
