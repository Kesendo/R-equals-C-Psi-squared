namespace MirrorWorld;

// What a number-conserving Hamiltonian does to the bare sectors: the superposition leaves the grid.
// Dephasing alone pins every coherence at an exact rung -2*gamma*k (bare counts 2^N*C(N,k)). Turn H
// on and it mixes rungs within each (p,q) block by even steps, so eigenmodes get fractional <n_XY>
// and many land OFF the integer grid. What stays on-grid is T1-derived (DEGENERACY_PALINDROME,
// verified N=2..7): the edges stay N+1 (k=0 = identity + N magnetization projectors, the kernel; its
// Pi-image k=N, the drain), and even N spikes the center k=N/2 (the Pi axis sits on the grid there).
// The on-grid folds below are ADOPTED (the proven d_total numbers), not recomputed.
public static class Redistribution
{
    // d_total(k) on-grid, H on (DEGENERACY_PALINDROME Result 3). null = outside the adopted table.
    public static int[]? OnGrid(int n) => n switch
    {
        2 => new[] { 3, 10, 3 },
        3 => new[] { 4, 14, 14, 4 },
        4 => new[] { 5, 20, 152, 20, 5 },
        _ => null,
    };

    // bare per-k = 2^N * C(N,k), the empty-world sector multiplicities.
    public static int[] Bare(int n)
    {
        var f = new int[n + 1];
        for (int k = 0; k <= n; k++) f[k] = (1 << n) * (int)Block.Binomial(n, k);
        return f;
    }
}
