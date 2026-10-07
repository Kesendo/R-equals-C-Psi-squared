namespace MirrorWorld;

// What a number-conserving Hamiltonian does to the bare sectors: the superposition leaves the grid.
// Dephasing alone pins every coherence at an exact rung -2*gamma*k (bare counts 2^N*C(N,k)). Turn H
// on and it mixes rungs within each (p,q) block by even steps, so eigenmodes get fractional <n_XY>
// and many land OFF the integer grid. What stays on-grid is T1-derived (DEGENERACY_PALINDROME,
// verified N=2..7): the edges stay N+1 (k=0 = identity + N magnetization projectors, the kernel; its
// Pi-image k=N, the drain), and even N spikes the center k=N/2 (the Pi axis sits on the grid there).
// The on-grid folds below are ADOPTED (the verified d_total numbers), not recomputed. They hold at
// generic coupling: N=2's 10 holds for gamma <= 2J (Pauli convention), and the N=3 numbers rise
// (to 4, 16, 16, 4) at F50's exceptional couplings, where real modes of the diagonal blocks join
// the rungs (PROOF_WEIGHT1_DEGENERACY, the count at exceptional couplings).
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
        // The largest N=16 entry is 843,448,320; at N=17 the central count already
        // exceeds Int32.MaxValue. Keep the adopted int-valued read within its domain.
        if (n < 0 || n > 16)
            throw new ArgumentOutOfRangeException(nameof(n), n, "bare multiplicities fit Int32 only for 0 <= N <= 16");
        var f = new int[n + 1];
        for (int k = 0; k <= n; k++) f[k] = checked((1 << n) * (int)Block.Binomial(n, k));
        return f;
    }
}
