using System;
using System.Collections.Generic;
using System.Linq;
using System.Numerics;

namespace RCPsiSquared.Core.Numerics;

/// <summary>q-parametric monodromy: as the parameter q traces a closed loop in the complex plane,
/// the d roots of an algebraic function λ(q) are permuted. That permutation group (over all loops) is
/// the monodromy group, which over a function field equals the Galois group. This is where the Galois
/// structure lives spectrally: NOT in the fixed-q geometry, but in how the eigenvalue strands braid as
/// q moves. A loop around a √-branch point (simple discriminant zero, a defective EP) swaps two roots;
/// a loop around a transversal crossing of analytic sheets (double discriminant zero, a diabolic point)
/// returns the identity. The tracker follows each root by nearest-neighbour continuity along a fine
/// discretisation of the loop, robust to the principal-branch jumps of the per-root closed forms (it
/// samples the root SET, which is continuous, not any single branch).</summary>
public static class Monodromy
{
    /// <summary>The permutation induced by one counter-clockwise loop of radius <paramref name="radius"/>
    /// about <paramref name="center"/>, discretised into <paramref name="steps"/> arcs. Returns perm where
    /// perm[i] = the final position (in the start ordering) of the root that started at index i. Identity
    /// ⟺ no braiding enclosed; a transposition ⟺ one √-branch point; a k-cycle ⟺ a k-fold branch.</summary>
    public static int[] Permutation(Func<Complex, Complex[]> rootsAt, Complex center, double radius, int steps)
    {
        Complex Q(double theta) => center + radius * new Complex(Math.Cos(theta), Math.Sin(theta));

        var r0 = rootsAt(Q(0));
        int d = r0.Length;
        var current = (Complex[])r0.Clone();
        var labels = Enumerable.Range(0, d).ToArray();      // labels[j] = start index of the root now at current[j]

        for (int s = 1; s <= steps; s++)
        {
            var next = rootsAt(Q(2.0 * Math.PI * s / steps));
            int[] assign = NearestBijection(current, next);  // assign[j] = index in current matched to next[j]
            var newLabels = new int[d];
            for (int j = 0; j < d; j++) newLabels[j] = labels[assign[j]];
            current = next;
            labels = newLabels;
        }

        // current ≈ r0 as a set; read off start-index → end-position-in-r0-ordering.
        var perm = new int[d];
        for (int j = 0; j < d; j++)
            perm[labels[j]] = NearestIndex(r0, current[j]);
        return perm;
    }

    // greedy global-nearest bijection: match each next root to a distinct current root, closest pairs first.
    private static int[] NearestBijection(Complex[] current, Complex[] next)
    {
        int d = current.Length;
        var pairs = new List<(double dist, int j, int m)>(d * d);
        for (int j = 0; j < d; j++)
            for (int m = 0; m < d; m++)
                pairs.Add(((next[j] - current[m]).Magnitude, j, m));
        pairs.Sort((a, b) => a.dist.CompareTo(b.dist));

        var assign = new int[d];
        Array.Fill(assign, -1);
        var usedCurrent = new bool[d];
        int done = 0;
        foreach (var (_, j, m) in pairs)
        {
            if (assign[j] != -1 || usedCurrent[m]) continue;
            assign[j] = m;
            usedCurrent[m] = true;
            if (++done == d) break;
        }
        return assign;
    }

    private static int NearestIndex(Complex[] set, Complex z)
    {
        int best = 0;
        double bd = double.PositiveInfinity;
        for (int m = 0; m < set.Length; m++)
        {
            double dd = (set[m] - z).Magnitude;
            if (dd < bd) { bd = dd; best = m; }
        }
        return best;
    }
}
