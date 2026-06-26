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
/// returns the identity. A LASSO (a path from a common base point out around a branch point and back)
/// expresses that branch point's transposition in the base point's labelling, so the transpositions of
/// all branch points live in ONE labelling and can be assembled into the monodromy = Galois group. The
/// tracker follows each root by nearest-neighbour continuity along a fine discretisation, robust to the
/// principal-branch jumps of the per-root closed forms (it samples the continuous root SET).</summary>
public static class Monodromy
{
    /// <summary>The permutation induced by following the roots along an arbitrary CLOSED path (path[0]
    /// ≈ path[^1]). Returns perm where perm[i] = the final position (in the start ordering) of the root
    /// that started at index i.</summary>
    public static int[] PermutationAlongPath(Func<Complex, Complex[]> rootsAt, IReadOnlyList<Complex> path)
    {
        var r0 = rootsAt(path[0]);
        int d = r0.Length;
        var current = (Complex[])r0.Clone();
        var labels = Enumerable.Range(0, d).ToArray();      // labels[j] = start index of the root now at current[j]

        for (int k = 1; k < path.Count; k++)
        {
            var next = rootsAt(path[k]);
            int[] assign = NearestBijection(current, next);  // assign[j] = index in current matched to next[j]
            var newLabels = new int[d];
            for (int j = 0; j < d; j++) newLabels[j] = labels[assign[j]];
            current = next;
            labels = newLabels;
        }

        var perm = new int[d];
        for (int j = 0; j < d; j++)
            perm[labels[j]] = NearestIndex(r0, current[j]);
        return perm;
    }

    /// <summary>Follow the d roots along an OPEN path, keeping each strand's start-label by
    /// nearest-neighbour continuity, and return the FULL per-step trajectory: traj[k][p] = the position
    /// of the strand labelled p (its index at path[0]) when q = path[k]. Unlike
    /// <see cref="PermutationAlongPath"/> (closed loop, returns only the net end-permutation), this
    /// exposes every intermediate position, so a strand can be watched flowing toward a coalescence
    /// (an EP / a diabolic point) along a real-q sweep. Same nearest-bijection matching, just emitting
    /// the trajectory instead of discarding it.</summary>
    public static Complex[][] TrajectoryAlongPath(Func<Complex, Complex[]> rootsAt, IReadOnlyList<Complex> path)
    {
        var r0 = rootsAt(path[0]);
        int d = r0.Length;
        var traj = new Complex[path.Count][];
        traj[0] = (Complex[])r0.Clone();                    // strand p starts as r0[p]
        var current = (Complex[])r0.Clone();
        var labels = Enumerable.Range(0, d).ToArray();      // labels[j] = start-label of the root now at current[j]

        for (int k = 1; k < path.Count; k++)
        {
            var next = rootsAt(path[k]);
            int[] assign = NearestBijection(current, next); // assign[j] = index in current matched to next[j]
            var newLabels = new int[d];
            var slot = new Complex[d];
            for (int j = 0; j < d; j++) newLabels[j] = labels[assign[j]];
            for (int j = 0; j < d; j++) slot[newLabels[j]] = next[j];   // place each root in its start-label slot
            traj[k] = slot;
            current = next;
            labels = newLabels;
        }
        return traj;
    }

    /// <summary>The permutation induced by one counter-clockwise circular loop of radius
    /// <paramref name="radius"/> about <paramref name="center"/>, discretised into
    /// <paramref name="steps"/> arcs. Identity ⟺ no braiding enclosed; a transposition ⟺ one √-branch
    /// point; a k-cycle ⟺ a k-fold branch.</summary>
    public static int[] Permutation(Func<Complex, Complex[]> rootsAt, Complex center, double radius, int steps)
    {
        var path = new Complex[steps + 1];
        for (int s = 0; s <= steps; s++)
        {
            double theta = 2.0 * Math.PI * s / steps;
            path[s] = center + radius * new Complex(Math.Cos(theta), Math.Sin(theta));
        }
        return PermutationAlongPath(rootsAt, path);
    }

    /// <summary>A lasso: the closed path baseQ → (straight in) → a full circle of radius
    /// <paramref name="radius"/> around <paramref name="branchPoint"/> → (straight back) → baseQ. Tracking
    /// the roots along it yields that branch point's monodromy expressed in baseQ's labelling, so lassos to
    /// different branch points share one labelling and their transpositions can be assembled.</summary>
    public static Complex[] Lasso(Complex baseQ, Complex branchPoint, double radius, int density = 240)
    {
        var dir = branchPoint - baseQ;
        double dist = dir.Magnitude;
        Complex u = dir / dist;                               // baseQ → branchPoint
        Complex enter = branchPoint - radius * u;             // point on the circle nearest baseQ
        var pts = new List<Complex>();

        int nLine = Math.Max(20, (int)(density * Math.Max(0.0, dist - radius)));
        for (int k = 0; k <= nLine; k++) pts.Add(baseQ + (enter - baseQ) * ((double)k / nLine));

        double th0 = Math.Atan2(enter.Imaginary - branchPoint.Imaginary, enter.Real - branchPoint.Real);
        int nCirc = Math.Max(120, density);
        for (int k = 1; k <= nCirc; k++)
        {
            double th = th0 + 2.0 * Math.PI * k / nCirc;
            pts.Add(branchPoint + radius * new Complex(Math.Cos(th), Math.Sin(th)));
        }

        for (int k = 1; k <= nLine; k++) pts.Add(enter + (baseQ - enter) * ((double)k / nLine));
        return pts.ToArray();
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
