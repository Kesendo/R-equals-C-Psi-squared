using System.Collections.Generic;
using System.Numerics;
using Xunit;

namespace RCPsiSquared.Core.Tests.TestHelpers;

/// <summary>Greedy nearest-neighbour multiset-equality assertion for Complex eigenvalue
/// arrays. Necessary because XY-chain spectra are highly degenerate and sort-then-zip
/// mis-pairs conjugate clusters across independent eigensolvers.</summary>
public static class MultisetAssert
{
    public static void NearestNeighbourEqual(
        IReadOnlyList<Complex> actual,
        IReadOnlyList<Complex> expected,
        double tolerance,
        string? context = null)
    {
        Assert.Equal(actual.Count, expected.Count);
        var taken = new bool[expected.Count];
        for (int i = 0; i < actual.Count; i++)
        {
            var x = actual[i];
            int bestIdx = -1;
            double bestDist = double.MaxValue;
            for (int j = 0; j < expected.Count; j++)
            {
                if (taken[j]) continue;
                double d = (x - expected[j]).Magnitude;
                if (d < bestDist) { bestDist = d; bestIdx = j; }
            }
            string ctxPrefix = context is null ? "" : $"[{context}] ";
            Assert.True(bestIdx >= 0 && bestDist < tolerance,
                $"{ctxPrefix}No multiset match within {tolerance}: actual[{i}]={x}; " +
                (bestIdx >= 0 ? $"closest expected[{bestIdx}]={expected[bestIdx]} at distance {bestDist:E3}" : "no candidate"));
            taken[bestIdx] = true;
        }
    }
}
