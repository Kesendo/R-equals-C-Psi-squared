using System;
using System.Linq;
using System.Numerics;
using RCPsiSquared.Core.Numerics;
using RCPsiSquared.Diagnostics.Foundation;
using Xunit;
using Xunit.Abstractions;

namespace RCPsiSquared.Diagnostics.Tests.Foundation;

public sealed class RouteBA2PositiveAnchorTests
{
    private const int K = 4;
    private const double WAnchor = 5.10083102;
    private const double LambdaAnchor = -4.79196037;
    private readonly ITestOutputHelper _out;

    public RouteBA2PositiveAnchorTests(ITestOutputHelper output) => _out = output;

    [Theory]
    [InlineData(+1)]
    [InlineData(-1)]
    [Trait("Category", "ROUTE_B_A2_CHARACTER")]
    public void PositiveA2Anchor_BothPhysicalQLifts_AreTwinScalarDiabolic(int sign)
    {
        // N=5 is path K=4. The unit-hop parameter is w=4q² in the physical C# pencil.
        // These decimal anchors are independent of the eigensolver being checked.
        var q = new Complex(sign * Math.Sqrt(WAnchor) / 2, 0);
        var lambda = new Complex(LambdaAnchor, 0);
        var distances = PathKMonodromyScout.AllRootsROdd(K, q)
            .Select(root => (root - lambda).Magnitude).OrderBy(distance => distance).ToArray();
        _out.WriteLine($"q={q}; nearest distances={string.Join(", ", distances.Take(3).Select(d => d.ToString("G17")))}");
        Assert.True(distances[1] < 1e-4, $"Second root is {distances[1]:E6} from the A2 anchor at q={q}");
        Assert.True(distances[2] > 100 * distances[1],
            $"The pair is not isolated: second={distances[1]:E6}, third={distances[2]:E6}");
        double radius = (distances[1] + distances[2]) / 2;
        Assert.Equal(2, distances.Count(distance => distance < radius));

        var reading = PathKMonodromyScout.CharacterizePencilAtROdd(K, q, lambda, radius);
        double relativeDeparture = reading.Character.Departure / reading.Character.CompressionNorm;
        _out.WriteLine($"radius={radius:G17}; kind={reading.Character.Kind}; alg={reading.Character.Algebraic}; geo={reading.Character.Geometric}; " +
            $"relative departure={relativeDeparture:G17}; dephase departure={reading.DephaseScalarDeparture:G17}; hop departure={reading.HopScalarDeparture:G17}");
        _out.WriteLine($"Dephase eigenvalues={string.Join(", ", reading.DephaseEigenvalues)}; hop eigenvalues={string.Join(", ", reading.HopEigenvalues)}");

        Assert.Equal(WAnchor, 4 * q.Real * q.Real, precision: 12);
        Assert.Equal(EpCharacter.EpKind.Diabolic, reading.Character.Kind);
        Assert.Equal(2, reading.Character.Algebraic);
        Assert.Equal(2, reading.Character.Geometric);
        Assert.True(relativeDeparture < 1e-6, $"Relative departure={relativeDeparture:E6}");
        Assert.True(reading.DephaseScalarDeparture < 1e-6, $"Dephase scalar departure={reading.DephaseScalarDeparture:E6}");
        // At real q the projected hop is anti-Hermitian and normal, so its Frobenius norm is
        // the l2 norm of its eigenvalues. A real repeated λ requires hop₂=0 for twin-scalarity;
        // self-normalized HopScalarDeparture is ill-conditioned as this zero test.
        double hopNorm = Math.Sqrt(reading.HopEigenvalues.Sum(value => value.Magnitude * value.Magnitude));
        double relativeHopNorm = hopNorm / Math.Max(1, reading.Character.CompressionNorm);
        _out.WriteLine($"Hop norm={hopNorm:G17}; compressed scale={reading.Character.CompressionNorm:G17}; relative hop norm={relativeHopNorm:G17}");
        Assert.True(relativeHopNorm < 1e-6, $"Relative hop norm={relativeHopNorm:E6}");
    }
}
