using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Diagnostics.Knowledge;
using Xunit;

namespace RCPsiSquared.Runtime.Tests.PolarityArchitecture;

/// <summary>Dedicated cumulative-cardinality assertion for the YParity-axis bucket
/// in <see cref="PolarityCubeMap"/>. This file is the single point of truth for the
/// expected total count of Z2Axis.YParity claims; each new regime (F102, F103, F105,
/// F106, F107, F109, F110, ...) bumps the count here rather than editing one
/// regime's test file. Per-regime <c>F87Z2CubedRefinement{NxKy}RegistrationTests</c>
/// files assert only their own claim's membership.</summary>
public class PolarityCubeMapYParityCardinalityTests
{
    [Fact]
    public void PolarityCubeMap_YParityClaims_TotalCount_MatchesExpected()
    {
        // F102 (YParityIndependenceAtK3) + F103 (F87Z2CubedRefinementN4K3) + F105
        // (F87Z2CubedRefinementN5K3) + F106 (F87Z2CubedRefinementN4K4) + F107
        // (TrulyYParityZeroPurity, first derived-not-empirical) + F109
        // (MotherSoftYParityOnePurity, second derived-not-empirical, fully
        // unconditional after F108 Part 1+2+3 closure 2026-05-25) + F110
        // (HardCellYInversionPattern, Tier1Candidate Aspect A closed-form + B/C
        // empirical) are the current Z2Axis.YParity members. F108 Part 1+3 are
        // BitB-axis; F108 Part 2 is BitA-axis (none are YParity, so not counted
        // here). When a new YParity regime lands (F111+ etc), bump this expected
        // value to match.
        const int expectedCount = 7;

        var registry = KnowledgeRegistryFactory.BuildDefault();
        var cubeMap = registry.Get<PolarityCubeMap>();
        Assert.NotNull(cubeMap);
        Assert.Equal(expectedCount, cubeMap.YParityClaims.Count);
    }
}
