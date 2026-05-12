using RCPsiSquared.Core.Symmetry;

namespace RCPsiSquared.Core.F89PathK;

/// <summary>F89PathK convention pointer: F89-side primitives use H = J·(XX+YY)
/// (NO 1/2 prefactor) and column-major vec(M)[b·d + a] = M[a, b]. The J convention
/// is locked at <see cref="F90F86C2BridgeIdentity.JConventionFactor"/> = 2.0
/// relative to F86's H = (J/2)·(XX+YY). Parallel C# port of
/// <c>simulations/_f89_pathk_lib.py</c>.</summary>
public static class F89PathKConvention
{
}
