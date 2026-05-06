namespace RCPsiSquared.Runtime;

// Public (not internal) so the smoke test in RCPsiSquared.Runtime.Tests can reference
// typeof(Placeholder) without an InternalsVisibleTo attribute. Removed once
// ClaimRegistry replaces it in Task 4.
public static class Placeholder
{
}
