namespace RCPsiSquared.Orchestration;

// Public (not internal) so the smoke test in RCPsiSquared.Orchestration.Tests can reference
// typeof(Placeholder) without an InternalsVisibleTo attribute. Removed once
// KnowledgeCli replaces it in Task 14.
public static class Placeholder
{
}
