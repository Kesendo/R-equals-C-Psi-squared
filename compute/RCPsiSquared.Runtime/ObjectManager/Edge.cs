namespace RCPsiSquared.Runtime.ObjectManager;

/// <summary>One typed inheritance edge: <paramref name="Child"/> depends on <paramref name="Parent"/>
/// because the child's factory called <c>b.Get&lt;Parent&gt;()</c>. <paramref name="Reason"/>
/// is a short human description provided by the registration site (defaults to a generic string
/// when not set explicitly).</summary>
public sealed record Edge(Type Parent, Type Child, string Reason);
