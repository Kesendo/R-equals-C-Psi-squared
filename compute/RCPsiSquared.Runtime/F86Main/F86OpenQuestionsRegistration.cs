using RCPsiSquared.Core.F86;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.F86Main;

/// <summary>Marker type for <see cref="OpenQuestionCollection{TFamilyMarker}"/> at the F86
/// family.</summary>
public sealed class F86Marker { }

/// <summary>Registers <see cref="F86OpenQuestions.Standard"/> as a single
/// <see cref="OpenQuestionCollection{TFamilyMarker}"/> at the OpenQuestion Tier. The
/// collection is registered alongside the explicit <see cref="FullIrreducibleSrpClassQuestion"/> and
/// inherits from it. That question in turn inherits from the secured
/// <see cref="ShiftedGeneratorSectorwisePClaim"/>.</summary>
public static class F86OpenQuestionsRegistration
{
    public static ClaimRegistryBuilder RegisterF86OpenQuestions(this ClaimRegistryBuilder builder) =>
        builder
            .Register<FullIrreducibleSrpClassQuestion>(b =>
            {
                _ = b.Get<ShiftedGeneratorSectorwisePClaim>();
                return new FullIrreducibleSrpClassQuestion();
            })
            .Register<OpenQuestionCollection<F86Marker>>(b =>
            {
                _ = b.Get<FullIrreducibleSrpClassQuestion>();
                return new OpenQuestionCollection<F86Marker>(
                    items: F86OpenQuestions.Standard,
                    anchor: "compute/RCPsiSquared.Core/F86/F86OpenQuestions.cs");
            });
}
