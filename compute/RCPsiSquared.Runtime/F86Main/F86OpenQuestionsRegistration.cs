using RCPsiSquared.Core.F86;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.F86Main;

/// <summary>Marker type for <see cref="OpenQuestionCollection{TFamilyMarker}"/> at the F86
/// family.</summary>
public sealed class F86Marker { }

/// <summary>Registers <see cref="F86OpenQuestions.Standard"/> as a single
/// <see cref="OpenQuestionCollection{TFamilyMarker}"/> at the OpenQuestion Tier, beside the
/// explicit <see cref="FullIrreducibleSrpClassQuestion"/>, which inherits from the secured
/// <see cref="ShiftedGeneratorSectorwisePClaim"/>. The collection's items are Q_peak-response
/// questions (the HWHM_left/Q_peak closed form, the 4-mode construction at c ≥ 3, the σ_0
/// asymptote), none of them about the SRP class, so the collection hangs under
/// <see cref="QEpLaw"/>, the Q-axis anchor those responses are measured against.</summary>
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
                _ = b.Get<QEpLaw>();
                return new OpenQuestionCollection<F86Marker>(
                    items: F86OpenQuestions.Standard,
                    anchor: "compute/RCPsiSquared.Core/F86/F86OpenQuestions.cs");
            });
}
