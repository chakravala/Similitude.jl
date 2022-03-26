using Similitude, Test

@test 1000normal(molarmass(SI2019)) ≈ normal(molarmass(EMU2019))
@test normal(luminousefficacy(SI2019)) ≈ 1e7*normal(luminousefficacy(EMU2019))
