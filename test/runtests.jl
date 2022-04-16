using Similitude, Test

@test 1000normal(molarmass(Metric)) ≈ normal(molarmass(Gauss))
@test normal(luminousefficacy(Metric)) ≈ 1e7*normal(luminousefficacy(Gauss))
