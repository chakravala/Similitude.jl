using Similitude, Test

@test 1000normal(molarmass(Metric)) == 1normal(molarmass(Gauss))
@test 1normal(luminousefficacy(Metric)) == 1e7*normal(luminousefficacy(Gauss))
