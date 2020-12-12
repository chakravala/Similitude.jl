using UnitfulSystems, Test

@test molarmass(SI2019) == molarmass(EMU2019)
@test luminousefficacy(SI2019) == luminousefficacy(EMU2019)
