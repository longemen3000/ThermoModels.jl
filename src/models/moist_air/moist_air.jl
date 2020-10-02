#Moist air, ported from  pjabardo/Psychro.jl 


"""
    MoistAir
available specs: (Pressure,Saturation=true)*
available specs: (Pressure,Temperature,Mass Fraction)*
Calculates properties of moist air, using ASHRAE correlations [1]. Water saturation curves are calculated via the WaterIF97 model.

there

"""
struct ASHRAEMoistAir <: GasModel end
