using DifferentialEquations,Plots,Parameters,LinearAlgebra,JLD,Statistics
include("functions.jl")
include("../../_GlobalFunctions/_includes.jl")



NGp = NextGen2PopParams()
NGp.η_0I = -10
NGp.κ=0.0203


tspan = (0.,15000.)
const AMP=0.0
const fr = 2.


u0 = complex(rand(10))

p = NGp

prob = ODEProblem(nextgen_stp_de,u0,tspan,p)

save_start = 9950
save_end = 10200
sol = solve(prob,maxiters=10e20,saveat=save_start:0.01:save_end)


rE = (FR.(sol[1,:],NGp.τE))
rI = (FR.(sol[2,:],NGp.τI))

plot(sol.t,rE)
plot!(sol.t,rI)


