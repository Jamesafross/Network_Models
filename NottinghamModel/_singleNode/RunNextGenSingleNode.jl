using DifferentialEquations,Plots,Parameters,LinearAlgebra,JLD,Statistics
include("functions.jl")
include("../_globalFunctions/_includes.jl")


NGp = NextGen2PopParams()

NGp.η_0I = 25.

NGp.κ=0.0203

tspan = (0.,50000.)
const AMP=0. # or 10.
const fr = 1.

const U0 = 0.1



u0 =rand(14)

#u0[1] = 2.257
#u0[2] = 0.036776
#u0[3] = 1.6943
#u0[4] = -0.98019
#u0[5] = 3.3855
#u0[6] = u0[5]
#u0[7] = 2.257
#u0[8] = u0[7]
#u0[9] = 0.073552
#u0[10] = u0[9]
#u0[11] = 0.11033
#u0[12] = u0[11]


#u0[39] = 100.

u0[13] = NGp.κSEI

p = NGp
probSS = SteadyStateProblem(nextgen_stp_de,u0,p)
prob = ODEProblem(nextgen_stp_de,u0,tspan,p)

save_start = 9950
save_end = 10200
sol = solve(prob,BS3(),maxiters=10e20,saveat=collect(save_start:0.01:save_end))



plot(sol.t[1:end],sol[1,1:end],label="rE")
plot!(sol.t[1:end],sol[2,1:end],label="rI")

