using DifferentialEquations
using Makie

β=0.5
function surfaceTensor!(du,u,p,t)
    du[1] = cos(u[3])
    du[2] = sin(u[3])
    du[3] = 2+β*u[2]-sin(u[3])/u[1]
end
u0 = [0.1;0.0;0.0]
tspan = (0.0,1.5)
prob = ODEProblem(surfaceTensor!,u0,tspan)
sol = solve(prob)


ps=hcat(sol.u...)
lines(ps[1:2,:])