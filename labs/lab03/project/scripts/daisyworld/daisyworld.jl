using DrWatson
@quickactivate "project"
using Agents
using DataFrames
using CairoMakie

include(srcdir("daisyworld.jl"))

model = daisyworld()

daisycolor(a::Daisy) = a.breed

plotkwargs = (
    agent_color = daisycolor, agent_size = 20, agent_marker = '✿',
    heatarray = :temperature,
    heatkwargs = (colorrange = (-20, 60),),
)

plt1, _ = abmplot(model; plotkwargs...)
plt1

step!(model, 5)
plt2, _ = abmplot(model; plotkwargs...)
plt2

step!(model, 40)
plt3, _ = abmplot(model; plotkwargs...)
plt3

save(plotsdir("daisy_step001.png"), plt1)
save(plotsdir("daisy_step005.png"), plt2)
save(plotsdir("daisy_step040.png"), plt3)
