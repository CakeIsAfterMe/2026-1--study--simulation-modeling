using DrWatson
@quickactivate "project"
@isdefined(SIRDES) || include(srcdir("sir_model.jl"))
using .SIRDES
using Random, StatsPlots, DataFrames, CSV

script_name = "sir_des"
mkpath(plotsdir())
mkpath(datadir("sims"))

tmax = 40.0
u0 = [990, 10, 0]      # S, I, R
p = [0.05, 10.0, 0.25] # β, c, γ

Random.seed!(1234)

des_model = MakeSIRModel(u0, p)
activate(des_model)
sir_run(des_model, tmax)
data_des = out(des_model)
println("Число событий: ", nrow(data_des) - 1)
println(sir_metrics(data_des))

filename = "sir_$(u0[1])_$(u0[2])_$(p[1])_$(p[2])_$(p[3]).csv"
CSV.write(datadir("sims", filename), data_des)
println("Таблица сохранена в data/sims/", filename)

plt_des = @df data_des plot(
    :t,
    [:S :I :R],
    labels = ["S" "I" "R"],
    xlab = "Время",
    ylab = "Численность",
    title = "Дискретно-событийная SIR модель",
    linewidth = 2,
)
savefig(plt_des, plotsdir("sir_des.png"))
plt_des
