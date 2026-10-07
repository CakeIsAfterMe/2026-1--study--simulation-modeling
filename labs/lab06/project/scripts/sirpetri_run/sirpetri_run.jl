using DrWatson
@quickactivate "project"
using Random
@isdefined(SIRPetri) || include(srcdir("SIRPetri.jl"))
using .SIRPetri
using DataFrames, CSV, Plots

script_name = "sirpetri_run"
mkpath(plotsdir())
mkpath(datadir())

β = 0.3
γ = 0.1
tmax = 100.0

net, u0, states = build_sir_network(β, γ)
u0

df_det = simulate_deterministic(net, u0, (0.0, tmax), saveat = 0.5, rates = [β, γ])
CSV.write(datadir("sir_det.csv"), df_det)
println("Пик I (детерминированная модель): ", round(maximum(df_det.I), digits = 2))
println("Время пика: ", df_det.time[argmax(df_det.I)])
println("Конечное R: ", round(df_det.R[end], digits = 2))

Random.seed!(123)
df_stoch = simulate_stochastic(net, u0, (0.0, tmax), rates = [β, γ])
CSV.write(datadir("sir_stoch.csv"), df_stoch)
println("Число событий: ", nrow(df_stoch) - 1)
println("Пик I (стохастическая модель): ", maximum(df_stoch.I))
println("Время последнего события: ", round(df_stoch.time[end], digits = 2))
println("Конечное R: ", df_stoch.R[end])

p_det = plot_sir(df_det)
savefig(p_det, plotsdir("sir_det_dynamics.png"))
p_det

p_stoch = plot_sir(df_stoch)
savefig(p_stoch, plotsdir("sir_stoch_dynamics.png"))
p_stoch

println("Базовый прогон завершён. Результаты в data/ и plots/")
