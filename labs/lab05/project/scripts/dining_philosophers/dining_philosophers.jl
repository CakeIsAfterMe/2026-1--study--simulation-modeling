using DrWatson
@quickactivate "project"
include(srcdir("DiningPhilosophers.jl"))
using .DiningPhilosophers
using DataFrames, CSV, Plots, Random

N = 5
tmax = 50.0
seed = 42

println("=== Классическая сеть (без арбитра) ===")
net_classic, u0_classic, _ = build_classical_network(N)
df_classic = simulate_stochastic(net_classic, u0_classic, tmax; rng = Xoshiro(seed))
CSV.write(datadir("dining_classic.csv"), df_classic)
dead = detect_deadlock(df_classic, net_classic)
println("Deadlock обнаружен: $dead")
println("Число событий: ", nrow(df_classic) - 1)
println("Время последнего события: ", round(df_classic.time[end], digits = 2))

plot_classic = plot_marking_evolution(df_classic, N)
savefig(plot_classic, plotsdir("classic_simulation.png"))
plot_classic

println("=== Сеть с арбитром ===")
net_arb, u0_arb, _ = build_arbiter_network(N)
df_arb = simulate_stochastic(net_arb, u0_arb, tmax; rng = Xoshiro(seed))
CSV.write(datadir("dining_arbiter.csv"), df_arb)
dead_arb = detect_deadlock(df_arb, net_arb)
println("Deadlock обнаружен: $dead_arb")
println("Число событий: ", nrow(df_arb) - 1)
println("Время последнего события: ", round(df_arb.time[end], digits = 2))

plot_arb = plot_marking_evolution(df_arb, N)
savefig(plot_arb, plotsdir("arbiter_simulation.png"))
plot_arb

df_ode = simulate_ode(net_classic, u0_classic, tmax)
CSV.write(datadir("dining_classic_ode.csv"), df_ode)
println("=== Классическая сеть, ODE ===")
println("Think_1 в конце: ", round(df_ode[end, "Think_1"], digits = 3))
println("Hungry_1 в конце: ", round(df_ode[end, "Hungry_1"], digits = 3))
println("Eat_1 в конце: ", round(df_ode[end, "Eat_1"], digits = 3))
println("Fork_1 в конце: ", round(df_ode[end, "Fork_1"], digits = 3))

plot_ode = plot_marking_evolution(df_ode, N)
savefig(plot_ode, plotsdir("classic_ode.png"))
plot_ode
