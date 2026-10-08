using DrWatson
@quickactivate "project"
@isdefined(SIRDES) || include(srcdir("sir_model.jl"))
using .SIRDES
using Random, Plots, DataFrames, CSV

script_name = "sir_des_demography"
mkpath(plotsdir())
mkpath(datadir())

tmax = 200.0
u0 = [990, 10, 0]
p = [0.05, 10.0, 0.25]
μ = 0.05

N0 = sum(u0)
S_star = N0 * (p[3] + μ) / (p[1] * p[2])
I_star = μ * (N0 - S_star) / (p[3] + μ)
R_star = N0 - S_star - I_star
println("Равновесие: S* = ", round(S_star, digits = 1), ", I* = ", round(I_star, digits = 1), ", R* = ", round(R_star, digits = 1))

Random.seed!(1234)
m_dem = MakeXModel(u0, p; μ = μ)
activate_x(m_dem)
run_x(m_dem, tmax)
df_dem = out_x(m_dem)
CSV.write(datadir("sir_des_demography.csv"), df_dem)
println("Число событий: ", nrow(df_dem) - 1)
println("Больных в конце: ", df_dem.I[end])

df_avg = DataFrame(
    group = ["S", "I", "R", "N"],
    des = [round(time_average(df_dem.t, df_dem[!, g], 100.0, tmax), digits = 1) for g in ["S", "I", "R", "N"]],
    theory = round.([S_star, I_star, R_star, N0], digits = 1),
)
println(df_avg)

df_dem_ode = sir_ode_extended(u0, p, tmax; μ = μ)

p1 = plot(df_dem.t, [df_dem.S df_dem.I df_dem.R], label = ["S" "I" "R"], xlabel = "Время", ylabel = "Численность", title = "Модель с демографией", linewidth = 2)
plot!(p1, df_dem_ode.t, [df_dem_ode.S df_dem_ode.I df_dem_ode.R], label = ["S, ОДУ" "I, ОДУ" "R, ОДУ"], linestyle = :dash)
p2 = plot(df_dem.t, df_dem.I, label = "I", xlabel = "Время", ylabel = "I", title = "Больные", linewidth = 2)
plot!(p2, df_dem_ode.t, df_dem_ode.I, label = "I, ОДУ", linestyle = :dash, linewidth = 2)
hline!(p2, [I_star], linestyle = :dot, label = "I*")
plt_dem = plot(p1, p2, layout = (1, 2), size = (1000, 420), left_margin = 5Plots.mm, bottom_margin = 5Plots.mm)
savefig(plt_dem, plotsdir("sir_des_demography.png"))
plt_dem
