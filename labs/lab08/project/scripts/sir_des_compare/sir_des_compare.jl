using DrWatson
@quickactivate "project"
@isdefined(SIRDES) || include(srcdir("sir_model.jl"))
using .SIRDES
using Random, Plots, DataFrames, CSV

script_name = "sir_des_compare"
mkpath(plotsdir())
mkpath(datadir())

tmax = 40.0
u0 = [990, 10, 0]
p = [0.05, 10.0, 0.25]

Random.seed!(1234)
m_exp = MakeSIRModel(u0, p)
activate(m_exp)
sir_run(m_exp, tmax)
df_exp = out(m_exp)

Random.seed!(1234)
m_fix = MakeSIRModel(u0, p; fixed = true)
activate(m_fix)
sir_run(m_fix, tmax)
df_fix = out(m_fix)

df_ode = sir_ode_solution(u0, p, tmax)

df_cmp = DataFrame([
    merge((version = "Экспоненциальная",), map(float, sir_metrics(df_exp))),
    merge((version = "Фиксированная",), map(float, sir_metrics(df_fix))),
    merge((version = "ОДУ",), map(x -> round(x, digits = 1), sir_metrics(df_ode))),
])
CSV.write(datadir("sir_des_compare.csv"), df_cmp)
println(df_cmp)

p1 = plot(df_exp.t, df_exp.I, label = "Экспоненциальная", xlabel = "Время", ylabel = "I", title = "Больные", linewidth = 2)
plot!(p1, df_fix.t, df_fix.I, label = "Фиксированная", linewidth = 2)
plot!(p1, df_ode.t, df_ode.I, label = "ОДУ", linestyle = :dash, linewidth = 2)
p2 = plot(df_exp.t, df_exp.R, label = "Экспоненциальная", xlabel = "Время", ylabel = "R", title = "Переболевшие", linewidth = 2, legend = :topleft)
plot!(p2, df_fix.t, df_fix.R, label = "Фиксированная", linewidth = 2)
plot!(p2, df_ode.t, df_ode.R, label = "ОДУ", linestyle = :dash, linewidth = 2)
plt_cmp = plot(p1, p2, layout = (1, 2), size = (1000, 420), left_margin = 5Plots.mm, bottom_margin = 5Plots.mm)
savefig(plt_cmp, plotsdir("sir_des_compare.png"))
plt_cmp
