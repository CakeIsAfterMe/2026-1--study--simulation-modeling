using DrWatson
@quickactivate "project"
@isdefined(RossModel) || include(srcdir("RossModel.jl"))
using .RossModel
using DataFrames, CSV, Plots, Statistics
using StableRNGs

script_name = "ross_run"
mkpath(plotsdir())
mkpath(datadir())

RUNS = 5
N = 10
S = 3
R = 1
LAMBDA = 100
MU = 1

ross_rng = StableRNG(42)
runs = [sim_repair(ross_rng; N = N, S = S, R = R, lambda = LAMBDA, mu = MU, verbose = true) for i = 1:RUNS]
crash_times = [r[1] for r in runs]
println("Average crash time: ", sum(crash_times) / RUNS)

T_theory = ross_analytic(N = N, S = S, R = R, lambda = LAMBDA, mu = MU)
println("Аналитическое среднее время до падения: ", round(T_theory, digits = 1))

df_monitor = DataFrame([
    merge((run = i, crash_time = round(runs[i][1], digits = 1)), map(x -> round(x, digits = 4), ross_monitor(runs[i][2], runs[i][1]; R = R))) for i = 1:RUNS
])
CSV.write(datadir("ross_run.csv"), df_monitor)
println(df_monitor)

timeline = ross_timeline(runs[1][2]; N = N, S = S)
T_first = runs[1][1]
p1 = plot(
    timeline.times,
    timeline.healthy,
    seriestype = :steppost,
    label = "Исправные машины",
    xlabel = "Time",
    ylabel = "Machines",
    title = "Весь прогон",
)
hline!(p1, [N], linestyle = :dash, label = "N")
p2 = plot(
    timeline.times,
    timeline.healthy,
    seriestype = :steppost,
    linewidth = 2,
    label = "Исправные машины",
    xlabel = "Time",
    ylabel = "Machines",
    title = "Перед падением",
    xlims = (T_first - 50, T_first + 2),
)
hline!(p2, [N], linestyle = :dash, label = "N")
plt_ross = plot(p1, p2, layout = (1, 2), size = (1000, 420), left_margin = 5Plots.mm, bottom_margin = 5Plots.mm)
savefig(plt_ross, plotsdir("ross_healthy.png"))
plt_ross

println("Базовый прогон модели Росса завершён. Результаты в data/ и plots/")
