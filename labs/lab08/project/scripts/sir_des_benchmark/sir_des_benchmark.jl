using DrWatson
@quickactivate "project"
@isdefined(SIRDES) || include(srcdir("sir_model.jl"))
using .SIRDES
using Random, Plots, DataFrames, CSV, BenchmarkTools

script_name = "sir_des_benchmark"
mkpath(plotsdir())
mkpath(datadir())

tmax = 40.0
p = [0.05, 10.0, 0.25]

function prepared_model(n, p)
    m = MakeSIRModel([n - 10, 10, 0], p)
    activate(m)
    return m
end

Random.seed!(1234)

bench = @benchmark sir_run(m, $tmax) setup = (m = prepared_model(10000, $p)) evals = 1 samples = 3 seconds = 600
display(bench)

function time_run(n, p, tmax)
    m = prepared_model(n, p)
    return @elapsed sir_run(m, tmax)
end

sizes = [1000, 2000, 5000, 10000]
times = [time_run(n, p, tmax) for n in sizes]
df_bench = DataFrame(N = sizes, time_s = round.(times, digits = 3))
CSV.write(datadir("sir_des_benchmark.csv"), df_bench)
println(df_bench)

plt_bench = plot(
    df_bench.N,
    df_bench.time_s,
    marker = :circle,
    linewidth = 2,
    label = "Время прогона",
    xlabel = "Размер популяции N",
    ylabel = "Время, с",
    title = "Производительность",
    legend = :topleft,
)
savefig(plt_bench, plotsdir("sir_des_benchmark.png"))
plt_bench
