using DrWatson
@quickactivate "project"
using DataFrames, Plots, CSV
using Statistics

datafile = datadir("beta_scan_all.csv")
isfile(datafile) || @warn "Файл не найден, сначала запустите sir_scan_beta.jl"
df = CSV.read(datafile, DataFrame)
first(df, 5)

grouped = combine(
    groupby(df, [:beta]),
    :peak => mean => :peak,
    :final_inf => mean => :final_inf,
    :final_rec => mean => :final_rec,
    :deaths => mean => :deaths,
)
println(grouped)

p1 = plot(grouped.beta, grouped.peak, label = "Пик", xlabel = "β",
    ylabel = "Доля инфицированных", marker = :circle, linewidth = 2)
plot!(p1, grouped.beta, grouped.final_inf, label = "Конечная", marker = :square)

p2 = plot(grouped.beta, grouped.deaths, label = false, xlabel = "β",
    ylabel = "Число умерших", marker = :circle, linewidth = 2, color = :red)

p3 = plot(grouped.beta, grouped.final_rec, label = false, xlabel = "β",
    ylabel = "Доля выздоровевших", marker = :circle, linewidth = 2, color = :green)

plt = plot(p1, p2, p3, layout = (3, 1), size = (800, 900))

savefig(plt, plotsdir("comprehensive_analysis.png"))
