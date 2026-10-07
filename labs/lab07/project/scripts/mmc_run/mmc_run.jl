using DrWatson
@quickactivate "project"
@isdefined(MMcModel) || include(srcdir("MMcModel.jl"))
using .MMcModel
using DataFrames, CSV, Plots, Statistics

script_name = "mmc_run"
mkpath(plotsdir())
mkpath(datadir())

num_customers = 10
num_servers = 2
mu = 1.0 / 2
lam = 0.9

df = simulate_mmc(lam = lam, mu = mu, c = num_servers, n = num_customers, seed = 123, verbose = true)
CSV.write(datadir("mmc_run.csv"), df)

println(df)
println("Среднее время ожидания: ", round(mean(df.wait), digits = 3))
println("Среднее время обслуживания: ", round(mean(df.service), digits = 3))
println("Доля заявок, которым пришлось ждать: ", mean(df.wait .> 1e-9))

p1 = plot(
    xlabel = "Time",
    ylabel = "Customer",
    title = "Ожидание и обслуживание",
    yticks = 1:num_customers,
    legend = :bottomright,
)
for row in eachrow(df)
    plot!(p1, [row.arrival, row.start], [row.id, row.id], linewidth = 6, color = :orange, label = row.id == 1 ? "Ожидание" : "")
    plot!(p1, [row.start, row.finish], [row.id, row.id], linewidth = 6, color = :steelblue, label = row.id == 1 ? "Обслуживание" : "")
end
savefig(p1, plotsdir("mmc_timeline.png"))
p1

size_data = system_size(df)
p2 = plot(
    size_data.times,
    size_data.counts,
    seriestype = :steppost,
    linewidth = 2,
    label = "Заявок в системе",
    xlabel = "Time",
    ylabel = "Customers",
    title = "Число заявок в системе",
)
hline!(p2, [num_servers], linestyle = :dash, label = "Число каналов")
savefig(p2, plotsdir("mmc_system_size.png"))
p2

println("Базовый прогон M/M/c завершён. Результаты в data/ и plots/")
