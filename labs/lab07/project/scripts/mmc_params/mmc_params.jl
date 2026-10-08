using DrWatson
@quickactivate "project"
@isdefined(MMcModel) || include(srcdir("MMcModel.jl"))
using .MMcModel
using DataFrames, CSV, Plots, Statistics

script_name = "mmc_params"
mkpath(plotsdir())
mkpath(datadir())

mmc_dict = Dict(
    :c => [2, 3, 4],
    :lam => [0.5, 0.7, 0.9],
    :mu => 0.5,
    :n => 20000,
    :seed => 123,
)

mmc_list = dict_list(mmc_dict)
length(mmc_list)

function run_mmc(p)
    res = simulate_mmc(lam = p[:lam], mu = p[:mu], c = p[:c], n = p[:n], seed = p[:seed])
    th = mmc_analytic(p[:lam], p[:mu], p[:c])
    return (
        c = p[:c],
        lam = p[:lam],
        rho = round(th.rho, digits = 3),
        Wq_sim = round(mean(res.wait), digits = 4),
        Wq_theory = round(th.Wq, digits = 4),
        Pwait_sim = round(mean(res.wait .> 1e-9), digits = 4),
        Pwait_theory = round(th.Pwait, digits = 4),
    )
end

df_mmc = DataFrame([run_mmc(p) for p in mmc_list])
sort!(df_mmc, [:c, :lam])
CSV.write(datadir("mmc_params.csv"), df_mmc)
println(df_mmc)

p1 = plot(xlabel = "λ", ylabel = "Wq", title = "Время ожидания", legend = :topleft)
p2 = plot(xlabel = "λ", ylabel = "Pwait", title = "Вероятность ожидания", legend = :topleft)
for c in sort(unique(df_mmc.c))
    part = df_mmc[df_mmc.c .== c, :]
    scatter!(p1, part.lam, part.Wq_sim, label = "c = $c, модель")
    plot!(p1, part.lam, part.Wq_theory, linestyle = :dash, label = "c = $c, теория")
    scatter!(p2, part.lam, part.Pwait_sim, label = "c = $c, модель")
    plot!(p2, part.lam, part.Pwait_theory, linestyle = :dash, label = "c = $c, теория")
end
plt_mmc = plot(p1, p2, layout = (1, 2), size = (1000, 420), left_margin = 5Plots.mm, bottom_margin = 5Plots.mm)

savefig(plt_mmc, plotsdir("mmc_params.png"))
println("Результаты сохранены в data/mmc_params.csv и plots/mmc_params.png")
