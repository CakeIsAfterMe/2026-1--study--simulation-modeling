using DrWatson
@quickactivate "project"
@isdefined(SIRDES) || include(srcdir("sir_model.jl"))
using .SIRDES
using Random, Plots, DataFrames, CSV

script_name = "sir_des_params"
mkpath(plotsdir())
mkpath(datadir())

des_dict = Dict(
    :β => [0.03, 0.05, 0.07],
    :c => [5.0, 10.0, 15.0],
    :γ => [0.25, 0.5],
    :u0 => [[990, 10, 0]],
    :tmax => 40.0,
    :seed => 1234,
)

des_list = dict_list(des_dict)
length(des_list)

function run_des(p)
    q = [p[:β], p[:c], p[:γ]]
    Random.seed!(p[:seed])
    m = MakeSIRModel(p[:u0], q)
    activate(m)
    sir_run(m, p[:tmax])
    des = sir_metrics(out(m))
    ode = sir_metrics(sir_ode_solution(p[:u0], q, p[:tmax]))
    return (
        β = p[:β],
        c = p[:c],
        γ = p[:γ],
        R0 = round(p[:β] * p[:c] / p[:γ], digits = 2),
        peak_I = des.peak_I,
        t_peak = round(des.t_peak, digits = 1),
        final_R = des.final_R,
        peak_I_ode = round(ode.peak_I, digits = 1),
        final_R_ode = round(ode.final_R, digits = 1),
    )
end

df_params = DataFrame([run_des(p) for p in des_list])
sort!(df_params, [:γ, :β, :c])
CSV.write(datadir("sir_des_params.csv"), df_params)
println(df_params)

p1 = scatter(df_params.R0, df_params.peak_I, label = "DES", xlabel = "R0", ylabel = "Peak I", title = "Пик эпидемии", legend = :topleft)
scatter!(p1, df_params.R0, df_params.peak_I_ode, marker = :xcross, label = "ОДУ")
p2 = scatter(df_params.R0, df_params.final_R, label = "DES", xlabel = "R0", ylabel = "R(tmax)", title = "Переболело", legend = :topleft)
scatter!(p2, df_params.R0, df_params.final_R_ode, marker = :xcross, label = "ОДУ")
vline!(p1, [1.0], linestyle = :dash, label = "R0 = 1")
vline!(p2, [1.0], linestyle = :dash, label = "R0 = 1")
plt_params = plot(p1, p2, layout = (1, 2), size = (1000, 420), left_margin = 5Plots.mm, bottom_margin = 5Plots.mm)

savefig(plt_params, plotsdir("sir_des_params.png"))
println("Результаты сохранены в data/sir_des_params.csv и plots/sir_des_params.png")
