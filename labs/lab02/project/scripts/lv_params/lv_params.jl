using DrWatson
@quickactivate "project"

using DifferentialEquations
using DataFrames
using CSV
using LaTeXStrings
using Plots
using Statistics

script_name = "lv_params"
mkpath(plotsdir(script_name))
mkpath(datadir(script_name))

function lotka_volterra!(du, u, p, t)
    x, y = u
    α, β, δ, γ = p
    @inbounds begin
        du[1] = α*x - β*x*y
        du[2] = δ*x*y - γ*y
    end
    nothing
end

param_grid = Dict(
    :α => [0.05, 0.1, 0.2, 0.3],
    :β => [0.02],
    :δ => [0.01],
    :γ => [0.3],
    :u0 => [[40.0, 9.0]],
    :tmax => [200.0],
    :save_dt => [0.1],
)

all_params = dict_list(param_grid)
println("Всего комбинаций параметров: ", length(all_params))

function find_peaks(signal, time)
    peaks = Float64[]
    for i in 2:length(signal)-1
        if signal[i] > signal[i-1] && signal[i] > signal[i+1]
            push!(peaks, time[i])
        end
    end
    return peaks
end

function run_lv(params::Dict)
    @unpack α, β, δ, γ, u0, tmax, save_dt = params
    prob = ODEProblem(lotka_volterra!, u0, (0.0, tmax), [α, β, δ, γ])
    sol = solve(prob, Tsit5(), reltol = 1e-8, abstol = 1e-10, saveat = save_dt)
    prey = [u[1] for u in sol.u]
    predator = [u[2] for u in sol.u]
    peaks = find_peaks(prey, sol.t)
    period = length(peaks) > 1 ? mean(diff(peaks)) : NaN
    return Dict(
        :t => sol.t,
        :prey => prey,
        :predator => predator,
        :x_star => γ / δ,
        :y_star => α / β,
        :period => period,
        :period_theory => 2π / sqrt(α * γ),
    )
end

results = []

plt_prey = plot(xlabel="Время", ylabel="Популяция жертв",
    title="Жертвы при разных α", legend=:topright, grid=true, size=(900, 500))

plt_predator = plot(xlabel="Время", ylabel="Популяция хищников",
    title="Хищники при разных α", legend=:topright, grid=true, size=(900, 500))

plt_phase = plot(xlabel="Популяция жертв (x)", ylabel="Популяция хищников (y)",
    title="Фазовые портреты при разных α", legend=:topright, grid=true, size=(800, 600))

for params in all_params
    res = run_lv(params)
    push!(results, (α = params[:α],
        y_star = res[:y_star],
        period = round(res[:period], digits=2),
        period_theory = round(res[:period_theory], digits=2),
        prey_min = round(minimum(res[:prey]), digits=2),
        prey_max = round(maximum(res[:prey]), digits=2),
        predator_min = round(minimum(res[:predator]), digits=2),
        predator_max = round(maximum(res[:predator]), digits=2)))
    lbl = "α = $(params[:α])"
    plot!(plt_prey, res[:t], res[:prey], label=lbl, linewidth=1.5)
    plot!(plt_predator, res[:t], res[:predator], label=lbl, linewidth=1.5)
    plot!(plt_phase, res[:prey], res[:predator], label=lbl, linewidth=1.5)
end

results_df = DataFrame(results)
println(results_df)

plt_prey

plt_predator

plt_phase

plt_period = plot(results_df.α, results_df.period,
    seriestype=:scatter, markersize=8, color=:red, label="Численное решение",
    xlabel="Скорость размножения жертв, α", ylabel="Период колебаний",
    title="Зависимость периода от α", legend=:topright, grid=true, size=(800, 400))

α_range = 0.05:0.005:0.3
plot!(plt_period, α_range, 2π ./ sqrt.(α_range .* 0.3),
    linestyle=:dash, color=:blue, linewidth=1.5, label="Теория")
plt_period

savefig(plt_prey, plotsdir(script_name, "lv_scan_prey.png"))
savefig(plt_predator, plotsdir(script_name, "lv_scan_predator.png"))
savefig(plt_phase, plotsdir(script_name, "lv_scan_phase.png"))
savefig(plt_period, plotsdir(script_name, "lv_scan_period.png"))
CSV.write(datadir(script_name, "lv_scan.csv"), results_df)
