# # Модель SIR: набор параметров
#
# Исследуем, как число контактов $c$ влияет на ход эпидемии.
# Остальные параметры фиксированы: $\beta = 0.05$, $\gamma = 0.25$.
# Уменьшение $c$ соответствует карантину.

# ## Активация проекта и загрузка пакетов

using DrWatson
@quickactivate "project"

using DifferentialEquations
using DataFrames
using CSV
using LaTeXStrings
using Plots

script_name = "sir_params"
mkpath(plotsdir(script_name))
mkpath(datadir(script_name))

# ## Определение модели

function sir_ode!(du, u, p, t)
    (S, I, R) = u
    (β, c, γ) = p
    N = S + I + R
    @inbounds begin
        du[1] = -β * c * I / N * S
        du[2] = β * c * I / N * S - γ * I
        du[3] = γ * I
    end
    nothing
end

# ## Сетка параметров
#
# Все параметры собраны в `Dict`. Меняется только число контактов $c$.

param_grid = Dict(
    :β => [0.05],
    :c => [5.0, 10.0, 15.0, 20.0],
    :γ => [0.25],
    :u0 => [[990.0, 10.0, 0.0]],
    :tmax => [100.0],
    :δt => [0.1],
)

# Функция `dict_list` строит все комбинации параметров:

all_params = dict_list(param_grid)
println("Всего комбинаций параметров: ", length(all_params))

# ## Функция одного эксперимента
#
# Функция решает задачу для одного набора параметров и возвращает
# решение и основные показатели.

function run_sir(params::Dict)
    @unpack β, c, γ, u0, tmax, δt = params
    prob = ODEProblem(sir_ode!, u0, (0.0, tmax), [β, c, γ])
    sol = solve(prob, Tsit5(), saveat = δt)
    S = [u[1] for u in sol.u]
    I = [u[2] for u in sol.u]
    R = [u[3] for u in sol.u]
    N = sum(u0)
    peak_idx = argmax(I)
    return Dict(
        :t => sol.t,
        :S => S,
        :I => I,
        :R => R,
        :R0 => c * β / γ,
        :peak_value => I[peak_idx],
        :peak_time => sol.t[peak_idx],
        :final_share => R[end] / N * 100,
    )
end

# ## Запуск всех экспериментов

results = []

plt_I = plot(xlabel="Время, дни", ylabel="Количество инфицированных",
    title="Число заразных при разных c", legend=:topright, grid=true, size=(800, 500))

plt_R = plot(xlabel="Время, дни", ylabel="Количество переболевших",
    title="Число переболевших при разных c", legend=:bottomright, grid=true, size=(800, 500))

for params in all_params
    res = run_sir(params)
    push!(results, (c = params[:c],
        R0 = res[:R0],
        peak_value = round(res[:peak_value], digits=1),
        peak_time = round(res[:peak_time], digits=1),
        final_share = round(res[:final_share], digits=1)))
    lbl = "c = $(params[:c]), R0 = $(round(res[:R0], digits=1))"
    plot!(plt_I, res[:t], res[:I], label=lbl, linewidth=2)
    plot!(plt_R, res[:t], res[:R], label=lbl, linewidth=2)
end

# ## Сводная таблица

results_df = DataFrame(results)
println(results_df)

# ## Графики
#
# Число заразных:

plt_I

# Число переболевших:

plt_R

# Зависимость высоты пика и итоговой доли переболевших от $R_0$:

plt_peak = plot(results_df.R0, results_df.peak_value,
    seriestype=:scatter, markersize=8, color=:red, label="Пик заразных",
    xlabel=L"R_0", ylabel="Количество людей",
    title="Высота пика эпидемии", legend=:topleft, grid=true, size=(800, 400))

# Пунктиром показан порог коллективного иммунитета $1 - 1/R_0$.

plt_share = plot(results_df.R0, results_df.final_share,
    seriestype=:scatter, markersize=8, color=:blue, label="Доля переболевших",
    xlabel=L"R_0", ylabel="Доля популяции, %",
    title="Итоговая доля переболевших", legend=:bottomright, grid=true, size=(800, 400))

R0_range = 1.0:0.05:4.0
plot!(plt_share, R0_range, (1 .- 1 ./ R0_range) .* 100,
    linestyle=:dash, color=:purple, linewidth=1.5, label="Порог коллективного иммунитета")
plt_share

# ## Сохранение результатов

savefig(plt_I, plotsdir(script_name, "sir_scan_infected.png"))
savefig(plt_R, plotsdir(script_name, "sir_scan_recovered.png"))
savefig(plt_peak, plotsdir(script_name, "sir_scan_peak.png"))
savefig(plt_share, plotsdir(script_name, "sir_scan_share.png"))
CSV.write(datadir(script_name, "sir_scan.csv"), results_df)
