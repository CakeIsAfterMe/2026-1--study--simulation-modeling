# # Модель SIR в виде сети Петри: набор параметров
#
# При $\beta$ от 0.1 до 0.8 заражается всё население сразу, потому что
# скорость заражения в модели равна $\beta S I$ и при $S \approx 1000$
# получается очень большой. Чтобы увидеть порог эпидемии, возьмём
# маленькие значения $\beta$ и несколько значений $\gamma$.
# Порог определяется числом $R_0 = \beta N / \gamma$, где $N = 1000$:
# эпидемия развивается только при $R_0 > 1$.

# ## Активация проекта и загрузка пакетов

using DrWatson
@quickactivate "project"
@isdefined(SIRPetri) || include(srcdir("SIRPetri.jl"))
using .SIRPetri
using DataFrames, CSV, Plots

script_name = "sirpetri_params"
mkpath(plotsdir())
mkpath(datadir())

# ## Набор параметров
#
# Параметры собраны в словарь. Функция `dict_list` строит все
# комбинации: 6 значений $\beta$ и 3 значения $\gamma$ --- всего 18 запусков.

param_dict = Dict(
    :β => [0.00005, 0.0001, 0.0002, 0.0003, 0.0005, 0.001],
    :γ => [0.05, 0.1, 0.2],
    :tmax => 300.0,
)

params_list = dict_list(param_dict)
length(params_list)

# ## Функция одного эксперимента
#
# Функция запускает детерминированную симуляцию и возвращает
# число $R_0$, пик числа инфицированных, время пика
# и конечное число выздоровевших.

function run_experiment(p)
    β = p[:β]
    γ = p[:γ]
    net, u0, _ = build_sir_network(β, γ)
    df = simulate_deterministic(net, u0, (0.0, p[:tmax]), saveat = 0.1, rates = [β, γ])
    k = argmax(df.I)
    return (
        β = β,
        γ = γ,
        R0 = round(β * sum(u0) / γ, digits = 2),
        peak_I = round(df.I[k], digits = 1),
        t_peak = df.time[k],
        final_R = round(df.R[end], digits = 1),
    )
end

# ## Запуск экспериментов

results = []
for p in params_list
    push!(results, run_experiment(p))
end

df_params = DataFrame(results)
sort!(df_params, [:γ, :β])
CSV.write(datadir("sir_params.csv"), df_params)
println(df_params)

# ## Графики
#
# Слева --- пик числа инфицированных, справа --- конечное число
# выздоровевших в зависимости от $\beta$ для каждого значения $\gamma$.

p1 = plot(xlabel = "β", ylabel = "Peak I", title = "Пик эпидемии", legend = :topleft)
p2 = plot(xlabel = "β", ylabel = "Final R", title = "Переболело", legend = :bottomright)
for γ in sort(unique(df_params.γ))
    sub = df_params[df_params.γ .== γ, :]
    plot!(p1, sub.β, sub.peak_I, marker = :circle, linewidth = 2, label = "γ = $γ")
    plot!(p2, sub.β, sub.final_R, marker = :circle, linewidth = 2, label = "γ = $γ")
end
plt = plot(p1, p2, layout = (1, 2), size = (900, 400), left_margin = 5Plots.mm, bottom_margin = 5Plots.mm)

# ## Сохранение графика

savefig(plt, plotsdir("sir_params.png"))
println("Результаты сохранены в data/sir_params.csv и plots/sir_params.png")
