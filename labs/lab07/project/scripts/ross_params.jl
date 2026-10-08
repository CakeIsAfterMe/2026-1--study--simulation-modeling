# # Модель Росса: набор параметров
#
# Время до падения системы случайное и меняется от прогона к прогону
# очень сильно, поэтому пяти прогонов для сравнения с теорией мало.
# Здесь моделирование повторяется по 50 раз для разного числа машин,
# разного резерва и одного или двух ремонтников.

# ## Активация проекта и загрузка пакетов

using DrWatson
@quickactivate "project"
@isdefined(RossModel) || include(srcdir("RossModel.jl"))
using .RossModel
using DataFrames, CSV, Plots, Statistics
using StableRNGs

script_name = "ross_params"
mkpath(plotsdir())
mkpath(datadir())

# ## Набор параметров
#
# Функция `dict_list` строит все комбинации: 3 значения числа машин,
# 3 значения резерва и 2 значения числа ремонтников --- всего 18 вариантов.

ross_dict = Dict(
    :N => [10, 15, 20],
    :S => [1, 2, 3],
    :R => [1, 2],
    :runs => 50,
    :seed => 42,
)

ross_list = dict_list(ross_dict)
length(ross_list)

# ## Функция одного эксперимента
#
# Функция выполняет заданное число прогонов и возвращает среднее время
# до падения из модели и из теории, а также среднюю загрузку
# ремонтников и среднюю длину очереди на ремонт.

function run_ross(p)
    rng = StableRNG(p[:seed])
    times = Float64[]
    utils = Float64[]
    queues = Float64[]
    for i = 1:p[:runs]
        T, events = sim_repair(rng; N = p[:N], S = p[:S], R = p[:R])
        mon = ross_monitor(events, T; R = p[:R])
        push!(times, T)
        push!(utils, mon.utilization)
        push!(queues, mon.mean_queue)
    end
    return (
        N = p[:N],
        S = p[:S],
        R = p[:R],
        T_sim = round(mean(times), digits = 1),
        T_theory = round(ross_analytic(N = p[:N], S = p[:S], R = p[:R]), digits = 1),
        utilization = round(mean(utils), digits = 4),
        mean_queue = round(mean(queues), digits = 4),
    )
end

# ## Запуск экспериментов

df_ross = DataFrame([run_ross(p) for p in ross_list])
sort!(df_ross, [:R, :N, :S])
CSV.write(datadir("ross_params.csv"), df_ross)
println(df_ross)

# ## Графики
#
# Точки --- результат моделирования, пунктир --- теория.
# Слева один ремонтник, справа два. Ось времени логарифмическая.

p1 = plot(xlabel = "N", ylabel = "E[T]", title = "Один ремонтник", yscale = :log10, legend = :topright)
p2 = plot(xlabel = "N", ylabel = "E[T]", title = "Два ремонтника", yscale = :log10, legend = :topright)
for s in sort(unique(df_ross.S))
    one = df_ross[(df_ross.S .== s) .& (df_ross.R .== 1), :]
    two = df_ross[(df_ross.S .== s) .& (df_ross.R .== 2), :]
    scatter!(p1, one.N, one.T_sim, label = "S = $s, модель")
    plot!(p1, one.N, one.T_theory, linestyle = :dash, label = "S = $s, теория")
    scatter!(p2, two.N, two.T_sim, label = "S = $s, модель")
    plot!(p2, two.N, two.T_theory, linestyle = :dash, label = "S = $s, теория")
end
plt_params = plot(p1, p2, layout = (1, 2), size = (1000, 420), left_margin = 5Plots.mm, bottom_margin = 5Plots.mm)

# ## Сохранение графика

savefig(plt_params, plotsdir("ross_params.png"))
println("Результаты сохранены в data/ross_params.csv и plots/ross_params.png")
