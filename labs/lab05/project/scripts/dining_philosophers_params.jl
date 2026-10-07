# # Обедающие философы: набор параметров
#
# Алгоритм Гиллеспи случайный: при каждом запуске получается своя
# траектория и своё время наступления блокировки. Поэтому повторим
# моделирование для набора параметров: разное число философов, две
# сети и по 20 зёрен генератора на каждый вариант.

# ## Активация проекта и загрузка пакетов

using DrWatson
@quickactivate "project"
@isdefined(DiningPhilosophers) || include(srcdir("DiningPhilosophers.jl"))
using .DiningPhilosophers
using DataFrames, CSV, Plots, Random
using Statistics

# ## Набор параметров
#
# Параметры собраны в словарь. Функция `dict_list` строит все
# комбинации: 3 значения $N$, 2 сети и 20 зёрен --- всего 120 запусков.

param_dict = Dict(
    :N => [3, 5, 7],
    :network => [:classic, :arbiter],
    :seed => collect(1:20),
    :tmax => 50.0,
)

params_list = dict_list(param_dict)
length(params_list)

# ## Функция одного эксперимента
#
# Функция строит сеть, запускает моделирование и возвращает:
# наступила ли блокировка, время последнего события и число приёмов
# пищи. Приём пищи считаем по моментам, когда число едящих уменьшается.

function run_experiment(p)
    N = p[:N]
    net, u0, _ = p[:network] == :classic ? build_classical_network(N) : build_arbiter_network(N)
    df = simulate_stochastic(net, u0, p[:tmax]; rng = Xoshiro(p[:seed]))
    eating = vec(sum(Matrix(df[:, ["Eat_$i" for i = 1:N]]), dims = 2))
    meals = count(diff(eating) .< 0)
    return (
        N = N,
        network = p[:network],
        seed = p[:seed],
        deadlock = detect_deadlock(df, net),
        t_end = df.time[end],
        meals = meals,
    )
end

# ## Запуск экспериментов

results = []
for p in params_list
    push!(results, run_experiment(p))
end

df_all = DataFrame(results)
CSV.write(datadir("dining_params_all.csv"), df_all)
first(df_all, 5)

# ## Усреднение по повторам
#
# Для каждого сочетания $N$ и сети считаем долю запусков с блокировкой,
# среднее время последнего события и среднее число приёмов пищи.

grouped = combine(
    groupby(df_all, [:network, :N]),
    :deadlock => mean => :deadlock_share,
    :t_end => mean => :mean_t_end,
    :meals => mean => :mean_meals,
)
sort!(grouped, [:network, :N])
println(grouped)

# ## Графики
#
# Слева среднее время до блокировки в классической сети, справа среднее
# число приёмов пищи в обеих сетях.

classic = filter(row -> row.network == :classic, grouped)
arbiter = filter(row -> row.network == :arbiter, grouped)

p1 = plot(
    classic.N,
    classic.mean_t_end,
    marker = :circle,
    linewidth = 2,
    xlabel = "Число философов N",
    ylabel = "Время до блокировки",
    label = "Классическая сеть",
    title = "Время до deadlock",
)
p2 = plot(
    classic.N,
    classic.mean_meals,
    marker = :circle,
    linewidth = 2,
    xlabel = "Число философов N",
    ylabel = "Приёмов пищи",
    label = "Классическая сеть",
    title = "Число приёмов пищи",
)
plot!(p2, arbiter.N, arbiter.mean_meals, marker = :square, linewidth = 2, label = "Сеть с арбитром")
plt = plot(p1, p2, layout = (1, 2), size = (900, 400))

# ## Сохранение графика

savefig(plt, plotsdir("params_scan.png"))
println("Результаты сохранены в data/dining_params_all.csv и plots/params_scan.png")
