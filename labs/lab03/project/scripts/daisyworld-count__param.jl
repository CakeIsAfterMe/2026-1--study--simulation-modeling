# # Daisyworld: динамика числа маргариток для набора параметров
#
# Построим график числа маргариток для нескольких наборов параметров.

# ## Активация проекта и загрузка пакетов

using DrWatson
@quickactivate "project"
using Agents
using DataFrames
using CairoMakie

include(srcdir("daisyworld.jl"))

# ## Сбор данных

black(a) = a.breed == :black
white(a) = a.breed == :white
adata = [(black, count), (white, count)]

# ## Параметры эксперимента
#
# Меняются максимальный возраст (25 и 40) и начальная доля белых
# маргариток (0.2 и 0.8).

param_dict = Dict(
    :griddims => (30, 30),
    :max_age => [25, 40],
    :init_white => [0.2, 0.8],
    :init_black => 0.2,
    :albedo_white => 0.75,
    :albedo_black => 0.25,
    :surface_albedo => 0.4,
    :solar_change => 0.005,
    :solar_luminosity => 1.0,
    :scenario => :default,
    :seed => 165,
)

# Список всех комбинаций:

params_list = dict_list(param_dict)
length(params_list)

# ## Запуск для всех наборов
#
# Для каждого набора выполняем 1000 шагов, строим график и сохраняем его.
# В конце печатаем число маргариток на последнем шаге.

for params in params_list

    model = daisyworld(; params...)

    agent_df, model_df = run!(model, 1000; adata)
    figure = Figure(size = (600, 400));
    ax = figure[1, 1] = Axis(figure, xlabel = "tick", ylabel = "daisy count")
    blackl = lines!(ax, agent_df[!, :time], agent_df[!, :count_black], color = :black)
    whitel = lines!(ax, agent_df[!, :time], agent_df[!, :count_white], color = :orange)
    Legend(figure[1, 2], [blackl, whitel], ["black", "white"], labelsize = 12)

    plt_name = savename("daisy-count", params) * ".png"

    save(plotsdir(plt_name), figure)

    println("max_age = ", params[:max_age], ", init_white = ", params[:init_white],
        ": black = ", agent_df[end, :count_black], ", white = ", agent_df[end, :count_white])

end
