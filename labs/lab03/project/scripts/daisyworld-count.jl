# # Daisyworld: динамика числа маргариток
#
# Построим график изменения числа чёрных и белых маргариток в
# зависимости от модельного времени.

# ## Активация проекта и загрузка пакетов

using DrWatson
@quickactivate "project"
using Agents
using DataFrames
using CairoMakie

include(srcdir("daisyworld.jl"))

# ## Сбор данных
#
# Функции `black` и `white` проверяют вид маргаритки. На каждом шаге
# считаем, сколько агентов каждого вида.

black(a) = a.breed == :black
white(a) = a.breed == :white
adata = [(black, count), (white, count)]

# ## Запуск модели
#
# Солнечная постоянная равна 1.0 и не меняется. Выполняем 1000 шагов.

model = daisyworld(; solar_luminosity = 1.0)

agent_df, model_df = run!(model, 1000; adata)
first(agent_df, 5)

# ## Построение графика

figure = Figure(size = (600, 400));
ax = figure[1, 1] = Axis(figure, xlabel = "tick", ylabel = "daisy count")
blackl = lines!(ax, agent_df[!, :time], agent_df[!, :count_black], color = :black)
whitel = lines!(ax, agent_df[!, :time], agent_df[!, :count_white], color = :orange)
Legend(figure[1, 2], [blackl, whitel], ["black", "white"], labelsize = 12)
figure

# ## Сохранение графика

save(plotsdir("daisy_count.png"), figure)
