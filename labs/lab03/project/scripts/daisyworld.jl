# # Daisyworld: базовая визуализация
#
# Модель «мир маргариток». На клеточной сетке растут чёрные и белые
# маргаритки. Чёрные поглощают свет и нагревают поверхность, белые
# отражают свет и охлаждают её. Маргаритки размножаются только при
# подходящей температуре.
#
# Построим тепловую карту температуры поверхности и покажем на ней
# маргаритки в начальный момент, через 5 шагов и ещё через 40 шагов.

# ## Активация проекта и загрузка пакетов

using DrWatson
@quickactivate "project"
using Agents
using DataFrames
using CairoMakie

# Код модели находится в файле `src/daisyworld.jl`:

include(srcdir("daisyworld.jl"))

# ## Создание модели
#
# Модель создаётся с параметрами по умолчанию: сетка 30 на 30, по 20%
# клеток занято белыми и чёрными маргаритками.

model = daisyworld()

# ## Настройка отображения
#
# Цвет маргаритки совпадает с её видом. Фоном служит температура клетки.

daisycolor(a::Daisy) = a.breed

plotkwargs = (
    agent_color = daisycolor, agent_size = 20, agent_marker = '✿',
    heatarray = :temperature,
    heatkwargs = (colorrange = (-20, 60),),
)

# ## Начальное состояние

plt1, _ = abmplot(model; plotkwargs...)
plt1

# ## Состояние через 5 шагов

step!(model, 5)
plt2, _ = abmplot(model; plotkwargs...)
plt2

# ## Состояние ещё через 40 шагов

step!(model, 40)
plt3, _ = abmplot(model; plotkwargs...)
plt3

# ## Сохранение графиков

save(plotsdir("daisy_step001.png"), plt1)
save(plotsdir("daisy_step005.png"), plt2)
save(plotsdir("daisy_step040.png"), plt3)
