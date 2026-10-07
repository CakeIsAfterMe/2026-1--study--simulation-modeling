# # Модель SIR в виде сети Петри: анимация
#
# Скрипт строит GIF-анимацию: на каждом кадре столбчатая диаграмма
# показывает, сколько людей находится в группах $S$, $I$ и $R$
# в данный момент времени.

# ## Активация проекта и загрузка пакетов

using DrWatson
@quickactivate "project"
@isdefined(SIRPetri) || include(srcdir("SIRPetri.jl"))
using .SIRPetri
using DataFrames, Plots

script_name = "sirpetri_animate"
mkpath(plotsdir())

# ## Параметры и симуляция
#
# Решение сохраняется с шагом 0.2, чтобы анимация была плавной.

β = 0.3
γ = 0.1
tmax = 100.0

net, u0, states = build_sir_network(β, γ)
df = simulate_deterministic(net, u0, (0.0, tmax), saveat = 0.2, rates = [β, γ])
nrow(df)

# ## Построение анимации
#
# В анимацию берётся каждая пятая точка решения, то есть один кадр
# на единицу времени.

frames = 1:5:nrow(df)
anim = @animate for k in frames
    bar(
        ["S", "I", "R"],
        [df.S[k], df.I[k], df.R[k]],
        ylim = (0, 1000),
        legend = false,
        xlabel = "Compartment",
        ylabel = "Population",
        title = "SIR dynamics at t = $(round(df.time[k], digits = 1))",
    )
end
println("Число кадров: ", length(frames))

# ## Сохранение анимации

gif(anim, plotsdir("sir_animation.gif"), fps = 10)
println("Анимация сохранена в plots/sir_animation.gif")
