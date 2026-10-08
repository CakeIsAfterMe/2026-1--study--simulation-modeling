# # Дискретно-событийная модель SEIR
#
# Добавим латентный период: заразившийся человек сначала попадает
# в группу $E$ и ещё не заражает других, а через случайное время
# со средним $1/\sigma$ становится больным. Сравним эту модель с SIR.

# ## Активация проекта и загрузка пакетов

using DrWatson
@quickactivate "project"
@isdefined(SIRDES) || include(srcdir("sir_model.jl"))
using .SIRDES
using Random, Plots, DataFrames, CSV

script_name = "sir_des_seir"
mkpath(plotsdir())
mkpath(datadir())

# ## Параметры модели
#
# Параметры заражения те же, латентный период в среднем длится
# $1/\sigma = 2$.

tmax = 60.0
u0 = [990, 10, 0]
p = [0.05, 10.0, 0.25]
σ = 0.5

# ## Модель SIR

Random.seed!(1234)
m_sir = MakeXModel(u0, p)
activate_x(m_sir)
run_x(m_sir, tmax)
df_sir = out_x(m_sir)

# ## Модель SEIR

Random.seed!(1234)
m_seir = MakeXModel(u0, p; latent = true, σ = σ)
activate_x(m_seir)
run_x(m_seir, tmax)
df_seir = out_x(m_seir)

# ## Детерминированные модели

df_sir_ode = sir_ode_extended(u0, p, tmax)
df_seir_ode = sir_ode_extended(u0, p, tmax; σ = σ)

# ## Сравнение показателей

df_seir_cmp = DataFrame([
    merge((model = "SIR, DES",), map(float, sir_metrics(df_sir))),
    merge((model = "SEIR, DES",), map(float, sir_metrics(df_seir))),
    merge((model = "SIR, ОДУ",), map(x -> round(x, digits = 1), sir_metrics(df_sir_ode))),
    merge((model = "SEIR, ОДУ",), map(x -> round(x, digits = 1), sir_metrics(df_seir_ode))),
])
CSV.write(datadir("sir_des_seir.csv"), df_seir_cmp)
println(df_seir_cmp)

# ## Графики
#
# Слева число больных в двух моделях, справа группы $E$ и $I$
# в модели SEIR. Пунктир --- детерминированные модели.

p1 = plot(df_sir.t, df_sir.I, label = "SIR", xlabel = "Время", ylabel = "I", title = "Больные", linewidth = 2)
plot!(p1, df_seir.t, df_seir.I, label = "SEIR", linewidth = 2)
plot!(p1, df_sir_ode.t, df_sir_ode.I, label = "SIR, ОДУ", linestyle = :dash)
plot!(p1, df_seir_ode.t, df_seir_ode.I, label = "SEIR, ОДУ", linestyle = :dash)
p2 = plot(df_seir.t, [df_seir.E df_seir.I], label = ["E" "I"], xlabel = "Время", ylabel = "Численность", title = "Модель SEIR", linewidth = 2)
plt_seir = plot(p1, p2, layout = (1, 2), size = (1000, 420), left_margin = 5Plots.mm, bottom_margin = 5Plots.mm)
savefig(plt_seir, plotsdir("sir_des_seir.png"))
plt_seir
