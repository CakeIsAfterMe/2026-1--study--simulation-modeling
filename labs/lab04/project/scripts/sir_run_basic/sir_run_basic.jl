using DrWatson
@quickactivate "project"

using Agents, DataFrames, Plots
using JLD2

include(srcdir("sir_model.jl"))

params = Dict(
    :Ns => [1000, 1000, 1000],
    :β_und => [0.5, 0.5, 0.5],
    :β_det => [0.05, 0.05, 0.05],
    :infection_period => 14,
    :detection_time => 7,
    :death_rate => 0.02,
    :reinfection_probability => 0.1,
    :Is => [0, 0, 1],
    :seed => 42,
    :n_steps => 100,
)

γ = 1 / params[:infection_period]
R0 = params[:β_und][1] / γ
println("R0 = β/γ = ", R0)

model = initialize_sir(; params...)

times = Int[]
S_vals = Int[]
I_vals = Int[]
R_vals = Int[]
total_vals = Int[]

for step = 1:params[:n_steps]
    Agents.step!(model, 1)

    push!(times, step)
    push!(S_vals, susceptible_count(model))
    push!(I_vals, infected_count(model))
    push!(R_vals, recovered_count(model))
    push!(total_vals, total_count(model))
end

agent_df = DataFrame(time = times, susceptible = S_vals, infected = I_vals, recovered = R_vals)
model_df = DataFrame(time = times, total = total_vals)
first(agent_df, 5)

peak_idx = argmax(I_vals)
println("Пик инфицированных: ", I_vals[peak_idx], " на ", times[peak_idx], " день")
println("В конце: S = ", S_vals[end], ", I = ", I_vals[end], ", R = ", R_vals[end])
println("Умерло: ", sum(params[:Ns]) - total_vals[end])

plt = plot(
    agent_df.time,
    agent_df.susceptible,
    label = "Восприимчивые",
    xlabel = "Дни",
    ylabel = "Количество",
    linewidth = 2,
)
plot!(plt, agent_df.time, agent_df.infected, label = "Инфицированные", linewidth = 2)
plot!(plt, agent_df.time, agent_df.recovered, label = "Выздоровевшие", linewidth = 2)
plot!(plt, agent_df.time, model_df.total, label = "Всего живых", linestyle = :dash)
plt

savefig(plt, plotsdir("sir_basic_dynamics.png"))
@save datadir("sir_basic_agent.jld2") agent_df
@save datadir("sir_basic_model.jld2") model_df
