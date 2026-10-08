module SIRDES

using ResumableFunctions, ConcurrentSim, Distributions, DataFrames, Random
using DifferentialEquations

export SIRPerson, SIRModel, MakeSIRModel, activate, sir_run, out
export sir_metrics, sir_ode_solution
export XPerson, XModel, MakeXModel, activate_x, run_x, out_x, sir_ode_extended, time_average

# Вспомогательные функции для обновления массивов состояния
function increment!(a::Array{Int64})
    push!(a, a[length(a)] + 1)
end
function decrement!(a::Array{Int64})
    push!(a, a[length(a)] - 1)
end
function carryover!(a::Array{Int64})
    push!(a, a[length(a)])
end

# Структуры данных
mutable struct SIRPerson
    id::Int64
    status::Symbol   # :S, :I, :R
end

mutable struct SIRModel
    sim::ConcurrentSim.Simulation    # Тип Simulation, не Environment
    β::Float64
    c::Float64
    γ::Float64
    fixed::Bool                      # true --- болезнь длится ровно 1/γ
    ta::Array{Float64}
    Sa::Array{Int64}
    Ia::Array{Int64}
    Ra::Array{Int64}
    allIndividuals::Array{SIRPerson}
end

# Функции обновления статистики при событиях
function infection_update!(sim::ConcurrentSim.Simulation, m::SIRModel)
    push!(m.ta, ConcurrentSim.now(sim))
    decrement!(m.Sa)
    increment!(m.Ia)
    carryover!(m.Ra)
end

function recovery_update!(sim::ConcurrentSim.Simulation, m::SIRModel)
    push!(m.ta, ConcurrentSim.now(sim))
    carryover!(m.Sa)
    decrement!(m.Ia)
    increment!(m.Ra)
end

# Основная логика жизни индивида
@resumable function live(env::ConcurrentSim.Simulation, individual::SIRPerson, m::SIRModel)
    while individual.status == :S
        @yield timeout(env, rand(Exponential(1 / m.c)))
        alter = individual
        while alter == individual
            N = length(m.allIndividuals)
            index = rand(DiscreteUniform(1, N))
            alter = m.allIndividuals[index]
        end
        if alter.status == :I
            if rand(Uniform(0, 1)) < m.β
                individual.status = :I
                infection_update!(env, m)
            end
        end
    end
    if individual.status == :I
        if m.fixed
            @yield timeout(env, 1 / m.γ)
        else
            @yield timeout(env, rand(Exponential(1 / m.γ)))
        end
        individual.status = :R
        recovery_update!(env, m)
    end
end

# Функции создания и запуска модели
function MakeSIRModel(u0, p; fixed = false)
    (S, I, R) = u0
    N = S + I + R
    (β, c, γ) = p
    sim = ConcurrentSim.Simulation()   # Создаём именно Simulation
    allIndividuals = SIRPerson[]
    for i = 1:S
        push!(allIndividuals, SIRPerson(i, :S))
    end
    for i = (S+1):(S+I)
        push!(allIndividuals, SIRPerson(i, :I))
    end
    for i = (S+I+1):N
        push!(allIndividuals, SIRPerson(i, :R))
    end
    ta = Float64[0.0]
    Sa = Int64[S]
    Ia = Int64[I]
    Ra = Int64[R]
    SIRModel(sim, β, c, γ, fixed, ta, Sa, Ia, Ra, allIndividuals)
end

function activate(m::SIRModel)
    [@process live(m.sim, individual, m) for individual in m.allIndividuals]
end

function sir_run(m::SIRModel, tf::Float64)
    ConcurrentSim.run(m.sim, tf)
end

function out(m::SIRModel)
    result = DataFrame()
    result[!, :t] = m.ta
    result[!, :S] = m.Sa
    result[!, :I] = m.Ia
    result[!, :R] = m.Ra
    return result
end

# Ключевые показатели прогона: пик I, время пика, итоговое R
function sir_metrics(df)
    k = argmax(df.I)
    return (peak_I = df.I[k], t_peak = df.t[k], final_R = df.R[end])
end

# Детерминированная модель SIR для сравнения:
# скорость заражения β c S I / N, скорость выздоровления γ I
function sir_ode_solution(u0, p, tmax; saveat = 0.1)
    N = sum(u0)
    function f!(du, u, q, t)
        β, c, γ = q
        S, I, R = u
        du[1] = -β * c * S * I / N
        du[2] = β * c * S * I / N - γ * I
        du[3] = γ * I
    end
    prob = ODEProblem(f!, Float64.(u0), (0.0, tmax), p)
    sol = solve(prob, Tsit5(), saveat = saveat)
    return DataFrame(t = sol.t, S = sol[1, :], I = sol[2, :], R = sol[3, :])
end

# ======================================================================
# Расширенная модель для дополнительных заданий:
# демография (смерть и рождение), вакцинация и латентный период (SEIR)
# ======================================================================

mutable struct XPerson
    id::Int64
    status::Symbol      # :S, :E, :I, :R, :D (умер)
    change_at::Float64  # момент перехода E -> I или I -> R
end

mutable struct XModel
    sim::ConcurrentSim.Simulation
    β::Float64
    c::Float64
    γ::Float64
    σ::Float64          # интенсивность перехода E -> I
    μ::Float64          # интенсивность смерти (и рождения)
    N0::Int64           # начальный размер популяции
    latent::Bool        # true --- модель SEIR
    ta::Array{Float64}
    Sa::Array{Int64}
    Ea::Array{Int64}
    Ia::Array{Int64}
    Ra::Array{Int64}
    people::Array{XPerson}
    infections::Int64   # сколько всего было заражений
    vaccinated::Int64   # сколько человек привито
end

# Запись изменения численности групп в текущий момент
function record!(env, m::XModel, dS, dE, dI, dR)
    push!(m.ta, ConcurrentSim.now(env))
    push!(m.Sa, m.Sa[end] + dS)
    push!(m.Ea, m.Ea[end] + dE)
    push!(m.Ia, m.Ia[end] + dI)
    push!(m.Ra, m.Ra[end] + dR)
end

# Случайный живой собеседник, отличный от самого человека
function pick_alter(m::XModel, ind::XPerson)
    alter = ind
    while alter === ind || alter.status == :D
        alter = m.people[rand(1:length(m.people))]
    end
    return alter
end

# Контакт здорового человека: возможно заражение
function contact!(env, m::XModel, ind::XPerson)
    alter = pick_alter(m, ind)
    if alter.status == :I && rand() < m.β
        m.infections += 1
        if m.latent
            ind.status = :E
            ind.change_at = ConcurrentSim.now(env) + rand(Exponential(1 / m.σ))
            record!(env, m, -1, 1, 0, 0)
        else
            ind.status = :I
            ind.change_at = ConcurrentSim.now(env) + rand(Exponential(1 / m.γ))
            record!(env, m, -1, 0, 1, 0)
        end
    end
end

# Переход E -> I или I -> R
function progress!(env, m::XModel, ind::XPerson)
    if ind.status == :E
        ind.status = :I
        ind.change_at = ConcurrentSim.now(env) + rand(Exponential(1 / m.γ))
        record!(env, m, 0, -1, 1, 0)
    elseif ind.status == :I
        ind.status = :R
        ind.change_at = Inf
        record!(env, m, 0, 0, -1, 1)
    end
end

# Смерть: человек убирается из своей группы
function die!(env, m::XModel, ind::XPerson)
    record!(
        env,
        m,
        ind.status == :S ? -1 : 0,
        ind.status == :E ? -1 : 0,
        ind.status == :I ? -1 : 0,
        ind.status == :R ? -1 : 0,
    )
    ind.status = :D
end

# Момент следующего собственного события человека (без учёта смерти)
function next_time(env, m::XModel, ind::XPerson)
    if ind.status == :S
        return ConcurrentSim.now(env) + rand(Exponential(1 / m.c))
    elseif ind.status == :E || ind.status == :I
        return ind.change_at
    else
        return Inf
    end
end

# Жизнь человека в расширенной модели
@resumable function live_x(env::ConcurrentSim.Simulation, ind::XPerson, m::XModel)
    death_at = Inf
    if m.μ > 0
        death_at = ConcurrentSim.now(env) + rand(Exponential(1 / m.μ))
    end
    done = false
    while !done
        t_next = next_time(env, m, ind)
        if death_at <= t_next
            if isinf(death_at)
                done = true
            else
                @yield timeout(env, death_at - ConcurrentSim.now(env))
                die!(env, m, ind)
                done = true
            end
        else
            @yield timeout(env, t_next - ConcurrentSim.now(env))
            if ind.status == :S
                contact!(env, m, ind)
            else
                progress!(env, m, ind)
            end
        end
    end
end

# Рождение нового здорового человека
function new_person!(env, m::XModel)
    ind = XPerson(length(m.people) + 1, :S, Inf)
    push!(m.people, ind)
    record!(env, m, 1, 0, 0, 0)
    return ind
end

# Процесс рождений с постоянной интенсивностью μ N0
@resumable function births(env::ConcurrentSim.Simulation, m::XModel)
    while true
        @yield timeout(env, rand(Exponential(1 / (m.μ * m.N0))))
        ind = new_person!(env, m)
        @process live_x(env, ind, m)
    end
end

# Мгновенная вакцинация доли fraction здоровых: они переходят в R
function vaccinate_now!(env, m::XModel, fraction)
    sus = filter(x -> x.status == :S, m.people)
    k = round(Int, fraction * length(sus))
    for x in shuffle(sus)[1:k]
        x.status = :R
        record!(env, m, -1, 0, 0, 1)
    end
    m.vaccinated = k
end

# Событие вакцинации по таймеру
@resumable function vaccinate(env::ConcurrentSim.Simulation, m::XModel, time::Float64, fraction::Float64)
    @yield timeout(env, time)
    vaccinate_now!(env, m, fraction)
end

function MakeXModel(u0, p; latent = false, σ = 0.5, μ = 0.0)
    (S, I, R) = u0
    (β, c, γ) = p
    sim = ConcurrentSim.Simulation()
    people = XPerson[]
    for i = 1:S
        push!(people, XPerson(i, :S, Inf))
    end
    for i = 1:I
        push!(people, XPerson(S + i, :I, rand(Exponential(1 / γ))))
    end
    for i = 1:R
        push!(people, XPerson(S + I + i, :R, Inf))
    end
    XModel(sim, β, c, γ, σ, μ, S + I + R, latent,
        Float64[0.0], Int64[S], Int64[0], Int64[I], Int64[R], people, 0, 0)
end

function activate_x(m::XModel; vaccination_time = nothing, fraction = 0.0)
    for ind in m.people
        @process live_x(m.sim, ind, m)
    end
    if m.μ > 0
        @process births(m.sim, m)
    end
    if vaccination_time !== nothing
        @process vaccinate(m.sim, m, Float64(vaccination_time), Float64(fraction))
    end
end

function run_x(m::XModel, tf::Float64)
    ConcurrentSim.run(m.sim, tf)
end

function out_x(m::XModel)
    df = DataFrame(t = m.ta, S = m.Sa, E = m.Ea, I = m.Ia, R = m.Ra)
    df.N = df.S .+ df.E .+ df.I .+ df.R
    return df
end

# Среднее по времени значение ступенчатого ряда x(t) на отрезке [t0, t1]
function time_average(t, x, t0, t1)
    total = 0.0
    for k = 1:(length(t)-1)
        a = max(t[k], t0)
        b = min(t[k+1], t1)
        if b > a
            total += x[k] * (b - a)
        end
    end
    if t[end] < t1
        total += x[end] * (t1 - max(t[end], t0))
    end
    return total / (t1 - t0)
end

# Детерминированная модель с латентным периодом и демографией
function sir_ode_extended(u0, p, tmax; σ = Inf, μ = 0.0, saveat = 0.1)
    N0 = sum(u0)
    latent = isfinite(σ)
    function f!(du, u, q, t)
        β, c, γ = q
        S, E, I, R = u
        N = S + E + I + R
        inf = β * c * S * I / N
        du[1] = μ * N0 - inf - μ * S
        du[2] = latent ? inf - σ * E - μ * E : 0.0
        du[3] = (latent ? σ * E : inf) - γ * I - μ * I
        du[4] = γ * I - μ * R
    end
    prob = ODEProblem(f!, Float64[u0[1], 0.0, u0[2], u0[3]], (0.0, tmax), p)
    sol = solve(prob, Tsit5(), saveat = saveat)
    return DataFrame(t = sol.t, S = sol[1, :], E = sol[2, :], I = sol[3, :], R = sol[4, :])
end

end # module
