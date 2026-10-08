module MMcModel

using ResumableFunctions
using ConcurrentSim
using Distributions
using StableRNGs
using Random
using DataFrames

export simulate_mmc, mmc_analytic, system_size

"""
    customer(env, server, id, t_a, d_s, rng, rec, verbose)

Поведение одной заявки: приход, ожидание свободного канала,
обслуживание и уход. Моменты времени записываются в матрицу `rec`.
"""
@resumable function customer(
    env::Environment,
    server::Resource,
    id::Integer,
    t_a::Float64,
    d_s::Distribution,
    rng::AbstractRNG,
    rec::Matrix{Float64},
    verbose::Bool,
)
    @yield timeout(env, t_a) # customer arrives
    rec[id, 1] = now(env)
    if verbose
        println("Customer $id arrived: ", now(env))
    end
    @yield request(server) # customer starts service
    rec[id, 2] = now(env)
    if verbose
        println("Customer $id entered service: ", now(env))
    end
    @yield timeout(env, rand(rng, d_s)) # server is busy
    @yield unlock(server) # customer exits service
    rec[id, 3] = now(env)
    if verbose
        println("Customer $id exited service: ", now(env))
    end
end

"""
    simulate_mmc(; lam, mu, c, n, seed=123, verbose=false)

Моделирует систему M/M/c: `lam` --- интенсивность входящего потока,
`mu` --- интенсивность обслуживания одним каналом, `c` --- число каналов,
`n` --- число заявок. Возвращает DataFrame с моментами прихода,
начала и конца обслуживания каждой заявки.
"""
function simulate_mmc(; lam, mu, c, n, seed = 123, verbose = false)
    rng = StableRNG(seed)
    arrival_dist = Exponential(1 / lam) # interarrival time distribution
    service_dist = Exponential(1 / mu) # service time distribution
    rec = zeros(n, 3)

    sim = Simulation() # initialize simulation environment
    server = Resource(sim, c) # initialize servers
    arrival_time = 0.0
    for i = 1:n # initialize customers
        arrival_time += rand(rng, arrival_dist)
        @process customer(sim, server, i, arrival_time, service_dist, rng, rec, verbose)
    end
    run(sim) # run simulation

    df = DataFrame(id = 1:n, arrival = rec[:, 1], start = rec[:, 2], finish = rec[:, 3])
    df.wait = df.start .- df.arrival
    df.service = df.finish .- df.start
    return df
end

"""
    mmc_analytic(lam, mu, c)

Аналитические характеристики системы M/M/c в стационарном режиме.
"""
function mmc_analytic(lam, mu, c)
    rho = lam / (c * mu)
    a = c * rho
    P0 = 1 / (sum(a^n / factorial(n) for n = 0:(c-1)) + a^c / (factorial(c) * (1 - rho)))
    Pwait = a^c / (factorial(c) * (1 - rho)) * P0
    Lq = rho / (1 - rho) * Pwait
    Wq = Lq / lam
    W = Wq + 1 / mu
    L = lam * W
    return (rho = rho, P0 = P0, Pwait = Pwait, Lq = Lq, Wq = Wq, W = W, L = L)
end

"""
    system_size(df)

По таблице заявок строит число заявок в системе как функцию времени.
"""
function system_size(df)
    ev = vcat([(t, 1) for t in df.arrival], [(t, -1) for t in df.finish])
    sort!(ev, by = first)
    times = [0.0]
    counts = [0]
    for (t, d) in ev
        push!(times, t)
        push!(counts, counts[end] + d)
    end
    return (times = times, counts = counts)
end

end # module
