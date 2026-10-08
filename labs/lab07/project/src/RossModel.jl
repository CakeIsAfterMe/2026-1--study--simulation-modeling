module RossModel

using ResumableFunctions
using ConcurrentSim
using Distributions
using StableRNGs
using Random

export sim_repair, ross_analytic, ross_monitor, ross_timeline

"""
    machine(env, repair_facility, spares, rng, F, G, events)

Жизненный цикл одной машины: ожидание в резерве, работа до отказа,
замена резервной машиной, ремонт и возврат в резерв.
В вектор `events` записываются моменты отказа, начала и конца ремонта.
"""
@resumable function machine(
    env::Environment,
    repair_facility::Resource,
    spares::Store{Process},
    rng::AbstractRNG,
    F::Distribution,
    G::Distribution,
    events::Vector{Tuple{Float64,Symbol}},
)
    while true
        try
            @yield timeout(env, Inf)
        catch
        end
        @yield timeout(env, rand(rng, F))
        push!(events, (now(env), :fail))
        get_spare = take!(spares)
        @yield get_spare | timeout(env)
        if state(get_spare) != ConcurrentSim.idle
            @yield interrupt(value(get_spare))
        else
            throw(StopSimulation("No more spares!"))
        end
        @yield request(repair_facility)
        push!(events, (now(env), :start))
        @yield timeout(env, rand(rng, G))
        push!(events, (now(env), :finish))
        @yield unlock(repair_facility)
        @yield put!(spares, active_process(env))
    end
end

"""
    start_sim(env, repair_facility, spares, n, s, rng, F, G, events)

Запускает `n` работающих машин и кладёт `s` машин в резерв.
"""
@resumable function start_sim(
    env::Environment,
    repair_facility::Resource,
    spares::Store{Process},
    n::Int,
    s::Int,
    rng::AbstractRNG,
    F::Distribution,
    G::Distribution,
    events::Vector{Tuple{Float64,Symbol}},
)
    for i = 1:n
        proc = @process machine(env, repair_facility, spares, rng, F, G, events)
        @yield interrupt(proc)
    end
    for i = 1:s
        proc = @process machine(env, repair_facility, spares, rng, F, G, events)
        @yield put!(spares, proc)
    end
end

"""
    sim_repair(rng; N=10, S=3, R=1, lambda=100, mu=1, verbose=false)

Один прогон модели Росса: `N` работающих машин, `S` резервных,
`R` ремонтников, `lambda` --- среднее время работы до отказа,
`mu` --- среднее время ремонта.
Возвращает время падения системы и журнал событий.
"""
function sim_repair(rng::AbstractRNG; N = 10, S = 3, R = 1, lambda = 100, mu = 1, verbose = false)
    F = Exponential(lambda)
    G = Exponential(mu)
    events = Tuple{Float64,Symbol}[]
    sim = Simulation()
    repair_facility = Resource(sim, R)
    spares = Store{Process}(sim)
    @process start_sim(sim, repair_facility, spares, N, S, rng, F, G, events)
    msg = run(sim)
    stop_time = now(sim)
    if verbose
        println("At time $stop_time: $msg")
    end
    return stop_time, events
end

"""
    ross_analytic(; N=10, S=3, R=1, lambda=100, mu=1)

Среднее время до падения системы. Состояние --- число сломанных машин
b = 0..S. Для средних времён T_b до падения решается система
линейных уравнений.
"""
function ross_analytic(; N = 10, S = 3, R = 1, lambda = 100, mu = 1)
    fail = N / lambda
    A = zeros(S + 1, S + 1)
    for b = 0:S
        rep = min(b, R) / mu
        A[b+1, b+1] = fail + rep
        if b < S
            A[b+1, b+2] = -fail
        end
        if b > 0
            A[b+1, b] = -rep
        end
    end
    T = A \ ones(S + 1)
    return T[1]
end

"""
    ross_monitor(events, T; R=1)

По журналу событий считает загрузку ремонтников, среднюю длину
очереди на ремонт и среднее число сломанных машин за время `T`.
"""
function ross_monitor(events, T; R = 1)
    queue = 0
    busy = 0
    broken = 0
    area_queue = 0.0
    area_busy = 0.0
    area_broken = 0.0
    t_prev = 0.0
    for (t, kind) in events
        dt = t - t_prev
        area_queue += queue * dt
        area_busy += busy * dt
        area_broken += broken * dt
        if kind == :fail
            broken += 1
            queue += 1
        elseif kind == :start
            queue -= 1
            busy += 1
        else
            busy -= 1
            broken -= 1
        end
        t_prev = t
    end
    return (
        utilization = area_busy / (R * T),
        mean_queue = area_queue / T,
        mean_broken = area_broken / T,
    )
end

"""
    ross_timeline(events; N=10, S=3)

Число исправных машин как функция времени.
"""
function ross_timeline(events; N = 10, S = 3)
    times = [0.0]
    healthy = [N + S]
    for (t, kind) in events
        if kind == :fail
            push!(times, t)
            push!(healthy, healthy[end] - 1)
        elseif kind == :finish
            push!(times, t)
            push!(healthy, healthy[end] + 1)
        end
    end
    return (times = times, healthy = healthy)
end

end # module
