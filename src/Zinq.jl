module Zinq

include("QuantumDynamics.jl")

using .QuantumDynamics, TimerOutputs

export julia_main

function get_timer_table()
    return join(split(sprint(print_timer, TimerOutputs.DEFAULT_TIMER), '\n')[5:end - 1] |> x -> (x[2] = uppercase(x[2]); x), '\n')
end

function print_header()
    println("ZINQ")
end

function julia_main()::Cint
    print_header()

    reset_timer!()

    @timeit "QUANTUM DYNAMICS" run_qd()

    println(get_timer_table())

    return 0
end

end # module Zinq
