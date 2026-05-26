module Zinq

include("QuantumDynamics.jl")

using .QuantumDynamics, TOML, TimerOutputs

export julia_main

function get_timer_table()
    return join(split(sprint(print_timer, TimerOutputs.DEFAULT_TIMER), '\n')[5:end - 1] |> x -> (x[2] = uppercase(x[2]); x), '\n')
end

function print_header()
    println("ZINQ")
end

function julia_main()::Cint
    print_header()

    input_file = isempty(ARGS) ? "input.toml" : ARGS[1]

    if !isfile(input_file)
        println(stderr, "INPUT FILE '$input_file' NOT FOUND"); return 1
    end

    reset_timer!()

    @timeit "QUANTUM DYNAMICS" run_qd(TOML.parsefile(input_file))

    println(get_timer_table())

    return 0
end

end # module Zinq
