import copy, json, numpy, os

input_1d = { "zinq" : [{
    "name" : "quantum_dynamics",
    "options" : {
        "grid" : {
            "limits" : [[0.1, 16]]
        },
        "initial_conditions" : {
            "gamma" : [2]
        },
        "log_intervals" : {
            "iteration" : 10
        },
        "potential" : {
            "file" : {
                "ndim" : 1
            }
        },
        "write" : {}
    }
}]}

input_2d = { "zinq" : [{
    "name" : "quantum_dynamics",
    "options" : {
        "grid" : {
            "limits" : [[-10, 16], [-8, 8]]
        },
        "initial_conditions" : {
            "gamma" : [2, 2]
        },
        "log_intervals" : {
            "iteration" : 10
        },
        "potential" : {
            "file" : {
                "ndim" : 2
            }
        },
        "write" : {}
    }
}]}

os.makedirs("result", exist_ok=True)

mass, x0 = 5674.35, 10

for e in list(numpy.arange(0.1, 1.0, 0.1)) + list(numpy.arange(1, 10, 1)) + list(numpy.arange(10, 100, 10)) + [100]:

    if e == 0: continue

    inp_1d = copy.deepcopy(input_1d)

    inp_1d["zinq"][0]["options"]["initial_conditions"]["momentum"] = [-numpy.sqrt(2 * mass * e)]
    inp_1d["zinq"][0]["options"]["initial_conditions"]["position"] = [x0]
    inp_1d["zinq"][0]["options"]["initial_conditions"]["mass"] = mass
    inp_1d["zinq"][0]["options"]["initial_conditions"]["state"] = 0

    inp_1d["zinq"][0]["options"]["grid"]["points"] = 1024

    inp_1d["zinq"][0]["options"]["iterations"] = 250
    inp_1d["zinq"][0]["options"]["time_step"] = 5

    inp_1d_2s = copy.deepcopy(inp_1d)

    inp_1d["zinq"][0]["options"]["write"]["population"] = f"POPULATION_EXACT_1D_E={e:.1f}.mat"
    inp_1d_2s["zinq"][0]["options"]["write"]["population"] = f"POPULATION_EXACT_1D_2S_E={e:.1f}.mat"

    inp_1d["zinq"][0]["options"]["potential"]["file"]["nstate"] = 3
    inp_1d_2s["zinq"][0]["options"]["potential"]["file"]["nstate"] = 2

    inp_1d["zinq"][0]["options"]["potential"]["file"]["path"] = "../data/U_1D.dat"
    inp_1d_2s["zinq"][0]["options"]["potential"]["file"]["path"] = "../data/U_1D_2S.mat"

    with open(f"result/input_1d_e={e:.1f}.json", "w") as input_1d_file:
        json.dump(inp_1d, input_1d_file, indent=4)

    with open(f"result/input_1d_2s_e={e:.1f}.json", "w") as input_1d_2s_file:
        json.dump(inp_1d_2s, input_1d_2s_file, indent=4)

    for b in range(0, 10, 1):

        inp_2d = copy.deepcopy(input_2d.copy())

        limits = input_2d["zinq"][0]["options"]["grid"]["limits"][1]

        inp_2d["zinq"][0]["options"]["initial_conditions"]["momentum"] = [-numpy.sqrt(2 * mass * e), 0]
        inp_2d["zinq"][0]["options"]["initial_conditions"]["position"] = [x0, b]
        inp_2d["zinq"][0]["options"]["initial_conditions"]["mass"] = mass
        inp_2d["zinq"][0]["options"]["initial_conditions"]["state"] = 0

        inp_2d["zinq"][0]["options"]["grid"]["limits"][1] = [limits[0] + b, limits[1] + b]
        inp_2d["zinq"][0]["options"]["grid"]["points"] = 1024

        inp_2d["zinq"][0]["options"]["iterations"] = 250
        inp_2d["zinq"][0]["options"]["time_step"] = 5

        inp_2d_2s = copy.deepcopy(inp_2d)

        inp_2d["zinq"][0]["options"]["write"]["population"] = f"POPULATION_EXACT_2D_B={b}_E={e:.1f}.mat"
        inp_2d_2s["zinq"][0]["options"]["write"]["population"] = f"POPULATION_EXACT_2D_2S_B={b}_E={e:.1f}.mat"

        inp_2d["zinq"][0]["options"]["potential"]["file"]["nstate"] = 3
        inp_2d_2s["zinq"][0]["options"]["potential"]["file"]["nstate"] = 2

        inp_2d["zinq"][0]["options"]["potential"]["file"]["path"] = "../data/U_2D.mat"
        inp_2d_2s["zinq"][0]["options"]["potential"]["file"]["path"] = "../data/U_2D_2S.mat"

        with open(f"result/input_2d_b={b}_e={e:.1f}.json", "w") as input_2d_file:
            json.dump(inp_2d, input_2d_file, indent=4)

        with open(f"result/input_2d_2s_b={b}_e={e:.1f}.json", "w") as input_2d_2s_file:
            json.dump(inp_2d_2s, input_2d_2s_file, indent=4)
