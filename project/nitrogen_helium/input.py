import copy, json, numpy, os

USE_SAMPLED_POINTS = True

SAMPLED_POINTS = [0.040, 0.048, 0.051, 0.054, 0.059, 0.061, 0.062, 0.064, 0.073, 0.078, 0.090, 0.099, 0.107, 0.118, 0.130, 0.144, 0.156, 0.173, 0.192, 0.215, 0.237, 0.261, 0.290, 0.312, 0.330, 0.361, 0.396, 0.431, 0.444, 0.493, 0.513, 0.581, 0.689, 0.822, 1.040, 1.370, 1.990, 2.800, 4.550, 10.000]

inputs_exc_1d_2s, inputs_exc_1d_3s = [], []
inputs_exc_2d_2s, inputs_exc_2d_3s = [], []

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
            "limits" : [[-12, 16], [-12, 12]]
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

mass, x0 = 5674.685752290152, 10

ALL_POINTS = list(numpy.arange(0.04, 1, 0.001)) + list(numpy.arange(1, 10.01, 0.01))

for e in (SAMPLED_POINTS if USE_SAMPLED_POINTS else ALL_POINTS):

    match e:
        case _ if e <=  0.5: points = 1024
        case _ if e <=  1.0: points = 2048
        case _ if e <= 10.0: points = 4096
        case _: raise ValueError("ENERGY OUT OF RANGE")

    time_step = 1

    if e == 0: continue

    inp_exc_1d_3s = copy.deepcopy(input_1d)

    inp_exc_1d_3s["zinq"][0]["options"]["initial_conditions"]["momentum"] = [-numpy.sqrt(2 * mass * e)]
    inp_exc_1d_3s["zinq"][0]["options"]["initial_conditions"]["position"] = [x0]
    inp_exc_1d_3s["zinq"][0]["options"]["initial_conditions"]["mass"] = mass
    inp_exc_1d_3s["zinq"][0]["options"]["initial_conditions"]["state"] = 0

    inp_exc_1d_3s["zinq"][0]["options"]["grid"]["points"] = points

    inp_exc_1d_3s["zinq"][0]["options"]["iterations"] = int(1500 / numpy.sqrt(e) / time_step)
    inp_exc_1d_3s["zinq"][0]["options"]["time_step"] = time_step

    inp_exc_1d_2s = copy.deepcopy(inp_exc_1d_3s)

    inp_exc_1d_3s["zinq"][0]["options"]["write"]["population"] = f"POPULATION_EXACT_EXC_1D_3S_E={e:.3f}.mat"
    inp_exc_1d_2s["zinq"][0]["options"]["write"]["population"] = f"POPULATION_EXACT_EXC_1D_2S_E={e:.3f}.mat"

    inp_exc_1d_3s["zinq"][0]["options"]["potential"]["file"]["nstate"] = 3
    inp_exc_1d_2s["zinq"][0]["options"]["potential"]["file"]["nstate"] = 2

    inp_exc_1d_3s["zinq"][0]["options"]["potential"]["file"]["path"] = "../data/U_1D_3S.dat"
    inp_exc_1d_2s["zinq"][0]["options"]["potential"]["file"]["path"] = "../data/U_1D_2S.mat"

    inputs_exc_1d_3s.append(inp_exc_1d_3s)
    inputs_exc_1d_2s.append(inp_exc_1d_2s)

    for b in list(numpy.arange(0, 3, 0.25)) + list(numpy.arange(3, 4, 0.1)) + [4, 7, 10]:

        inp_exc_2d_3s = copy.deepcopy(input_2d)

        limits = input_2d["zinq"][0]["options"]["grid"]["limits"][1]

        inp_exc_2d_3s["zinq"][0]["options"]["initial_conditions"]["momentum"] = [-numpy.sqrt(2 * mass * e), 0]
        inp_exc_2d_3s["zinq"][0]["options"]["initial_conditions"]["position"] = [x0, b]
        inp_exc_2d_3s["zinq"][0]["options"]["initial_conditions"]["mass"] = mass
        inp_exc_2d_3s["zinq"][0]["options"]["initial_conditions"]["state"] = 0

        inp_exc_2d_3s["zinq"][0]["options"]["grid"]["limits"][1] = [limits[0] + b, limits[1] + b]
        inp_exc_2d_3s["zinq"][0]["options"]["grid"]["points"] = points

        inp_exc_2d_3s["zinq"][0]["options"]["iterations"] = int(1500 / numpy.sqrt(e) / time_step)
        inp_exc_2d_3s["zinq"][0]["options"]["time_step"] = time_step

        inp_exc_2d_2s = copy.deepcopy(inp_exc_2d_3s)

        inp_exc_2d_3s["zinq"][0]["options"]["write"]["population"] = f"POPULATION_EXACT_EXC_2D_3S_B={b:.2f}_E={e:.3f}.mat"
        inp_exc_2d_2s["zinq"][0]["options"]["write"]["population"] = f"POPULATION_EXACT_EXC_2D_2S_B={b:.2f}_E={e:.3f}.mat"

        inp_exc_2d_3s["zinq"][0]["options"]["potential"]["file"]["nstate"] = 3
        inp_exc_2d_2s["zinq"][0]["options"]["potential"]["file"]["nstate"] = 2

        inp_exc_2d_3s["zinq"][0]["options"]["potential"]["file"]["path"] = "../data/U_2D_3S.mat"
        inp_exc_2d_2s["zinq"][0]["options"]["potential"]["file"]["path"] = "../data/U_2D_2S.mat"

        inputs_exc_2d_3s.append(inp_exc_2d_3s)
        inputs_exc_2d_2s.append(inp_exc_2d_2s)

input_exc_1d_2s_all = { "zinq" : []}
input_exc_1d_3s_all = { "zinq" : []}
input_exc_2d_2s_all = { "zinq" : []}
input_exc_2d_3s_all = { "zinq" : []}

for inp in inputs_exc_1d_2s: input_exc_1d_2s_all["zinq"].append(inp["zinq"][0])
for inp in inputs_exc_1d_3s: input_exc_1d_3s_all["zinq"].append(inp["zinq"][0])
for inp in inputs_exc_2d_2s: input_exc_2d_2s_all["zinq"].append(inp["zinq"][0])
for inp in inputs_exc_2d_3s: input_exc_2d_3s_all["zinq"].append(inp["zinq"][0])

with open(f"result/input_exact_exc_1d_2s_all.json", "w") as input_exc_1d_2s_all_file:
    json.dump(input_exc_1d_2s_all, input_exc_1d_2s_all_file, indent=4)

for inp in inputs_exc_2d_2s:
    
    filename = inp["zinq"][0]["options"]["write"]["population"].lower().replace("population", "input").replace("mat", "json")

    with open(f"result/{filename}", "w") as input_exc_2d_2s_file:
        json.dump(inp, input_exc_2d_2s_file, indent=4)
