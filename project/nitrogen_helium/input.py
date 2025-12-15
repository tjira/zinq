import copy, json, numpy, os

USE_SAMPLED_POINTS = True

SAMPLED_POINTS = [0.040, 0.045, 0.048, 0.049, 0.050, 0.051, 0.052, 0.053, 0.054, 0.055, 0.058, 0.059, 0.060, 0.061, 0.062, 0.063, 0.064, 0.065, 0.069, 0.071, 0.073, 0.075, 0.076, 0.078, 0.084, 0.088, 0.090, 0.092, 0.097, 0.099, 0.101, 0.105, 0.107, 0.108, 0.110, 0.115, 0.117, 0.118, 0.119, 0.121, 0.125, 0.127, 0.129, 0.130, 0.131, 0.133, 0.135, 0.139, 0.141, 0.143, 0.144, 0.145, 0.147, 0.153, 0.156, 0.158, 0.160, 0.168, 0.171, 0.173, 0.177, 0.184, 0.187, 0.190, 0.192, 0.195, 0.198, 0.201, 0.207, 0.210, 0.213, 0.215, 0.217, 0.220, 0.229, 0.232, 0.235, 0.237, 0.239, 0.242, 0.255, 0.258, 0.261, 0.265, 0.269, 0.283, 0.290, 0.296, 0.307, 0.312, 0.317, 0.326, 0.330, 0.333, 0.337, 0.341, 0.353, 0.357, 0.361, 0.365, 0.369, 0.382, 0.387, 0.392, 0.396, 0.401, 0.407, 0.424, 0.431, 0.437, 0.444, 0.451, 0.459, 0.482, 0.491, 0.499, 0.507, 0.515, 0.524, 0.562, 0.572, 0.581, 0.600, 0.651, 0.665, 0.678, 0.691, 0.705, 0.771, 0.790, 0.809, 0.835, 0.865, 0.899, 0.976, 1.000, 1.030, 1.080, 1.280, 1.320, 1.370, 1.410, 1.450, 1.670, 1.790, 1.900, 2.000, 2.120, 2.530, 2.790, 3.050, 4.070, 4.570, 5.340, 6.430, 7.310, 9.550, 18.800, 26.200, 32.000]

inputs_exc_1d_2s, inputs_exc_1d_3s = [], []
inputs_exc_2d_2s, inputs_exc_2d_3s = [], []

input_1d = { "zinq" : [{
    "name" : "quantum_dynamics",
    "options" : {
        "grid" : {
            "limits" : [[0.1, 18]]
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
            "limits" : [[-10, 18], [-8, 8]]
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

ALL_POINTS = list(numpy.arange(0.04, 1, 0.001)) + list(numpy.arange(1, 10, 0.01)) + list(numpy.arange(10, 32.1, 0.1))

for e in (SAMPLED_POINTS if USE_SAMPLED_POINTS else ALL_POINTS):

    match e:
        case _ if e <=   0.5: points = 512
        case _ if e <=   1.0: points = 1024
        case _ if e <=  10.0: points = 2048
        case _ if e <=  50.0: points = 4096
        case _ if e <= 100.0: points = 8192
        case _: raise ValueError("ENERGY OUT OF RANGE")

    time_step = 1

    if e == 0: continue

    inp_exc_1d_3s = copy.deepcopy(input_1d)

    inp_exc_1d_3s["zinq"][0]["options"]["initial_conditions"]["momentum"] = [-numpy.sqrt(2 * mass * e)]
    inp_exc_1d_3s["zinq"][0]["options"]["initial_conditions"]["position"] = [x0]
    inp_exc_1d_3s["zinq"][0]["options"]["initial_conditions"]["mass"] = mass
    inp_exc_1d_3s["zinq"][0]["options"]["initial_conditions"]["state"] = 0

    inp_exc_1d_3s["zinq"][0]["options"]["grid"]["points"] = points

    inp_exc_1d_3s["zinq"][0]["options"]["iterations"] = int(1100 / numpy.sqrt(e) / time_step)
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

    for b in range(0, 10, 1):

        inp_exc_2d_3s = copy.deepcopy(input_2d)

        limits = input_2d["zinq"][0]["options"]["grid"]["limits"][1]

        inp_exc_2d_3s["zinq"][0]["options"]["initial_conditions"]["momentum"] = [-numpy.sqrt(2 * mass * e), 0]
        inp_exc_2d_3s["zinq"][0]["options"]["initial_conditions"]["position"] = [x0, b]
        inp_exc_2d_3s["zinq"][0]["options"]["initial_conditions"]["mass"] = mass
        inp_exc_2d_3s["zinq"][0]["options"]["initial_conditions"]["state"] = 0

        inp_exc_2d_3s["zinq"][0]["options"]["grid"]["limits"][1] = [limits[0] + b, limits[1] + b]
        inp_exc_2d_3s["zinq"][0]["options"]["grid"]["points"] = 1024

        inp_exc_2d_3s["zinq"][0]["options"]["iterations"] = 250
        inp_exc_2d_3s["zinq"][0]["options"]["time_step"] = 5

        inp_exc_2d_2s = copy.deepcopy(inp_exc_2d_3s)

        inp_exc_2d_3s["zinq"][0]["options"]["write"]["population"] = f"POPULATION_EXACT_EXC_2D_3S_B={b}_E={e:.3f}.mat"
        inp_exc_2d_2s["zinq"][0]["options"]["write"]["population"] = f"POPULATION_EXACT_EXC_2D_2S_B={b}_E={e:.3f}.mat"

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

with open(f"result/input_exc_1d_2s_all.json", "w") as input_exc_1d_2s_all_file:
    json.dump(input_exc_1d_2s_all, input_exc_1d_2s_all_file, indent=4)

with open(f"result/input_exc_1d_3s_all.json", "w") as input_exc_1d_3s_all_file:
    json.dump(input_exc_1d_3s_all, input_exc_1d_3s_all_file, indent=4)

with open(f"result/input_exc_2d_2s_all.json", "w") as input_exc_2d_2s_all_file:
    json.dump(input_exc_2d_2s_all, input_exc_2d_2s_all_file, indent=4)

with open(f"result/input_exc_2d_3s_all.json", "w") as input_exc_2d_3s_all_file:
    json.dump(input_exc_2d_3s_all, input_exc_2d_3s_all_file, indent=4)
