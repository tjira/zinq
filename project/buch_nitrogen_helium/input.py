import copy, json, numpy

with open("input_1d.json", "r") as input_1d_file:
    input_1d = json.loads(input_1d_file.read())

with open("input_2d.json", "r") as input_2d_file:
    input_2d = json.loads(input_2d_file.read())

mass_1d = input_1d["zinq"][0]["options"]["initial_conditions"]["mass"]
mass_2d = input_2d["zinq"][0]["options"]["initial_conditions"]["mass"]

if mass_1d != mass_2d: raise ValueError("MASS IN 1D AND 2D INPUT FILES DO NOT MATCH")

x0_1d = input_1d["zinq"][0]["options"]["initial_conditions"]["position"][0]
x0_2d = input_2d["zinq"][0]["options"]["initial_conditions"]["position"][0]

if x0_1d != x0_2d: raise ValueError("POSITION IN 1D AND 2D INPUT FILES DO NOT MATCH")

for e in [1] + list(range(0, 100, 10)) + [100]:

    if e == 0: continue

    inp_1d = copy.deepcopy(input_1d)

    pop_file_1d = f"POPULATION_EXACT_1D_E={e}.mat"

    inp_1d["zinq"][0]["options"]["initial_conditions"]["momentum"][0] = -numpy.sqrt(2 * mass_1d * e)
    inp_1d["zinq"][0]["options"]["initial_conditions"]["position"][0] = x0_1d

    inp_1d["zinq"][0]["options"]["write"]["population"] = pop_file_1d

    with open(f"input_1d_e={e}.json", "w") as input_1d_file:
        json.dump(inp_1d, input_1d_file, indent=4)

    for b in range(0, 10, 1):

        inp_2d = copy.deepcopy(input_2d.copy())

        pop_file_2d = f"POPULATION_EXACT_2D_B={b}_E={e}.mat"

        limits = input_2d["zinq"][0]["options"]["grid"]["limits"][1]

        inp_2d["zinq"][0]["options"]["initial_conditions"]["momentum"] = [-numpy.sqrt(2 * mass_2d * e), 0]
        inp_2d["zinq"][0]["options"]["initial_conditions"]["position"] = [x0_2d, b]

        inp_2d["zinq"][0]["options"]["grid"]["limits"][1] = [limits[0] + b, limits[1] + b]

        inp_2d["zinq"][0]["options"]["write"]["population"] = pop_file_2d

        with open(f"input_2d_b={b}_e={e}.json", "w") as input_2d_file:
            json.dump(inp_2d, input_2d_file, indent=4)
