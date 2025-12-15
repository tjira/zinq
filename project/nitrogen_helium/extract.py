import numpy, os

EXC_1D_2S, EXC_1D_2S_PROB = [], []

for filename in os.listdir("result"):
    if filename.startswith("POPULATION_EXACT_EXC_1D_2S_E=") and filename.endswith(".mat"):
        EXC_1D_2S.append(filename)

for filename in EXC_1D_2S:

    e = filename.split("E=")[1].split(".mat")[0]

    data = numpy.loadtxt(os.path.join("result", filename), skiprows=1)

    EXC_1D_2S_PROB.append([float(e), data[-1, 2]])

EXC_1D_2S_PROB = numpy.array(EXC_1D_2S_PROB)[numpy.argsort(numpy.array(EXC_1D_2S_PROB)[:, 0])]

numpy.savetxt("PROB_EXC_1D_2S.mat", EXC_1D_2S_PROB, header=f"{EXC_1D_2S_PROB.shape[0]} {EXC_1D_2S_PROB.shape[1]}", comments="")
