import numpy, os, pandas, scipy

EXC_1D_2S, EXC_1D_2S_PROB = [], []
EXC_2D_2S, EXC_2D_2S_PROB = [], []

def findPlateau(y, window, percentile):
    slope = numpy.gradient(scipy.ndimage.uniform_filter1d(y, size=window))

    is_flat = numpy.abs(slope) < numpy.nanpercentile(numpy.abs(slope), percentile)

    return numpy.where(is_flat)[0][-1]

def extractProbability(column):
    start = findPlateau(column, 100, 5)

    return column[start:min(start + 10, len(column))].mean()

for filename in os.listdir("result"):
    if filename.startswith("POPULATION_EXACT_EXC_1D_2S_E=") and filename.endswith(".mat"):
        EXC_1D_2S.append(filename)
    if filename.startswith("POPULATION_EXACT_EXC_2D_2S_B=") and filename.endswith(".mat"):
        EXC_2D_2S.append(filename)

for filename in EXC_1D_2S:

    e = filename.split("E=")[1].split(".mat")[0]

    data = numpy.loadtxt(os.path.join("result", filename), skiprows=1)

    EXC_1D_2S_PROB.append([float(e), extractProbability(data[:, 2])])

for filename in EXC_2D_2S:

    e, b = filename.split("E=")[1].split(".mat")[0], filename.split("_B=")[1].split("_E=")[0]

    data = numpy.loadtxt(os.path.join("result", filename), skiprows=1)

    EXC_2D_2S_PROB.append([float(b), float(e), extractProbability(data[:, 2])])

EXC_1D_2S_PROB = numpy.array(EXC_1D_2S_PROB)[numpy.argsort(numpy.array(EXC_1D_2S_PROB)[:, 0])] if EXC_1D_2S_PROB else numpy.array([])
EXC_2D_2S_PROB = numpy.array(EXC_2D_2S_PROB)[numpy.argsort(numpy.array(EXC_2D_2S_PROB)[:, 0])] if EXC_2D_2S_PROB else numpy.array([])

if EXC_1D_2S_PROB.size != 0: numpy.savetxt("PROB_EXC_1D_2S.mat", EXC_1D_2S_PROB, header=f"{EXC_1D_2S_PROB.shape[0]} {EXC_1D_2S_PROB.shape[1]}", comments="")

for b in numpy.unique(EXC_2D_2S_PROB[:, 0]):

    subset = EXC_2D_2S_PROB[EXC_2D_2S_PROB[:, 0] == b][:, 1:]

    numpy.savetxt(f"PROB_EXC_2D_2S_B={b:.2f}.mat", subset[numpy.argsort(subset[:, 0])], header=f"{subset.shape[0]} {subset.shape[1]}", comments="")

pandas.DataFrame(EXC_2D_2S_PROB, columns=["b", "E", "P"]).pivot(index="b", columns="E", values="P").rename_axis(index="b\\E", columns=None).to_csv("PROB_EXC_2D_2S.csv", index=True)
