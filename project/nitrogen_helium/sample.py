import scipy.interpolate, matplotlib.pyplot, numpy

FILENAME, POINTS = "PROB_EXC_1D_2S.mat", 40

def sample(x, y, n):
    selected = [0, len(x) - 1]

    while len(selected) < n:

        current = sorted(selected)

        interpf = scipy.interpolate.interp1d(x[current], y[current], kind="linear", fill_value="extrapolate")

        errors = numpy.abs(y - interpf(x))

        errors[current] = 0

        selected.append(numpy.argmax(errors))

    return sorted(selected)

data = numpy.loadtxt(FILENAME, skiprows=1)

indices = sample(data[:, 0], data[:, 1], POINTS)

print(f"SAMPLED_POINTS = [{', '.join(f'{x:.3f}' for x in data[:, 0][indices])}]")

matplotlib.pyplot.plot(data[:, 0], data[:, 1], label="Original", linewidth=3, alpha=0.5)

matplotlib.pyplot.plot(data[:, 0][indices], data[:, 1][indices], label="Sampled Interpolation")

matplotlib.pyplot.legend(); matplotlib.pyplot.tight_layout(); matplotlib.pyplot.show()
