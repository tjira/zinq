#!/usr/bin/env python3
import argparse, math, numpy as np, sympy as sp


def evaluate(expr, r):
    symbols = [sp.Symbol(f"r{i + 1}") for i in range(len(r))]

    expr = sp.sympify(expr).subs({sp.Symbol("xyz"[i]): symbols[i] for i in range(min(len(r), 3))})

    return np.broadcast_to(sp.lambdify(symbols, expr, "numpy")(*r), r[0].shape).astype(float)


def getGrid(bounds, npoints):
    ndim = len(bounds)

    grids_1d = [np.linspace(bounds[i][0], bounds[i][1], npoints[i], endpoint=True) for i in range(ndim)]

    return np.column_stack([m.ravel() for m in np.meshgrid(*grids_1d, indexing="ij")])


def getData(functions, r_cols, ncol):
    evaluated = [evaluate(f, r_cols) for f in functions]

    idx = getIndices(len(functions))

    return np.array(evaluated)[idx].transpose(2, 0, 1).reshape(ncol, -1)


def getIndices(nfunc):
    N = (math.isqrt(1 + 8 * nfunc) - 1) // 2

    r, c = np.minimum(*np.ogrid[:N, :N]), np.maximum(*np.ogrid[:N, :N])

    return r * N - r * (r + 1) // 2 + c


def getLimits(grid, npoint):
    bounds = [(grid[2 * i], grid[2 * i + 1]) for i in range(len(grid) // 2)]

    return bounds, [npoint[0]] * (len(grid) // 2) if len(npoint) == 1 else npoint


def main():
    parser = argparse.ArgumentParser(description="Potential matrix generator")

    parser.add_argument("-o", "--output", type=str, default="POT.mat", help="Output file path")
    parser.add_argument("-n", "--npoint", nargs="+", type=int, default=[1024], help="Grid points per dimension")
    parser.add_argument("-g", "--grid", nargs="+", type=float, required=True, help="Grid ranges")
    parser.add_argument("funcs", nargs="+", type=str, help="Potential matrix element functions")

    args = parser.parse_args()

    grid = getGrid(*getLimits(args.grid, args.npoint))

    data = getData(args.funcs, [grid[:, i] for i in range(len(args.grid) // 2)], grid.shape[0])

    header = f"{grid.shape[0]} {grid.shape[1] + data.shape[1]}"

    np.savetxt(args.output, np.column_stack((grid, data)), fmt="%20.14f", header=header, comments="")


if __name__ == "__main__":
    main()
