import numpy as np

U_1D = np.loadtxt("data/U_1D.dat", skiprows=1)

grid_2D = np.stack([G.ravel() for G in np.meshgrid(*[np.linspace(-U_1D[-1, 0], U_1D[-1, 0], 1024)] * 2, indexing="ij")], axis=-1).reshape(-1, 2)

V00_2D = np.interp(np.linalg.norm(grid_2D, axis=1), U_1D[:, 0], U_1D[:, 1])
V01_2D = np.interp(np.linalg.norm(grid_2D, axis=1), U_1D[:, 0], U_1D[:, 2])
V02_2D = np.interp(np.linalg.norm(grid_2D, axis=1), U_1D[:, 0], U_1D[:, 3])
V10_2D = np.interp(np.linalg.norm(grid_2D, axis=1), U_1D[:, 0], U_1D[:, 4])
V11_2D = np.interp(np.linalg.norm(grid_2D, axis=1), U_1D[:, 0], U_1D[:, 5])
V12_2D = np.interp(np.linalg.norm(grid_2D, axis=1), U_1D[:, 0], U_1D[:, 6])
V20_2D = np.interp(np.linalg.norm(grid_2D, axis=1), U_1D[:, 0], U_1D[:, 7])
V21_2D = np.interp(np.linalg.norm(grid_2D, axis=1), U_1D[:, 0], U_1D[:, 8])
V22_2D = np.interp(np.linalg.norm(grid_2D, axis=1), U_1D[:, 0], U_1D[:, 9])

np.savetxt("data/U_1D_2S.mat", np.column_stack((U_1D[:, 0], U_1D[:, 1], U_1D[:, 2], U_1D[:, 4], U_1D[:, 5])), fmt="%20.14f", header=f"{U_1D.shape[0]} 5", comments="")

np.savetxt("data/U_2D_2S.mat", np.column_stack((grid_2D, V00_2D, V01_2D,         V10_2D, V11_2D                                )), fmt="%20.14f", header=f"{grid_2D.shape[0]} 6",  comments="")
np.savetxt("data/U_2D.mat",    np.column_stack((grid_2D, V00_2D, V01_2D, V02_2D, V10_2D, V11_2D, V12_2D, V20_2D, V21_2D, V22_2D)), fmt="%20.14f", header=f"{grid_2D.shape[0]} 11", comments="")
