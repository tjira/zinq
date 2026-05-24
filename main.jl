using Pkg; Pkg.activate(@__DIR__)

using Zinq, BenchmarkTools

ic = InitialConditions((-10.0,), (15.0,), (2.0,), 1, true)
grid = Grid(((-24.0, 32.0),), 2048)
m = 2000.0
iters = 3000
dt = 1.0

r = gen_grid_r(grid)
k = gen_grid_k(grid)

W = gen_wfn(ic, r, grid.npoint, 2)
V = tully_1(r..., 0.01)
A, U = get_pot_eigen(V)
K = get_prop_k(grid, m, k, dt)
R = get_prop_r(A, U, dt / 2)

normalize!(W, grid)
W_k = to_kspace(W)

for i in 1:iters + 1
    if i > 1
        propagate!(W, R, K)
    end

    pop = calc_pop(W, grid)
    println("Iteration: $i, Population: $pop")
end

# norm = calc_norm(W, grid)
# pe = calc_pe(W, V, grid)
# ke = calc_ke(to_kspace(W), grid, m, k)
# pop = calc_pop(W, grid)
# pos = calc_pos(W, grid, r)
# mom = calc_mom(to_kspace(W), grid, k)
#
# println("Final Norm: $norm")
# println("Final PE: $pe")
# println("Final KE: $ke")
# println("Final Total Energy: $(ke + pe)")
# println("Final Pop: $pop")
# println("Final Pos: $pos")
# println("Final Mom: $mom")
