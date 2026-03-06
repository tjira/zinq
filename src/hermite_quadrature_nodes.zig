/// Precomputed nodes and weights for Gaussian-Hermite quadrature. This file is auto-generated.

const std = @import("std");

const eigenproblem_solver = @import("eigenproblem_solver.zig");
const global_variables = @import("global_variables.zig");
const real_matrix = @import("real_matrix.zig");

const RealMatrix = real_matrix.RealMatrix;

const eigensystemHermitian = eigenproblem_solver.eigensystemHermitian;

const MAX_HERMITE_QUADRATURE_POINTS = global_variables.MAX_HERMITE_QUADRATURE_POINTS;

/// Gets the nodes for Gaussian-Hermite quadrature for a given n.
pub fn getNodes(comptime T: type, n: usize) ![]const T {
    return (try getPrecalculatedNodesAndWeights(T, n)).nodes;
}

/// Gets the precalculated nodes and weights for Gaussian-Hermite quadrature for a given n.
pub fn getPrecalculatedNodesAndWeights(comptime T: type, n: usize) !struct {nodes: []const T, weights: []const T} {

    if (n > MAX_HERMITE_QUADRATURE_POINTS or n == 0) return error.InvalidQuadraturePoints;

    const HERMITE_NODES_AND_WEIGHTS = comptime blk: {
        @setEvalBranchQuota(10000000); break :blk getNodesAndWeightsTo(T, MAX_HERMITE_QUADRATURE_POINTS);
    };

    return .{.nodes = HERMITE_NODES_AND_WEIGHTS.nodes[n], .weights = HERMITE_NODES_AND_WEIGHTS.weights[n]};
}

/// Gets the weights for Gaussian-Hermite quadrature for a given n.
pub fn getWeights(comptime T: type, n: usize) ![]const T {
    return (try getPrecalculatedNodesAndWeights(T, n)).weights;
}

/// Computes the nodes and weights for Gaussian-Hermite quadrature using the Newton-Raphson method.
pub fn getNodesAndWeights(comptime T: type, comptime n: usize) struct {nodes: [n]T, weights: [n]T} {
    var nodes: [n]T = undefined; var weights: [n]T = undefined; const m = (n + 1) / 2;
        
    var z: T = undefined; var dp: T = undefined; const nf: T = @floatFromInt(n); const iters = 100;

    inline for (0..m) |i| {

        if (i == 0) {
            z = std.math.sqrt(2.0 * nf + 1.0) - 1.85575 * std.math.pow(T, 2.0 * nf + 1.0, -0.16667);
        } else if (i == 1) {
            z -= 1.14 * std.math.pow(T, nf, 0.426) / z;
        } else if (i == 2) {
            z = 1.86 * z - 0.86 * nodes[n - 1];
        } else if (i == 3) {
            z = 1.91 * z - 0.91 * nodes[n - 2];
        } else {
            z = 2.0 * z - nodes[n - i + 1];
        }

        inline for (0..iters) |j| {

            var p1: T = 1 / std.math.pow(T, std.math.pi, 0.25); var p2: T = 0;

            for (1..n + 1) |k| {

                const kf: T = @floatFromInt(k); const p3 = p2; p2 = p1; 
                
                p1 = z * std.math.sqrt(2.0 / kf) * p2 - std.math.sqrt((kf - 1.0) / kf) * p3;
            }

            dp = std.math.sqrt(2.0 * nf) * p2; const z_old = z; z = z - p1 / dp;

            if (@abs(z - z_old) < 10.0 * std.math.floatEps(T)) break;

            if (j == iters - 1) @compileError("FAILED TO CONVERGE TO A SOLUTION FOR HERMITE QUADRATURE NODES");
        }

        nodes[i] = -z; nodes[n - i - 1] = z; weights[i] = 2.0 / (dp * dp); weights[n - i - 1] = weights[i];
    }

    return .{.nodes = nodes, .weights = weights};
}

/// Precomputes the nodes and weights for Gaussian-Hermite quadrature for all n from 1 to the specified maximum n, and stores them in static arrays.
pub fn getNodesAndWeightsTo(comptime T: type, comptime n: usize) struct {nodes: [n + 1][]const T, weights: [n + 1][]const T} {
    var nodes: [n + 1][]const T = undefined; var weights: [n + 1][]const T = undefined;

    nodes[0] = &.{}; weights[0] = &.{};

    inline for (1..n + 1) |i| {

        const nodes_and_weights = getNodesAndWeights(T, i);

        const static_stoage = struct {
            const static_nodes = nodes_and_weights.nodes;
            const static_weights = nodes_and_weights.weights;
        };

        nodes[i] = &static_stoage.static_nodes; weights[i] = &static_stoage.static_weights;
    }

    return .{.nodes = nodes, .weights = weights};
}
