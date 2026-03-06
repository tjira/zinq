/// Precomputed nodes and weights for Gaussian-Hermite quadrature. This file is auto-generated.

const std = @import("std");

const eigenproblem_solver = @import("eigenproblem_solver.zig");
const error_handling = @import("error_handling.zig");
const global_variables = @import("global_variables.zig");
const real_matrix = @import("real_matrix.zig");

const RealMatrix = real_matrix.RealMatrix;

const eigensystemHermitian = eigenproblem_solver.eigensystemHermitian;
const throw = error_handling.throw;

const MAX_HERMITE_QUADRATURE_POINTS = global_variables.MAX_HERMITE_QUADRATURE_POINTS;

/// Gets the nodes for Gaussian-Hermite quadrature for a given n.
pub fn getNodes(comptime T: type, n: usize) ![]const T {
    return (try getPrecalculatedNodesAndWeights(T, n)).nodes;
}

/// Gets the precalculated nodes and weights for Gaussian-Hermite quadrature for a given n.
pub fn getPrecalculatedNodesAndWeights(comptime T: type, n: usize) !struct {nodes: []const T, weights: []const T} {

    if (n > MAX_HERMITE_QUADRATURE_POINTS) {

        const return_type = @typeInfo(@TypeOf(getPrecalculatedNodesAndWeights)).@"fn".return_type.?;

        return throw(return_type, "REQUESTED NUMBER OF QUADRATURE POINTS EXCEEDS THE MAXIMUM SUPPORTED", .{});
    }

    const HERMITE_NODES_AND_WEIGHTS = comptime naw: {
        break :naw getNodesAndWeightsTo(f64, MAX_HERMITE_QUADRATURE_POINTS) catch {
            @compileError("FAILED TO COMPUTE HERMITE QUADRATURE NODES AND WEIGHTS");
        };
    };

    return .{.nodes = HERMITE_NODES_AND_WEIGHTS.nodes[n], .weights = HERMITE_NODES_AND_WEIGHTS.weights[n]};
}

/// Gets the weights for Gaussian-Hermite quadrature for a given n.
pub fn getWeights(comptime T: type, n: usize) ![]const T {
    return (try getPrecalculatedNodesAndWeights(T, n)).weights;
}

/// Computes the nodes and weights for Gaussian-Hermite quadrature using the Golub-Welsch algorithm.
pub fn getNodesAndWeights(comptime T: type, comptime n: usize) !struct {nodes: [n]T, weights: [n]T} {
    var nodes: [n]T = undefined; var weights: [n]T = undefined;

    var workspace: [3 * n * n]T = undefined; for (0..workspace.len) |i| workspace[i] = 0;

    var J  = RealMatrix(T){.rows = n, .cols = n, .data = workspace[0 * n * n..1 * n * n]};
    var JJ = RealMatrix(T){.rows = n, .cols = n, .data = workspace[1 * n * n..2 * n * n]};
    var JA = RealMatrix(T){.rows = n, .cols = n, .data = workspace[2 * n * n..3 * n * n]};

    inline for (1..n) |i| {
        J.ptr(i, i - 1).* = std.math.sqrt(0.5 * @as(T, @floatFromInt(i))); J.ptr(i - 1, i).* = J.at(i, i - 1);
    }

    try eigensystemHermitian(T, &JJ, &JA, J);

    inline for (0..n) |i| {
        nodes[i] = JJ.at(i, i); weights[i] = std.math.sqrt(std.math.pi) * JA.at(0, i) * JA.at(0, i);
    }

    inline for (0..n / 2) |i| {

        const j = n - 1 - i;

        const average_node = (@abs(nodes[i]) + @abs(nodes[j])) / 2;

        nodes[i] = -average_node; nodes[j] = average_node;

        const average_weight = (weights[i] + weights[j]) / 2;

        weights[i] = average_weight; weights[j] = average_weight;
    }

    return .{.nodes = nodes, .weights = weights};
}

/// Precomputes the nodes and weights for Gaussian-Hermite quadrature for all n from 1 to the specified maximum n, and stores them in static arrays.
pub fn getNodesAndWeightsTo(comptime T: type, comptime n: usize) !struct {nodes: [n + 1][]const T, weights: [n + 1][]const T} {
    @setEvalBranchQuota(100_000_000);

    var nodes: [n + 1][]const T = undefined; var weights: [n + 1][]const T = undefined;

    nodes[0] = &.{}; weights[0] = &.{};

    inline for (1..n + 1) |i| {

        const nodes_and_weights = try getNodesAndWeights(T, i);

        const static_stoage = struct {
            const static_nodes = nodes_and_weights.nodes;
            const static_weights = nodes_and_weights.weights;
        };

        nodes[i] = &static_stoage.static_nodes; weights[i] = &static_stoage.static_weights;
    }

    return .{.nodes = nodes, .weights = weights};
}
