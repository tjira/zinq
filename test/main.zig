const std = @import("std");
const zinq = @import("zinq");

test {
    _ = @import("quantum_dynamics.zig");
    _ = @import("classical_dynamics.zig");
    _ = @import("hartree_fock.zig");
    _ = @import("moller_plesset.zig");
    _ = @import("configuration_interaction.zig");
}
