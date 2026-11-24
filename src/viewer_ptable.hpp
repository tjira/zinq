#pragma once

#include "viewer_math.hpp"

#include <string>
#include <unordered_map>

/// Physical and visual characteristics of an element including covalent and van der Waals radii, CPK color, and mass.
struct AtomInfo {
    float radius; float covalent; vec3 color; double mass;
};

extern const std::unordered_map<std::string, AtomInfo> ptable;

/// Retrieves element properties by chemical symbol with a fallback for unknown elements.
const AtomInfo& get_atom_info(const std::string& symbol);
