#pragma once

#include "TriangleMesh.hpp"

namespace LayoutOptimizer {

// Get mapping of vertices to new, better positions
std::vector<uint32_t> optimize_layout(const TriangleMesh& mesh, const std::vector<uint32_t>& clusters);

} // namespace