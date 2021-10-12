#pragma once

#include <Eigen/Dense>
#include <unordered_map>
#include <vector>

namespace LayoutMaker {

enum class MultiLevel {
	eVertexLaplacian
};

std::unordered_map<uint32_t, uint32_t> get_mapping_optimized_layout(
	const std::vector<Eigen::Array3i>& faces,
	uint32_t num_vertices,
	const MultiLevel multi_level_approach
);
}