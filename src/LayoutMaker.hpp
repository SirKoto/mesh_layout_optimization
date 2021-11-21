#pragma once

#include <Eigen/Dense>
#include <unordered_map>
#include <vector>

namespace LayoutMaker {

enum class MultiLevel {
	eVertexLaplacian
};

std::vector<uint32_t> get_mapping_optimized_layout(
	const std::vector<Eigen::Array3i>& faces,
	uint32_t num_vertices,
	const MultiLevel multi_level_approach,
	const uint32_t max_depth,
	const uint32_t max_cluster_size,
	const uint32_t max_number_interations_eigen,
	const float eigen_error
);
}