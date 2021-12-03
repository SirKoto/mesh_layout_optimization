#pragma once

#include <Eigen/Dense>
#include <unordered_map>
#include <vector>
#include "TriangleMesh.hpp"

namespace LayoutMaker {

std::vector<uint32_t> get_mapping_optimized_layout(
	const std::shared_ptr<TriangleMesh> p_mesh,
	const uint32_t max_depth,
	const uint32_t max_cluster_size,
	const uint32_t max_spectral_size,
	const uint32_t max_number_interations_eigen,
	const float eigen_error
);
}