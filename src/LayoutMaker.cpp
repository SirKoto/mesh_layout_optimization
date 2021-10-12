#include "LayoutMaker.hpp"

#include <Eigen/SparseCore>
#include <Spectra/SymEigsSolver.h>
#include <Spectra/MatOp/SparseSymMatProd.h>
#include <unordered_set>
#include <unordered_map>
#include <iostream>

namespace LayoutMaker {

std::unordered_map<uint32_t, uint32_t>
vertex_laplacian_layout(const std::vector<Eigen::Array3i>& faces,
	const std::unordered_multimap<uint32_t, uint32_t>& vert2face,
	const std::unordered_set<uint32_t>& vertices_indices,
	const uint32_t num_max_vertices) {



	std::vector<uint32_t> vert_degree(num_max_vertices, 0);
	std::vector< Eigen::Triplet<float>> triplet_list;
	triplet_list.reserve(3 * vertices_indices.size());
	// Fill connectivity
	for (uint32_t i : vertices_indices) {

		auto range = vert2face.equal_range(i);
		std::for_each(range.first, range.second,
			[&](const std::pair<const uint32_t, uint32_t>& f_id) {
				const Eigen::Array3i& face = faces.at(f_id.second);
				for (uint32_t j = 0; j < 3; ++j) {
					uint32_t v = face[j];
					if (v != i) {
						if (vertices_indices.count(v) != 0) {
							vert_degree[i] += 1;
							triplet_list.push_back(Eigen::Triplet<float>(i, v, -1.f));
						}
					}
				}
			});
	}

	// fill degree triplets
	for (uint32_t i : vertices_indices) {
		triplet_list.push_back(Eigen::Triplet<float>(i, i, (float)vert_degree[i]));
	}
	

	// Square sparse matrix
	Eigen::SparseMatrix<float> laplacian(num_max_vertices, num_max_vertices);
	laplacian.setFromTriplets(triplet_list.begin(), triplet_list.end());
	laplacian.makeCompressed();

	Spectra::SparseSymMatProd<float> op(laplacian);
	// Compute 2 eigenvectors
	Spectra::SymEigsSolver<Spectra::SparseSymMatProd<float>> eigs(op, 2, 4);
	eigs.init();

	eigs.compute(Spectra::SortRule::LargestAlge);

	// Get results
	if (eigs.info() == Spectra::CompInfo::Successful) {
		std::cout << "Eigenvalues\n" <<eigs.eigenvalues() << std::endl;
	}
	else {
		std::cout << "No eigenvalues" << std::endl;
	}

	return {};
}


std::unordered_map<uint32_t, uint32_t>
get_mapping_optimized_layout(
	const std::vector<Eigen::Array3i>& faces,
	uint32_t num_vertices,
	const MultiLevel multi_level_approach)

{
	std::unordered_set<uint32_t> vert_indices;
	vert_indices.reserve(num_vertices);
	for (uint32_t i = 0; i < num_vertices; ++i) {
		vert_indices.insert(i);
	}

	std::unordered_multimap<uint32_t, uint32_t> vert2face;
	vert2face.reserve(num_vertices);
	for (uint32_t f = 0; f < (uint32_t)faces.size(); ++f) {
		for (uint32_t j = 0; j < 3; ++j) {
			vert2face.insert({faces[f][j], f});
		}
	}

	if (multi_level_approach == MultiLevel::eVertexLaplacian) {
		return vertex_laplacian_layout(faces, vert2face, vert_indices, num_vertices);
	}

	return {};
}

}