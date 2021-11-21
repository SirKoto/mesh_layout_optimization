#include "LayoutMaker.hpp"

#include <Eigen/SparseCore>
#include <Spectra/SymEigsSolver.h>
#include <Spectra/MatOp/SparseSymMatProd.h>
#include <unordered_set>
#include <unordered_map>
#include <iostream>

namespace LayoutMaker {

struct LayoutContext {

	uint32_t next_id;
	const uint32_t max_depth;
	const uint32_t max_cluster_size;
	std::vector<uint32_t> final_cluster;
	const uint32_t max_iterations_eigen;
	const float error_eigen;

	LayoutContext(uint32_t max_depth,
		uint32_t max_cluster_size,
		uint32_t num_vertices,
		uint32_t max_iterations_eigen,
		float error_eigen) :
		next_id(0), max_depth(max_depth), max_cluster_size(max_cluster_size),
		max_iterations_eigen(max_iterations_eigen),
		error_eigen(error_eigen)
	{
		final_cluster.resize(num_vertices, 0);
	}
};

void
vertex_laplacian_layout(
	LayoutContext& context,
	const uint32_t depth,
	const std::vector<Eigen::Array3i>& faces,
	const std::unordered_multimap<uint32_t, uint32_t>& vert2face,
	const std::unordered_set<uint32_t>& vertices_indices) {

	assert(!vertices_indices.empty());
	std::unordered_map<uint32_t, uint32_t> old2new_vert;
	std::vector<uint32_t> new2old_vert(vertices_indices.size());
	old2new_vert.reserve(vertices_indices.size());
	{
		std::unordered_set<uint32_t>::const_iterator it = vertices_indices.begin();
		for (uint32_t i = 0; i < (uint32_t)vertices_indices.size(); ++i, ++it) {
			old2new_vert.insert({ *it, i });
			new2old_vert[i] = *it;
		}
	}

	// Compute second smallest eigenvector
	std::vector<uint32_t> vert_degree(vertices_indices.size(), 0);
	std::vector< Eigen::Triplet<float>> triplet_list;
	triplet_list.reserve(3 * vertices_indices.size());
	// Fill connectivity
	std::unordered_set<uint32_t>::const_iterator it = vertices_indices.begin();
	for (uint32_t v_new = 0; v_new < (uint32_t)vertices_indices.size(); ++v_new, ++it) {
		uint32_t v_old = *it;
		auto range = vert2face.equal_range(v_old);
		std::for_each(range.first, range.second,
			[&](const std::pair<const uint32_t, uint32_t>& f_id) {
				const Eigen::Array3i& face = faces.at(f_id.second);
				for (uint32_t j = 0; j < 3; ++j) {
					uint32_t v2_old = face[j];
					if (v2_old != v_old) {
						if (vertices_indices.count(v2_old) != 0) {
							vert_degree[v_new] += 1;
							uint32_t v2_new = old2new_vert.at(v2_old);
							triplet_list.push_back(Eigen::Triplet<float>(v_new, v2_new, -1.f));
						}
					}
				}
			});
	}

	// fill degree triplets
	for (uint32_t i = 0; i < (uint32_t)vertices_indices.size(); ++i) {
		triplet_list.push_back(Eigen::Triplet<float>(i, i, (float)vert_degree[i]));
	}


	// Square sparse matrix
	Eigen::SparseMatrix<float> laplacian(vertices_indices.size(), vertices_indices.size());
	laplacian.setFromTriplets(triplet_list.begin(), triplet_list.end());
	laplacian.makeCompressed();

	Spectra::SparseSymMatProd<float> op(laplacian);
	// Get Fiedler vector
	// Compute second smallest eigenvector
	Spectra::SymEigsSolver<Spectra::SparseSymMatProd<float>> eigs(op, 2, 4);
	eigs.init();

	Eigen::Index num_values = eigs.compute(Spectra::SortRule::SmallestAlge,
		context.max_iterations_eigen, context.error_eigen,
		Spectra::SortRule::LargestAlge);

	if (num_values != 2) {
		std::cerr << "Error: num eigenvalues computed is " << num_values << std::endl;
		return;
	}
	// Get results
	if (eigs.info() != Spectra::CompInfo::Successful) {
		std::cout << "Error: No eigenvalues. Computation not successful!" << std::endl;
		return;
	}

	float eigenvalue = eigs.eigenvalues()[0];
	if (eigenvalue <= 0) {
		std::cerr << "Error: Fiedler eigenvalue is less than 0. Not a connected graph!!" << std::endl;
		std::cerr << "Computed eigenvalues " << eigenvalue << std::endl;

		return;
	}

	const Eigen::VectorXf eigenvectors = eigs.eigenvectors().col(0);


	// Output or continue
	{
		uint32_t size_cluster_0 = 0;
		for (uint32_t i = 0; i < (uint32_t)eigenvectors.size(); ++i) {
			if (eigenvectors[i] < 0.0f) {
				size_cluster_0 += 1;
			}
		}
		std::unordered_set<uint32_t> indices_0, indices_1;
		indices_0.reserve(size_cluster_0);
		indices_1.reserve((uint32_t)eigenvectors.size() - size_cluster_0);

		for (uint32_t i = 0; i < (uint32_t)eigenvectors.size(); ++i) {
			if (eigenvectors[i] < 0.0f) {
				indices_0.insert(new2old_vert[i]);
			}
			else {
				indices_1.insert(new2old_vert[i]);
			}
		}

		if (depth < context.max_depth &&indices_0.size() > context.max_cluster_size) {
			vertex_laplacian_layout(context, depth + 1, faces, vert2face,
				indices_0);
		}
		else {
			const uint32_t id = context.next_id++;
			for (uint32_t idx : indices_0) {
				context.final_cluster[idx] = id;
			}
		}

		if (depth < context.max_depth && indices_1.size() > context.max_cluster_size) {
			vertex_laplacian_layout(context, depth + 1, faces, vert2face,
				indices_1);
		}
		else {
			const uint32_t id = context.next_id++;
			for (uint32_t idx : indices_1) {
				context.final_cluster[idx] = id;
			}
		}
	}

}


std::vector<uint32_t>
get_mapping_optimized_layout(
	const std::vector<Eigen::Array3i>& faces,
	uint32_t num_vertices,
	const MultiLevel multi_level_approach,
	const uint32_t max_depth,
	const uint32_t max_cluster_size,
	const uint32_t max_number_interations_eigen,
	const float eigen_error)

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

	LayoutContext context(max_depth, max_cluster_size, (uint32_t)vert2face.size(), max_number_interations_eigen, eigen_error);
	if (multi_level_approach == MultiLevel::eVertexLaplacian) {
		{
			vertex_laplacian_layout(context, 0, faces, vert2face, vert_indices);
		}
		return context.final_cluster;
	}

	return {};
}

}