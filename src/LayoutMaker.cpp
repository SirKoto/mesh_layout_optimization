#include "LayoutMaker.hpp"

#include <Eigen/SparseCore>
#include <Spectra/SymEigsSolver.h>
#include <Spectra/MatOp/SparseSymMatProd.h>
#include <unordered_set>
#include <unordered_map>
#include <iostream>
#include <numeric>
#include <stack>

namespace LayoutMaker {

struct LayoutContext {

	uint32_t next_id;
	const std::shared_ptr<TriangleMesh> p_mesh;
	const uint32_t max_depth;
	const uint32_t max_cluster_size;
	const uint32_t max_spectral_size;
	std::vector<uint32_t> final_cluster;
	const uint32_t max_iterations_eigen;
	const float error_eigen;

	LayoutContext(
		const std::shared_ptr<TriangleMesh> p_mesh,
		uint32_t max_depth,
		uint32_t max_cluster_size,
		uint32_t max_spectral_size,
		uint32_t max_iterations_eigen,
		float error_eigen) :
		next_id(0), p_mesh(p_mesh), max_depth(max_depth),
		max_cluster_size(max_cluster_size),
		max_spectral_size(max_spectral_size),
		max_iterations_eigen(max_iterations_eigen),
		error_eigen(error_eigen)
	{
		final_cluster.resize(p_mesh->get_vertices().size(), 0);
	}
};

void vertex_laplacian_layout(
	LayoutContext& context,
	const uint32_t depth,
	const std::unordered_multimap<uint32_t, uint32_t>& vert2face,
	const std::unordered_set<uint32_t>& vertices_indices) {

	
	// Termination if conditions fulfilled
	if (vertices_indices.empty()) {
		return;
	}
	if (depth >= context.max_depth || vertices_indices.size() <= context.max_cluster_size) {
		const uint32_t id = context.next_id++;
		for (uint32_t idx : vertices_indices) {
			context.final_cluster[idx] = id;
		}
		return;
	}

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
				const Eigen::Array3i& face = context.p_mesh->get_faces().at(f_id.second);
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

		vertex_laplacian_layout(context, depth + 1, vert2face,
				indices_0);
		
		vertex_laplacian_layout(context, depth + 1, vert2face,
				indices_1);
		
	}

}



void vertex_clustering_layout(
	LayoutContext& context,
	const std::unordered_multimap<uint32_t, uint32_t>& vert2face) {
	if (context.p_mesh->get_vertices().empty()) {
		return;
	}

	const std::vector<Eigen::Vector3f>& vertices_mesh = context.p_mesh->get_vertices();

	struct OctNodeTask {
		std::vector<uint32_t> vertices;
		uint32_t depth;
		Eigen::Vector3f mid_coord;
	};

	Eigen::Vector3f minBBox = Eigen::Vector3f::Constant( std::numeric_limits<float>::infinity());
	Eigen::Vector3f maxBBox = Eigen::Vector3f::Constant(-std::numeric_limits<float>::infinity());

	for (const Eigen::Vector3f& v : context.p_mesh->get_vertices()) {
		minBBox = minBBox.cwiseMin(v);
		maxBBox = maxBBox.cwiseMax(v);
	}

	const float octree_size = (maxBBox - minBBox).maxCoeff();

	std::stack<OctNodeTask> tasks;
	std::array<std::vector<uint32_t>, 8> child_verts;

	std::unordered_set<uint32_t> vert_indices_spectral;
	vert_indices_spectral.reserve(vertices_mesh.size() * 3 / 2);
	

	// Create root node
	{
		OctNodeTask root;
		root.vertices.resize(context.p_mesh->get_vertices().size());
		std::iota(root.vertices.begin(), root.vertices.end(), 0);
		root.depth = 0;
		root.mid_coord = (maxBBox + minBBox) * 0.5f;
		tasks.push(std::move(root));
	}

	// process octree
	while (!tasks.empty()) {
		const OctNodeTask task = std::move(tasks.top());
		tasks.pop();

		const float size_node = octree_size / static_cast<float>(1 << task.depth);

		// Clear childs
		for (auto& v : child_verts)  v.clear();

		// Classify vertices into the 8 childs
		for (const uint32_t& i : task.vertices) {
			const Eigen::Vector3f dir = vertices_mesh[i] - task.mid_coord;
			uint32_t k =
				((dir.x() >= 0.f ? 1 : 0) << 0) +
				((dir.y() >= 0.f ? 1 : 0) << 1) +
				((dir.z() >= 0.f ? 1 : 0) << 2);

			child_verts[k].push_back(i);
		}

		// Keep generating tasks or cluster
		for (uint32_t k = 0; k < 8; ++k) {
			if (child_verts[k].size() < context.max_spectral_size) {
				// Reuse set
				vert_indices_spectral.clear();
				vert_indices_spectral.insert(child_verts[k].begin(), child_verts[k].end());
				// Spectral classification
				vertex_laplacian_layout(
					context,
					0, // depth
					vert2face,
					vert_indices_spectral
				);
			}
			else {
				Eigen::Vector3f dir = { k & 0b1 ? 1.f : -1.f, k & 0b10 ? 1.f : -1.f, k & 0b100 ? 1.f : -1.f };
				OctNodeTask newT;
				newT.depth = task.depth + 1;
				newT.mid_coord = task.mid_coord + 0.25f * size_node * dir;
				newT.vertices = child_verts[k];
				tasks.push(std::move(newT));
			}
		}
		
	}

}

std::vector<uint32_t>
get_mapping_optimized_layout(
	const std::shared_ptr<TriangleMesh> p_mesh,
	const uint32_t max_depth,
	const uint32_t max_cluster_size,
	const uint32_t max_spectral_size,
	const uint32_t max_number_interations_eigen,
	const float eigen_error)

{

	uint32_t num_vertices = (uint32_t) p_mesh->get_vertices().size();

	std::unordered_set<uint32_t> vert_indices;
	vert_indices.reserve(num_vertices);
	for (uint32_t i = 0; i < num_vertices; ++i) {
		vert_indices.insert(i);
	}

	std::unordered_multimap<uint32_t, uint32_t> vert2face;
	vert2face.reserve(num_vertices);
	for (uint32_t f = 0; f < (uint32_t)p_mesh->get_faces().size(); ++f) {
		for (uint32_t j = 0; j < 3; ++j) {
			vert2face.insert({ p_mesh->get_faces()[f][j], f});
		}
	}

	LayoutContext context(p_mesh, max_depth, max_cluster_size, max_spectral_size,
		max_number_interations_eigen, eigen_error);
		
	vertex_clustering_layout(context, vert2face);
		
	return context.final_cluster;
}

}