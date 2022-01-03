#include "LayoutOptimizer.hpp"

#include <unordered_set>

namespace LayoutOptimizer {


struct Edge {
	Edge(uint32_t v0, uint32_t v1) {
		d = std::minmax(v0, v1);
	}

	bool operator==(const Edge& o) const {
		return d.first == o.d.first && d.second == o.d.second;
	}

	std::pair<uint32_t, uint32_t> d;
};

struct EdgeHash {
	size_t operator()(const Edge& e) const
	{
		size_t hash1 = std::hash<uint32_t>{}(e.d.first);
		size_t hash2 = std::hash<uint32_t>{}(e.d.second);
		return hash1 ^ hash2;
	}
};

std::vector<uint32_t> optimize_layout(const TriangleMesh& mesh, const std::vector<uint32_t>& clusters)
{
	assert(!clusters.empty());
	assert(mesh.get_vertices().size() == clusters.size());

	std::vector<uint32_t> new_layout(clusters.size(), std::numeric_limits<uint32_t>::max());

	const int32_t num_clusters = 1 + *std::max_element(clusters.begin(), clusters.end());

	// Cluster_to_vert has for each cluster a sorted list of all the vertices in the cluster
	std::vector<std::vector<uint32_t>> cluster_to_vert(num_clusters);
	for (uint32_t i = 0; i < (uint32_t)clusters.size(); ++i) {
		cluster_to_vert[clusters[i]].push_back(i);
	}

	// Create set with all the edges
	std::unordered_set<Edge, EdgeHash> edges_set;
	edges_set.reserve((mesh.get_faces().size() * 3) / 2);
	for (const Eigen::Array3i& face : mesh.get_faces()) {
		for (uint32_t i = 0; i < 3; ++i) {
			Edge edge(face[i], face[(i + 1) % 3]);
			edges_set.insert(edge);
		}
	}

	std::vector<uint32_t> offsets(num_clusters, 0);
	for (uint32_t i = 1; i < (uint32_t)num_clusters; ++i) {
		offsets[i] = (uint32_t)cluster_to_vert[i - 1].size() + offsets[i - 1];
	}

	std::vector<uint32_t> tmp;

#pragma omp parallel for schedule(dynamic) firstprivate(tmp)
	for (int32_t c = 0; c < num_clusters; ++c) {
		std::vector<uint32_t>& cluster = cluster_to_vert[c];
		const uint32_t cluster_size = (uint32_t)cluster.size();
		tmp.clear();

		int32_t edge_span = std::numeric_limits<int32_t>::max();
		do {
			int32_t new_edge_span = 0;
			for (uint32_t i = 0; i < cluster_size; ++i) {
				for (uint32_t j = i + 1; j < cluster_size; ++j) {
					// if edge is contained
					if (edges_set.count(Edge(cluster[i], cluster[j]))) {
						new_edge_span += j - i;
					}
				}
			}

			if (new_edge_span < edge_span) {
				edge_span = new_edge_span;
				tmp = cluster;
			}

		} while (std::next_permutation(cluster.begin(), cluster.end()));

		std::copy(tmp.begin(), tmp.end(), new_layout.data() + offsets[c]);

		assert((uint32_t)tmp.size() == cluster_size);
	}

	return new_layout;
}

} // namespace