#pragma once

#include <Eigen/Dense>
#include <tinyply.h>
#include <vector>
#include <cstdint>


class TriangleMesh {
public:
	TriangleMesh(const char* path);

	void print_debug_info() const;

	void write_mesh_ply(const char* fileName, const std::vector<Eigen::Array3<uint8_t>>& colors = {}) const;

	const std::vector<Eigen::Vector3f>& get_vertices() const {
		return m_vertices;
	}

	const std::vector<Eigen::Array3i>& get_faces() const {
		return m_faces;
	}

	void rearrange_vertices(const std::vector<uint32_t>& old2new);

private:

	void parse_ply(const char* path);

	// Variables
	std::vector<Eigen::Vector3f> m_vertices;
	std::vector<Eigen::Array3i> m_faces;


};

