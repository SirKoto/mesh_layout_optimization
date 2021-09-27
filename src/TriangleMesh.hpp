#pragma once

#include <Eigen/Dense>
#include <tinyply.h>
#include <vector>
#include <cstdint>


class TriangleMesh {
public:
	TriangleMesh(const char* path);

	void print_debug_info() const;

	void write_mesh_ply(const char* fileName) const;

private:

	void parse_ply(const char* path);

	// Variables
	std::vector<Eigen::Vector3f> m_vertices;
	std::vector<Eigen::Array3i> m_faces;


};

