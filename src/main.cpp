#include <iostream>
#include "Args.hpp"
#include "TriangleMesh.hpp"
#include "LayoutMaker.hpp"
#include <chrono>

void print_usage() {
    std::cout << 
        "./mesh_layout_opt [options=?]\n"
        "\t-in=input mesh path (.ply) [mandatory]\n"
        "\t-out=output mesh path\n"
        "\t-max_iterations=int"
        "\t-h or --help to see this information\n"
        << std::endl;
}

int main(int argc, char** argv) {
    Args args(argc, argv);

    if (args.has("h") || args.has("-help")) {
        print_usage();
        return 0;
    }

    std::string in;
    if (!args.has("in")) {
        std::cerr << "Missing *in* parameter" << std::endl;
        print_usage();
        return 1;
    }
    else {
        in = args.get("in");
    }
    std::string out;
    if (!args.has("out")) {
        out = "out.ply";
    }
    else {
        out = args.get("out");
    }
    uint32_t max_number_interations_eigen = 100000;
    if (args.has("max_iterations")) {
        max_number_interations_eigen = (uint16_t)std::stoi(args.get("max_iterations"));
    }
    
    TriangleMesh mesh(in.c_str());

    mesh.print_debug_info();

    const auto ini_timer = std::chrono::high_resolution_clock::now();

    std::vector<uint32_t> clusters =
    LayoutMaker::get_mapping_optimized_layout(
        mesh.get_faces(),
        (uint32_t)mesh.get_vertices().size(),
        LayoutMaker::MultiLevel::eVertexLaplacian,
        max_number_interations_eigen);

    const auto end_timer = std::chrono::high_resolution_clock::now();
    const std::chrono::duration<double> duration = end_timer - ini_timer;

    std::cout << "Clustering took " << duration.count() << " s." << std::endl;

    // find max id
    uint32_t num_colors = 0;
    for (uint32_t i = 0; i < (uint32_t)clusters.size(); ++i) {
        num_colors = std::max(num_colors, clusters[i]);
    }
    
    // Random colors
    std::vector<Eigen::Array3<uint8_t>> color_map(num_colors + 1);
    for (Eigen::Array3<uint8_t>& color : color_map) {
        color.setRandom();
    }

    std::vector<Eigen::Array3<uint8_t>> colors(clusters.size());
    for (uint32_t i = 0; i < (uint32_t)colors.size(); ++i) {
        colors[i] = color_map[clusters[i]];
    }
    mesh.write_mesh_ply("out.ply", colors);
}