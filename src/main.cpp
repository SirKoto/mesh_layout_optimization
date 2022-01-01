#undef NDEBUG

#include <iostream>
#include "Args.hpp"
#include "TriangleMesh.hpp"
#include "LayoutMaker.hpp"
#include <chrono>

void print_usage() {
    std::cout << 
        "./mesh_layout_opt [options=?]\n"
        "\t-in=input mesh path (.ply) [mandatory]\n"
        "\t-mode=int [default=0]\n"
        "\t\t0: generate mesh with patches\n"
        "\t\t1: optimise mesh layout\n"
        "\t-out=output mesh path\n"
        "\t-max_iterations=int [default=100000]\n"
        "\t-error=float [default=1.0e-5]\n"
        "\t-max_deph=int [default=10]\n"
        "\t-max_cluster_size=int [default=100]\n"
        "\t-max_spectral_size=int [default=100000]\n"
        "\t-h or --help to see this information\n"
        << std::endl;
}

void assess_clustering_quality(const TriangleMesh &mesh, const std::vector<uint32_t>& clusters) {
    std::map<uint32_t, uint32_t> cluster_sizes;
    std::map<uint32_t, std::vector<uint32_t>> cluster_to_vertices;

    for (size_t i = 0; i < clusters.size(); ++i) {
        uint32_t cluster = clusters[i];
        auto it = cluster_sizes.find(cluster);
        if (it == cluster_sizes.end()) {
            cluster_sizes.insert({ cluster, 1 });
        }
        else {
            it->second += 1;
        }

        cluster_to_vertices[cluster].push_back((uint32_t)i);
    }



    double_t average_cluster_size = 0.0;
    uint32_t min_cluster_size = std::numeric_limits<uint32_t>::max();
    uint32_t max_cluster_size = 0;

    for (const auto& it : cluster_sizes) {
        assert(it.second == (uint32_t)cluster_to_vertices[it.first].size());
        average_cluster_size += (double_t)it.second;
        min_cluster_size = std::min(min_cluster_size, it.second);
        max_cluster_size = std::max(max_cluster_size, it.second);
    }
    average_cluster_size /= (double_t)cluster_sizes.size();

    std::cout << "Num clusters: " << cluster_sizes.size() <<
        "\nClusters quality:\n"
        "\tAverage size: " << average_cluster_size << "\n"
        "\tMin Size: " << min_cluster_size << "\n"
        "\tMax Size: " << max_cluster_size << std::endl;
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
        max_number_interations_eigen = (uint32_t)std::stoi(args.get("max_iterations"));
    }

    uint32_t max_depth = 10;
    if (args.has("max_deph")) {
        max_depth = (uint32_t)std::stoi(args.get("max_deph"));
    }

    uint32_t max_cluster_size = 100;
    if (args.has("max_cluster_size")) {
        max_cluster_size = (uint32_t)std::stoi(args.get("max_cluster_size"));
    }
    uint32_t max_spectral_size = 10000;
    if (args.has("max_spectral_size")) {
        max_spectral_size = (uint32_t)std::stoi(args.get("max_spectral_size"));
    }

    int32_t mode = 0;
    if (args.has("mode")) {
        mode = std::stoi(args.get("mode"));
        if (mode < 0 || mode > 1) {
            print_usage();
            return 1;
        }
    }

    float error = 1.0e-5f;
    if (args.has("error")) {
        error = std::stof(args.get("error"));
    }
    
    std::shared_ptr<TriangleMesh> mesh = std::make_shared<TriangleMesh>(in.c_str());

    mesh->print_debug_info();

    std::cout << "Starting clustering:\n"
        "\tMax depth: " << max_depth << "\n"
        "\tMax Cluster size: " << max_cluster_size << "\n"
        "\tMax Spectral size: " << max_spectral_size << "\n"
        "\tMax iterations Eigen: " << max_number_interations_eigen << "\n"
        "\tError: " << error << std::endl;

    const auto ini_timer = std::chrono::high_resolution_clock::now();

    std::vector<uint32_t> clusters =
    LayoutMaker::get_mapping_optimized_layout(
        mesh,
        max_depth, max_cluster_size, max_spectral_size,
        max_number_interations_eigen, error);

    const auto end_timer = std::chrono::high_resolution_clock::now();
    const std::chrono::duration<double> duration = end_timer - ini_timer;

    std::cout << "Clustering took " << duration.count() << " s." << std::endl;

    assert(clusters.size() == mesh->get_vertices().size());

    assess_clustering_quality(*mesh, clusters);


    if (mode == 0) {
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
        mesh->write_mesh_ply("out.ply", colors);
    }
}