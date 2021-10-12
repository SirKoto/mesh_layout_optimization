#include <iostream>
#include "Args.hpp"
#include "TriangleMesh.hpp"
#include "LayoutMaker.hpp"

void print_usage() {
    std::cout << 
        "./mesh_layout_opt [options=?]\n"
        "\t-in=input mesh path (.ply) [mandatory]\n"
        "\t-out=output mesh path\n"
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
    
    TriangleMesh mesh(in.c_str());

    mesh.print_debug_info();

    LayoutMaker::get_mapping_optimized_layout(
        mesh.get_faces(),
        (uint32_t)mesh.get_vertices().size(),
        LayoutMaker::MultiLevel::eVertexLaplacian);

    mesh.write_mesh_ply("out.ply");
}