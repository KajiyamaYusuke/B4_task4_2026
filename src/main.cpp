#include <iostream>
#include "FEMSolver.h"
#include "Mesh.h"

int main() {
    Mesh mesh;
    if (!mesh.load("../input/node.dat", "../input/element.dat")) {
        std::cerr << "Failed to load mesh data from input/ directory." << std::endl;
        return 1;
    }

    FEMSolver solver(mesh);
    solver.solve();
    solver.write_result("../output/result.dat");

    std::cout << "Laplace FEM solve completed." << std::endl;
    std::cout << "Nodes: " << mesh.num_nodes() << ", Elements: " << mesh.num_elements() << std::endl;
    std::cout << "Result written to output/result.dat" << std::endl;

    return 0;
}
