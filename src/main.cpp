#include <iostream>
#include <stdexcept>

#include "FEMSolver.h"
#include "Mesh.h"

int main() {
    try {
        Mesh mesh;
        mesh.loadFromFiles("../input/node.dat", "../input/element.dat");

        FEMSolver solver(mesh);
        solver.solve();
        solver.writeSolution("../output/result.dat");

        std::cout << "Laplace FEM solve completed.\n";
        std::cout << "Nodes: " << mesh.nodeCount() << "\n";
        std::cout << "Elements: " << mesh.elementCount() << "\n";
        std::cout << "Result written to ../output/result.dat\n";
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << '\n';
        return 1;
    }

    return 0;
}
