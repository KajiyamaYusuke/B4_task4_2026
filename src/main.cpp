#include <iostream>
#include "FEM.h"

int main()
{
    FEM fem;

    std::cout << "Read mesh" << std::endl;
    fem.readMesh("../input/node.dat", "../input/element.dat");

    std::cout << "Assemble matrix" << std::endl;
    fem.assemble();

    std::cout << "Apply boundary condition" << std::endl;
    fem.applyBoundaryCondition();

    std::cout << "Solve" << std::endl;
    fem.solve();

    std::cout << "Write result" << std::endl;
    fem.writeResult("result.dat");

    std::cout << "Finish" << std::endl;

    return 0;
}