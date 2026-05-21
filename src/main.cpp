#include "Mesh.h"
#include "FEMSolver.h"

int main()
{
    Mesh mesh;

    mesh.loadNode("../input/node.dat");

    mesh.loadElement("../input/element.dat");

    FEMSolver solver(mesh);

    solver.assemble();

    solver.applyBoundaryCondition();

    solver.solve();

    solver.output();

    return 0;
}