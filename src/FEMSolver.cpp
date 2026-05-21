#include "FEMSolver.h"

#include <fstream>
#include <iostream>

FEMSolver::FEMSolver(const Mesh& m)
    : mesh(m)
{
    int n = mesh.nodeCount;

    K = Eigen::MatrixXd::Zero(n, n);

    F = Eigen::VectorXd::Zero(n);

    U = Eigen::VectorXd::Zero(n);
}

void FEMSolver::assemble()
{
    std::cout << "assemble" << std::endl;
}

void FEMSolver::applyBoundaryCondition()
{
    std::cout << "boundary condition" << std::endl;
}

void FEMSolver::solve()
{
    U = K.colPivHouseholderQr().solve(F);

    std::cout << "solve" << std::endl;
}

void FEMSolver::output()
{
    std::ofstream ofs("../output/result.dat");

    for (int i = 0; i < U.size(); i++) {

        ofs << U(i) << std::endl;
    }

    std::cout << "output" << std::endl;
}