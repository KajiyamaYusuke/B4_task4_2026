#pragma once

#include "Mesh.h"

#include <Eigen/Dense>

class FEMSolver {
public:
    explicit FEMSolver(const Mesh& m);

    void assemble();
    void applyBoundaryCondition();
    void solve();
    void output();

private:
    const Mesh& mesh;
    Eigen::MatrixXd K;
    Eigen::VectorXd F;
    Eigen::VectorXd U;
};
