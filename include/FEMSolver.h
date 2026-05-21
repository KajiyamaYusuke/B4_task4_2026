#ifndef FEM_SOLVER_H
#define FEM_SOLVER_H

#include "Mesh.h"

#include <Eigen/Dense>

#include <string>
#include <vector>

class FEMSolver {
public:
    explicit FEMSolver(const Mesh& mesh);

    void solve();
    void writeSolution(const std::string& filename) const;
    const Eigen::VectorXd& solution() const { return solution_; }

private:
    const Mesh& mesh_;
    Eigen::VectorXd solution_;
    std::vector<bool> isDirichlet_;
    std::vector<double> dirichletValue_;

    static constexpr double LEFT_VALUE = 1.0;
    static constexpr double RIGHT_VALUE = 0.0;
    static constexpr double BOUNDARY_TOL = 1.0e-8;

    void setupBoundaryConditions();
    void assembleGlobalStiffness(Eigen::MatrixXd& K) const;
    Eigen::Matrix4d assembleElementStiffness(int elementIndex) const;
    void applyDirichletConditions(Eigen::MatrixXd& K, Eigen::VectorXd& b) const;
};

#endif
