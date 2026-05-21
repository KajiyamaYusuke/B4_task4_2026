#ifndef FEM_SOLVER_H
#define FEM_SOLVER_H

#include "Mesh.h"

#include <string>
#include <vector>

class FEMSolver {
public:
    explicit FEMSolver(const Mesh& mesh);

    void solve();
    void writeSolution(const std::string& filename) const;
    const std::vector<double>& solution() const { return solution_; }

private:
    const Mesh& mesh_;
    std::vector<double> solution_;
    std::vector<bool> isDirichlet_;
    std::vector<double> dirichletValue_;

    static constexpr double LEFT_VALUE = 1.0;
    static constexpr double RIGHT_VALUE = 0.0;
    static constexpr double BOUNDARY_TOL = 1.0e-8;

    void setupBoundaryConditions();
    void assembleGlobalStiffness(std::vector<std::vector<double>>& K) const;
    void assembleElementStiffness(int elementIndex, double Ke[4][4]) const;
    void applyDirichletConditions(std::vector<std::vector<double>>& K, std::vector<double>& b) const;
};

#endif
