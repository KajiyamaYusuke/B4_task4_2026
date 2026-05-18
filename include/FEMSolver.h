#ifndef FEM_SOLVER_H
#define FEM_SOLVER_H

#include <string>
#include <vector>
#include "Mesh.h"

class FEMSolver {
public:
    explicit FEMSolver(const Mesh& mesh);

    void solve();
    void write_result(const std::string& path) const;
    const std::vector<double>& solution() const { return phi_; }

private:
    const Mesh& mesh_;
    std::vector<double> phi_;
    std::vector<std::vector<double>> stiffness_;
    std::vector<double> load_;

    void assemble();
    void apply_boundary_conditions();
    void solve_linear_system();
    void assemble_element(const Element& elem, double ke[4][4]) const;
    bool is_bottom(int node) const;
    bool is_top(int node) const;
    static void shape_derivatives(int i, double xi, double eta, double& dN_dxi, double& dN_deta);
};

#endif
