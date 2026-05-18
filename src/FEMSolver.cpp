#include "FEMSolver.h"

#include <cmath>
#include <fstream>
#include <iostream>

namespace {

constexpr double kGauss = 0.5773502691896262;  // 1/sqrt(3)
constexpr double kTol = 1.0e-9;

}  // namespace

FEMSolver::FEMSolver(const Mesh& mesh) : mesh_(mesh) {}

void FEMSolver::solve() {
    const int n = mesh_.num_nodes();
    stiffness_.assign(n, std::vector<double>(n, 0.0));
    load_.assign(n, 0.0);
    phi_.assign(n, 0.0);

    assemble();
    apply_boundary_conditions();
    solve_linear_system();
}

void FEMSolver::assemble() {
    const int ne = mesh_.num_elements();
    for (int e = 0; e < ne; ++e) {
        double ke[4][4];
        assemble_element(mesh_.element(e), ke);

        const Element& elem = mesh_.element(e);
        for (int i = 0; i < 4; ++i) {
            for (int j = 0; j < 4; ++j) {
                stiffness_[elem.n[i]][elem.n[j]] += ke[i][j];
            }
        }
    }
}

void FEMSolver::assemble_element(const Element& elem, double ke[4][4]) const {
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            ke[i][j] = 0.0;
        }
    }

    const double gp[2] = {-kGauss, kGauss};
    for (int a = 0; a < 2; ++a) {
        for (int b = 0; b < 2; ++b) {
            const double xi = gp[a];
            const double eta = gp[b];

            double dN_dxi[4];
            double dN_deta[4];
            for (int i = 0; i < 4; ++i) {
                shape_derivatives(i, xi, eta, dN_dxi[i], dN_deta[i]);
            }

            double J11 = 0.0, J12 = 0.0, J21 = 0.0, J22 = 0.0;
            for (int i = 0; i < 4; ++i) {
                const Node& nd = mesh_.node(elem.n[i]);
                J11 += dN_dxi[i] * nd.x;
                J12 += dN_dxi[i] * nd.y;
                J21 += dN_deta[i] * nd.x;
                J22 += dN_deta[i] * nd.y;
            }

            const double detJ = J11 * J22 - J12 * J21;
            const double inv_det = 1.0 / detJ;
            const double dxi_dx = J22 * inv_det;
            const double dxi_dy = -J12 * inv_det;
            const double deta_dx = -J21 * inv_det;
            const double deta_dy = J11 * inv_det;

            double dN_dx[4];
            double dN_dy[4];
            for (int i = 0; i < 4; ++i) {
                dN_dx[i] = dN_dxi[i] * dxi_dx + dN_deta[i] * deta_dx;
                dN_dy[i] = dN_dxi[i] * dxi_dy + dN_deta[i] * deta_dy;
            }

            const double weight = std::abs(detJ);
            for (int i = 0; i < 4; ++i) {
                for (int j = 0; j < 4; ++j) {
                    ke[i][j] += weight * (dN_dx[i] * dN_dx[j] + dN_dy[i] * dN_dy[j]);
                }
            }
        }
    }
}

bool FEMSolver::is_bottom(int node) const {
    return std::abs(mesh_.node(node).y) < kTol;
}

bool FEMSolver::is_top(int node) const {
    return std::abs(mesh_.node(node).y - 1.0) < kTol;
}

void FEMSolver::apply_boundary_conditions() {
    const int n = mesh_.num_nodes();
    for (int i = 0; i < n; ++i) {
        double value = 0.0;
        bool fixed = false;

        if (is_bottom(i)) {
            value = 0.0;
            fixed = true;
        } else if (is_top(i)) {
            value = 1.0;
            fixed = true;
        }

        if (!fixed) {
            continue;
        }

        for (int j = 0; j < n; ++j) {
            if (j != i) {
                load_[j] -= stiffness_[j][i] * value;
                stiffness_[j][i] = 0.0;
                stiffness_[i][j] = 0.0;
            }
        }
        stiffness_[i][i] = 1.0;
        load_[i] = value;
    }
}

void FEMSolver::solve_linear_system() {
    const int n = mesh_.num_nodes();
    std::vector<std::vector<double>> a = stiffness_;
    std::vector<double> b = load_;

    for (int col = 0; col < n; ++col) {
        int pivot = col;
        double max_val = std::abs(a[col][col]);
        for (int row = col + 1; row < n; ++row) {
            const double v = std::abs(a[row][col]);
            if (v > max_val) {
                max_val = v;
                pivot = row;
            }
        }

        if (pivot != col) {
            std::swap(a[pivot], a[col]);
            std::swap(b[pivot], b[col]);
        }

        const double diag = a[col][col];
        if (std::abs(diag) < 1.0e-14) {
            continue;
        }

        for (int row = col + 1; row < n; ++row) {
            const double factor = a[row][col] / diag;
            for (int k = col; k < n; ++k) {
                a[row][k] -= factor * a[col][k];
            }
            b[row] -= factor * b[col];
        }
    }

    for (int i = n - 1; i >= 0; --i) {
        double sum = b[i];
        for (int j = i + 1; j < n; ++j) {
            sum -= a[i][j] * phi_[j];
        }
        phi_[i] = (std::abs(a[i][i]) > 1.0e-14) ? sum / a[i][i] : sum;
    }
}

void FEMSolver::write_result(const std::string& path) const {
    std::ofstream out(path);
    if (!out) {
        std::cerr << "Failed to open output file: " << path << std::endl;
        return;
    }

    out << "# x y phi\n";
    for (int i = 0; i < mesh_.num_nodes(); ++i) {
        const Node& nd = mesh_.node(i);
        out << nd.x << " " << nd.y << " " << phi_[i] << "\n";
    }
}

void FEMSolver::shape_derivatives(int i, double xi, double eta, double& dN_dxi, double& dN_deta) {
    switch (i) {
        case 0:
            dN_dxi = -0.25 * (1.0 - eta);
            dN_deta = -0.25 * (1.0 - xi);
            break;
        case 1:
            dN_dxi = 0.25 * (1.0 - eta);
            dN_deta = -0.25 * (1.0 + xi);
            break;
        case 2:
            dN_dxi = 0.25 * (1.0 + eta);
            dN_deta = 0.25 * (1.0 + xi);
            break;
        case 3:
            dN_dxi = -0.25 * (1.0 + eta);
            dN_deta = 0.25 * (1.0 - xi);
            break;
        default:
            dN_dxi = 0.0;
            dN_deta = 0.0;
            break;
    }
}
