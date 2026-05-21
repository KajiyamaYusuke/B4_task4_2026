#include "FEMSolver.h"

#include "LinearSolver.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <stdexcept>

namespace {

const double GAUSS_POINT = 1.0 / std::sqrt(3.0);

void shapeFunctionDerivatives(double xi, double eta, double dNdXi[4], double dNdEta[4]) {
    dNdXi[0] = -(1.0 - eta) / 4.0;
    dNdXi[1] = (1.0 - eta) / 4.0;
    dNdXi[2] = (1.0 + eta) / 4.0;
    dNdXi[3] = -(1.0 + eta) / 4.0;

    dNdEta[0] = -(1.0 - xi) / 4.0;
    dNdEta[1] = -(1.0 + xi) / 4.0;
    dNdEta[2] = (1.0 + xi) / 4.0;
    dNdEta[3] = (1.0 - xi) / 4.0;
}

}  // namespace

FEMSolver::FEMSolver(const Mesh& mesh)
    : mesh_(mesh),
      solution_(Eigen::VectorXd::Zero(mesh.nodeCount())),
      isDirichlet_(mesh.nodeCount(), false),
      dirichletValue_(mesh.nodeCount(), 0.0) {
    setupBoundaryConditions();
}

void FEMSolver::setupBoundaryConditions() {
    double xMin = mesh_.node(0)[0];
    double xMax = mesh_.node(0)[0];
    const int nodeCount = mesh_.nodeCount();
    for (int i = 0; i < nodeCount; ++i) {
        const double x = mesh_.node(i)[0];
        xMin = std::min(xMin, x);
        xMax = std::max(xMax, x);
    }

    for (int i = 0; i < nodeCount; ++i) {
        const double x = mesh_.node(i)[0];
        if (std::abs(x - xMin) < BOUNDARY_TOL) {
            isDirichlet_[i] = true;
            dirichletValue_[i] = LEFT_VALUE;
        } else if (std::abs(x - xMax) < BOUNDARY_TOL) {
            isDirichlet_[i] = true;
            dirichletValue_[i] = RIGHT_VALUE;
        }
    }
}

Eigen::Matrix4d FEMSolver::assembleElementStiffness(int elementIndex) const {
    Eigen::Matrix4d Ke = Eigen::Matrix4d::Zero();

    const std::array<int, 4>& elem = mesh_.element(elementIndex);
    const double gaussPoints[2] = {-GAUSS_POINT, GAUSS_POINT};

    for (double xi : gaussPoints) {
        for (double eta : gaussPoints) {
            double dNdXi[4] = {};
            double dNdEta[4] = {};
            shapeFunctionDerivatives(xi, eta, dNdXi, dNdEta);

            double dxDxi = 0.0;
            double dyDxi = 0.0;
            double dxDeta = 0.0;
            double dyDeta = 0.0;
            for (int a = 0; a < 4; ++a) {
                const double x = mesh_.node(elem[a])[0];
                const double y = mesh_.node(elem[a])[1];
                dxDxi += dNdXi[a] * x;
                dyDxi += dNdXi[a] * y;
                dxDeta += dNdEta[a] * x;
                dyDeta += dNdEta[a] * y;
            }

            const double detJ = dxDxi * dyDeta - dyDxi * dxDeta;
            if (std::abs(detJ) < 1.0e-14) {
                throw std::runtime_error("Non-positive Jacobian determinant.");
            }

            const double invDetJ = 1.0 / detJ;
            Eigen::Vector4d dNdx;
            Eigen::Vector4d dNdy;
            for (int a = 0; a < 4; ++a) {
                dNdx[a] = invDetJ * (dyDeta * dNdXi[a] - dyDxi * dNdEta[a]);
                dNdy[a] = invDetJ * (-dxDeta * dNdXi[a] + dxDxi * dNdEta[a]);
            }

            Ke += (dNdx * dNdx.transpose() + dNdy * dNdy.transpose()) * detJ;
        }
    }

    return Ke;
}

void FEMSolver::assembleGlobalStiffness(Eigen::MatrixXd& K) const {
    const int nodeCount = mesh_.nodeCount();
    K.setZero(nodeCount, nodeCount);

    const int elementCount = mesh_.elementCount();
    for (int e = 0; e < elementCount; ++e) {
        const Eigen::Matrix4d Ke = assembleElementStiffness(e);
        const std::array<int, 4>& elem = mesh_.element(e);

        for (int i = 0; i < 4; ++i) {
            const int globalI = elem[i];
            for (int j = 0; j < 4; ++j) {
                const int globalJ = elem[j];
                K(globalI, globalJ) += Ke(i, j);
            }
        }
    }
}

void FEMSolver::applyDirichletConditions(Eigen::MatrixXd& K, Eigen::VectorXd& b) const {
    const int nodeCount = mesh_.nodeCount();
    for (int i = 0; i < nodeCount; ++i) {
        if (!isDirichlet_[i]) {
            continue;
        }

        const double value = dirichletValue_[i];
        for (int k = 0; k < nodeCount; ++k) {
            if (!isDirichlet_[k]) {
                b[k] -= K(k, i) * value;
            }
        }

        K.row(i).setZero();
        K.col(i).setZero();
        K(i, i) = 1.0;
        b[i] = value;
    }
}

void FEMSolver::solve() {
    const int nodeCount = mesh_.nodeCount();
    Eigen::MatrixXd K;
    Eigen::VectorXd b = Eigen::VectorXd::Zero(nodeCount);

    assembleGlobalStiffness(K);
    applyDirichletConditions(K, b);
    solution_ = LinearSolver::solve(K, b);
}

void FEMSolver::writeSolution(const std::string& filename) const {
    std::ofstream ofs(filename);
    if (!ofs) {
        throw std::runtime_error("Failed to open output file: " + filename);
    }

    ofs << std::scientific << std::setprecision(8);
    const int nodeCount = mesh_.nodeCount();
    for (int i = 0; i < nodeCount; ++i) {
        ofs << i << " "
            << mesh_.node(i)[0] << " "
            << mesh_.node(i)[1] << " "
            << solution_[i] << '\n';
    }
}
