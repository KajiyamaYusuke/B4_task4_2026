#include "FEMSolver.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <stdexcept>

namespace {

constexpr int NodesPerElement = 4;

std::array<double, NodesPerElement> dNdxi(double, double eta)
{
    return {
        -0.25 * (1.0 - eta),
        0.25 * (1.0 - eta),
        0.25 * (1.0 + eta),
        -0.25 * (1.0 + eta),
    };
}

std::array<double, NodesPerElement> dNdeta(double xi, double)
{
    return {
        -0.25 * (1.0 - xi),
        -0.25 * (1.0 + xi),
        0.25 * (1.0 + xi),
        0.25 * (1.0 - xi),
    };
}

} // namespace

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
    const int n = mesh.nodeCount;
    K = Eigen::MatrixXd::Zero(n, n);
    F = Eigen::VectorXd::Zero(n);

    const double a = 1.0 / std::sqrt(3.0);
    const std::array<double, 2> gaussPoints = { -a, a };

    for (const auto& element : mesh.elements) {
        if (element.size() != NodesPerElement) {
            throw std::runtime_error("element must have 4 nodes");
        }

        Eigen::Matrix4d Ke = Eigen::Matrix4d::Zero();

        for (double xi : gaussPoints) {
            for (double eta : gaussPoints) {
                const auto dxi = dNdxi(xi, eta);
                const auto deta = dNdeta(xi, eta);

                Eigen::Matrix2d J = Eigen::Matrix2d::Zero();

                for (int i = 0; i < NodesPerElement; i++) {
                    const int nodeId = element[i];
                    const double x = mesh.nodes[nodeId][0];
                    const double y = mesh.nodes[nodeId][1];

                    J(0, 0) += dxi[i] * x;
                    J(0, 1) += deta[i] * x;
                    J(1, 0) += dxi[i] * y;
                    J(1, 1) += deta[i] * y;
                }

                const double detJ = J.determinant();

                if (detJ <= 0.0) {
                    throw std::runtime_error("element has non-positive Jacobian");
                }

                const Eigen::Matrix2d invJT = J.inverse().transpose();
                std::array<Eigen::Vector2d, NodesPerElement> gradN;

                for (int i = 0; i < NodesPerElement; i++) {
                    gradN[i] = invJT * Eigen::Vector2d(dxi[i], deta[i]);
                }

                for (int i = 0; i < NodesPerElement; i++) {
                    for (int j = 0; j < NodesPerElement; j++) {
                        Ke(i, j) += gradN[i].dot(gradN[j]) * detJ;
                    }
                }
            }
        }

        for (int i = 0; i < NodesPerElement; i++) {
            const int row = element[i];

            for (int j = 0; j < NodesPerElement; j++) {
                const int col = element[j];
                K(row, col) += Ke(i, j);
            }
        }
    }

    std::cout << "assemble" << std::endl;
}

void FEMSolver::applyBoundaryCondition()
{
    double yMin = std::numeric_limits<double>::max();
    double yMax = std::numeric_limits<double>::lowest();

    for (const auto& node : mesh.nodes) {
        yMin = std::min(yMin, node[1]);
        yMax = std::max(yMax, node[1]);
    }

    const double tolerance = std::max(1.0e-12, (yMax - yMin) * 1.0e-8);

    for (int nodeId = 0; nodeId < mesh.nodeCount; nodeId++) {
        const double y = mesh.nodes[nodeId][1];
        double value = 0.0;
        bool isDirichlet = false;

        if (std::abs(y - yMin) <= tolerance) {
            value = 0.0;
            isDirichlet = true;
        } else if (std::abs(y - yMax) <= tolerance) {
            value = 1.0;
            isDirichlet = true;
        }

        if (!isDirichlet) {
            continue;
        }

        F -= K.col(nodeId) * value;
        K.row(nodeId).setZero();
        K.col(nodeId).setZero();
        K(nodeId, nodeId) = 1.0;
        F(nodeId) = value;
    }

    std::cout << "boundary condition" << std::endl;
}

void FEMSolver::solve()
{
    U = K.colPivHouseholderQr().solve(F);

    std::cout << "solve" << std::endl;
}

void FEMSolver::output()
{
    const std::filesystem::path outputPath = "../output/result.dat";
    std::filesystem::create_directories(outputPath.parent_path());

    std::ofstream ofs(outputPath);

    if (!ofs) {
        std::cerr << "output file open error: " << outputPath << std::endl;
        return;
    }

for (int i = 0; i < U.size(); i++) {

    double x = mesh.nodes[i][0];
    double y = mesh.nodes[i][1];

    ofs << x << " "
        << y << " "
        << U(i) << std::endl;
}

    std::cout << "output" << std::endl;
}
