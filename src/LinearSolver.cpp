#include "LinearSolver.h"

#include <stdexcept>

Eigen::VectorXd LinearSolver::solve(const Eigen::MatrixXd& A, const Eigen::VectorXd& b) {
    if (A.rows() != A.cols() || A.rows() != b.size()) {
        throw std::runtime_error("Linear system size mismatch.");
    }

    const Eigen::PartialPivLU<Eigen::MatrixXd> lu(A);
    return lu.solve(b);
}
