#ifndef LINEAR_SOLVER_H
#define LINEAR_SOLVER_H

#include <Eigen/Dense>

class LinearSolver {
public:
    static Eigen::VectorXd solve(const Eigen::MatrixXd& A, const Eigen::VectorXd& b);
};

#endif
