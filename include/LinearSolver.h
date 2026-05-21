#ifndef LINEAR_SOLVER_H
#define LINEAR_SOLVER_H

#include <vector>

class LinearSolver {
public:
    static std::vector<double> solve(std::vector<std::vector<double>> A, std::vector<double> b);
};

#endif
