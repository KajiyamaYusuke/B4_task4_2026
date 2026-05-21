#include "LinearSolver.h"

#include <cmath>
#include <stdexcept>

std::vector<double> LinearSolver::solve(std::vector<std::vector<double>> A, std::vector<double> b) {
    const int n = static_cast<int>(b.size());
    if (n == 0 || static_cast<int>(A.size()) != n) {
        throw std::runtime_error("Linear system size mismatch.");
    }

    for (int k = 0; k < n; ++k) {
        int pivot = k;
        double maxAbs = std::abs(A[k][k]);
        for (int i = k + 1; i < n; ++i) {
            const double value = std::abs(A[i][k]);
            if (value > maxAbs) {
                maxAbs = value;
                pivot = i;
            }
        }

        if (maxAbs < 1.0e-14) {
            throw std::runtime_error("Singular matrix encountered during solve.");
        }

        if (pivot != k) {
            std::swap(A[k], A[pivot]);
            std::swap(b[k], b[pivot]);
        }

        for (int i = k + 1; i < n; ++i) {
            const double factor = A[i][k] / A[k][k];
            for (int j = k; j < n; ++j) {
                A[i][j] -= factor * A[k][j];
            }
            b[i] -= factor * b[k];
        }
    }

    std::vector<double> x(n, 0.0);
    for (int i = n - 1; i >= 0; --i) {
        double sum = b[i];
        for (int j = i + 1; j < n; ++j) {
            sum -= A[i][j] * x[j];
        }
        x[i] = sum / A[i][i];
    }

    return x;
}
