#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <array>
#include <cmath>
#include <cstdlib>
#include <algorithm>
#include <string>

class LaplaceFEM2D {
private:
    struct Node {
        double x;
        double y;
    };

    struct Triangle {
        std::array<int, 3> node;
    };

    std::vector<Node> nodes;
    std::vector<Triangle> triangles;

    std::vector<std::vector<double>> K;
    std::vector<double> F;
    std::vector<double> phi;

public:
    void solve(const std::string& nodeFile,
               const std::string& elementFile)
    {
        readNodes(nodeFile);
        readElements(elementFile);

        int n = static_cast<int>(nodes.size());

        K.assign(n, std::vector<double>(n, 0.0));
        F.assign(n, 0.0);
        phi.assign(n, 0.0);

        assembleStiffnessMatrix();
        applyDirichletBoundaryCondition();
        solveLinearSystem();
    }

    void writeResult(const std::string& filename) const
    {
        std::ofstream fout(filename);

        if (!fout) {
            std::cerr << "Error: cannot open " << filename << std::endl;
            std::exit(1);
        }

        for (size_t i = 0; i < nodes.size(); i++) {
            fout << nodes[i].x << " "
                 << nodes[i].y << " "
                 << phi[i] << "\n";
        }
    }

private:
    std::vector<double> parseNumbers(const std::string& line) const
    {
        std::stringstream ss(line);
        std::vector<double> values;
        double v;

        while (ss >> v) {
            values.push_back(v);
        }

        return values;
    }

    void readNodes(const std::string& filename)
    {
        std::ifstream fin(filename);

        if (!fin) {
            std::cerr << "Error: cannot open " << filename << std::endl;
            std::exit(1);
        }

        std::string line;

        while (std::getline(fin, line)) {
            if (line.empty()) continue;
            if (line[0] == '#') continue;

            std::vector<double> values = parseNumbers(line);

            if (values.size() < 2) continue;

            Node nd;

            // node.dat が
            // x y
            // または
            // id x y
            // のどちらでも読めるように、最後の2列を座標として使う
            nd.x = values[values.size() - 2];
            nd.y = values[values.size() - 1];

            nodes.push_back(nd);
        }

        if (nodes.empty()) {
            std::cerr << "Error: no nodes were read." << std::endl;
            std::exit(1);
        }

        std::cout << "Number of nodes = " << nodes.size() << std::endl;
    }

    void readElements(const std::string& filename)
    {
        std::ifstream fin(filename);

        if (!fin) {
            std::cerr << "Error: cannot open " << filename << std::endl;
            std::exit(1);
        }

        std::string line;
        std::vector<std::array<int, 4>> quads;

        while (std::getline(fin, line)) {
            if (line.empty()) continue;
            if (line[0] == '#') continue;

            std::vector<double> valuesDouble = parseNumbers(line);

            if (valuesDouble.size() < 4) continue;

            std::vector<int> values;
            for (double v : valuesDouble) {
                values.push_back(static_cast<int>(v));
            }

            std::array<int, 4> q;

            // element.dat が
            // n1 n2 n3 n4
            // または
            // id n1 n2 n3 n4
            // のどちらでも読めるように、最後の4列を節点番号として使う
            int m = static_cast<int>(values.size());
            q[0] = values[m - 4];
            q[1] = values[m - 3];
            q[2] = values[m - 2];
            q[3] = values[m - 1];

            quads.push_back(q);
        }

        if (quads.empty()) {
            std::cerr << "Error: no elements were read." << std::endl;
            std::exit(1);
        }

        // 節点番号が 1 始まりなら 0 始まりに直す
        int minIndex = quads[0][0];
        for (const auto& q : quads) {
            for (int i = 0; i < 4; i++) {
                minIndex = std::min(minIndex, q[i]);
            }
        }

        if (minIndex == 1) {
            for (auto& q : quads) {
                for (int i = 0; i < 4; i++) {
                    q[i]--;
                }
            }
        }

        // 四角形要素を三角形2つに分割
        for (const auto& q : quads) {
            Triangle t1;
            Triangle t2;

            t1.node = {q[0], q[1], q[2]};
            t2.node = {q[0], q[2], q[3]};

            triangles.push_back(t1);
            triangles.push_back(t2);
        }

        std::cout << "Number of quadrilateral elements = "
                  << quads.size() << std::endl;
        std::cout << "Number of triangular elements = "
                  << triangles.size() << std::endl;
    }

    void assembleStiffnessMatrix()
    {
        for (const auto& tri : triangles) {
            int n1 = tri.node[0];
            int n2 = tri.node[1];
            int n3 = tri.node[2];

            if (n1 < 0 || n2 < 0 || n3 < 0 ||
                n1 >= static_cast<int>(nodes.size()) ||
                n2 >= static_cast<int>(nodes.size()) ||
                n3 >= static_cast<int>(nodes.size())) {
                std::cerr << "Error: element has invalid node index."
                          << std::endl;
                std::exit(1);
            }

            double x1 = nodes[n1].x;
            double y1 = nodes[n1].y;
            double x2 = nodes[n2].x;
            double y2 = nodes[n2].y;
            double x3 = nodes[n3].x;
            double y3 = nodes[n3].y;

            double area = 0.5 * std::fabs(
                (x2 - x1) * (y3 - y1)
              - (x3 - x1) * (y2 - y1)
            );

            if (area < 1.0e-14) {
                std::cerr << "Error: element area is zero."
                          << std::endl;
                std::exit(1);
            }

            std::array<double, 3> b = {
                y2 - y3,
                y3 - y1,
                y1 - y2
            };

            std::array<double, 3> c = {
                x3 - x2,
                x1 - x3,
                x2 - x1
            };

            std::array<int, 3> id = {n1, n2, n3};

            for (int i = 0; i < 3; i++) {
                for (int j = 0; j < 3; j++) {
                    double ke =
                        (b[i] * b[j] + c[i] * c[j])
                        / (4.0 * area);

                    K[id[i]][id[j]] += ke;
                }
            }
        }
    }

    void applyDirichletBoundaryCondition()
    {
        double ymin = nodes[0].y;
        double ymax = nodes[0].y;

        for (const auto& nd : nodes) {
            ymin = std::min(ymin, nd.y);
            ymax = std::max(ymax, nd.y);
        }

        double eps = 1.0e-10;

        for (int i = 0; i < static_cast<int>(nodes.size()); i++) {
            bool isBoundary = false;
            double value = 0.0;

            // 下底 phi = 0
            if (std::fabs(nodes[i].y - ymin) < eps) {
                isBoundary = true;
                value = 0.0;
            }

            // 上底 phi = 1
            if (std::fabs(nodes[i].y - ymax) < eps) {
                isBoundary = true;
                value = 1.0;
            }

            if (isBoundary) {
                // 対称性を保つため、既知値の寄与を右辺へ移す
                for (int r = 0; r < static_cast<int>(nodes.size()); r++) {
                    F[r] -= K[r][i] * value;
                }

                for (int j = 0; j < static_cast<int>(nodes.size()); j++) {
                    K[i][j] = 0.0;
                    K[j][i] = 0.0;
                }

                K[i][i] = 1.0;
                F[i] = value;
            }
        }
    }

    void solveLinearSystem()
    {
        int n = static_cast<int>(nodes.size());

        for (int k = 0; k < n; k++) {
            int pivotRow = k;
            double maxAbs = std::fabs(K[k][k]);

            for (int i = k + 1; i < n; i++) {
                if (std::fabs(K[i][k]) > maxAbs) {
                    maxAbs = std::fabs(K[i][k]);
                    pivotRow = i;
                }
            }

            if (maxAbs < 1.0e-14) {
                std::cerr << "Error: zero pivot at row "
                          << k << std::endl;
                std::exit(1);
            }

            if (pivotRow != k) {
                std::swap(K[k], K[pivotRow]);
                std::swap(F[k], F[pivotRow]);
            }

            double pivot = K[k][k];

            for (int j = k; j < n; j++) {
                K[k][j] /= pivot;
            }
            F[k] /= pivot;

            for (int i = 0; i < n; i++) {
                if (i == k) continue;

                double factor = K[i][k];

                for (int j = k; j < n; j++) {
                    K[i][j] -= factor * K[k][j];
                }

                F[i] -= factor * F[k];
            }
        }

        phi = F;
    }
};

int main()
{
    LaplaceFEM2D solver;

    solver.solve("../input/node.dat", "../input/element.dat");
    solver.writeResult("../result/laplace_result.dat");

    std::cout << "Calculation finished." << std::endl;
    std::cout << "Result is written to ../result/laplace_result.dat"
              << std::endl;

    return 0;
}
