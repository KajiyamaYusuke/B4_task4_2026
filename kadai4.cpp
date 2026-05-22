#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>
#include <stdexcept>
#include <sys/stat.h>
#if __has_include(<Eigen/Dense>)
#include <Eigen/Dense>
#elif __has_include(<eigen3/Eigen/Dense>)
#include <eigen3/Eigen/Dense>
#else
#error "Eigen is not found. Install Eigen or compile with -I/usr/include/eigen3."
#endif

using namespace std;
using namespace Eigen;

class FEMSolver {
private:
    vector<Vector2d> nodes;
    vector<Vector4i> elements;

    MatrixXd K;
    VectorXd F;
    VectorXd phi;

    int numNodes;
    int numElements;
    string baseDir;

    bool fileExists(const string& filename) {
        ifstream file(filename);
        return file.good();
    }

    void makeDirectory(const string& dirname) {
        mkdir(dirname.c_str(), 0755);
    }

    void applyDirichlet(int nodeIndex, double value) {
        F -= K.col(nodeIndex) * value;

        K.row(nodeIndex).setZero();
        K.col(nodeIndex).setZero();

        K(nodeIndex, nodeIndex) = 1.0;
        F(nodeIndex) = value;
    }

public:
    FEMSolver() : numNodes(0), numElements(0), baseDir("") {}

    void readNodes(const string& filename) {
        ifstream file(filename);
        if (!file) {
            throw runtime_error("Cannot open node file: " + filename);
        }

        double x, y;

        while (file >> x >> y) {
            nodes.push_back(Vector2d(x, y));
        }

        numNodes = nodes.size();
        if (numNodes == 0) {
            throw runtime_error("No node data was read from: " + filename);
        }
    }

    void readElements(const string& filename) {
        ifstream file(filename);
        if (!file) {
            throw runtime_error("Cannot open element file: " + filename);
        }

        int n1, n2, n3, n4;

        while (file >> n1 >> n2 >> n3 >> n4) {
            elements.push_back(Vector4i(n1, n2, n3, n4));
        }

        numElements = elements.size();
        if (numElements == 0) {
            throw runtime_error("No element data was read from: " + filename);
        }
    }

    Vector4d shapeFunction(double s, double t) {
        Vector4d N;

        N(0) = 0.25 * (1.0 - s) * (1.0 - t);
        N(1) = 0.25 * (1.0 + s) * (1.0 - t);
        N(2) = 0.25 * (1.0 + s) * (1.0 + t);
        N(3) = 0.25 * (1.0 - s) * (1.0 + t);

        return N;
    }

    Matrix<double, 2, 4> dNdST(double s, double t) {
        Matrix<double, 2, 4> dN;

        dN(0,0) = -0.25 * (1.0 - t);
        dN(0,1) =  0.25 * (1.0 - t);
        dN(0,2) =  0.25 * (1.0 + t);
        dN(0,3) = -0.25 * (1.0 + t);

        dN(1,0) = -0.25 * (1.0 - s);
        dN(1,1) = -0.25 * (1.0 + s);
        dN(1,2) =  0.25 * (1.0 + s);
        dN(1,3) =  0.25 * (1.0 - s);

        return dN;
    }

    Matrix4d elementStiffness(Vector4i elem) {

        Matrix4d ke = Matrix4d::Zero();

        double gauss[2] = {
            -1.0 / sqrt(3.0),
             1.0 / sqrt(3.0)
        };

        for (int i = 0; i < 2; i++) {
            for (int j = 0; j < 2; j++) {

                double s = gauss[i];
                double t = gauss[j];

                Matrix<double, 2, 4> dNst = dNdST(s, t);

                Matrix2d J = Matrix2d::Zero();

                for (int a = 0; a < 4; a++) {

                    double x = nodes[elem(a)](0);
                    double y = nodes[elem(a)](1);

                    J(0,0) += dNst(0,a) * x;
                    J(0,1) += dNst(0,a) * y;
                    J(1,0) += dNst(1,a) * x;
                    J(1,1) += dNst(1,a) * y;
                }

                double detJ = J.determinant();

                Matrix2d invJ = J.inverse();

                Matrix<double, 2, 4> dNxy;

                dNxy = invJ * dNst;

                ke += (dNxy.transpose() * dNxy) * detJ;
            }
        }

        return ke;
    }

    void assemble() {

        K = MatrixXd::Zero(numNodes, numNodes);
        F = VectorXd::Zero(numNodes);

        for (int e = 0; e < numElements; e++) {

            Vector4i elem = elements[e];

            Matrix4d ke = elementStiffness(elem);

            for (int i = 0; i < 4; i++) {
                for (int j = 0; j < 4; j++) {

                    int I = elem(i);
                    int J = elem(j);

                    if (I < 0 || I >= numNodes || J < 0 || J >= numNodes) {
                        throw runtime_error("Element contains an invalid node index.");
                    }

                    K(I, J) += ke(i, j);
                }
            }
        }
    }

    void applyBoundaryCondition() {

        for (int i = 0; i < numNodes; i++) {

            double x = nodes[i](0);

            // 左端 x=0 → phi=1
            if (fabs(x - 0.0) < 1e-6) {
                applyDirichlet(i, 1.0);
            }

            // 右端 x=5 → phi=0
            if (fabs(x - 5.0) < 1e-6) {
                applyDirichlet(i, 0.0);
            }
        }
    }

    void solve() {

        phi = K.colPivHouseholderQr().solve(F);
    }

    void output(const string& filename) {

        ofstream file(filename);
        if (!file) {
            throw runtime_error("Cannot open output file: " + filename);
        }

        for (int i = 0; i < numNodes; i++) {

            file
                << nodes[i](0) << " "
                << nodes[i](1) << " "
                << phi(i) << endl;
        }

        file.close();
    }

    void run() {
        if (fileExists("input/node.dat") && fileExists("input/element.dat")) {
            baseDir = "";
        } else if (fileExists("../input/node.dat") && fileExists("../input/element.dat")) {
            baseDir = "../";
        } else {
            throw runtime_error("Cannot find input/node.dat and input/element.dat.");
        }

        makeDirectory(baseDir + "output");

        readNodes(baseDir + "input/node.dat");
        readElements(baseDir + "input/element.dat");

        assemble();

        applyBoundaryCondition();

        solve();

        output(baseDir + "output/result.dat");

        cout << "Finished FEM Laplace Solve" << endl;
    }
};

int main() {

    try {
        FEMSolver fem;

        fem.run();
    } catch (const exception& e) {
        cerr << "Error: " << e.what() << endl;
        return 1;
    }

    return 0;
}
