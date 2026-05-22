#ifndef FEM_H
#define FEM_H

#include <Eigen/Dense>
#include <array>
#include <string>
#include <vector>

class FEM {
public:
    FEM();

    void readMesh(const std::string& nodeFile,
                  const std::string& elementFile);

    void assemble();

    void applyBoundaryCondition();

    void solve();

    void writeResult(const std::string& filename) const;

private:
    // 節点座標 nodes[i] = (x_i, y_i)
    std::vector<Eigen::Vector2d> nodes_;

    // 四角形要素 elements_[e] = {n1, n2, n3, n4}
    std::vector<std::array<int, 4>> elements_;

    // 全体剛性行列
    Eigen::MatrixXd K_;

    // 右辺ベクトル
    Eigen::VectorXd f_;

    // 未知量 phi
    Eigen::VectorXd phi_;

    // 境界判定
    double ymin_;
    double ymax_;

private:
    void shapeFunctionDerivative(double s,
                                 double t,
                                 Eigen::Vector4d& dN_ds,
                                 Eigen::Vector4d& dN_dt) const;
};

#endif