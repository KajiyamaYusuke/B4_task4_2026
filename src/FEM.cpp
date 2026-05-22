#include "FEM.h"

#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <limits>

FEM::FEM()
    : ymin_(0.0), ymax_(0.0)
{
}

void FEM::readMesh(const std::string& nodeFile, const std::string& elementFile)
{
    nodes_.clear();
    elements_.clear();

    std::ifstream nodeIn(nodeFile.c_str());

    if(!nodeIn){
        std::cerr << "Error: cannot open " << nodeFile << std::endl;
        std::exit(1);
    }

    double x, y;

    ymin_ =  std::numeric_limits<double>::max();
    ymax_ = -std::numeric_limits<double>::max();

    while(nodeIn >> x >> y){
        nodes_.push_back(Eigen::Vector2d(x, y));

        if(y < ymin_) ymin_ = y;
        if(y > ymax_) ymax_ = y;
    }

    std::ifstream elemIn(elementFile.c_str());

    if(!elemIn){
        std::cerr << "Error: cannot open " << elementFile << std::endl;
        std::exit(1);
    }

    int n1, n2, n3, n4;

    while(elemIn >> n1 >> n2 >> n3 >> n4){
        elements_.push_back(std::array<int, 4>{{n1, n2, n3, n4}});
    }

    const int numberOfNodes = static_cast<int>(nodes_.size());

    K_   = Eigen::MatrixXd::Zero(numberOfNodes, numberOfNodes);
    f_   = Eigen::VectorXd::Zero(numberOfNodes);
    phi_ = Eigen::VectorXd::Zero(numberOfNodes);

    std::cout << "number of nodes    = " << nodes_.size() << std::endl;
    std::cout << "number of elements = " << elements_.size() << std::endl;
    std::cout << "ymin = " << ymin_ << std::endl;
    std::cout << "ymax = " << ymax_ << std::endl;
}

void FEM::shapeFunctionDerivative(double s,
                                  double t,
                                  Eigen::Vector4d& dN_ds,
                                  Eigen::Vector4d& dN_dt) const
{
    // 形状関数
    // N1 = 1/4 (1-s)(1-t)
    // N2 = 1/4 (1+s)(1-t)
    // N3 = 1/4 (1+s)(1+t)
    // N4 = 1/4 (1-s)(1+t)

    dN_ds(0) = -0.25 * (1.0 - t);
    dN_ds(1) =  0.25 * (1.0 - t);
    dN_ds(2) =  0.25 * (1.0 + t);
    dN_ds(3) = -0.25 * (1.0 + t);

    dN_dt(0) = -0.25 * (1.0 - s);
    dN_dt(1) = -0.25 * (1.0 + s);
    dN_dt(2) =  0.25 * (1.0 + s);
    dN_dt(3) =  0.25 * (1.0 - s);
}

void FEM::assemble()
{
    const double gaussPoint[2] = {
        -1.0 / std::sqrt(3.0),
         1.0 / std::sqrt(3.0)
    };

    const double weight[2] = {1.0, 1.0};

    for(size_t e = 0; e < elements_.size(); e++){

        const std::array<int, 4>& element = elements_[e];

        Eigen::Matrix4d Ke = Eigen::Matrix4d::Zero();

        for(int i = 0; i < 2; i++){
            for(int j = 0; j < 2; j++){

                const double s = gaussPoint[i];
                const double t = gaussPoint[j];

                Eigen::Vector4d dN_ds;
                Eigen::Vector4d dN_dt;

                shapeFunctionDerivative(s, t, dN_ds, dN_dt);

                // ヤコビ行列 J
                Eigen::Matrix2d J = Eigen::Matrix2d::Zero();

                for(int a = 0; a < 4; a++){
                    const int nodeID = element[a];

                    const double x = nodes_[nodeID](0);
                    const double y = nodes_[nodeID](1);

                    J(0,0) += dN_ds(a) * x; // dx/ds
                    J(0,1) += dN_ds(a) * y; // dy/ds
                    J(1,0) += dN_dt(a) * x; // dx/dt
                    J(1,1) += dN_dt(a) * y; // dy/dt
                }

                const double detJ = J.determinant();

                if(detJ <= 0.0){
                    std::cerr << "Error: detJ <= 0 at element " << e << std::endl;
                    std::cerr << "Check element node order." << std::endl;
                    std::exit(1);
                }

                const Eigen::Matrix2d invJ = J.inverse();

                // 1行目：dN/dx
                // 2行目：dN/dy
                Eigen::Matrix<double, 2, 4> B;

                for(int a = 0; a < 4; a++){
                    const Eigen::Vector2d dN_global
                        = invJ * Eigen::Vector2d(dN_ds(a), dN_dt(a));

                    B(0,a) = dN_global(0);
                    B(1,a) = dN_global(1);
                }

                // 要素剛性行列
                // Ke = ∫ B^T B detJ ds dt
                Ke += B.transpose() * B * detJ * weight[i] * weight[j];
            }
        }

        // 全体行列へ
        for(int a = 0; a < 4; a++){
            for(int b = 0; b < 4; b++){
                const int A = element[a];
                const int Bnode = element[b];

                K_(A, Bnode) += Ke(a, b);
            }
        }
    }
}

void FEM::applyBoundaryCondition()
{
    // 下底で phi = 0
    // 上底で phi = 1
    // 左右ノイマン0

    const double eps = 1.0e-8;

    for(int i = 0; i < static_cast<int>(nodes_.size()); i++){

        const double y = nodes_[i](1);

        bool isDirichlet = false;
        double value = 0.0;

        if(std::abs(y - ymin_) < eps){
            isDirichlet = true;
            value = 0.0;
        }

        if(std::abs(y - ymax_) < eps){
            isDirichlet = true;
            value = 1.0;
        }

        if(isDirichlet){
            f_ -= K_.col(i) * value;

            K_.row(i).setZero();
            K_.col(i).setZero();
            K_(i,i) = 1.0;

            f_(i) = value;
        }
    }
}

void FEM::solve()
{
    phi_ = K_.colPivHouseholderQr().solve(f_);
}

void FEM::writeResult(const std::string& filename) const
{
    std::ofstream out(filename.c_str());

    if(!out){
        std::cerr << "Error: cannot open " << filename << std::endl;
        std::exit(1);
    }

    for(int i = 0; i < static_cast<int>(nodes_.size()); i++){
        out << nodes_[i](0) << " "
            << nodes_[i](1) << " "
            << phi_(i) << std::endl;
    }
}
