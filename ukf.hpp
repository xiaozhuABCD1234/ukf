#pragma once

#include <Eigen/Dense>
#include <iostream>

class UKF
{
private:
    // 状态维度（位置3+速度3）
    static constexpr int n_x = 6;
    // 观测维度
    static constexpr int n_z = 3;
    // UKF参数
    double alpha_, beta_, kappa_,lambda_;
    // 状态向量 [px,py,pz,vx,vy,vz]
    Eigen::VectorXd x_;
    // 状态协方差矩阵
    Eigen::MatrixXd P_;
    // 过程噪声协方差矩阵
    Eigen::MatrixXd Q_;
    // 测量噪声协方差矩阵
    Eigen::MatrixXd R_;
    // sigma点权重
    Eigen::VectorXd weights_m_;
    Eigen::VectorXd weights_c_;
    // 辅助函数
    Eigen::MatrixXd generateSigmaPoints() const;
    Eigen::VectorXd motionModel(const Eigen::VectorXd &x, double dt) const;
    Eigen::Vector3d observationModel(const Eigen::VectorXd &x) const;

public:
    // 初始化状态
    UKF();
    // 根据时间步长dt推演状态
    void predict(double dt);
    // 根据测量值z更新状态
    void update(const Eigen::Vector3d &z);
    // 获取当前估计的位置和速度
    Eigen::Vector3d getPosition() { return x_.head<3>(); }
    Eigen::Vector3d getVelocity() { return x_.tail<3>(); }
};