#include "ukf.hpp"
#include <iostream>
#include <random>
#include <cmath>
#include <fstream>

int main() {
    UKF ukf;
    double dt = 0.02; // 50 Hz

    // 设置真实轨迹（匀速直线）
    Eigen::VectorXd true_state(6);
    true_state << 0.0, 0.0, 1.0, 0.5, 0.3, 0.0; // start at (0,0,1), move in xy-plane

    // 随机数生成器（模拟相机噪声）
    std::random_device rd;
    std::mt19937 gen(rd());
    std::normal_distribution<double> noise_xy(0.0, 0.005); // ±5mm
    std::normal_distribution<double> noise_z(0.0, 0.010);  // ±10mm

    std::cout << "Step | True Pos (m)       | Est Pos (m)        | Error (mm)\n";
    std::cout << "---------------------------------------------------------------\n";

    std::ofstream csv("/home/ubuntu/code/ukf/output/ukf_results.csv");
    csv << "step,true_x,true_y,true_z,est_x,est_y,est_z,error_mm\n";

    for (int step = 0; step < 100; ++step) {
        // 真实状态演化
        true_state.segment<3>(0) += true_state.segment<3>(3) * dt;

        // 模拟相机测量（加噪声）
        Eigen::Vector3d measurement = true_state.head<3>();
        measurement.x() += noise_xy(gen);
        measurement.y() += noise_xy(gen);
        measurement.z() += noise_z(gen);

        // UKF 步骤
        ukf.predict(dt);
        ukf.update(measurement);

        // 输出结果
        Eigen::Vector3d est = ukf.getPosition();
        Eigen::Vector3d true_pos = true_state.head<3>();
        double error_mm = (est - true_pos).norm() * 1000.0;

        printf("%4d | [%5.2f,%5.2f,%5.2f] | [%5.2f,%5.2f,%5.2f] | %7.2f\n",
               step,
               true_pos.x(), true_pos.y(), true_pos.z(),
               est.x(), est.y(), est.z(),
               error_mm);

        csv << step << ","
            << true_pos.x() << "," << true_pos.y() << "," << true_pos.z() << ","
            << est.x() << "," << est.y() << "," << est.z() << ","
            << error_mm << "\n";
    }

    return 0;
}