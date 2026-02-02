#include "ukf.hpp"

UKF::UKF() : alpha_(1e-3), beta_(2.0), kappa_(0.0) // 控制sigma点分布，高斯分布最优为2
{
    // 计算 lambda
    lambda_ = alpha_ * alpha_ * (n_x + kappa_) - n_x;
    // 初始化权重
    int n_sigma = 2 * n_x + 1;
    weights_m_ = Eigen::VectorXd::Zero(n_sigma);
    weights_c_ = Eigen::VectorXd::Zero(n_sigma);
    // 第0个点的权重
    weights_m_(0) = lambda_ / (lambda_ + n_x);
    weights_c_(0) = weights_m_(0) + (1 - alpha_ * alpha_ + beta_);
    // 其余点的权重
    double weight_rest = 1.0 / (2 * (lambda_ + n_x));
    for (int i = 1; i < n_sigma; ++i)
    {
        weights_m_(i) = weight_rest;
        weights_c_(i) = weight_rest;
    }
    // 初始化状态
    x_ = Eigen::VectorXd::Zero(n_x);
    // 初始化协方差矩阵p：表示初始不确定性
    // 位置不确定 ，速度不确定
    P_ = Eigen::MatrixXd::Zero(n_x, n_x);
    P_.block<3, 3>(0, 0) = Eigen::Matrix3d::Identity() * 1.0; // 位置初始不确定性
    P_.block<3, 3>(3, 3) = Eigen::Matrix3d::Identity() * 1.0; // 速度初始不确定性
    // 过程噪声协方差矩阵q
    Q_ = Eigen::MatrixXd::Zero(n_x, n_x);
    double q_pos = 0.01; // 位置过程噪声(调高会使得预测更发散)
    double q_vel = 0.1;  // 速度过程噪声
    Q_.block<3, 3>(0, 0) = Eigen::Matrix3d::Identity() * q_pos;
    Q_.block<3, 3>(3, 3) = Eigen::Matrix3d::Identity() * q_vel;
    // 测量噪声协方差矩阵r
    // 例如：±5mm，方差=(0.005m)^2=2.5e-5
    R_ = Eigen::MatrixXd::Zero(n_z, n_z);
    R_(0, 0) = 2.5e-5; // x方向测量噪声
    R_(1, 1) = 2.5e-5; // y方向测量噪声
    R_(2, 2) = 2.5e-5; // z方向测量噪声
}
Eigen::MatrixXd UKF::generateSigmaPoints() const
{
    Eigen::MatrixXd Xsig(n_x, 2 * n_x + 1);
    Xsig.col(0) = x_; // 第0个sigma点是均值
    // 计算（n_x+lambda）*p的cholesky分解
    Eigen::MatrixXd P_scaled = (n_x + lambda_) * P_;
    Eigen::LLT<Eigen::MatrixXd> llt(P_scaled);
    if (llt.info() != Eigen::Success)
    {
        std::cerr << "LLT decomposition failed!" << std::endl;
        return Xsig;
    }
    Eigen::MatrixXd L = llt.matrixL(); // 下三角矩阵
    // 生成其余2*n_x个点
    for (int i=0; i < n_x; ++i)
    {
        Xsig.col(i + 1) = x_ + L.col(i);         // +方向
        Xsig.col(i + 1 + n_x) = x_ - L.col(i); //- 方向
    }
    return Xsig;
}
Eigen::VectorXd UKF::motionModel(const Eigen::VectorXd &x,double dt)const
{
    Eigen::VectorXd x_new = x;
    //位置更新：p=p+v*dt
    x_new.segment<3>(0)+=x.segment<3>(3)*dt;
    return x_new;
}
Eigen::Vector3d UKF::observationModel(const Eigen::VectorXd &x)const
{
    return x.head<3>();
}
void UKF::predict(double dt)
{
    //生成当前sigma点
    Eigen::MatrixXd Xsig = generateSigmaPoints();
    int n_sigma = 2*n_x + 1;
    //对每个sigma点应用运动模型
    Eigen::MatrixXd Xsig_pred(n_x, n_sigma);
    for(int i=0;i<n_sigma;++i)
    {
        Xsig_pred.col(i)=motionModel(Xsig.col(i),dt);
    }
    //计算预测均值
    x_.setZero();
    for(int i=0;i<n_sigma;++i)
    {
        x_+=weights_m_(i)*Xsig_pred.col(i);
    }

    //计算预测协方差
    P_.setZero();
    for(int i=0;i<n_sigma;++i)
    {
        Eigen::VectorXd diff=Xsig_pred.col(i)-x_;
        P_+=weights_c_(i)*(diff*diff.transpose());
    }
    P_ += Q_; // 加上过程噪声
}
void UKF::update(const Eigen::Vector3d &z)
{
    //用当前状态生成sigma点（注意：此时x_和p_是预测后的）
    Eigen::MatrixXd Xsig=generateSigmaPoints();
    int n_sigma=2*n_x+1;
    //通过观测模型得到预测测量值的sigma点
    Eigen::MatrixXd Zsig(n_z,n_sigma);
    for(int i=0;i<n_sigma;++i)
    {
        Zsig.col(i)=observationModel(Xsig.col(i));
    } 
    //计算预测观测均值
    Eigen::Vector3d z_pred=Eigen::Vector3d::Zero();
    for(int i=0;i<n_sigma;++i)
    {
        z_pred+=weights_m_(i)*Zsig.col(i);
    }
    //计算观测协方差
    Eigen::Matrix3d S=Eigen::Matrix3d::Zero();
    for(int i=0;i<n_sigma;++i)
    {
        Eigen::Vector3d diff=Zsig.col(i)-z_pred;
        S+=weights_c_(i)*(diff*diff.transpose());
    }
    S+=R_; //加上观测噪声
    //计算状态与观测的协方差
    Eigen::MatrixXd P_xz=Eigen::MatrixXd::Zero(n_x,n_z);
    for(int i=0;i<n_sigma;++i)
    {
        Eigen::VectorXd dx=Xsig.col(i)-x_;
        Eigen::Vector3d dz=Zsig.col(i)-z_pred;
        P_xz+=weights_c_(i)*(dx*dz.transpose());
    }
    //计算卡尔曼增益
    Eigen::MatrixXd K=P_xz*S.inverse();
    //更新状态和协方差
    Eigen::Vector3d z_diff=z - z_pred;
    x_+=K*z_diff;
    P_-=K*S*K.transpose();
    //强制对称
    P_=(P_+P_.transpose())*0.5;
}

