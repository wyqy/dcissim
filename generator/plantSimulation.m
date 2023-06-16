function ret_struct = plantSimulation(plant_info, excitation, rs)
%PLANTSIMULATION 针对噪声的Monte-Carlo实验

    % 获取系统参数
    x_size = size(plant_info.A, 1); y_size = size(plant_info.C, 1); u_size = size(plant_info.B, 2);
    xyu_size = x_size + y_size + u_size;
    signal_sample = size(excitation, 2);
    covariance = plant_info.cov;
    % 构造噪声变换矩阵
    nozero_locs = diag(covariance) ~= 0;
    scale_transfer = zeros(xyu_size, xyu_size);
    scale_transfer(nozero_locs, nozero_locs) = ctranspose(chol(covariance(nozero_locs, nozero_locs)));

    % 生成噪声信号
    noise = zeros(xyu_size, signal_sample);
    for iter = 1:xyu_size, noise(iter, :) = randn(rs{iter}, 1, signal_sample); end
    noise = scale_transfer*noise;

    % 生成仿真数据
    x_init = zeros(x_size, 1);  % 零初始化
    [xk, yk, uk] = plantModel(plant_info, x_init, excitation, noise);
    
    % 返回值
    ret_struct = struct('xk', xk, 'yk', yk, 'uk', uk, 'noise', noise);

end

