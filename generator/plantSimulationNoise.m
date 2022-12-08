function ret_struct = plantSimulationNoise(plant_info, excitation, rs)
%PLANTSIMULATIONNOISE 相同系统, 不同噪声

    % 获取系统参数
    x_size = size(plant_info.A, 1); y_size = size(plant_info.C, 1); u_size = size(plant_info.B, 2);
    xyu_size = x_size + y_size + u_size;
    signal_sample = size(excitation, 2);
    covariance = plant_info.cov;

    % 生成噪声信号
    noise = zeros(xyu_size, signal_sample);
    for iter = 1:xyu_size, noise(iter, :) = randn(rs{iter}, 1, signal_sample); end
    nozero_selection = diag(covariance) ~= 0;
    covariance_rational = zeros(xyu_size, xyu_size);
    covariance_rational(nozero_selection, nozero_selection) = covariance(nozero_selection, nozero_selection);
    scale_transfer = zeros(xyu_size, xyu_size);
    scale_transfer(nozero_selection, nozero_selection) = sqrtm(covariance_rational(nozero_selection, nozero_selection));
    noise = scale_transfer*noise;

    % 生成仿真数据
    x_init = zeros(x_size, 1);  % 零初始化
    [xn, yn, un] = plantModel(plant_info, x_init, excitation, noise);
    
    % 返回值
    ret_struct = struct('xn', xn, 'yn', yn, 'un', un, 'noise', noise);

end

