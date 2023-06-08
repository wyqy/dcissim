function [xk_corrupted, yk_corrupted, uk_corrupted] = plantModel(plant_info, x_init, uk, noise)
%PLANTMODEL 模型仿真主函数
% x(k+1) = Ax(k) + Bu(k);
% y(k) = Cx(k) + Du(k);

    % 参数提取
    plantss_a = plant_info.A;
    plantss_b = plant_info.B;
    plantss_c = plant_info.C;
    plantss_d = plant_info.D;

    % 参数计算
    x_size = size(plantss_a, 1);
    y_size = size(plantss_c, 1);
    u_size = size(plantss_d, 2);
    signal_sample = size(uk, 2);

    % 返回值初始化
    xk_corrupted = zeros(x_size, signal_sample);
    yk_corrupted = zeros(y_size, signal_sample);
    uk_corrupted = zeros(u_size, signal_sample);
    
    % 切分噪声
    x_noise = noise(1:x_size, :);
    y_noise = noise(x_size+1:x_size+y_size, :);
    u_noise = noise(x_size+y_size+1:end, :);

    % 系统仿真
    xk_corrupted(:, 1) = x_init;
    for iter = 1:signal_sample
        if iter < signal_sample, xk_corrupted(:, iter+1) = plantss_a*xk_corrupted(:, iter) + plantss_b*uk(:, iter) + x_noise(:, iter); end
        yk_corrupted(:, iter) = plantss_c*xk_corrupted(:, iter) + plantss_d*uk(:, iter) + y_noise(:, iter);
        uk_corrupted(:, iter) = uk(:, iter) + u_noise(:, iter);
    end

end

