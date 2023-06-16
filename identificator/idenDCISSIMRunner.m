function ret_struct = idenDCISSIMRunner(iden_struct, yk, uk)
%IDENDCISSIMRUNNER discrete-cISSIM系统辨识 - 运行时程序
%   此处显示详细说明

    % 保存步进
    persistent is_inited k
    persistent x_size_est yv_est uv_est xv_est
    persistent mat_a_est mat_b_est mat_c_est mat_d_est
    persistent noise_para_est cov_est kalman_est
    
    % 提取参数
    mat_s = iden_struct.mat_s;
    mat_v = iden_struct.mat_v;  % 单周期内的v
    regressor = iden_struct.regressor;
    v0 = regressor.v0;
    freq_list = regressor.freq_list;
    phi_list = regressor.phi_list;

    prior = iden_struct.prior;
    x_size_prior = prior.xsize;
    aorder_prior = prior.aorder;

    y_size = iden_struct.y_size;
    u_size = iden_struct.u_size;
    x_size_upbound = iden_struct.x_size_upbound;
    samples_period = iden_struct.samples_period;

    algo_type = iden_struct.algo_type;
    xsize_est_type = iden_struct.xsize_est_type;
    aorder_est_type = iden_struct.aorder_est_type;
    als_est_type = iden_struct.als_est_type;

    plant_d_type = iden_struct.plant_d_type;
    cov_cross_type = iden_struct.cov_cross_type;

    % 初始化参数
    if isempty(is_inited)
        v_size = size(mat_s, 1);
        k = 0; x_size_est = x_size_prior;
        yv_est = zeros(y_size, v_size); uv_est = zeros(u_size, v_size); xv_est = zeros(x_size_est, v_size); x0_est = zeros(x_size_est, 1);
        mat_a_est = zeros(x_size_est); mat_b_est = zeros(x_size_est, u_size); mat_c_est = zeros(y_size, 1); mat_d_est = zeros(y_size, u_size);
        cov_est = zeros(x_size_est+y_size+u_size); kalman_est = zeros(x_size_est, y_size);
        noise_para_est = struct('cov_all', cov_est, 'kalman', kalman_est);
        is_inited = 1;
    end

    % 离线or在线
    switch algo_type
        case 'offline'  % 离线辨识
            % 去掉第一个周期和最末的不完整周期
            yk = yk(:, samples_period+1:end);
            uk = uk(:, samples_period+1:end);
            signal_sample = samples_period*fix(size(yk, 2)/samples_period);
            yk = yk(:, 1:signal_sample);
            uk = uk(:, 1:signal_sample);
            % ISIM辨识
            [yv_est, uv_est] = idenISIM(yk, uk, mat_v, samples_period, 'ordinary');
            % ynoise = yn - yv_est*vn; figure; plot(ynoise(1, :));
            % unoise = un - uv_est*vn; figure; plot(unoise(1, :));
            % SIM辨识
            [mat_a_est, mat_c_est, x_size_est] = idenACN(yv_est, uv_est, mat_s, x_size_upbound, x_size_prior, xsize_est_type);  % 询问阶数
            [mat_b_est, mat_d_est, xv_est] = idenBDX(yv_est, uv_est, mat_s, mat_a_est, mat_c_est, x_size_est, plant_d_type);
            % 方差辨识
            if ~strcmp(als_est_type, 'none')
                para_est = struct('v0', v0, 'X', xv_est, 'Y', yv_est, 'U', uv_est, 'A', mat_a_est, 'B', mat_b_est, 'C', mat_c_est, 'D', mat_d_est);
                noise_para_est = idenCovariance(yk, uk, mat_v, para_est, als_est_type, cov_cross_type, aorder_prior, aorder_est_type);  % (可能)询问阶数
            end
        case 'online'  % 在线辨识
            % 记录迭代步
            k = k + 1;
            % ISIM辨识
            vk = sin((freq_list*k) + phi_list);
            [yv_est, uv_est] = idenISIM(yk, uk, vk, samples_period, 'recursive');
            % 整数周期
            if mod(k, samples_period) == 0
                % SIM辨识
                [mat_a_est, mat_c_est, x_size_est] = idenACN(yv_est, uv_est, mat_s, x_size_upbound, x_size_prior, xsize_est_type);
                [mat_b_est, mat_d_est, xv_est] = idenBDX(yv_est, uv_est, mat_s, mat_a_est, mat_c_est, x_size_est, plant_d_type);
                % 无方差辨识
            end
    end
    
    % 整合方差矩阵
    cov_est = noise_para_est.cov_all;
    kalman_est = noise_para_est.kalman;
    
    % 返回值
    % 构造状态空间
    ret_sys = idss(mat_a_est, mat_b_est, mat_c_est, mat_d_est, kalman_est);
    ret_struct = struct('ss', ret_sys, 'v0', v0, 'freqs', freq_list, 'X', xv_est, 'Y', yv_est, 'U', uv_est, ...
        'x_size', x_size_est, 'A', mat_a_est, 'B', mat_b_est, 'C', mat_c_est, 'D', mat_d_est, ...
        'K', kalman_est, 'cov', cov_est, 'noise_para', noise_para_est);

end

