function ret_struct = idenDCISSIMRunner(iden_struct, yn, un)
%IDENDCISSIMRUNNER discrete-cISSIM系统辨识 - 运行时程序
%   此处显示详细说明

    % 保存步进
    persistent is_inited n 
    persistent x_size_sim y_isim u_isim mat_x_sim mat_a_sim mat_b_sim mat_c_sim mat_d_sim
    persistent covariance_struct cov kalman_gain
    
    % 提取参数
    y_size = iden_struct.y_size;
    u_size = iden_struct.u_size;
    period_samples = iden_struct.period_samples;
    cutted_periods = iden_struct.cutted_periods;
    dcissim_type = iden_struct.dcissim_type;
    isim_excitation_type = iden_struct.isim_excitation_type;
    x_size_upbound = iden_struct.x_size_upbound;
    sim_ss_bdx_type = iden_struct.sim_ss_bdx_type;
    sim_ss_d_type = iden_struct.sim_ss_d_type;
    cov_cross_type = iden_struct.cov_cross_type;
    cov_est_type = iden_struct.cov_est_type;
    % 初始化器增加的参数
    frequencies = iden_struct.frequencies;
    mat_s = iden_struct.mat_s;
    % 在线辨识时需要的参数
    regressor = iden_struct.regressor;
    regressor_frequencies = regressor.excitation_frequencies;
    regressor_phi = regressor.excitation_phi;
    % 自动化辨识需要的参数
    sim_x_size_type = iden_struct.sim_x_size_type;
    cov_order_type = iden_struct.cov_order_type;
    sim_x_size = iden_struct.sim_x_size;
    cov_order = iden_struct.cov_order;
    % 其它参数
    v0 = regressor.v0;
    % 方差估计参数

    % 初始化参数
    if isempty(is_inited)
        v_size = size(mat_s, 1);
        n = 0;
        x_size_sim = sim_x_size; y_isim = zeros(y_size, v_size); u_isim = zeros(u_size, v_size);
        mat_x_sim = zeros(x_size_sim, v_size); mat_a_sim = zeros(x_size_sim); mat_b_sim = zeros(x_size_sim, u_size); mat_c_sim = zeros(y_size, 1); mat_d_sim = zeros(y_size, u_size);
        cov = zeros(x_size_sim+y_size+u_size); kalman_gain = zeros(x_size_sim, y_size); covariance_struct = struct('cov', cov, 'kalman', kalman_gain);
        is_inited = 1;
    end

    % 离线or在线
    switch dcissim_type
        case 'offline'  % 离线辨识
            % 参数计算与预处理
            vn = idenRegressor(period_samples, frequencies, size(yn, 2), 'ordinary').vn;  % 直接计算
            cutted_yn = yn(:, cutted_periods*period_samples+1:end);
            cutted_un = un(:, cutted_periods*period_samples+1:end);
            cutted_vn = vn(:, cutted_periods*period_samples+1:end);
            switch isim_excitation_type
                case 'full'
                case 'reduced'
                    signal_sample = period_samples*fix(size(cutted_yn, 2)/period_samples);
                    cutted_yn = cutted_yn(:, 1:signal_sample);
                    cutted_un = cutted_un(:, 1:signal_sample);
                    cutted_vn = cutted_vn(:, 1:signal_sample);
                otherwise, cutted_yn = 0; cutted_un = 0;
            end
            % ISIM辨识
            [y_isim, u_isim] = idenISIM(cutted_yn, cutted_un, cutted_vn, 'ordinary');
            % yy = cutted_yn - y_isim*cutted_vn; figure; plot(yy(1, :));
            % uu = cutted_un - u_isim*cutted_vn; figure; plot(uu(1, :));
            % SIM辨识
            [mat_a_sim, mat_c_sim, x_size_sim] = idenACN(y_isim, u_isim, mat_s, x_size_upbound, sim_x_size, sim_x_size_type);  % 询问阶数
            [mat_b_sim, mat_d_sim, mat_x_sim] = idenBDX(y_isim, u_isim, mat_s, mat_a_sim, mat_c_sim, x_size_sim, sim_ss_bdx_type, sim_ss_d_type);
            % 方差辨识
            mat_b_sim = real(mat_b_sim); mat_d_sim = real(mat_d_sim);
            covariance_struct = idenCovariance(yn, un, vn, u_isim, mat_a_sim, mat_b_sim, mat_c_sim, mat_d_sim, cov_est_type, cov_cross_type, cov_order, cov_order_type, cutted_periods*period_samples);  % (可能)询问阶数
        case 'online'  % 在线辨识
            % 记录迭代步
            n = n + 1;
            % ISIM辨识
            vn = sin((regressor_frequencies*n) + regressor_phi);
            [y_isim, u_isim] = idenISIM(yn, un, vn, 'recursive');
            % SIM辨识
            [mat_a_sim, mat_c_sim, x_size_sim] = idenACN(y_isim, u_isim, mat_s, x_size_upbound, sim_x_size, sim_x_size_type);
            [mat_b_sim, mat_d_sim, mat_x_sim] = idenBDX(y_isim, u_isim, mat_s, mat_a_sim, mat_c_sim, x_size_sim, sim_ss_bdx_type, sim_ss_d_type);
            % 节省计算时间
            % if mod(n, period_samples) == 0, end
            % 无方差辨识
    end
    
    % 防止可能的虚参数
    mat_b_sim = real(mat_b_sim);
    mat_d_sim = real(mat_d_sim);
    mat_x_sim = real(mat_x_sim);
    % 整合方差矩阵
    cov = covariance_struct.cov;
    kalman_gain = covariance_struct.kalman;
    % 减小存储空间
    mat_s = sparse(mat_s);
    % 返回值
    ret_struct = struct('x_size', x_size_sim, 'S', mat_s, 'X', mat_x_sim, 'Y', y_isim, 'U', u_isim, 'v0', v0, 'freqs', frequencies, ...
        'A', mat_a_sim, 'B', mat_b_sim, 'C', mat_c_sim, 'D', mat_d_sim, ...
        'K', kalman_gain, 'cov', cov);

end

