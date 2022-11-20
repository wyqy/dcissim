function ret_struct = idenDCISSIMRunner(iden_struct, yn, un)
%IDENDCISSIMRUNNER discrete-cISSIM系统辨识 - 运行时程序
%   此处显示详细说明

    % 保存步进
    persistent is_inited n 
    persistent x_size_sim y_isim u_isim mat_x_sim mat_a_sim mat_b_sim mat_c_sim mat_d_sim
    persistent cov_zr_innovation cov_zrt_innovation cov_zrt_real innovation_mat_k innovation_mat_covariance
    
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
    % 初始化器增加的参数
    frequencies = iden_struct.frequencies;
    mat_s = iden_struct.mat_s;
    % 在线辨识时需要的参数
    online_regressor = iden_struct.online_regressor;
    online_regressor_frequencies = online_regressor.excitation_frequencies;
    online_regressor_phi = online_regressor.excitation_phi;
    online_sim_x_size_type = iden_struct.online_sim_x_size_type;
    online_cov_order_type = iden_struct.online_cov_order_type;
    online_sim_x_size = iden_struct.online_sim_x_size;
    online_cov_order = iden_struct.online_cov_order;

    % 初始化参数
    if isempty(is_inited)
        v_size = size(mat_s, 1);
        n = 0;
        x_size_sim = 1; y_isim = zeros(y_size, v_size); u_isim = zeros(u_size, v_size);
        mat_x_sim = zeros(1, v_size); mat_a_sim = 0; mat_b_sim = zeros(1, u_size); mat_c_sim = zeros(y_size, 1); mat_d_sim = zeros(y_size, u_size);
        cov_zr_innovation = zeros(1+y_size, 1+y_size); cov_zrt_innovation = zeros(1+y_size+u_size, u_size); cov_zrt_real = zeros(1+y_size+u_size, u_size); innovation_mat_k = zeros(1, y_size); innovation_mat_covariance = zeros(y_size, y_size);
        is_inited = 1;
    end

    % 离线or在线
    switch dcissim_type
        case 'offline'  % 离线辨识
            % 参数计算与预处理
            switch isim_excitation_type
                case 'full'
                    yn = yn(:, cutted_periods*period_samples+1:end);
                    un = un(:, cutted_periods*period_samples+1:end);
                    signal_sample = size(yn, 2);
                case 'reduced'
                    yn = yn(:, cutted_periods*period_samples+1:end);
                    un = un(:, cutted_periods*period_samples+1:end);
                    signal_sample = period_samples*fix(size(yn, 2)/period_samples);
                    yn = yn(:, 1:signal_sample);
                    un = un(:, 1:signal_sample);
                otherwise, yn = 0; un = 0; signal_sample = 1;
            end
            % ISIM辨识
            regressor_struct = idenRegressor(period_samples, frequencies, signal_sample, 'ordinary');  % 直接计算
            vn = regressor_struct.vn;
            [y_isim, u_isim] = idenISIM(yn, un, vn, 'ordinary');
            % SIM辨识
            [mat_a_sim, mat_c_sim, x_size_sim] = idenACN(y_isim, u_isim, mat_s, x_size_upbound, online_sim_x_size, online_sim_x_size_type);  % 询问阶数
            [mat_b_sim, mat_d_sim, mat_x_sim] = idenBDX(y_isim, u_isim, mat_s, mat_a_sim, mat_c_sim, x_size_sim, sim_ss_bdx_type, sim_ss_d_type);
            % 方差辨识
            [cov_zr_innovation, cov_zrt_innovation, cov_zrt_real, innovation_mat_k, innovation_mat_covariance] = idenCovariance(yn, un, vn, y_isim, u_isim, mat_a_sim, mat_c_sim, online_cov_order, online_cov_order_type);  % (可能)询问阶数
        case 'online'  % 在线辨识
            % 记录迭代步
            n = n + 1;
            % ISIM辨识
            vn = sin((online_regressor_frequencies*n) + online_regressor_phi);
            [y_isim, u_isim] = idenISIM(yn, un, vn, 'recursive');
            % 节省计算时间
            if mod(n, period_samples) == 0
                % SIM辨识
                [mat_a_sim, mat_c_sim, x_size_sim] = idenACN(y_isim, u_isim, mat_s, x_size_upbound, online_sim_x_size, online_sim_x_size_type);
                [mat_b_sim, mat_d_sim, mat_x_sim] = idenBDX(y_isim, u_isim, mat_s, mat_a_sim, mat_c_sim, x_size_sim, sim_ss_bdx_type, sim_ss_d_type);
                % 方差辨识
                [cov_zr_innovation, cov_zrt_innovation, cov_zrt_real, innovation_mat_k, innovation_mat_covariance] = idenCovariance(yn, un, vn, y_isim, u_isim, mat_a_sim, mat_c_sim, online_cov_order, online_cov_order_type);  % (可能)询问阶数
            end
        case 'offline-test'  % 离线-间隔n次进行 (避免速度过慢, 下同)
            % 参数计算与预处理
            yn = yn(:, cutted_periods*period_samples+1:end);
            un = un(:, cutted_periods*period_samples+1:end);
            signal_sample = size(yn, 2);
            % ISIM辨识
            regressor_struct = idenRegressor(period_samples, frequencies, signal_sample, 'ordinary');  % 直接计算
            vn = regressor_struct.vn;
            [y_isim, u_isim] = idenISIM(yn, un, vn, 'ordinary');
            % SIM辨识
            [mat_a_sim, mat_c_sim, x_size_sim] = idenACN(y_isim, u_isim, mat_s, x_size_upbound, online_sim_x_size, online_sim_x_size_type);  % 询问阶数
            [mat_b_sim, mat_d_sim, mat_x_sim] = idenBDX(y_isim, u_isim, mat_s, mat_a_sim, mat_c_sim, x_size_sim, sim_ss_bdx_type, sim_ss_d_type);
            % 方差辨识
            [cov_zr_innovation, cov_zrt_innovation, cov_zrt_real, innovation_mat_k, innovation_mat_covariance] = idenCovariance(yn, un, vn, y_isim, u_isim, mat_a_sim, mat_c_sim, online_cov_order, online_cov_order_type);  % (可能)询问阶数
    end
    
    % 防止可能的虚参数
    mat_b_sim = real(mat_b_sim);
    mat_d_sim = real(mat_d_sim);
    mat_x_sim = real(mat_x_sim);
    % 整合方差矩阵
    cov_ztrt_innovation = cov_zrt_innovation(1:x_size_sim+y_size, :); cov_ztrt_real = cov_zrt_real(1:x_size_sim+y_size, :);
    cov_tt_innovation = cov_zrt_innovation(x_size_sim+y_size+1:end, :); cov_tt_real = cov_zrt_real(x_size_sim+y_size+1:end, :);
    mat_cov_innovation = [cov_zr_innovation cov_ztrt_innovation; cov_ztrt_innovation.' cov_tt_innovation];
    mat_cov_halfreal = [cov_zr_innovation cov_ztrt_real; cov_ztrt_real.' cov_tt_real];

    % 返回值
    ret_struct = struct('x_size', x_size_sim, 'X', mat_x_sim, 'Y', y_isim, 'U', u_isim, ...
        'A', mat_a_sim, 'B', mat_b_sim, 'C', mat_c_sim, 'D', mat_d_sim, 'K', innovation_mat_k, ...
        'cov_innovation', innovation_mat_covariance, 'cov_io_innovation', mat_cov_innovation, 'cov_io_halfreal', mat_cov_halfreal);

end

