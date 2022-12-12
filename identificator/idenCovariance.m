function covariance_struct = idenCovariance(yn, un, vn, u_isim, mat_a_sim, mat_b_sim, mat_c_sim, mat_d_sim, cov_est_type, cov_cross_type, cov_order, cov_order_type, cutted_samples)
%IDENCOVARIANCE 辨识方差参数
% 针对存在互协方差与否提供不同的估计
% 存在互协方差: 返回cov_innovation_all (zrt), cov_innovation_kalman
% 不存在互协方差: 返回cov_real_all (zrt)

    % 计算各阶矩
    % 参数计算
    x_size = size(mat_a_sim, 1); y_size = size(yn, 1); u_size = size(un, 1);
    signal_sample = size(yn, 2);
    % 扰动量计算, 忽略暂态过程
    res_model = struct('A', mat_a_sim, 'B', mat_b_sim, 'C', mat_c_sim, 'D', mat_d_sim);
    res_est_xinit = zeros(x_size, 1); res_ss_noise = zeros(x_size+y_size+u_size, signal_sample);
    res_est_un = u_isim*vn;
    [~, res_est_yn, ~] = plantModel(res_model, res_est_xinit, un, res_ss_noise);
    % z0 = -real(mat_x_sim)*vn(:, 1); rstate = zrTransientState(z0, mat_a_sim, mat_c_sim, signal_sample);
    rn = yn - res_est_yn; rn = rn(:, cutted_samples+1:end);
    tn = un - res_est_un; tn = tn(:, cutted_samples+1:end);

    % 阶数估算
    % 估算阶数上限, 按矩阵谱半径估计
    test_order_max = ceil(log(1e-5)/log(max(abs(eig(mat_a_sim)))));
    test_rr_moment_norm = zeros(test_order_max+1, 1);
    % 询问或者直接估计/固定值
    switch cov_order_type
        case 'ask'
            % 计算
            for iter_test = 0:test_order_max, test_rr_moment_norm(iter_test+1) = norm(estimationCaller(rn, rn, iter_test, 0), 'fro'); end
            % 绘图-询问
            disp('The Frobenius norms of R''s auto-covariance are:')
            fig_ask = figure;
            plot(0:test_order_max, test_rr_moment_norm);
            uiwait(fig_ask);
            cov_order_estimated = input('Enter the estimated sizes of states: ');
            % 无效输入, 估计合适的值
            if ~(isscalar(cov_order_estimated) && isnumeric(cov_order_estimated)), cov_order_estimated = orderCaller(test_rr_moment_norm); end
        case 'estimate'
            for iter_test = 0:test_order_max, test_rr_moment_norm(iter_test+1) = norm(estimationCaller(rn, rn, iter_test, 0), 'fro'); end
            cov_order_estimated = orderCaller(test_rr_moment_norm);
        case 'fixed'
            cov_order_estimated = cov_order;
        case 'null'
            cov_all = zeros(x_size+y_size+u_size);
            kalman_gain = zeros(x_size, y_size);
            covariance_struct = struct('cov', cov_all, 'kalman', kalman_gain);
            return;
        otherwise, cov_order_estimated = 1;
    end

    % 各阶矩计算
    % 准备
    est_rr_moment = cell(cov_order_estimated+1, 1);
    est_rt_moment = cell(cov_order_estimated+1, 1);
    % 计算
    for iter_moment = 0:cov_order_estimated
        est_rr_moment{iter_moment+1} = estimationCaller(rn, rn, iter_moment, 0);
        est_rt_moment{iter_moment+1} = estimationCaller(rn, tn, iter_moment, 0);
    end
    est_tt_moment = estimationCaller(tn, tn, 0, 0);

    switch cov_cross_type
        case 'valid'
            % 估计zr的协方差矩阵(innovation form)
            [cov_zr, kalman_gain] = zrInnovationMomentEstimation(est_rr_moment, mat_a_sim, mat_c_sim);
            % 估计zrt的协方差矩阵(innovation form)
            cov_innovation_zrt = zrtInnovationMomentEstimation(est_rt_moment, est_tt_moment, mat_a_sim, mat_c_sim, kalman_gain);
            % 参数处理
            cov_zrt_1 = cov_innovation_zrt(1:x_size+y_size, :); cov_zrt_2 = cov_innovation_zrt(x_size+y_size+1:end, :);
            cov_all = [cov_zr cov_zrt_1; cov_zrt_1.' cov_zrt_2];
        case 'null'
            % 估计zr的协方差矩阵(icm form)
            switch cov_est_type
                case 'simple', cov_zr = zrRealMomentEstimation(est_rr_moment, mat_a_sim, mat_c_sim);
                case 'classical', cov_zr = zrRealMomentEstimation2(rn, cov_order_estimated, mat_a_sim, mat_c_sim, cutted_samples);
            end
            % 估计zrt的协方差矩阵(original form)
            cov_zrt = zrtRealMomentEstimation(est_rt_moment, est_tt_moment, mat_a_sim, mat_c_sim);
            % 参数处理
            cov_zrt_2 = cov_zrt(x_size+y_size+1:end, :);
            cov_all = [cov_zr zeros(x_size+y_size, u_size); zeros(u_size, x_size+y_size) cov_zrt_2];
            % 其它填0
            kalman_gain = zeros(x_size, y_size);
        otherwise, cov_all = zeros(x_size+y_size+u_size); kalman_gain = zeros(x_size, y_size);
    end

    % 返回数据
    if min(eig(cov_zr)) < 0 || min(eig(cov_zrt_2)) < 0
        cov_all = zeros(x_size+y_size+u_size); kalman_gain = zeros(x_size, y_size);
    end
    covariance_struct = struct('cov', cov_all, 'kalman', kalman_gain);
    
end

function moment_estimated = estimationCaller(an, bn, order, offset)
    % 计算Expectation{a_{n+order}*b_n}, offset为开头截断的长度
    moment_estimation_samples = size(an, 2);
    moment_estimated = (1/(moment_estimation_samples-order-offset)).*(an(:, offset+order+1:end)*(bn(:, offset+1:end-order).'));
end
function order_estimated = orderCaller(moment_norm)
    % 估算order
    order_max = length(moment_norm)-1;
    for iter_order = 0:order_max
        if moment_norm(iter_order+1) < 1e-4*moment_norm(2)
            order_estimated = iter_order; break;
        end
    end
    if iter_order == order_max, order_estimated = order_max; end
end

function [cov_zr, cov_kalman] = zrInnovationMomentEstimation(est_rr_moment, mat_a_sim, mat_c_sim)
% 输出-状态协方差 - 存在互协方差
% 此处仅计算forward innovaiton form对应的矩

    % 参数计算
    w_size = size(mat_a_sim, 1);
    r_size = size(mat_c_sim, 1);
    order = length(est_rr_moment)-1;
    est_rr0 = est_rr_moment{1};

    % 估计wr_1
    % 参数计算
    est_wr1_mat_q = zeros(order*r_size, w_size);
    location_base = 0;
    for iter_wr1 = 1:order
        est_wr1_mat_q(location_base+1:location_base+r_size, :) = mat_c_sim*mpower(mat_a_sim, iter_wr1-1);
        location_base = location_base + r_size;
    end
    est_wr1_mat_a = kron(eye(r_size), est_wr1_mat_q);
    est_wr1_vec_b = zeros(order*r_size, r_size);
    location_base = 0;
    for iter_wr1 = 1:order
        est_wr1_vec_b(location_base+1:location_base+r_size, :) = est_rr_moment{iter_wr1+1};
        location_base = location_base + r_size;
    end
    est_wr1_vec_b = reshape(est_wr1_vec_b, [], 1);
    % 最小二乘估计
    est_wr1 = lsqminnorm(est_wr1_mat_a, est_wr1_vec_b);
    est_wr1 = reshape(est_wr1, [w_size r_size]);

    % 计算ww0和k
    [est_ww0, cov_kalman, ~, ~] = idare(mat_a_sim.', mat_c_sim.', zeros(w_size), -est_rr0, -est_wr1, []);
    cov_kalman = cov_kalman.';

    % 计算qsr
    est_cov_q = est_ww0 - mat_a_sim*est_ww0*(mat_a_sim.');
    est_cov_s = est_wr1 - mat_a_sim*est_ww0*(mat_c_sim.');
    est_cov_r = est_rr0 - mat_c_sim*est_ww0*(mat_c_sim.');
    % 返回值
    cov_zr = [est_cov_q est_cov_s; est_cov_s.' est_cov_r];

end

function cov_zr = zrRealMomentEstimation(est_rr_moment, mat_a_sim, mat_c_sim)
% 输出-状态协方差 - 不存在互协方差
% 此处计算真实矩

    % 参数计算
    x_size = size(mat_a_sim, 1); 
    y_size = size(mat_c_sim, 1);
    order = length(est_rr_moment);
    est_rr0 = est_rr_moment{1};

    % 估计P的一种估计值
    % 参数计算
    est_para_m = zeros((order-1)*y_size, x_size);
    location_base = 0;
    for iter_moment = 1:order-1
        est_para_m(location_base+1:location_base+y_size, :) = mat_c_sim*mpower(mat_a_sim, iter_moment);
        location_base = location_base + y_size;
    end
    est_para_m = kron(mat_c_sim, est_para_m);
    est_para_b = zeros((order-1)*y_size, y_size);
    location_base = 0;
    for iter_moment = 1:order-1
        est_para_b(location_base+1:location_base+y_size, :) = est_rr_moment{iter_moment+1};
        location_base = location_base + y_size;
    end
    est_para_b = reshape(est_para_b, [(order-1)*y_size*y_size 1]);
    % LS估计
    % est_para_p = lsqminnorm(est_para_m, est_para_b);
    % est_mat_p = reshape(est_para_p, [x_size x_size]);
    % SDP估计
    cvx_begin sdp quiet;
        variable est_mat_p(x_size, x_size) symmetric;
        minimize(norm(est_para_m * vec(est_mat_p) - est_para_b)); %, Inf));
        subject to;
            est_mat_p == semidefinite(x_size);
            est_rr0 - mat_c_sim * est_mat_p * mat_c_sim.' == semidefinite(y_size);
            est_mat_p - mat_a_sim * est_mat_p * mat_a_sim.' == semidefinite(x_size);
    cvx_end;

    % 估计cov_y : R
    est_cov_y = est_rr0 - mat_c_sim * est_mat_p * mat_c_sim.';
    % 估计cov_x : Q
    est_cov_x = est_mat_p - mat_a_sim * est_mat_p * mat_a_sim.';

    % 返回值
    cov_zr = blkdiag(est_cov_x, est_cov_y);

end

function cov_zr = zrRealMomentEstimation2(rn, order, mat_a_sim, mat_c_sim, cutted_samples) 
% 输出-状态协方差 - 不存在互协方差
% 进行状态估计

    % 参数计算
    x_size = size(mat_a_sim, 1); 
    y_size = size(mat_c_sim, 1);
    sample_size = size(rn, 2);
    % K, L, F矩阵计算
    riccati_q = 1e-5*eye(x_size); riccati_r = 1e-5*eye(y_size);
    mat_k_sim = idare(mat_a_sim.', mat_c_sim.', riccati_q, riccati_r, [], []);
    mat_k_sim = mat_k_sim.'; mat_l_sim = mat_a_sim*mat_k_sim;
    mat_f_sim = mat_a_sim - mat_l_sim*mat_c_sim;

    % 随机状态估计 (零状态初始化)
    en = zeros(y_size, sample_size);
    xn = zeros(x_size, 1);
    for iter_est = 1:sample_size
        en(:, iter_est) = rn(:, iter_est) - mat_c_sim*xn;
        xn = mat_a_sim*xn + mat_l_sim*en(:, iter_est);
    end
    en = en(:, cutted_samples+1:end);
    % innovation矩估计
    en_moment = cell(order, 1);
    for iter_order = 0:order
        en_moment{iter_order+1} = estimationCaller(en, en, iter_order, 0);
    end

    % 估计P的一种估计值
    % 参数计算
    % L, N计算
    est_para_l = zeros(order*y_size, y_size);
    location_base = 0;
    for iter_order = 1:order
        est_para_l(location_base+1:location_base+y_size, :) = mat_c_sim*mpower(mat_f_sim, iter_order-1);
        location_base = location_base + y_size;
    end
    est_para_n = zeros(order*y_size, y_size);
    location_base = 0;
    est_para_n(location_base+1:location_base+y_size, :) = eye(y_size);
    location_base = location_base + y_size;
    for iter_order = 2:order
        est_para_n(location_base+1:location_base+y_size, :) = mat_c_sim*mpower(mat_f_sim, iter_order-2)*mat_l_sim;
        location_base = location_base + y_size;
    end
    % E, G计算
    est_para_e = kron(mat_c_sim, est_para_l) / (eye(x_size^2) - kron(mat_f_sim, mat_f_sim));
    est_para_g = est_para_e*kron(mat_l_sim, mat_l_sim) + kron(eye(y_size), est_para_n);
    % A, b计算
    est_para_b = zeros(order*y_size, y_size);
    location_base = 0;
    for iter_order = 1:order
        est_para_b(location_base+1:location_base+y_size, :) = en_moment{iter_order};
        location_base = location_base + y_size;
    end
    est_para_b = reshape(est_para_b, [order*y_size*y_size 1]);
    % LS估计
    % est_para_q = lsqminnorm([est_para_e est_para_g], est_para_b);
    % est_mat_q = reshape(est_para_q(1:x_size*x_size), [x_size x_size]);
    % est_mat_r = reshape(est_para_q(x_size*x_size+1:end), [y_size y_size]);
    % SDP估计
    cvx_begin sdp quiet;
        variable est_mat_q(x_size, x_size) symmetric;
        variable est_mat_r(y_size, y_size) symmetric;
        minimize(norm(est_para_e*vec(est_mat_q) + est_para_g*vec(est_mat_r) - est_para_b)); %, Inf));
        subject to;
            est_mat_q == semidefinite(x_size);
            est_mat_r == semidefinite(y_size);
    cvx_end;

    % 返回值
    cov_zr = blkdiag(est_mat_q, est_mat_r);

end

function cov_zrt = zrtInnovationMomentEstimation(est_rt_moment, est_tt_moment, mat_a_sim, mat_c_sim, innovation_mat_k)
% 输入协方差 - 存在互协方差
% 由于w的初值未知, 故使用解方程的方法求解
    
    % 参数计算
    w_size = size(mat_a_sim, 1);
    r_size = size(mat_c_sim, 1);
    t_size = size(est_tt_moment, 1);
    order = length(est_rt_moment)-1;
    
    % 估计wt
    % 参数计算
    est_wt_mat_q = zeros(order*r_size, w_size);
    location_base = 0;
    for iter_wt = 1:order
        est_wt_mat_q(location_base+1:location_base+r_size, :) = mat_c_sim*mpower(mat_a_sim, iter_wt-1)*innovation_mat_k;  % 注意此处!
        location_base = location_base + r_size;
    end
    est_wt_mat_a = kron(eye(t_size), est_wt_mat_q);
    est_wt_vec_b = zeros(order*r_size, t_size);
    location_base = 0;
    for iter_wt = 1:order
        est_wt_vec_b(location_base+1:location_base+r_size, :) = est_rt_moment{iter_wt+1};
        location_base = location_base + r_size;
    end
    est_wt_vec_b = reshape(est_wt_vec_b, [], 1);
    % 最小二乘估计
    est_wt = lsqminnorm(est_wt_mat_a, est_wt_vec_b);
    est_wt = reshape(est_wt, [w_size t_size]);

    % 估计tt
    est_tt = est_tt_moment;
    % 估计rt
    est_rt = est_rt_moment{1};
    % 返回值
    cov_zrt = [est_wt; est_rt; est_tt];

end

function cov_zrt = zrtRealMomentEstimation(~, est_tt_moment, mat_a_sim, mat_c_sim)
% 输入协方差 - 不存在互协方差
% 直接得到
    % 参数计算
    w_size = size(mat_a_sim, 1);
    r_size = size(mat_c_sim, 1);
    t_size = size(est_tt_moment, 1);
    
    % 不估计wt和rt
    est_wt = zeros(w_size, t_size); est_rt = zeros(r_size, t_size);
    % 估计tt
    est_tt = est_tt_moment;
    % 返回值
    cov_zrt = [est_wt; est_rt; est_tt];

end

%#ok<*DEFNU> 
%#ok<*EQEFF>
%#ok<*VUNUS>

