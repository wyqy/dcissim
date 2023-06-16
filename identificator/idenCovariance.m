function covariance_struct = idenCovariance(yk, uk, mat_v, para_est, als_est_type, cov_cross_type, cov_order, cov_order_type)
%IDENCOVARIANCE 辨识方差参数
% 针对存在互协方差与否提供不同的估计
% 存在互协方差: 返回cov_innovation_all (zrt), cov_innovation_kalman
% 不存在互协方差: 返回cov_real_all (zrt)

    % 计算各阶矩
    % 参数计算
    x_size = size(para_est.A, 1); y_size = size(yk, 1); u_size = size(uk, 1);
    N = size(yk, 2); T = size(mat_v, 2);

    % 激励计算
    vk = repmat(mat_v, [1, floor(N/T)]);
    % 输入激励
    ut_est = para_est.U*vk;
    % 输出激励
    % 1. 状态方程噪声
    % x0_est = para_est.X*para_est.v0;
    % [~, yt_est, ~] = plantModel(para_est, x0_est, uk, zeros(x_size+y_size+u_size, sample_size));
    % 2. 直接计算
    yt_est = para_est.Y*vk;
    % 误差计算
    rk = yk - yt_est;
    tk = uk - ut_est;

    % 阶数估算
    % 估算阶数上限, 按矩阵谱半径估计
    test_order_max = ceil(log(1e-5)/log(max(abs(eig(para_est.A)))));
    test_rr_moment_norm = zeros(test_order_max+1, 1);
    % 询问或者直接估计/固定值
    switch cov_order_type
        case 'ask'
            % 计算
            for iter_test = 0:test_order_max, test_rr_moment_norm(iter_test+1) = norm(estimationCaller(rk, rk, iter_test, 0), 'fro'); end
            % 绘图-询问
            disp('The Frobenius norms of R''s auto-covariance are:')
            fig_ask = figure;
            plot(0:test_order_max, test_rr_moment_norm);
            uiwait(fig_ask);
            moment_order = input('Enter the estimated sizes of states: ');
            % 无效输入, 估计合适的值
            if ~(isscalar(moment_order) && isnumeric(moment_order)), moment_order = orderCaller(test_rr_moment_norm); end
        case 'estimate'
            for iter_test = 0:test_order_max, test_rr_moment_norm(iter_test+1) = norm(estimationCaller(rk, rk, iter_test, 0), 'fro'); end
            moment_order = orderCaller(test_rr_moment_norm);
        case 'fixed'
            moment_order = cov_order;
        case 'null'
            cov_all = zeros(x_size+y_size+u_size);
            kalman_gain = zeros(x_size, y_size);
            covariance_struct = struct('cov', cov_all, 'kalman', kalman_gain);
            return;
        otherwise, moment_order = 1;
    end

    % 各阶矩计算
    % 准备
    moment_rr_est = cell(moment_order+1, 1);
    moment_rt_est = cell(moment_order+1, 1);
    % 计算
    for iter_moment = 0:moment_order
        moment_rr_est{iter_moment+1} = estimationCaller(rk, rk, iter_moment, 0);
        moment_rt_est{iter_moment+1} = estimationCaller(rk, tk, iter_moment, 0);
    end
    moment_tt_est = estimationCaller(tk, tk, 0, 0);

    switch cov_cross_type
        case 'valid'
            % 估计名义协方差矩阵(innovation form)
            [inno_zr, inno_zrt, kalman_gain] = innovationEstimation(moment_rr_est, moment_rt_est, moment_tt_est, para_est.A, para_est.C);
            % 参数处理
            cov_zrt_1 = inno_zrt(1:x_size+y_size, :);
            cov_zrt_2 = inno_zrt(x_size+y_size+1:end, :);
            cov_all = [inno_zr cov_zrt_1; cov_zrt_1.' cov_zrt_2];
        case 'null'
            % 估计zr的协方差矩阵(real form)
            switch als_est_type
                case 'simple', cov_zr = zrRealMomentEstimation_nocross_simple(moment_rr_est, para_est.A, para_est.C);
                case 'classical', cov_zr = zrRealMomentEstimation_nocross_classical(rk, moment_order, para_est.A, para_est.C);
            end
            % 估计zrt的协方差矩阵(real form)
            cov_zrt_2 = zrtRealMomentEstimation_nocross(moment_rt_est, moment_tt_est, para_est.A, para_est.C);
            % 参数处理
            cov_zrt_2 = cov_zrt_2(x_size+y_size+1:end, :);
            cov_all = [cov_zr zeros(x_size+y_size, u_size); zeros(u_size, x_size+y_size) cov_zrt_2];
            % 其它填0
            kalman_gain = zeros(x_size, y_size);
        otherwise, cov_all = zeros(x_size+y_size+u_size); kalman_gain = zeros(x_size, y_size);
    end

    % 返回数据
    if min(eig(cov_zr)) < 0 || min(eig(cov_zrt_2)) < 0
        cov_all = zeros(x_size+y_size+u_size); kalman_gain = zeros(x_size, y_size);
    end
    covariance_struct = struct('cov_all', cov_all, 'kalman', kalman_gain, ...
        'cov_xx', cov_all(1:x_size, 1:x_size), ...
        'cov_xy', cov_all(1:x_size, x_size+1:x_size+y_size), ...
        'cov_yy', cov_all(x_size+1:x_size+y_size, x_size+1:x_size+y_size), ...
        'cov_xu', cov_all(1:x_size, x_size+y_size+1:end), ...
        'cov_yu', cov_all(x_size+1:x_size+y_size, x_size+y_size+1:end), ...
        'cov_uu', cov_all(x_size+y_size+1:end, x_size+y_size+1:end));
    
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

function moment_estimated = estimationCaller(an, bn, order, offset)
    % 计算Expectation{a_{n+order}*b_n}, offset为开头截断的长度
    moment_estimation_samples = size(an, 2);
    moment_estimated = (1/(moment_estimation_samples-order-offset)).*(an(:, offset+order+1:end)*(bn(:, offset+1:end-order).'));
end

function [cov_zr, cov_zrt, cov_kalman] = innovationEstimation(moment_rr_est, moment_rt_est, moment_tt_est, mat_a_est, mat_c_est)
% 名义(新息)协方差估计 - 存在互协方差
% 输入-状态协方差: forward innovaiton form对应的矩
% 输入协方差: 由于w的初值未知, 故使用解方程的方法求解

    % 参数计算
    w_size = size(mat_a_est, 1);
    r_size = size(mat_c_est, 1);
    t_size = size(moment_tt_est, 1);
    order = length(moment_rr_est)-1;
    est_rr0 = moment_rr_est{1};

    % 估计wr_1
    % 参数计算
    est_wr1_mat_q = zeros(order*r_size, w_size);
    loc_base = 0;
    for iter_wr1 = 1:order
        est_wr1_mat_q(loc_base+1:loc_base+r_size, :) = mat_c_est*mpower(mat_a_est, iter_wr1-1);
        loc_base = loc_base + r_size;
    end
    est_wr1_mat_a = kron(eye(r_size), est_wr1_mat_q);
    est_wr1_vec_b = zeros(order*r_size, r_size);
    loc_base = 0;
    for iter_wr1 = 1:order
        est_wr1_vec_b(loc_base+1:loc_base+r_size, :) = moment_rr_est{iter_wr1+1};
        loc_base = loc_base + r_size;
    end
    est_wr1_vec_b = reshape(est_wr1_vec_b, [], 1);
    % 最小二乘估计
    est_wr1 = lsqminnorm(est_wr1_mat_a, est_wr1_vec_b);
    est_wr1 = reshape(est_wr1, [w_size r_size]);

    % 计算ww的自协方差和kalman增益
    [est_ww0, cov_kalman, ~, ~] = idare(mat_a_est.', mat_c_est.', zeros(w_size), -est_rr0, -est_wr1, []);
    cov_kalman = cov_kalman.';
    
    % 估计wt
    % 参数计算
    est_wt_mat_q = zeros(order*r_size, w_size);
    loc_base = 0;
    for iter_wt = 1:order
        est_wt_mat_q(loc_base+1:loc_base+r_size, :) = mat_c_est*mpower(mat_a_est, iter_wt-1)*cov_kalman;  % 注意此处!
        loc_base = loc_base + r_size;
    end
    est_wt_mat_a = kron(eye(t_size), est_wt_mat_q);
    est_wt_vec_b = zeros(order*r_size, t_size);
    loc_base = 0;
    for iter_wt = 1:order
        est_wt_vec_b(loc_base+1:loc_base+r_size, :) = moment_rt_est{iter_wt+1};
        loc_base = loc_base + r_size;
    end
    est_wt_vec_b = reshape(est_wt_vec_b, [], 1);
    % 最小二乘估计ww的名义方差
    est_wt = lsqminnorm(est_wt_mat_a, est_wt_vec_b);
    est_wt = reshape(est_wt, [w_size t_size]);

    % 估计ττ的真实方差
    est_tt = moment_tt_est;
    % 估计rτ的名义方差
    est_rt = moment_rt_est{1};
    
    % 计算qsr
    est_cov_q = est_ww0 - mat_a_est*est_ww0*(mat_a_est.');
    est_cov_s = est_wr1 - mat_a_est*est_ww0*(mat_c_est.');
    est_cov_r = est_rr0 - mat_c_est*est_ww0*(mat_c_est.');
    % 返回值
    cov_zr = [est_cov_q est_cov_s; est_cov_s.' est_cov_r];
    cov_zrt = [est_wt; est_rt; est_tt];

end

function cov_zr = zrRealMomentEstimation_nocross_simple(est_rr_moment, mat_a_sim, mat_c_sim)
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
    loc_base = 0;
    for iter_moment = 1:order-1
        est_para_m(loc_base+1:loc_base+y_size, :) = mat_c_sim*mpower(mat_a_sim, iter_moment);
        loc_base = loc_base + y_size;
    end
    est_para_m = kron(mat_c_sim, est_para_m);
    est_para_b = zeros((order-1)*y_size, y_size);
    loc_base = 0;
    for iter_moment = 1:order-1
        est_para_b(loc_base+1:loc_base+y_size, :) = est_rr_moment{iter_moment+1};
        loc_base = loc_base + y_size;
    end
    est_para_b = reshape(est_para_b, [(order-1)*y_size*y_size 1]);
    % LS估计
    % est_para_p = lsqminnorm(est_para_m, est_para_b);
    % est_mat_p = reshape(est_para_p, [x_size x_size]);
    % SDP估计
    % cvx_begin sdp quiet;
    %     variable est_mat_p(x_size, x_size) symmetric;
    %     variable err_mat_m((order-1)*y_size*y_size, x_size*x_size)
    %     variable err_vec_b((order-1)*y_size*y_size, 1)
    %     minimize(norm([err_mat_m err_vec_b], "fro"))
    %     subject to;
    %         (est_para_m + err_mat_m) * vec(est_mat_p) == est_para_b + err_vec_b;
    %         est_mat_p == semidefinite(x_size);
    %         est_rr0 - mat_c_sim * est_mat_p * mat_c_sim.' == semidefinite(y_size);
    %         est_mat_p - mat_a_sim * est_mat_p * mat_a_sim.' == semidefinite(x_size);
    % cvx_end;
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

function cov_zr = zrRealMomentEstimation_nocross_classical(rk, order, mat_a_est, mat_c_est) 
% 输出-状态协方差 - 不存在互协方差
% 进行状态估计

    % 参数计算
    x_size = size(mat_a_est, 1); 
    y_size = size(mat_c_est, 1);
    sample_size = size(rk, 2);
    % K, L, F矩阵计算
    riccati_q = 1e-5*eye(x_size); riccati_r = 1e-5*eye(y_size);
    mat_k_sim = idare(mat_a_est.', mat_c_est.', riccati_q, riccati_r, [], []);
    mat_k_sim = mat_k_sim.';
    mat_l_sim = mat_a_est*mat_k_sim;
    mat_f_sim = mat_a_est - mat_l_sim*mat_c_est;

    % 随机状态估计 (零状态初始化)
    en = zeros(y_size, sample_size);
    xn = zeros(x_size, 1);
    for iter_est = 1:sample_size
        en(:, iter_est) = rk(:, iter_est) - mat_c_est*xn;
        xn = mat_a_est*xn + mat_l_sim*en(:, iter_est);
    end
    % innovation矩估计
    moment_en = cell(order, 1);
    for iter_order = 0:order
        moment_en{iter_order+1} = estimationCaller(en, en, iter_order, 0);
    end

    % 估计P的一种估计值
    % 参数计算
    % L, N计算
    est_para_l = zeros(order*y_size, y_size);
    loc_base = 0;
    for iter_order = 1:order
        est_para_l(loc_base+1:loc_base+y_size, :) = mat_c_est*mpower(mat_f_sim, iter_order-1);
        loc_base = loc_base + y_size;
    end
    est_para_n = zeros(order*y_size, y_size);
    loc_base = 0;
    est_para_n(loc_base+1:loc_base+y_size, :) = eye(y_size);
    loc_base = loc_base + y_size;
    for iter_order = 2:order
        est_para_n(loc_base+1:loc_base+y_size, :) = mat_c_est*mpower(mat_f_sim, iter_order-2)*mat_l_sim;
        loc_base = loc_base + y_size;
    end
    % E, G计算
    est_para_e = kron(mat_c_est, est_para_l) / (eye(x_size^2) - kron(mat_f_sim, mat_f_sim));
    est_para_g = est_para_e*kron(mat_l_sim, mat_l_sim) + kron(eye(y_size), est_para_n);
    % A, b计算
    est_para_b = zeros(order*y_size, y_size);
    loc_base = 0;
    for iter_order = 1:order
        est_para_b(loc_base+1:loc_base+y_size, :) = moment_en{iter_order};
        loc_base = loc_base + y_size;
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

function cov_zrt = zrtRealMomentEstimation_nocross(~, est_tt_moment, mat_a_est, mat_c_est)
% 输入协方差 - 不存在互协方差
% 直接得到
    % 参数计算
    w_size = size(mat_a_est, 1);
    r_size = size(mat_c_est, 1);
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

