function covariance_struct = idenCovariance(yn, un, vn, u_isim, mat_a_sim, mat_b_sim, mat_c_sim, mat_d_sim, cov_cross_type, cov_order, cov_order_type)
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
    rn = yn - res_est_yn;
    tn = un - res_est_un;

    % 阶数估算
    % 准备计算0~0.5*samples的阶数
    test_order_max = fix(signal_sample/2);
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
            cov_zr = zrRealMomentEstimation(est_rr_moment, mat_a_sim, mat_c_sim);
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

    % 嵌套函数
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

function cov_zrt = zrtRealMomentEstimation(est_rt_moment, est_tt_moment, mat_a_sim, mat_c_sim)
% 输入协方差 - 不存在互协方差
% 解和上个函数类似的方程

    % 参数计算
    w_size = size(mat_a_sim, 1);
    r_size = size(mat_c_sim, 1);
    t_size = size(est_tt_moment, 1);
    % order = length(est_rt_moment)-1;
    
    % % 估计wt
    % % 参数计算
    % est_wt_mat_q = zeros(order*r_size, w_size);
    % location_base = 0;
    % for iter_wt = 1:order
    %     est_wt_mat_q(location_base+1:location_base+r_size, :) = mat_c_sim*mpower(mat_a_sim, iter_wt-1);
    %     location_base = location_base + r_size;
    % end
    % est_wt_mat_a = kron(eye(t_size), est_wt_mat_q);
    % est_wt_vec_b = zeros(order*r_size, t_size);
    % location_base = 0;
    % for iter_wt = 1:order
    %     est_wt_vec_b(location_base+1:location_base+r_size, :) = est_rt_moment{iter_wt+1};
    %     location_base = location_base + r_size;
    % end
    % est_wt_vec_b = reshape(est_wt_vec_b, [], 1);
    % % 最小二乘估计
    % est_wt = lsqminnorm(est_wt_mat_a, est_wt_vec_b);
    % est_wt = reshape(est_wt, [w_size t_size]);
    % % 估计rt
    % est_rt = est_rt_moment{1};
    
    % 不估计wt和rt
    est_wt = zeros(w_size, t_size); est_rt = zeros(r_size, t_size);
    % 估计tt
    est_tt = est_tt_moment;
    % 返回值
    cov_zrt = [est_wt; est_rt; est_tt];

end

%#ok<*EQEFF>
%#ok<*VUNUS>

