function [cov_zr_innovation, cov_zrt_innovation, cov_zrt_real, innovation_mat_k, innovation_mat_covariance] = idenCovariance(yn, un, vn, y_isim, u_isim, mat_a_sim, mat_c_sim, cov_order, cov_order_type)
%IDENCOVARIANCE 辨识方差参数

    % 计算各阶矩
    % 参数计算
    x_size = size(mat_a_sim, 1);
    y_size = size(yn, 1);
    u_size = size(un, 1);
    samples = size(yn, 2);
    % 扰动量计算
    rn = yn - y_isim*vn;
    tn = un - u_isim*vn;

    % 阶数估算
    % 准备计算0~0.5*samples的阶数
    test_order_max = fix(samples/2);
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
            cov_zr_innovation = zeros(x_size+y_size, x_size+y_size);
            cov_zrt_innovation = zeros(x_size+y_size+u_size, u_size);
            cov_zrt_real = zeros(x_size+y_size+u_size, u_size);
            innovation_mat_k = zeros(x_size, y_size);
            innovation_mat_covariance = zeros(y_size, y_size);
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

    % 估计zr的协方差矩阵(innovation form)
    [cov_zr_innovation, innovation_mat_k] = zrInnovationMomentEstimation(est_rr_moment, mat_a_sim, mat_c_sim);
    innovation_mat_covariance = cov_zr_innovation(x_size+1:end, x_size+1:end);
    % 估计zrt的协方差矩阵(innovation form)
    cov_zrt_innovation = zrtInnovationMomentEstimation(est_rt_moment, est_tt_moment, mat_a_sim, mat_c_sim, innovation_mat_k);
    % 估计zrt的协方差矩阵(original form)
    cov_zrt_real = zrtRealMomentEstimation(est_rt_moment, est_tt_moment, mat_a_sim, mat_c_sim);

    % 嵌套函数
    function moment_estimated = estimationCaller(an, bn, order, offset)
        % 计算Expectation{a_{n+order}*b_n}, offset为开头截断的长度
        moment_estimation_samples = size(an, 2);
        moment_estimated = (1/(moment_estimation_samples-order-offset)).*(an(:, offset+order+1:end)*(bn(:, offset+1:end-order).'));
    end
    function order_estimated = orderCaller(moment)
        % 估算order
        order_max = length(moment)-1;
        for iter_order = 0:order_max
            if moment(iter_order+1) < 1e-4*moment(2)
                order_estimated = iter_order; break;
            end
        end
    end

end

function [cov_zr_innovation, est_mat_k] = zrInnovationMomentEstimation(est_rr_moment, mat_a_sim, mat_c_sim)
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
    [est_ww0, est_mat_k, ~, ~] = idare(mat_a_sim.', mat_c_sim.', zeros(w_size), -est_rr0, -est_wr1, []);
    est_mat_k = est_mat_k.';

    % 计算qsr
    est_q = est_ww0 - mat_a_sim*est_ww0*(mat_a_sim.');
    est_s = est_wr1 - mat_a_sim*est_ww0*(mat_c_sim.');
    est_r = est_rr0 - mat_c_sim*est_ww0*(mat_c_sim.');
    % 返回值
    cov_zr_innovation = [est_q est_s; est_s.' est_r];

end

function cov_zrt_innovation = zrtInnovationMomentEstimation(est_rt_moment, est_tt_moment, mat_a_sim, mat_c_sim, innovation_mat_k)
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
    cov_zrt_innovation = [est_wt; est_rt; est_tt];

end

function cov_zrt_real = zrtRealMomentEstimation(est_rt_moment, est_tt_moment, mat_a_sim, mat_c_sim)
% 解和上个函数类似的方程

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
        est_wt_mat_q(location_base+1:location_base+r_size, :) = mat_c_sim*mpower(mat_a_sim, iter_wt-1);
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
    cov_zrt_real = [est_wt; est_rt; est_tt];

end


