%% 生成测试信号
% 基本参数
x_size = 2; y_size = 1; xy_size = x_size + y_size;
scale_mu = 3; scale_cov = 10;
sim_size = 2e5;
% 生成种子, 保证可复现
rng('shuffle'); seed = randi(intmax('uint32'), 'uint32');
cov_seed = RandStream('dsfmt19937', 'Seed', seed);
init_seed = RandStream('dsfmt19937', 'Seed', seed);
noise_seed = RandStream.create('mrg32k3a', 'NumStreams', xy_size, 'Seed', seed, 'CellOutput', true);
% 随机生成系统参数, 并保证系统渐近稳定
is_stable = false;
while ~is_stable
    temp_sys = drss(x_size, y_size, 0);
    is_stable = max(abs(eig(temp_sys.A, 'vector'))) ~= 1;
end
mat_a = temp_sys.A; mat_c = temp_sys.C;
% 随机生成状态初值期望和方差
mu_x0 = scale_mu*rand(cov_seed, x_size, 1);
cov_x0 = scale_cov*semidefMatrixBuilder(cov_seed, x_size);
% 随机生成噪声方差
cov_x = scale_cov*semidefMatrixBuilder(cov_seed, x_size);
cov_y = scale_cov*semidefMatrixBuilder(cov_seed, y_size);
cov_xy = blkdiag(cov_x, cov_y);
% 生成噪声, t = 0~N
noise = zeros(xy_size, sim_size+1);
for iter = 1:xy_size, noise(iter, :) = randn(noise_seed{iter}, 1, sim_size+1); end
noise = sqrtm(cov_xy)*noise;
wn = noise(1:x_size, :);
vn = noise(x_size+1:end, :);
% 生成初值, t = 0
x0 = zeros(x_size, 1);
for iter = 1:x_size, x0(iter) = randn(init_seed); end
x0 = mu_x0 + sqrtm(cov_x0)*x0;
% 生成DT ARMA数据, t = 1~N
xn = zeros(x_size, sim_size); xn(:, 1) = mat_a*x0 + wn(:, 1);
yn = zeros(y_size, sim_size);
for iter = 1:sim_size
    if iter < sim_size, xn(:, iter+1) = mat_a*xn(:, iter) + wn(:, iter+1); end
    yn(:, iter) = mat_c*xn(:, iter) + vn(:, iter+1);
end

% 稳态数据
% 计算稳态大致的起点, 以初值趋于零为依据
temp_init_response = x0;
moment_norm0 = norm(temp_init_response, 2);
for iter_loc = 1:sim_size
    temp_init_response = mat_a * temp_init_response;
    if norm(temp_init_response, 2) < 1e-6*moment_norm0, break; end
end
% 从下一个点开始采样
moment_yn = yn(:, iter_loc+1:end);
moment_yn_size = size(moment_yn, 2);


%% 使用EM算法求解
% 参数选取
em_yn_size = 1e4;  % 由于速度较慢, 使用较少的采样点
turns = 100;  % 执行轮数
% 调用函数
[em_mu_x0, em_cov_x0, em_cov_x, em_cov_y, em_xn, em_residual] = emEstimation(yn(:, 1:em_yn_size), mat_a, mat_c, 'end_condition', 'iteration', 'turn', turns);
em_cov_xy = blkdiag(em_cov_x, em_cov_y);


%% 使用ALS直接求解 (不论有无唯一解, 用伪逆给出最小方差解)
% 参数选取
moment_size = 1e3;
% 计算二阶矩
moment_yy = zeros(y_size, y_size, moment_size);
for iter_i = 0:moment_size-1
    moment_yy(:, :, iter_i+1) = momentEstimation(moment_yn, moment_yn, iter_i, 0);
end
% 调用函数
als_cov_xy = alsEstimation(moment_yy, mat_a, mat_c);


%% 结果分析&可视化
% 稳态的噪声自协方差和输出自协方差计算
verify_size = 50;
[original_moments, original_hnorms] = anaCov(cov_xy, mat_a, mat_c, verify_size);
[em_moments, em_hnorms] = anaCov(em_cov_xy, mat_a, mat_c, verify_size);
[als_moments, als_hnorms] = anaCov(als_cov_xy, mat_a, mat_c, verify_size);
disp(['Relatively MSE of cov(r): ', mat2str(mean((em_hnorms-original_hnorms).^2)/mean((original_hnorms).^2), 4), ', ', mat2str(mean((als_hnorms-original_hnorms).^2)/mean((original_hnorms).^2), 4)]);
disp(['Relatively MSE of cov(w-v): ', mat2str(mean((em_cov_xy-cov_xy).^2, 'all')/mean((cov_xy).^2, 'all'), 4), ', ', mat2str(mean((als_cov_xy-cov_xy).^2, 'all')/mean((cov_xy).^2, 'all'), 4)]);

% 绘图
% 状态估计相对误差
% figure; semilogy((abs((em_xn(:, 2:end)-xn)./xn)).'); legend;
% 迭代误差
figure; plot(em_residual);
% 稳态的噪声自协方差和输出自协方差
figure;
subplot(1, 2, 1);
semilogx(cell2mat(original_moments), '-.o'); hold on;
semilogx(cell2mat(em_moments), '-.x'); hold on;
legend;
subplot(1, 2, 2);
semilogx(cell2mat(original_moments), '-.o'); hold on;
semilogx(cell2mat(als_moments), '-.x'); hold on;
legend;


%% 针对LTI模型的EM算法
function [est_mu_x0, est_cov_x0, est_cov_x, est_cov_y, est_xn, est_residual] = emEstimation(yn, mat_a, mat_c, varargin)
% 隐参数为各个状态量
% 返回值同时还包括隐参数, 即(最后一轮的)状态变量的估计值
    
    % 迭代终止条件解析
    parser = inputParser;
    % 输入定义
    addParameter(parser, 'end_condition', 'residual', @(i)(ischar(i)));
    addParameter(parser, 'residual', 1e-4, @(i)(isnumeric(i)&&isscalar(i)));
    addParameter(parser, 'turn', 1e5, @(i)(isnumeric(i)&&isscalar(i)));
    parse(parser, varargin{:});
    end_type = parser.Results.end_condition;  % 迭代终止条件
    end_residual = parser.Results.residual;  % 终止残差
    end_turn = parser.Results.turn;  % 终止轮次

    % 参数准备
    turn = 0;
    end_turn = min(1e5, end_turn);  % 最多迭代次数
    x_size = size(mat_a, 1); y_size = size(mat_c, 1);
    sample_size = size(yn, 2);

    % 随机数种子
    rng('shuffle'); seed = randi(intmax('uint32'), 'uint32');
    init_seed = RandStream('dsfmt19937', 'Seed', seed);

    % 生成一组初始(迭代)解
    est_mu_x0 = zeros(x_size, 1);
    est_cov_x0 = semidefMatrixBuilder(init_seed, x_size);
    est_cov_x = semidefMatrixBuilder(init_seed, x_size);
    est_cov_y = semidefMatrixBuilder(init_seed, y_size);

    % 迭代初始化, t = 0~N
    yn_pad = [zeros(y_size, 1) yn];  % 为了后面程序格式规范
    conexp_xn = zeros(x_size, sample_size+1);
    conexp_cov_xn_now_now = zeros(x_size, x_size, sample_size+1);
    conexp_cov_xn_now_past = zeros(x_size, x_size, sample_size+1);  % 3rd-index看的是now而不是past
    kalman_xn_now_now = zeros(x_size, sample_size+1);
    kalman_cov_xn_now_now = zeros(x_size, x_size, sample_size+1);
    kalman_cov_xn_now_past = zeros(x_size, x_size, sample_size+1);  % 3rd-index看的是now而不是past

    % 返回残差曲线
    est_residual = zeros(end_turn, 1);

    % EM算法迭代
    while true

        %%% expectation:  在显参数和观测值条件下, 关于隐参数的期望
        % 隐参数的条件期望
        % 最优隐参数估计
        % Kalman初始化, t = 0~N
        kalman_xn_now_now(:, 1) = est_mu_x0;
        kalman_cov_xn_now_now(:, :, 1) = est_cov_x0;
        for iter_k = 2:sample_size+1  % t = 1~N
            % Kalman 状态估计
            % 计算
            kalman_xn_now_past = mat_a * kalman_xn_now_now(:, iter_k-1);
            kalman_cov_xn_now_past(:, :, iter_k) = mat_a * kalman_cov_xn_now_now(:, :, iter_k-1) * mat_a.' + est_cov_x;
            kalman_gain_now = kalman_cov_xn_now_past(:, :, iter_k) * mat_c.' / (mat_c * kalman_cov_xn_now_past(:, :, iter_k) * mat_c.' + est_cov_y);
            kalman_xn_now_now(:, iter_k) = kalman_xn_now_past + kalman_gain_now * (yn_pad(:, iter_k) - mat_c * kalman_xn_now_past);
            kalman_cov_xn_now_now(:, :, iter_k) = kalman_cov_xn_now_past(:, :, iter_k) - kalman_gain_now * mat_c * kalman_cov_xn_now_past(:, :, iter_k);
        end

        % 隐参数计算
        % 隐参数初始化
        conexp_xn(:, sample_size+1) = kalman_xn_now_now(:, sample_size+1);
        conexp_cov_xn_now_now(:, :, sample_size+1) = kalman_cov_xn_now_now(:, :, sample_size+1);
        conexp_cov_xn_now_past(:, :, sample_size+1) = (eye(x_size) - kalman_gain_now * mat_c) * mat_a * kalman_cov_xn_now_now(:, :, sample_size);
        % 反向迭代
        for iter_k = sample_size+1:-1:2  % t = 1~N
            % 辅助矩阵
            conexp_j_past = kalman_cov_xn_now_now(:, :, iter_k-1) * mat_a / kalman_cov_xn_now_past(:, :, iter_k);
            if iter_k >= 3, conexp_j_past_past = kalman_cov_xn_now_now(:, :, iter_k-2) * mat_a / kalman_cov_xn_now_past(:, :, iter_k-1); end
            % 计算
            conexp_xn(:, iter_k-1) = kalman_xn_now_now(:, iter_k-1) + conexp_j_past * (conexp_xn(:, iter_k) - mat_a * kalman_xn_now_now(:, iter_k-1));
            conexp_cov_xn_now_now(:, :, iter_k-1) = kalman_cov_xn_now_now(:, :, iter_k-1) + conexp_j_past * (conexp_cov_xn_now_now(:, :, iter_k) - kalman_cov_xn_now_past(:, :, iter_k)) * conexp_j_past.';
            if iter_k >= 3, conexp_cov_xn_now_past(:, :, iter_k-1) = kalman_cov_xn_now_now(:, :, iter_k-1) * conexp_j_past_past  + conexp_j_past * (conexp_cov_xn_now_past(:, :, iter_k) - mat_a * kalman_cov_xn_now_now(:, :, iter_k-1)) * conexp_j_past_past.'; end
        end
        % 条件期望计算
        conexp_para_cov_x_mat1 = sum(conexp_cov_xn_now_now(:, :, 1:sample_size), 3) + conexp_xn(:, 1:sample_size) * conexp_xn(:, 1:sample_size).';
        conexp_para_cov_x_mat2 = sum(conexp_cov_xn_now_past(:, :, 2:sample_size+1), 3) + conexp_xn(:, 2:sample_size+1) * conexp_xn(:, 1:sample_size).';
        conexp_para_cov_x_mat3 = sum(conexp_cov_xn_now_now(:, :, 2:sample_size+1), 3) + conexp_xn(:, 2:sample_size+1) * conexp_xn(:, 2:sample_size+1).';
        conexp_para_cov_x = conexp_para_cov_x_mat3 - conexp_para_cov_x_mat2 * mat_a.' - mat_a * conexp_para_cov_x_mat2.' + mat_a * conexp_para_cov_x_mat1 * mat_a.';
        conexp_para_cov_y_diff = yn_pad - mat_c * conexp_xn;
        conexp_para_cov_y = mat_c * sum(conexp_cov_xn_now_now(:, :, 2:sample_size+1), 3) * mat_c.' + conexp_para_cov_y_diff(:, 2:sample_size+1) * conexp_para_cov_y_diff(:, 2:sample_size+1).';

        %%% maximization: 对显参数的最大化
        argmax_mu_x0 = conexp_xn(:, 1);
        argmax_cov_x0 = conexp_cov_xn_now_now(:, :, 1);
        argmax_cov_x = conexp_para_cov_x./sample_size;
        argmax_cov_y = conexp_para_cov_y./sample_size;

        %%% 迭代条件检查
        turn = turn + 1;
        % 计算均方根误差, 并存储
        residual_1 = sum((est_mu_x0 - argmax_mu_x0).^2,'all');
        residual_2 = sum((est_cov_x0 - argmax_cov_x0).^2,'all');
        residual_3 = sum((est_cov_x - argmax_cov_x).^2,'all');
        residual_4 = sum((est_cov_y - argmax_cov_y).^2,'all');
        est_residual(turn) = sqrt((residual_1 + residual_2 + residual_3 + residual_4)/(x_size + x_size*x_size*(sample_size+1) + y_size*y_size*sample_size));
        % 判断终止条件 (注意: 最后一次的迭代值会被舍弃)
        switch end_type
            case 'residual'
                if (est_residual(turn) < end_residual) || isnan(est_residual(turn)) || (turn >= end_turn)
                    disp(['第', num2str(turn), '次EM算法参数估计完成']);
                    if isnan(est_residual(turn)) || (turn >= max_turn), disp('模型不收敛');
                    else, disp(['模型收敛, 参数均方根误差为: ', num2str(est_residual(turn))]); end
                    est_xn = conexp_xn;
                    break;
                end
            case 'iteration'
                if turn >= end_turn
                    disp(['第', num2str(turn), '次EM算法参数估计完成']);
                    disp(['参数均方根误差为: ', num2str(est_residual(turn))]);
                    est_xn = conexp_xn;
                    break;
                end
        end
        % 不满足跳出条件, 继续循环
        est_mu_x0 = argmax_mu_x0;
        est_cov_x0 = argmax_cov_x0;
        est_cov_x = argmax_cov_x;
        est_cov_y = argmax_cov_y;

    end

end

%% 简化ALS估计法
function cov_real_zr = alsEstimation(moment_yy, mat_a, mat_c)
% 输出-状态协方差 - 不存在互协方差
% 此处计算真实矩

    % 参数计算
    x_size = size(mat_a, 1);
    y_size = size(mat_c, 1);
    order = size(moment_yy, 3);

    % 估计ww_vv
    % 参数计算 - M, N
    est_mat_m = zeros(order*y_size, x_size);
    location_base = 0;
    for iter_moment = 1:order
        est_mat_m(location_base+1:location_base+y_size, :) = mat_c*mpower(mat_a, iter_moment-1);
        location_base = location_base + y_size;
    end
    est_mat_n = zeros(order*y_size, y_size);
    est_mat_n(1:y_size, 1:y_size) = eye(y_size);
    % 参数计算 - E, F
    est_mat_e = kron(mat_c, est_mat_m) / (eye(x_size^2) - kron(mat_a, mat_a));
    est_mat_f = kron(eye(y_size), est_mat_n);
    % 参数计算 - b
    est_vec_b = zeros(order*y_size, y_size);
    location_base = 0;
    for iter_moment = 1:order
        est_vec_b(location_base+1:location_base+y_size, :) = moment_yy(:, :, iter_moment);
        location_base = location_base + y_size;
    end
    est_mat_a = [est_mat_e est_mat_f];
    est_vec_b = reshape(est_vec_b, [], 1);
    % 最小二乘估计
    est_para = lsqminnorm(est_mat_a, est_vec_b);
    est_cov_zz = est_para(1:x_size*x_size); est_cov_zz = reshape(est_cov_zz, [x_size, x_size]);
    est_cov_vv = est_para(x_size*x_size+1:end); est_cov_vv = reshape(est_cov_vv, [y_size, y_size]);
    % SDP估计 (需要cvx工具箱)
    % cvx_begin quiet
    %     variable est_cov_zz(z_size, z_size) symmetric
    %     variable est_cov_vv(r_size, r_size) symmetric
    %     minimize(norm(est_mat_e*vec(est_cov_zz)+est_mat_f*vec(est_cov_vv)-est_vec_b)) %, Inf))
    %     subject to
    %         est_cov_zz == semidefinite(z_size)
    %         est_cov_vv == semidefinite(r_size)
    % cvx_end

    % 返回值
    cov_real_zr = [est_cov_zz zeros(x_size, y_size); zeros(y_size, x_size) est_cov_vv];

end

%% 正定矩阵生成(0~1)
function ret_mat = semidefMatrixBuilder(seed, n)
    ret_mat = rand(seed, n, n);
    ret_mat = 0.5.*(ret_mat+ret_mat.');
    ret_mat = ret_mat + n.*eye(n);
end


%% 矩估计
function moment_estimated = momentEstimation(an, bn, order, offset)
    % 计算Expectation{a_{n+order}*b_n}, offset为开头截断的长度
    moment_estimation_samples = size(an, 2);
    moment_estimated = (1/(moment_estimation_samples-order-offset)).*(an(:, offset+order+1:end)*(bn(:, offset+1:end-order).'));
end


%% 对比稳态矩
function [moments, hnorms] = anaCov(covariance, mat_a, mat_c, out_size)
%ANACOV 根据协方差和(渐近稳定)系统矩阵计算稳态输出矩(输入的矩仅包含状态和输出)
% 同时返回计算结果和Hnorm-2对比
    
    % 参数计算
    z_size = size(mat_a, 1);
    cov_zz = covariance(1:z_size, 1:z_size);
    cov_zr = covariance(1:z_size, z_size+1:end);
    cov_rr = covariance(z_size+1:end, z_size+1:end);

    % 输出准备
    moments = cell(out_size, 1);
    hnorms = zeros(out_size, 1);

    % 计算
    stable_zz = dlyap(mat_a, cov_zz);
    stable_zr1 = mat_a*stable_zz*(mat_c.') + cov_zr;
    % 零阶矩
    moments{1} = mat_c*stable_zz*(mat_c.') + cov_rr;
    [hnorms(1), ~] = anaNorm(moments{1});
    % 非零阶矩
    for iter = 1:out_size-1
        moments{iter+1} = mat_c*mpower(mat_a, iter-1)*stable_zr1;
        [hnorms(iter+1), ~] = anaNorm(moments{iter+1});
    end

end


%% 参考内容
% https://spaces.ac.cn/archives/5239
% https://zhuanlan.zhihu.com/p/438941218
% https://doi.org/10.1111/j.1467-9892.1982.tb00349.x
