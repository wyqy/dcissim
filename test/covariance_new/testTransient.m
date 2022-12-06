%% 原理验证
% x_size = 4; y_size = 1;
% rng(142857, 'simdTwister');
% is_stable = false;
% while ~is_stable
%     rand_sys = drss(x_size, y_size, 0);
%     is_stable = max(abs(eig(rand_sys.A, 'vector'))) <= 0.95;
% end
% A = rand_sys.A; C = rand_sys.C;
% mat_obsv = obsv(A, C);
% % 构造矩阵 - 1
% cov_mat_1 = kron(C, mat_obsv*A);
% % disp(['rank(mat_1) = ' num2str(rank(cov_mat_1)) ', with size: ' mat2str(size(cov_mat_1))])
% % 构造矩阵 - 2
% cov_mat_2 = kron(mat_obsv*A, mat_obsv*mpower(A, x_size+1));
% % disp(['rank(mat_2) = ' num2str(rank(cov_mat_2)) ', with size: ' mat2str(size(cov_mat_2))])

%% 实验验证
%% 生成测试信号
% 基本参数
x_size = 4; y_size = 1; xy_size = x_size + y_size;
scale_cov = 10; sim_size = 1e2; exp_size = 1e4;
% 生成种子, 保证可复现
seed = uint32(271828); rng(seed, 'simdTwister');
% rng('shuffle'); seed = randi(intmax('uint32'), 'uint32');
cov_seed = RandStream('dsfmt19937', 'Seed', seed);
noise_seed = RandStream.create('mrg32k3a', 'NumStreams', xy_size*exp_size, 'Seed', seed, 'CellOutput', true);
noise_seed = reshape(noise_seed, [xy_size exp_size]);
% 随机生成系统参数, 并保证系统渐近稳定
is_stable = false;
while ~is_stable
    rand_sys = drss(x_size, y_size, 0);
    is_stable = max(abs(eig(rand_sys.A, 'vector'))) <= 0.95;
end
% rand_sys = drss(x_size, y_size, 0);
% rand_sys.A = cumsum(diag(0.8.*ones(x_size, 1)), 2);  % 人为指定A矩阵
% rand_sys.C = cumsum(ones(y_size, x_size), 2);
mat_a = rand_sys.A; mat_c = rand_sys.C;
% 随机生成噪声方差
% cov_xy = scale_cov*semidefMatrixBuilder(cov_seed, xy_size);
cov_x = scale_cov*semidefMatrixBuilder(cov_seed, x_size);
cov_y = scale_cov*semidefMatrixBuilder(cov_seed, y_size);
cov_xy = blkdiag(cov_x, cov_y);
% 生成状态和状态噪声初值, t = 0
x0 = zeros(x_size, 1); x0 = repmat(x0, [1 exp_size]);
w0 = zeros(x_size, exp_size);
for iter_exp = 1:exp_size
    for iter_x = 1:x_size
        w0(iter_x, iter_exp) = randn(noise_seed{iter_x, iter_exp});
    end
end
% 生成噪声, t = 1~N, exp = 1~E
noise = zeros(xy_size, sim_size, exp_size);
for iter_exp = 1:exp_size
    for iter_xy = 1:xy_size
        noise(iter_xy, :, iter_exp) = randn(noise_seed{iter_xy, iter_exp}, 1, sim_size);
    end
end
noise = pagemtimes(sqrtm(cov_xy), noise);
wn = noise(1:x_size, :, :);
vn = noise(x_size+1:end, :, :);
% 生成DT ARMA数据, t = 1~N
xn = zeros(x_size, sim_size, exp_size);
yn = zeros(y_size, sim_size, exp_size);
covdim = @(x) reshape(x, [size(x, 1) size(x, 3)]);
for iter_sim = 1:sim_size
    if iter_sim == 1, xn(:, 1, :) = mat_a*x0 + w0;
    else, xn(:, iter_sim, :) = mat_a*covdim(xn(:, iter_sim-1, :)) + covdim(wn(:, iter_sim-1, :)); end
    yn(:, iter_sim, :) = mat_c*covdim(xn(:, iter_sim, :)) + covdim(vn(:, iter_sim, :));
end

%% 验证
diff_norm = testEstimation(yn, cov_xy, mat_a, mat_c);

%% 计算 & 验证
est_cov_xy = tranEstimation(yn, mat_a, mat_c);

% 稳态的噪声自协方差和输出自协方差计算
verify_size = 50;
[original_moments, original_hnorms] = anaCov(cov_xy, mat_a, mat_c, verify_size);
[est_moments, est_hnorms] = anaCov(est_cov_xy, mat_a, mat_c, verify_size);
disp(['abs(eig(A)): ', mat2str(abs(eig(rand_sys.A, 'vector')), 2)]);
disp(['Relatively MSE of cov: ', mat2str(mean((est_hnorms-original_hnorms).^2)/mean((original_hnorms).^2), 4)]);

% 绘图
% 稳态的噪声自协方差和输出自协方差
figure;
semilogx(cell2mat(original_moments), '-.o'); hold on;
semilogx(cell2mat(est_moments), '-.x'); hold on;
legend('original', 'identified');

%% 暂态估计法
function est_cov_xy = tranEstimation(yn, mat_a, mat_c)
% 输出-状态协方差 - 不存在互协方差
% 此处计算真实矩[Q S; S^T R], 而P_k为状态协方差(非稳态)
%#ok<*EQEFF>

    % 参数计算
    offset = 0;  % 用以建模暂态测量一开始的不稳定点(废弃)
    x_size = size(mat_a, 1); y_size = size(mat_c, 1);
    i_size = x_size; k_size = 2*x_size;
    l_size = i_size+k_size+2;

    % 样本方差估计量, 设符号表示为E(x_{k+i}x_{k}\tr),
    % 则至少要取x*(x+y)/y个量
    % 此处直接取2*x^2个 @ (y_size, y_size, i_size, k_size)
    % 同时, 对于R的估计, 直接取2*x^2+2个点(即前面所用的点).
    est_mom_zzik = zeros(y_size, y_size, i_size, k_size);
    for iter_k = 1:k_size
        for iter_i = 1:i_size
            est_mom_zzik(:, :, iter_i, iter_k) = momentEstimation(yn, iter_i, iter_k, offset);
        end
    end
    est_mom_zzl = zeros(y_size, y_size, l_size);
    for iter_l = 1:l_size
        est_mom_zzl(:, :, iter_l) = momentEstimation(yn, 0, iter_l, offset);
    end

    % A矩阵级数计算, 0~mom_offset+i_size+k_size
    series_mat_a = zeros(x_size, x_size, offset+i_size+k_size+2);
    for iter_series = 1:size(series_mat_a, 3)
        series_mat_a(:, :, iter_series) = mpower(mat_a, iter_series-1);
    end

    % 估计Q
    % A参数计算(不唯一!)
    est_para_a = zeros(i_size*k_size*(y_size^2), x_size^2); loc_base = 0;
    for iter_k = 1:k_size
        for iter_i = 1:i_size
            temp_a = zeros(y_size^2, x_size^2);
            for iter_j = 0:iter_k-1
                temp_a = temp_a + kron(mat_c*series_mat_a(:, :, offset+iter_j+1), mat_c*series_mat_a(:, :, offset+iter_j+iter_i+1));
            end
            est_para_a(loc_base+1:loc_base+y_size^2, :) = temp_a;
            loc_base = loc_base + y_size^2;
        end
    end
    % b参数计算(要对应!)
    est_para_b = zeros(i_size*k_size*(y_size^2), 1); loc_base = 0;
    for iter_k = 1:k_size
        for iter_i = 1:i_size
            est_para_b(loc_base+1:loc_base+y_size^2, :) = reshape(est_mom_zzik(:, :, iter_i, iter_k), [], 1);
            loc_base = loc_base + y_size^2;
        end
    end
    % SDP估计
    cvx_begin sdp;% quiet;
        variable est_cov_q(x_size, x_size) symmetric;
        minimize(norm(est_para_a * vec(est_cov_q) - est_para_b)); %, Inf));
        subject to
        est_cov_q == semidefinite(x_size);
    cvx_end;
    % OLS估计
    % est_cov_q = lsqr(est_para_a, est_para_b, 1e-7, 100);
    % est_cov_q = reshape(est_cov_q, [x_size x_size]);

    % 估计S(未编写)
    est_cov_s = zeros(x_size, y_size);
    
    % 估计R(多个值求平均)
    % 参数计算
    est_para_r = zeros(y_size, y_size, l_size);
    for iter_l = 1:l_size
        temp_r = zeros(y_size, y_size);
        for iter_j = 0:iter_l-1
            temp_r = temp_r + mat_c*series_mat_a(:, :, offset+iter_j+1)*est_cov_q*series_mat_a(:, :, offset+iter_j+1).'*mat_c.';
        end
        est_para_r(:, :, iter_l) = temp_r;
    end
    % 求平均估计
    est_cov_r = est_mom_zzl - est_para_r;
    est_cov_r = sum(est_cov_r, 3)./l_size;

    % 返回值
    est_cov_xy = [est_cov_q est_cov_s; est_cov_s.' est_cov_r];

end

%% 模型验证
function diff_norm = testEstimation(yn, noise_cov, mat_a, mat_c)
% 输出-状态协方差 - 不存在互协方差
% 此处计算真实矩[Q S; S^T R], 而P_k为状态协方差(非稳态)
%#ok<*EQEFF>

    % 参数计算
    x_size = size(mat_a, 1); y_size = size(mat_c, 1);
    i_size = 10; k_offset = 12;
    % 参数提取
    noise_q = noise_cov(1:x_size, 1:x_size);
    noise_s = noise_cov(1:x_size, x_size+1:end);
    noise_r = noise_cov(x_size+1:end, x_size+1:end);

    % 样本方差估计量, 设符号表示为E(x_{k+i}x_{k}\tr),
    % 则至少要取x*(x+y)/y个量
    % 此处直接取2*x^2个 @ (y_size, y_size, i_size, k_size)
    % 同时, 对于R的估计, 直接取2*x^2+2个点(即前面所用的点).
    est_mom_zz0 = momentEstimation(yn, 0, k_offset, 0);
    est_mom_zzi = zeros(y_size, y_size, i_size);
    for iter_i = 1:i_size
        est_mom_zzi(:, :, iter_i) = momentEstimation(yn, iter_i, k_offset, 0);
    end

    % 根据已知数据直接计算方差
    series_transient = zeros(x_size, x_size);
    for iter_temp = 0:k_offset-1
        series_transient = series_transient + mpower(mat_a, iter_temp)*noise_q*(mpower(mat_a, iter_temp).');
    end
    ana_mom_zz0 = mat_c*series_transient*(mat_c.') + noise_r;
    ana_mom_zzi = zeros(y_size, y_size, i_size);
    for iter_i = 1:i_size
        ana_mom_zzi(:, :, iter_i) = mat_c*mpower(mat_a, iter_i-1)*(mat_a*series_transient*(mat_c.') + noise_s);
    end

    % 对比
    diff_norm = zeros(i_size+1, 1);
    [diff_norm(1), ~] = anaNorm(est_mom_zz0-ana_mom_zz0);
    for iter_i = 1:i_size
        [diff_norm(iter_i+1), ~] = anaNorm(est_mom_zzi(:, :, iter_i)-ana_mom_zzi(:, :, iter_i));
    end

end


%% 正定矩阵生成(0~1)
function ret_mat = semidefMatrixBuilder(seed, n)
    ret_mat = rand(seed, n, n);
    ret_mat = 0.5.*(ret_mat+ret_mat.');
    ret_mat = ret_mat + n.*eye(n);
end

%% 矩估计
function ret = momentEstimation(xn, i, k, offset)
% 计算E(xn_{k+i}xn_{k}\tr), offset为开头截断的长度
    exp_size = size(xn, 3);
    covdim = @(x) reshape(x, [size(x, 1) size(x, 3)]);
    ret = (1/exp_size).*(covdim(xn(:, offset+k+i, :))*(covdim(xn(:, offset+k, :)).'));
end

%% 解析稳态矩
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
