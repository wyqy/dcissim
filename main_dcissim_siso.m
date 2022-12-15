%#ok<*REPMAT>
%#ok<*EXIST>
%#ok<*UNRCH>
%% 确定参数
% 仿真系统
para_xsize = 4; para_ysize = 1; para_usize = 1;
% para_xsize = 4; para_ysize = 1; para_usize = 1;

% 时间参数
para_sim_duration = 600;  % 总时间
para_sim_step = 1/1e2;  % 时间步长
para_sim_period = 30;  % 输入信号周期时间

% 扰动参数
para_sim_x_snr = 20;  % 状态信号基础信噪比
para_sim_y_snr = 20;  % 输出信号基础信噪比
para_sim_u_snr = 20;  % 输入信号基础信噪比
para_sim_xyu_snr = [repmat(para_sim_x_snr, [para_xsize 1]); repmat(para_sim_y_snr, [para_ysize 1]); repmat(para_sim_u_snr, [para_usize 1])];
% 算法参数
para_iden_cutted_period = 1;  % 切除最开始未进入中心流形的部分
para_iden_x_size_bound = 6;   % X的变量数上界
para_iden_sim_ss_bdx_type = 'analytical';  % BDX辨识方法: 'analytical', 'optimize'
para_iden_sim_ss_d_type = 'null';  % 是否假定D矩阵为0
para_iden_cov_cross_type = 'null';  % 是否假定存在互协方差

% 随机模型个数
if ~exist('para_experiment_count'), para_experiment_count = 1; end
result_original_cell = cell(para_experiment_count, 1);
result_full_cell = cell(para_experiment_count, 1); time_full = zeros(para_experiment_count, 1);
result_reduce_cell = cell(para_experiment_count, 1); time_reduce = zeros(para_experiment_count, 1);
result_sim_cell = cell(para_experiment_count, 1); time_sim = zeros(para_experiment_count, 1);

% 随机数种子
% rng('shuffle', 'twister'); seed = uint32(randi(intmax('uint32'), 1)); % 非固定参数
seed = uint32(271828);  % 固定参数, 271828 for covariance estimation
rng(seed, 'simdTwister');

%% 对每个随机模型运行一次
for iter_exp = 1:para_experiment_count
    %% 数据生成
    % 激励信号生成
    sim_excitation_struct = genExcitationSignal('type', 'random', 'signal_size', para_usize, ...
        'duration', para_sim_duration, 'step', para_sim_step, 'period', para_sim_period);
    % 导出信息
    sim_excitation_signal = sim_excitation_struct.signal;
    sim_excitation_samples = sim_excitation_struct.samples;
    sim_period_samples = sim_excitation_struct.period_samples;  % 单周期仿真采样数
    % 仿真状态生成
    sim_result_struct = plantSimulation('type', 'generated', 'x_size', para_xsize, 'y_size', para_ysize, 'u_size', para_usize, ...
        'period_samples', sim_period_samples, 'excitation', sim_excitation_signal, 'snr', para_sim_xyu_snr);
    un = sim_result_struct.un; yn = sim_result_struct.yn;
    result_original_cell{iter_exp} = sim_result_struct.plant_info;
    
    %% 辨识系统参数 - dcISSIM (full)
    tic;
    % 初始化
    iden_full_init = idenDCISSIMLaucher(un, 'y_size', size(yn, 1), 'u_size', size(un, 1), 'period_samples', sim_period_samples, 'cutted_periods', para_iden_cutted_period, ...
        'dcissim_type', 'offline', 'isim_excitation_type', 'full', 'x_size_upbound', para_iden_x_size_bound, 'sim_ss_bdx_type', para_iden_sim_ss_bdx_type, 'sim_ss_d_type', para_iden_sim_ss_d_type, 'cov_cross_type', para_iden_cov_cross_type, ...
        'sim_x_size_type', 'fixed', 'cov_order_type', 'estimate', 'sim_x_size', para_xsize);
    % 辨识
    result_full_cell{iter_exp} = idenDCISSIMRunner(iden_full_init, yn, un);
    % 计时
    time_full(iter_exp) = toc;

    %% 辨识系统参数 - dcISSIM (reduce)
    tic;
    % 初始化
    iden_reduce_init = idenDCISSIMLaucher(un, 'y_size', size(yn, 1), 'u_size', size(un, 1), 'period_samples', sim_period_samples, 'cutted_periods', para_iden_cutted_period, ...
        'dcissim_type', 'offline', 'isim_excitation_type', 'reduced', 'x_size_upbound', para_iden_x_size_bound, 'sim_ss_bdx_type', para_iden_sim_ss_bdx_type, 'sim_ss_d_type', para_iden_sim_ss_d_type, 'cov_cross_type', para_iden_cov_cross_type, ...
        'sim_x_size_type', 'fixed', 'cov_order_type', 'estimate', 'sim_x_size', para_xsize);
    % 辨识
    result_reduce_cell{iter_exp} = idenDCISSIMRunner(iden_reduce_init, yn, un);
    % 计时
    time_reduce(iter_exp) = toc;

    %% 辨识系统参数 - SIM
    tic;
    % 初始化
    iden_sim_data = iddata(yn.', un.', para_sim_step);
    iden_sim_data.Period = repmat(sim_period_samples, [para_usize, 1]);
    iden_sim_option = n4sidOptions('InitialState', 'zero', 'Focus', 'simulation', 'EstimateCovariance', true, 'Display', 'off');
    % 辨识
    result_sim_cell{iter_exp} = n4sid(iden_sim_data, para_xsize, 'Feedthrough', 0, iden_sim_option);
    % 计时
    time_sim(iter_exp) = toc;

end

%% 验证确定性参数辨识结果
% 准备参数, 第一行为原系统, 第二行为辨识系统
analysis_xsize = zeros(4, para_experiment_count);
analysis_error_h2 = zeros(3, para_experiment_count);
analysis_error_hinf = zeros(3, para_experiment_count);
% 计算
for iter_exp = 1:para_experiment_count
    % 准备
    temp_struct = result_original_cell{iter_exp}; analysis_original_ss = ss(temp_struct.A, temp_struct.B, temp_struct.C, temp_struct.D, para_sim_step);
    temp_struct = result_full_cell{iter_exp}; analysis_full_idss = idss(temp_struct.A, temp_struct.B, temp_struct.C, temp_struct.D, temp_struct.K, zeros(temp_struct.x_size, 1), para_sim_step, 'NoiseVariance', temp_struct.cov(temp_struct.x_size+1:temp_struct.x_size+para_ysize, temp_struct.x_size+1:temp_struct.x_size+para_ysize));
    temp_struct = result_reduce_cell{iter_exp}; analysis_reduce_idss = idss(temp_struct.A, temp_struct.B, temp_struct.C, temp_struct.D, temp_struct.K, zeros(temp_struct.x_size, 1), para_sim_step, 'NoiseVariance', temp_struct.cov(temp_struct.x_size+1:temp_struct.x_size+para_ysize, temp_struct.x_size+1:temp_struct.x_size+para_ysize));
    analysis_sim_idss = result_sim_cell{iter_exp};
    % xsize
    analysis_xsize(1, iter_exp) = para_xsize;
    analysis_xsize(2, iter_exp) = size(analysis_full_idss.A, 1);
    analysis_xsize(3, iter_exp) = size(analysis_reduce_idss.A, 1);
    analysis_xsize(4, iter_exp) = size(analysis_sim_idss.A, 1);
    % 误差系统的H范数 (仅关注第一个输入到第一个输出)
    [temp_z, temp_p] = ss2tf(analysis_original_ss.A, analysis_original_ss.B, analysis_original_ss.C, analysis_original_ss.D, 1); temp_z = temp_z(1, :); temp_p = temp_p(1, :); temp_ss_ori = tf(temp_z, temp_p, para_sim_step);
    [temp_z, temp_p] = ss2tf(analysis_full_idss.A, analysis_full_idss.B, analysis_full_idss.C, analysis_full_idss.D, 1); temp_z = temp_z(1, :); temp_p = temp_p(1, :); temp_ss_full = tf(temp_z, temp_p, para_sim_step);
    [temp_z, temp_p] = ss2tf(analysis_reduce_idss.A, analysis_reduce_idss.B, analysis_reduce_idss.C, analysis_reduce_idss.D, 1); temp_z = temp_z(1, :); temp_p = temp_p(1, :); temp_ss_reduce = tf(temp_z, temp_p, para_sim_step);    
    [temp_z, temp_p] = ss2tf(analysis_sim_idss.A, analysis_sim_idss.B, analysis_sim_idss.C, analysis_sim_idss.D, 1); temp_z = temp_z(1, :); temp_p = temp_p(1, :); temp_ss_sim = tf(temp_z, temp_p, para_sim_step);
    [analysis_error_h2(1, iter_exp), analysis_error_hinf(1, iter_exp)] = sys_compare(temp_ss_ori, temp_ss_full);
    [analysis_error_h2(2, iter_exp), analysis_error_hinf(2, iter_exp)] = sys_compare(temp_ss_ori, temp_ss_reduce);
    [analysis_error_h2(3, iter_exp), analysis_error_hinf(3, iter_exp)] = sys_compare(temp_ss_ori, temp_ss_sim);
end
% 求平均误差
analysis_criterion_h2 = mean(analysis_error_h2, 2); analysis_criterion_hinf = mean(analysis_error_hinf, 2);
disp(['H2 norm differences of determined paras: ', mat2str(analysis_criterion_h2(1), 2), ' / ', mat2str(analysis_criterion_h2(2), 2), ' / ', mat2str(analysis_criterion_h2(3), 2)]);
disp(['Hinf norm differences of determined paras: ', mat2str(analysis_criterion_hinf(1), 2), ' / ', mat2str(analysis_criterion_hinf(2), 2), ' / ', mat2str(analysis_criterion_hinf(3), 2)]);

%% 验证随机参数辨识结果
% 准备参数
analysis_cov_count = 10;  % 验证阶数
analysis_incov_error_norm = zeros(1, para_experiment_count);
analysis_outcov_error_norm = zeros(3, para_experiment_count);
analysis_outcov_series = zeros(4, analysis_cov_count, para_experiment_count);

% 计算方差
% 输入方差
for iter_exp = 1:para_experiment_count
    temp_ori_cov = result_original_cell{iter_exp}.cov(para_xsize+para_ysize+1:end, para_xsize+para_ysize+1:end);
    temp_full_cov = result_full_cell{iter_exp}.cov(analysis_xsize(2, iter_exp)+para_ysize+1:end, analysis_xsize(2, iter_exp)+para_ysize+1:end);
    [analysis_incov_error_norm(iter_exp), ~] = sys_compare(temp_ori_cov, temp_full_cov);
end
% 输出方差
for iter_exp = 1:para_experiment_count
    % 可视化 + 数据计算
    temp_ori_cov = result_original_cell{iter_exp}.cov(1:para_xsize+para_ysize, 1:para_xsize+para_ysize);
    temp_full_cov = result_full_cell{iter_exp}.cov(1:analysis_xsize(2, iter_exp)+para_ysize, 1:analysis_xsize(2, iter_exp)+para_ysize);
    temp_reduce_cov = result_reduce_cell{iter_exp}.cov(1:analysis_xsize(3, iter_exp)+para_ysize, 1:analysis_xsize(3, iter_exp)+para_ysize);
    temp_sim_cov = anaInnoCov(result_sim_cell{iter_exp}.K, result_sim_cell{iter_exp}.NoiseVariance);
    [~, ~, temp_ori_moments, analysis_outcov_series(1, :, iter_exp)] = anaCov(temp_ori_cov, result_original_cell{iter_exp}, analysis_cov_count);
    [~, ~, temp_full_moments, analysis_outcov_series(2, :, iter_exp)] = anaCov(temp_full_cov, result_full_cell{iter_exp}, analysis_cov_count);
    [~, ~, temp_reduce_moments, analysis_outcov_series(3, :, iter_exp)] = anaCov(temp_reduce_cov, result_reduce_cell{iter_exp}, analysis_cov_count);
    [~, ~, temp_sim_moments, analysis_outcov_series(4, :, iter_exp)] = anaCov(temp_sim_cov, result_sim_cell{iter_exp}, analysis_cov_count);
    % 用于数据计算
    [temp_outerror_1, ~] = sys_compare(temp_ori_moments{1}, temp_full_moments{1}); [temp_outerror_2, ~] = sys_compare(temp_ori_moments{2}, temp_full_moments{2}); analysis_outcov_error_norm(1, iter_exp) = temp_outerror_1 + temp_outerror_2;
    [temp_outerror_1, ~] = sys_compare(temp_ori_moments{1}, temp_reduce_moments{1}); [temp_outerror_2, ~] = sys_compare(temp_ori_moments{2}, temp_reduce_moments{2}); analysis_outcov_error_norm(2, iter_exp) = temp_outerror_1 + temp_outerror_2;
    [temp_outerror_1, ~] = sys_compare(temp_ori_moments{1}, temp_sim_moments{1}); [temp_outerror_2, ~] = sys_compare(temp_ori_moments{2}, temp_sim_moments{2}); analysis_outcov_error_norm(3, iter_exp) = temp_outerror_1 + temp_outerror_2;
end

% 求平均误差
analysis_criterion_incov = mean(analysis_incov_error_norm);
analysis_criterion_outcov = mean(analysis_outcov_error_norm, 2);
disp(['Differences of input covariances as norm: ', mat2str(analysis_criterion_incov, 2)]);
disp(['Differences of output covariances as norm: ', mat2str(analysis_criterion_outcov(1), 2), ' / ', mat2str(analysis_criterion_outcov(2), 2), ' / ', mat2str(analysis_criterion_outcov(3), 2)]);

%% H2差异最大, 中位数和最小的Bode图
% [~, analysis_hinf_sort] = sort(analysis_error_h2(1, :), 'descend');
% analysis_bode_location = analysis_hinf_sort(1);
% analysis_bode_location = analysis_hinf_sort(ceil(para_experiment_count/2));
analysis_bode_location = 93;
analysis_bode_cell = {result_original_cell{analysis_bode_location}, result_full_cell{analysis_bode_location}, result_reduce_cell{analysis_bode_location}, result_sim_cell{analysis_bode_location}};
% fig = anaPlotBode(analysis_bode_cell, 'sample', para_sim_step); sgtitle(fig, 'Bode plot for minimum d\_{Hinf}');

%% 输出方差差异最小的矩估计图
% [~, analysis_cov_sort] = sort(analysis_criterion_outcov(1, :), 'descend');
analysis_cov_location = 93;
% ax = anaPlotCov(analysis_outcov_series(:, :, analysis_cov_location)); title(ax, 'Coavariance plot for minimum d\_{Hinf}');

%% 输出时间
disp(['Time used: ', mat2str(mean(time_full), 2), ' / ', mat2str(mean(time_reduce), 2), ' / ', mat2str(mean(time_sim), 2)]);

%% 保存参数
if ~exist('save_valid'), save_valid = false; end
save_postfix = '_siso';
if save_valid
    save(['results\mat_auto' save_postfix '.mat'], '-v7.3');
end

%% 辅助命令
% mean(analysis_error_h2, 2)
% std(analysis_error_h2, 0, 2)
% %%
% mean(analysis_error_hinf, 2)
% std(analysis_error_hinf, 0, 2)
% %%
% mean(analysis_incov_error_norm)
% std(analysis_incov_error_norm, 0, 2)
% %%
% mean(analysis_outcov_error_norm, 2)
% std(analysis_outcov_error_norm, 0, 2)
% %%
% mean(time_full)
% mean(time_reduce)
% mean(time_sim)
% %%
% std(time_reduce)
% std(time_full)
% std(time_sim)

%% 辅助函数
function [compare_h2, compare_hinf] = sys_compare(sys_ori, sys_ind)
    [num_h2, num_hinf] = anaNorm(sys_ind - sys_ori);
    [den_h2, den_hinf] = anaNorm(sys_ori);
    compare_h2 = abs(num_h2)/den_h2;
    compare_hinf = abs(num_hinf)/den_hinf;
end

