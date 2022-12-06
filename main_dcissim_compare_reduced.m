%#ok<*REPMAT>
%% 确定参数
% 仿真系统
para_xsize = 4; para_ysize = 1; para_usize = 1;

% 是否全量辨识
para_isfull = false; %#ok<*UNRCH>
if para_isfull
    para_iden_isim_excitation_type = 'full';
    % 时间参数
    para_sim_duration = 500;  % 总时间
    para_sim_step = 1/1e2;  % 时间步长
    para_sim_period = 20;  % 输入信号周期时间
else
    para_iden_isim_excitation_type = 'reduced';
    % 时间参数
    para_sim_duration = 1000;  % 总时间
    para_sim_step = 1/1e2;  % 时间步长
    para_sim_period = 40;  % 输入信号周期时间
end

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
para_experiment_count = 1;
result_original_cell = cell(para_experiment_count, 1);
result_dcissim_cell = cell(para_experiment_count, 1); time_dcissim = zeros(para_experiment_count, 1);
result_sim_cell = cell(para_experiment_count, 1); time_sim = zeros(para_experiment_count, 1);

% 随机数种子
% rng('shuffle', 'twister'); seed = randi(intmax('uint32'), 'uint32'); % 非固定参数
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
    
    %% 辨识系统参数 - dcISSIM
    tic;
    % 初始化
    iden_dcissim_init = idenDCISSIMLaucher(un, 'y_size', size(yn, 1), 'u_size', size(un, 1), 'period_samples', sim_period_samples, 'cutted_periods', para_iden_cutted_period, ...
        'dcissim_type', 'offline', 'isim_excitation_type', para_iden_isim_excitation_type, 'x_size_upbound', para_iden_x_size_bound, 'sim_ss_bdx_type', para_iden_sim_ss_bdx_type, 'sim_ss_d_type', para_iden_sim_ss_d_type, 'cov_cross_type', para_iden_cov_cross_type, ...
        'sim_x_size_type', 'fixed', 'cov_order_type', 'estimate', 'sim_x_size', para_xsize);
    % 辨识
    result_dcissim_cell{iter_exp} = idenDCISSIMRunner(iden_dcissim_init, yn, un);
    % 计时
    time_dcissim(iter_exp) = toc;

    %% 辨识系统参数 - SIM
    tic;
    % 初始化
    iden_sim_data = iddata(yn.', un.', 1);
    iden_sim_data.Period = repmat(sim_period_samples, [para_usize, 1]);
    iden_sim_option = n4sidOptions('InitialState', 'zero', 'Focus', 'simulation', 'EstimateCovariance', true, 'Display', 'off');
    % 辨识
    result_sim_cell{iter_exp} = n4sid(iden_sim_data, para_xsize, 'Feedthrough', 0, iden_sim_option);
    % 计时
    time_sim(iter_exp) = toc;

end

%% 验证确定性参数辨识结果
% 准备参数, 第一行为原系统, 第二行为辨识系统
analysis_xsize = zeros(3, para_experiment_count);
analysis_eig = zeros(3, para_xsize, para_experiment_count);
analysis_error_h2 = zeros(2, para_experiment_count);
analysis_error_hinf = zeros(2, para_experiment_count);
% 计算
for iter_exp = 1:para_experiment_count
    % 准备
    temp_struct = result_original_cell{iter_exp};
    analysis_original_ss = ss(temp_struct.A, temp_struct.B, temp_struct.C, temp_struct.D, para_sim_step);
    temp_struct = result_dcissim_cell{iter_exp};
    analysis_dcissim_idss = idss(temp_struct.A, temp_struct.B, temp_struct.C, temp_struct.D, temp_struct.K, ...
        zeros(temp_struct.x_size, 1), para_sim_step, 'NoiseVariance', temp_struct.cov(temp_struct.x_size+1:temp_struct.x_size+para_ysize, temp_struct.x_size+1:temp_struct.x_size+para_ysize));
    analysis_sim_idss = result_sim_cell{iter_exp};
    % xsize
    analysis_xsize(1, iter_exp) = para_xsize;
    analysis_xsize(2, iter_exp) = result_dcissim_cell{iter_exp}.x_size;
    analysis_xsize(3, iter_exp) = size(analysis_sim_idss.A, 1);
    % 特征值
    analysis_eig(1, :, iter_exp) = sort(eig(analysis_original_ss.A), 'descend', 'ComparisonMethod', 'abs');
    temp_eig = sort(eig(analysis_dcissim_idss.A), 'descend', 'ComparisonMethod', 'abs');
    if length(temp_eig) > para_xsize, temp_eig = temp_eig(1:para_xsize);
    else, temp_eig = [temp_eig(1:length(temp_eig)); zeros(para_xsize-length(temp_eig), 1)]; end
    analysis_eig(2, :, iter_exp) = temp_eig;
    temp_eig = sort(eig(analysis_sim_idss.A), 'descend', 'ComparisonMethod', 'abs');
    if length(temp_eig) > para_xsize, temp_eig = temp_eig(1:para_xsize);
    else, temp_eig = [temp_eig(1:length(temp_eig)); zeros(para_xsize-length(temp_eig), 1)]; end
    analysis_eig(3, :, iter_exp) = temp_eig;
    % 误差系统的H范数 (仅关注第一个输入到第一个输出)
    [temp_z, temp_p] = ss2tf(analysis_original_ss.A, analysis_original_ss.B, analysis_original_ss.C, analysis_original_ss.D, 1); temp_z = temp_z(1, :); temp_p = temp_p(1, :); temp_ss_ori = tf(temp_z, temp_p, para_sim_step);
    [temp_z, temp_p] = ss2tf(analysis_dcissim_idss.A, analysis_dcissim_idss.B, analysis_dcissim_idss.C, analysis_dcissim_idss.D, 1); temp_z = temp_z(1, :); temp_p = temp_p(1, :); temp_ss_dcissim = tf(temp_z, temp_p, para_sim_step);
    [temp_z, temp_p] = ss2tf(analysis_sim_idss.A, analysis_sim_idss.B, analysis_sim_idss.C, analysis_sim_idss.D, 1); temp_z = temp_z(1, :); temp_p = temp_p(1, :); temp_ss_sim = tf(temp_z, temp_p, para_sim_step);
    [analysis_error_h2(1, iter_exp), analysis_error_hinf(1, iter_exp)] = sys_compare(temp_ss_ori, temp_ss_dcissim);
    [analysis_error_h2(2, iter_exp), analysis_error_hinf(2, iter_exp)] = sys_compare(temp_ss_ori, temp_ss_sim);
end
% 求平均误差
analysis_criterion_h2 = mean(analysis_error_h2, 2); analysis_criterion_hinf = mean(analysis_error_hinf, 2);
disp(['R2 of H2 norm differences of determined paras: ', mat2str(analysis_criterion_h2(1), 2), ' / ', mat2str(analysis_criterion_h2(2), 2)]);
disp(['R2 of Hinf norm differences of determined paras: ', mat2str(analysis_criterion_hinf(1), 2), ' / ', mat2str(analysis_criterion_hinf(2), 2)]);

%% 验证随机参数辨识结果
% 准备参数
analysis_cov_count = 10;  % 验证阶数
analysis_incov_error_norm = zeros(1, para_experiment_count); analysis_outcov_error_norm = zeros(analysis_cov_count, para_experiment_count);
analysis_outcov_norm = zeros(2, analysis_cov_count, para_experiment_count);
% 计算方差
% 输入方差
if para_isfull
    for iter_exp = 1:para_experiment_count
        temp_ori_cov = result_original_cell{iter_exp}.cov(para_xsize+para_ysize+1:end, para_xsize+para_ysize+1:end);
        temp_ind_cov = result_dcissim_cell{iter_exp}.cov(result_dcissim_cell{iter_exp}.x_size+para_ysize+1:end, result_dcissim_cell{iter_exp}.x_size+para_ysize+1:end);
        [analysis_incov_error_norm(iter_exp), ~] = sys_compare(temp_ori_cov, temp_ind_cov);
    end
end
% 输出方差
for iter_exp = 1:para_experiment_count
    temp_ori_cov = result_original_cell{iter_exp}.cov(1:para_xsize+para_ysize, 1:para_xsize+para_ysize);
    temp_ind_cov = result_dcissim_cell{iter_exp}.cov(1:result_dcissim_cell{iter_exp}.x_size+para_ysize, 1:result_dcissim_cell{iter_exp}.x_size+para_ysize);
    [temp_ori_cell, analysis_outcov_norm(1, :, iter_exp)] = anaCov(temp_ori_cov, result_original_cell{iter_exp}, analysis_cov_count);
    [temp_ind_cell, analysis_outcov_norm(2, :, iter_exp)] = anaCov(temp_ind_cov, result_dcissim_cell{iter_exp}, analysis_cov_count);
    for iter_cov = 1:analysis_cov_count
        analysis_outcov_error_norm(iter_cov, iter_exp) = sys_compare(temp_ori_cell{iter_cov}, temp_ind_cell{iter_cov});
    end
end
% 求平均误差
if para_isfull, analysis_criterion_incov = mean(analysis_incov_error_norm); end
analysis_criterion_outcov = mean(analysis_outcov_error_norm, 'all');
if para_isfull, disp(['R2 of dCov(input) as norm: ', mat2str(analysis_criterion_incov, 2)]); end
disp(['R2 of dCov(output) as norm: ', mat2str(analysis_criterion_outcov, 2)]);

%% H2差异最大, 中位数和最小的Bode图
[~, analysis_hinf_sort] = sort(analysis_error_h2(1, :), 'descend');

% analysis_bode_location = analysis_hinf_sort(1); fig = anaPlotBode(result_original_cell{analysis_bode_location}, result_dcissim_cell{analysis_bode_location}, 'sample', para_sim_step); sgtitle(fig, 'Bode plot for maximum d\_{Hinf}');
% analysis_bode_location = analysis_hinf_sort(ceil(para_experiment_count/2)); fig = anaPlotBode(result_original_cell{analysis_bode_location}, result_dcissim_cell{analysis_bode_location}, 'sample', para_sim_step); sgtitle(fig, 'Bode plot for medium d\_{Hinf}');
% analysis_bode_location = analysis_hinf_sort(end); fig = anaPlotBode(result_original_cell{analysis_bode_location}, result_dcissim_cell{analysis_bode_location}, 'sample', para_sim_step); sgtitle(fig, 'Bode plot for minimum d\_{Hinf}');

%% 输出方差差异最小的矩估计图
[~, analysis_cov_sort] = sort(mean(analysis_outcov_error_norm, 1), 'descend');

analysis_cov_location = analysis_cov_sort(end);
temp_cov_ori = shiftdim(analysis_outcov_norm(1, :, analysis_cov_location));
temp_cov_ind = shiftdim(analysis_outcov_norm(2, :, analysis_cov_location));
% ax = anaPlotCov(temp_cov_ori, temp_cov_ind); title(ax, 'Coavariance plot for minimum d\_{Hinf}');

%% 输出时间
disp(['Time used: ', mat2str(mean(time_dcissim), 2), ' / ', mat2str(mean(time_sim), 2)]);

%% 保存参数
save('results\reduced_many_auto.mat', '-v7.3');

%% 残差估计
% res_t = 0:para_sim_step:(sim_excitation_samples-1)*para_sim_step;
% res_model = result_identified_cell{end};
% % ISIM估计
% res_isim_vn = idenRegressor(sim_period_samples, res_model.freqs, sim_excitation_samples, 'ordinary').vn;
% res_isim_xn = res_model.X*res_isim_vn; figure; plot(res_t, sim_result_struct.xn(1, :), res_t, res_isim_xn(1, :));
% res_isim_yn = res_model.Y*res_isim_vn; figure; plot(res_t, sim_result_struct.yn(1, :), res_t, res_isim_yn(1, :));
% % cISSIM估计
% res_ss_xinit = zeros(res_model.x_size, 1);
% res_ss_noise = zeros(res_model.x_size+para_ysize+para_usize, sim_excitation_samples);
% [res_ss_xn, res_ss_yn, res_ss_un] = plantModel(res_model, res_ss_xinit, un, res_ss_noise);
% res_ss_dyn = yn - res_ss_yn; % figure; plot(res_ss_dyn(1, :));
% % 输出噪声估计
% res_isim_vn = idenRegressor(sim_period_samples, res_model.freqs, sim_excitation_samples, 'ordinary').vn;
% res_dun = un - res_model.U*res_isim_vn;
% % 对比功率谱密度 (取一个通道)
% figure;
% subplot(2, 2, 1); periodogram(res_ss_dyn(1, :), rectwin(sim_excitation_samples), sim_period_samples);
% subplot(2, 2, 2); periodogram(sim_result_struct.noise(para_xsize+1, :), rectwin(sim_excitation_samples), sim_period_samples);
% subplot(2, 2, 3); periodogram(res_dun(1, :), rectwin(sim_excitation_samples), sim_period_samples);
% subplot(2, 2, 4); periodogram(sim_result_struct.noise(para_xsize+para_ysize+1, :), rectwin(sim_excitation_samples), sim_period_samples);

%% 辅助函数
function [compare_h2, compare_hinf] = sys_compare(sys_ori, sys_ind)
    [num_h2, num_hinf] = anaNorm(sys_ind - sys_ori);
    [den_h2, den_hinf] = anaNorm(sys_ori);
    compare_h2 = abs(num_h2)/den_h2;
    compare_hinf = abs(num_hinf)/den_hinf;
end

