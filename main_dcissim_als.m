%#ok<*REPMAT>
%#ok<*EXIST>
%#ok<*UNRCH>
%% 确定参数
% 仿真系统
para_xsize = 2; para_ysize = para_xsize; para_usize = 1;
para_xyusize = para_xsize + para_ysize + para_usize;
% para_xsize = 4; para_ysize = 1; para_usize = 1;
% 时间参数
para_sim_duration = 600;  % 总时间
para_sim_step = 1/1e2;  % 时间步长
para_sim_period = 20;  % 输入信号周期时间
% 扰动参数
para_sim_x_snr = 20;  % 状态信号基础信噪比
para_sim_y_snr = 20;  % 输出信号基础信噪比
para_sim_u_snr = 20;  % 输入信号基础信噪比
para_sim_xyu_snr = [repmat(para_sim_x_snr, [para_xsize 1]); repmat(para_sim_y_snr, [para_ysize 1]); repmat(para_sim_u_snr, [para_usize 1])];
% 算法参数
para_iden_cutted_period = 1;  % 切除最开始未进入中心流形的部分
para_iden_x_size_bound = 4;   % X的变量数上界
para_iden_sim_ss_bdx_type = 'analytical';  % BDX辨识方法: 'analytical', 'optimize'
para_iden_sim_ss_d_type = 'null';  % 是否假定D矩阵为0
para_iden_cov_cross_type = 'null';  % 是否假定存在互协方差

% 随机模型个数
if ~exist('para_experiment_count'), para_experiment_count = 1; end
result_simple_cell = cell(para_experiment_count, 1); time_full = zeros(para_experiment_count, 1);
result_classic_cell = cell(para_experiment_count, 1); time_reduce = zeros(para_experiment_count, 1);

% 随机数种子
% rng('shuffle', 'twister'); seed = uint32(randi(intmax('uint32'), 1)); % 非固定参数
seed = uint32(3231441189);  % 固定参数
rng(seed, 'simdTwister');

%% 生成固定系统
% 激励信号生成
sim_excitation_struct = genExcitationSignal('type', 'random', 'signal_size', para_usize, ...
    'duration', para_sim_duration, 'step', para_sim_step, 'period', para_sim_period);
% 导出信息
sim_excitation_signal = sim_excitation_struct.signal;
sim_excitation_samples = sim_excitation_struct.samples;
sim_period_samples = sim_excitation_struct.period_samples;  % 单周期仿真采样数
% 仿真系统生成
sim_result_struct = plantSimulation('type', 'fixed', 'x_size', para_xsize, 'y_size', para_ysize, 'u_size', para_usize, ...
    'period_samples', sim_period_samples, 'excitation', sim_excitation_signal, 'snr', para_sim_xyu_snr);
result_original = sim_result_struct.plant_info;
% 噪声种子生成
sim_noise_rs = RandStream.create('mrg32k3a', 'NumStreams', para_xyusize, 'Seed', seed, 'CellOutput', true);

% 初始化算法
iden_simple_init = idenDCISSIMLaucher(sim_result_struct.un, 'y_size', para_ysize, 'u_size', para_usize, 'period_samples', sim_period_samples, 'cutted_periods', para_iden_cutted_period, ...
    'dcissim_type', 'offline', 'isim_excitation_type', 'full', 'x_size_upbound', para_iden_x_size_bound, 'sim_ss_bdx_type', para_iden_sim_ss_bdx_type, 'sim_ss_d_type', para_iden_sim_ss_d_type, 'cov_cross_type', para_iden_cov_cross_type, ...
    'sim_x_size_type', 'fixed', 'cov_order_type', 'estimate', 'sim_x_size', para_xsize, 'cov_est_type', 'simple');
iden_classic_init = idenDCISSIMLaucher(sim_result_struct.un, 'y_size', para_ysize, 'u_size', para_usize, 'period_samples', sim_period_samples, 'cutted_periods', para_iden_cutted_period, ...
    'dcissim_type', 'offline', 'isim_excitation_type', 'full', 'x_size_upbound', para_iden_x_size_bound, 'sim_ss_bdx_type', para_iden_sim_ss_bdx_type, 'sim_ss_d_type', para_iden_sim_ss_d_type, 'cov_cross_type', para_iden_cov_cross_type, ...
    'sim_x_size_type', 'fixed', 'cov_order_type', 'estimate', 'sim_x_size', para_xsize, 'cov_est_type', 'classical');

%% 对随机模型运行多次
for iter_exp = 1:para_experiment_count
    %% 重新生成扰动数据
    sim_result_struct = plantSimulationNoise(result_original, sim_excitation_signal, sim_noise_rs);
    un = sim_result_struct.un; yn = sim_result_struct.yn;
    
    %% 辨识方差参数 - dcISSIM (full)
    tic;
    result_simple_cell{iter_exp} = idenDCISSIMRunner(iden_simple_init, yn, un);
    % 计时
    time_full(iter_exp) = toc;
    
    %% 辨识方差参数 - dcISSIM (reduced)
    tic;
    result_classic_cell{iter_exp} = idenDCISSIMRunner(iden_classic_init, yn, un);
    % 计时
    time_reduce(iter_exp) = toc;

end

%% 变换辨识结果
for iter_exp = 1:para_experiment_count
    result_simple_cell{iter_exp} = anaSimilarTrans(result_original, result_simple_cell{iter_exp});
    result_classic_cell{iter_exp} = anaSimilarTrans(result_original, result_classic_cell{iter_exp});
end

%% 验证确定性参数辨识结果
% 准备参数, 第一行为原系统, 第二行为辨识系统
analysis_xsize = zeros(3, para_experiment_count);
analysis_error_h2 = zeros(2, para_experiment_count);
analysis_error_hinf = zeros(2, para_experiment_count);
% 计算
for iter_exp = 1:para_experiment_count
    % 准备
    temp_struct = result_original;
    analysis_original_ss = ss(temp_struct.A, temp_struct.B, temp_struct.C, temp_struct.D, para_sim_step);
    temp_struct = result_simple_cell{iter_exp};
    analysis_simple_idss = idss(temp_struct.A, temp_struct.B, temp_struct.C, temp_struct.D, temp_struct.K, ...
        zeros(temp_struct.x_size, 1), para_sim_step);%, 'NoiseVariance', temp_struct.cov(temp_struct.x_size+1:temp_struct.x_size+para_ysize, temp_struct.x_size+1:temp_struct.x_size+para_ysize));
    temp_struct = result_classic_cell{iter_exp};
    analysis_classic_idss = idss(temp_struct.A, temp_struct.B, temp_struct.C, temp_struct.D, temp_struct.K, ...
        zeros(temp_struct.x_size, 1), para_sim_step);%, 'NoiseVariance', temp_struct.cov(temp_struct.x_size+1:temp_struct.x_size+para_ysize, temp_struct.x_size+1:temp_struct.x_size+para_ysize));
    % xsize
    analysis_xsize(1, iter_exp) = para_xsize;
    analysis_xsize(2, iter_exp) = result_simple_cell{iter_exp}.x_size;
    analysis_xsize(3, iter_exp) = result_classic_cell{iter_exp}.x_size;
    % 误差系统的H范数 (仅关注第一个输入到第一个输出)
    [temp_z, temp_p] = ss2tf(analysis_original_ss.A, analysis_original_ss.B, analysis_original_ss.C, analysis_original_ss.D, 1); temp_z = temp_z(1, :); temp_p = temp_p(1, :); temp_ss_ori = tf(temp_z, temp_p, para_sim_step);
    [temp_z, temp_p] = ss2tf(analysis_simple_idss.A, analysis_simple_idss.B, analysis_simple_idss.C, analysis_simple_idss.D, 1); temp_z = temp_z(1, :); temp_p = temp_p(1, :); temp_ss_full = tf(temp_z, temp_p, para_sim_step);
    [temp_z, temp_p] = ss2tf(analysis_classic_idss.A, analysis_classic_idss.B, analysis_classic_idss.C, analysis_classic_idss.D, 1); temp_z = temp_z(1, :); temp_p = temp_p(1, :); temp_ss_reduce = tf(temp_z, temp_p, para_sim_step);
    [analysis_error_h2(1, iter_exp), analysis_error_hinf(1, iter_exp)] = sys_compare(temp_ss_ori, temp_ss_full);
    [analysis_error_h2(2, iter_exp), analysis_error_hinf(2, iter_exp)] = sys_compare(temp_ss_ori, temp_ss_reduce);
end
% 求平均误差
analysis_criterion_h2 = mean(analysis_error_h2, 2); analysis_criterion_hinf = mean(analysis_error_hinf, 2);
disp(['R2 of H2 norm differences of determined paras: ', mat2str(analysis_criterion_h2(1), 2), ' / ', mat2str(analysis_criterion_h2(2), 2)]);
disp(['R2 of Hinf norm differences of determined paras: ', mat2str(analysis_criterion_hinf(1), 2), ' / ', mat2str(analysis_criterion_hinf(2), 2)]);

%% 验证随机参数辨识结果
% 准备参数
analysis_cov_count = 10;  % 验证阶数
analysis_incov_error_norm = zeros(2, 1, para_experiment_count);
analysis_outcov_error_norm = zeros(2, 2, para_experiment_count);
analysis_inoutcov_error = zeros(2, para_xsize*para_xsize+para_ysize*para_ysize+para_usize*para_usize, para_experiment_count);
analysis_outcov_series = zeros(3, analysis_cov_count, para_experiment_count);

% 计算方差
for iter_exp = 1:para_experiment_count  %(默认para_ysize >= para_xsize)
    % 输入数据计算
    temp_ori_cov_tt = result_original.cov(para_xsize+para_ysize+1:end, para_xsize+para_ysize+1:end);
    temp_simple_cov_tt = result_simple_cell{iter_exp}.cov(analysis_xsize(2, iter_exp)+para_ysize+1:end, analysis_xsize(2, iter_exp)+para_ysize+1:end);
    temp_classic_cov_tt = result_classic_cell{iter_exp}.cov(analysis_xsize(2, iter_exp)+para_ysize+1:end, analysis_xsize(2, iter_exp)+para_ysize+1:end);
    analysis_inoutcov_error(1, para_xsize*para_xsize+para_ysize*para_ysize+1:end, iter_exp) = reshape(temp_simple_cov_tt-temp_ori_cov_tt, [], 1);
    analysis_inoutcov_error(2, para_xsize*para_xsize+para_ysize*para_ysize+1:end, iter_exp) = reshape(temp_classic_cov_tt-temp_ori_cov_tt, [], 1);
    [analysis_incov_error_norm(1, 1, iter_exp), ~] = sys_compare(temp_ori_cov_tt, temp_simple_cov_tt);
    [analysis_incov_error_norm(2, 1, iter_exp), ~] = sys_compare(temp_ori_cov_tt, temp_classic_cov_tt);
    % 输出数据计算
    temp_ori_cov_ww = result_original.cov(1:para_xsize, 1:para_xsize); temp_ori_cov_vv = result_original.cov(para_xsize+1:para_xsize+para_ysize, para_xsize+1:para_xsize+para_ysize);
    temp_simple_cov_ww = result_simple_cell{iter_exp}.cov(1:analysis_xsize(2, iter_exp), 1:analysis_xsize(2, iter_exp)); temp_simple_cov_vv = result_simple_cell{iter_exp}.cov(analysis_xsize(2, iter_exp)+1:analysis_xsize(2, iter_exp)+para_ysize, analysis_xsize(2, iter_exp)+1:analysis_xsize(2, iter_exp)+para_ysize);
    temp_classic_cov_ww = result_classic_cell{iter_exp}.cov(1:analysis_xsize(3, iter_exp), 1:analysis_xsize(3, iter_exp)); temp_classic_cov_vv = result_classic_cell{iter_exp}.cov(analysis_xsize(3, iter_exp)+1:analysis_xsize(3, iter_exp)+para_ysize, analysis_xsize(3, iter_exp)+1:analysis_xsize(3, iter_exp)+para_ysize);
    analysis_inoutcov_error(1, 1:para_xsize*para_xsize, iter_exp) = reshape(temp_simple_cov_ww-temp_ori_cov_ww, [], 1);
    analysis_inoutcov_error(2, 1:para_xsize*para_xsize, iter_exp) = reshape(temp_classic_cov_ww-temp_ori_cov_ww, [], 1);
    analysis_inoutcov_error(1, para_xsize*para_xsize+1:para_xsize*para_xsize+para_ysize*para_ysize, iter_exp) = reshape(temp_simple_cov_vv-temp_ori_cov_vv, [], 1);
    analysis_inoutcov_error(2, para_xsize*para_xsize+1:para_xsize*para_xsize+para_ysize*para_ysize, iter_exp) = reshape(temp_classic_cov_vv-temp_ori_cov_vv, [], 1);
    
    [~, ~, ~, analysis_outcov_series(1, :, iter_exp)] = anaCov(blkdiag(temp_ori_cov_ww, temp_ori_cov_vv), result_original, analysis_cov_count);
    [~, ~, ~, analysis_outcov_series(2, :, iter_exp)] = anaCov(blkdiag(temp_simple_cov_ww, temp_simple_cov_vv), result_simple_cell{iter_exp}, analysis_cov_count);
    [~, ~, ~, analysis_outcov_series(3, :, iter_exp)] = anaCov(blkdiag(temp_classic_cov_ww, temp_classic_cov_vv), result_classic_cell{iter_exp}, analysis_cov_count);
    
    [temp_simple_cov_ww, ~] = sys_compare(temp_ori_cov_ww, temp_simple_cov_ww); [temp_simple_cov_vv, ~] = sys_compare(temp_ori_cov_vv, temp_simple_cov_vv);
    [temp_classic_cov_ww, ~] = sys_compare(temp_ori_cov_ww, temp_classic_cov_ww); [temp_classic_cov_vv, ~] = sys_compare(temp_ori_cov_vv, temp_classic_cov_vv);
    analysis_outcov_error_norm(1, 1, iter_exp) = temp_simple_cov_ww ; analysis_outcov_error_norm(1, 2, iter_exp) = temp_simple_cov_vv ;
    analysis_outcov_error_norm(2, 1, iter_exp) = temp_classic_cov_ww; analysis_outcov_error_norm(2, 2, iter_exp) = temp_classic_cov_vv;
    
end

% 求平均误差
analysis_criterion_incov = mean(analysis_incov_error_norm, 3);
analysis_criterion_outcov = mean(analysis_outcov_error_norm, 3);
disp(['R2 of dww(output) as norm: ', mat2str(analysis_criterion_outcov(1, 1), 2), ' / ', mat2str(analysis_criterion_outcov(2, 1), 2)]);
disp(['R2 of dvv(output) as norm: ', mat2str(analysis_criterion_outcov(1, 2), 2), ' / ', mat2str(analysis_criterion_outcov(2, 2), 2)]);
disp(['R2 of dtt(input) as norm: ', mat2str(analysis_criterion_incov(1), 2)]);

%% H2差异最大, 中位数和最小的Bode图
[~, analysis_hinf_sort] = sort(analysis_error_h2(1, :), 'descend');
% analysis_bode_location = analysis_hinf_sort(1);
% analysis_bode_location = analysis_hinf_sort(ceil(para_experiment_count/2));
analysis_bode_location = analysis_hinf_sort(end);
analysis_bode_cell = {result_original, result_simple_cell{analysis_bode_location}, result_classic_cell{analysis_bode_location}};
% fig = anaPlotBode(analysis_bode_cell, 'sample', para_sim_step); % sgtitle(fig, 'Bode plot for minimum d\_{Hinf}');

%% 输出矩估计分布图
% [~, analysis_cov_sort] = sort(analysis_criterion_outcov(1, :), 'descend');
analysis_cov_location = 1;
% ax = anaPlotCov(analysis_outcov_series(:, :, analysis_cov_location)); % title(ax, 'Coavariance plot for minimum d\_{Hinf}');

%% 输出方差分布图
analysis_spot_inout_3d = cat(2, analysis_outcov_error_norm, analysis_incov_error_norm);
% analysis_spot_out_2d = analysis_outcov_error_norm;
% ax = anaPlotSpace(analysis_spot_inout_3d);
% ax = anaPlotSpace(analysis_spot_out_2d);
% ax = anaPlotSpace(analysis_inoutcov_error);

%% 输出时间
disp(['Time used: ', mat2str(mean(time_full), 2), ' / ', mat2str(mean(time_reduce), 2)]);

%% 保存参数
if ~exist('save_valid'), save_valid = false; end
save_postfix = '_als';
if save_valid
    save(['results\mat_auto' save_postfix '.mat'], '-v7.3');
end

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

%% 辅助命令
% mean(analysis_error_h2, 2)
% std(analysis_error_h2, 0, 2)
% 
% mean(analysis_error_hinf, 2)
% std(analysis_error_hinf, 0, 2)
% 
% mean(analysis_incov_error_norm, 3)
% std(analysis_incov_error_norm, 0, 3)
% 
% mean(analysis_outcov_error_norm, 3)
% std(analysis_outcov_error_norm, 0, 3)
% 
% mean(time_full)
% std(time_full)
% mean(time_reduce)
% std(time_reduce)

%% 辅助函数
function [compare_h2, compare_hinf] = sys_compare(sys_ori, sys_ind)
    [num_h2, num_hinf] = anaNorm(sys_ind - sys_ori);
    [den_h2, den_hinf] = anaNorm(sys_ori);
    compare_h2 = abs(num_h2)/den_h2;
    compare_hinf = abs(num_hinf)/den_hinf;
end

