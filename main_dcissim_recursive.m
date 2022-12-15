%#ok<*REPMAT>
%#ok<*EXIST>
%#ok<*UNRCH>
%% 确定参数
% 仿真系统
para_xsize = 2; para_ysize = para_xsize; para_usize = 1;
para_xyusize = para_xsize + para_ysize + para_usize;
% 时间参数
para_sim_duration = 600;  % 总时间
para_sim_step = 1/1e1;  % 时间步长
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
un = sim_result_struct.un;
yn = sim_result_struct.yn;
result_original = sim_result_struct.plant_info;

% 随机模型个数
para_experiment_period = fix(0.1*sim_period_samples);
para_experiment_count = fix(sim_excitation_samples/para_experiment_period);
result_full_cell = cell(para_experiment_count, 1);
result_reduce_cell = cell(para_experiment_count, 1);

%% 对每个时间节点运行一次, 不计算方差
%% 辨识方差参数 - dcISSIM (reduced)
clear idenDCISSIMRunner idenISIM; iden_location = 1;
iden_reduce_init = idenDCISSIMLaucher(un, 'y_size', para_ysize, 'u_size', para_usize, 'period_samples', sim_period_samples, 'cutted_periods', para_iden_cutted_period, ...
    'dcissim_type', 'online', 'isim_excitation_type', 'reduced', 'x_size_upbound', para_iden_x_size_bound, 'sim_ss_bdx_type', para_iden_sim_ss_bdx_type, 'sim_ss_d_type', para_iden_sim_ss_d_type, ...
    'sim_x_size_type', 'fixed', 'sim_x_size', para_xsize);
for iter_sample = 1:sim_excitation_samples
    if mod(iter_sample, para_experiment_period) == 0
        result_reduce_cell{iden_location} = idenDCISSIMRunner(iden_reduce_init, yn(:, iter_sample), un(:, iter_sample));
        iden_location = iden_location + 1;
    else
        idenDCISSIMRunner(iden_reduce_init, yn(:, iter_sample), un(:, iter_sample));
    end
end

%% 辨识方差参数 - dcISSIM (full)
clear idenDCISSIMRunner idenISIM; iden_location = 1;
iden_full_init = idenDCISSIMLaucher(un, 'y_size', para_ysize, 'u_size', para_usize, 'period_samples', sim_period_samples, 'cutted_periods', para_iden_cutted_period, ...
    'dcissim_type', 'online', 'isim_excitation_type', 'full', 'x_size_upbound', para_iden_x_size_bound, 'sim_ss_bdx_type', para_iden_sim_ss_bdx_type, 'sim_ss_d_type', para_iden_sim_ss_d_type, ...
    'sim_x_size_type', 'fixed', 'sim_x_size', para_xsize);
for iter_sample = 1:sim_excitation_samples
    
    if mod(iter_sample, para_experiment_period) == 0
        result_full_cell{iden_location} = idenDCISSIMRunner(iden_full_init, yn(:, iter_sample), un(:, iter_sample));
        iden_location = iden_location + 1;
    else
        idenDCISSIMRunner(iden_full_init, yn(:, iter_sample), un(:, iter_sample));
    end
end


%% 验证确定性参数辨识结果
% 准备参数, 第一行为原系统, 第二行为辨识系统
analysis_xsize = zeros(3, para_experiment_count);
analysis_h2 = zeros(3, para_experiment_count);
analysis_hinf = zeros(3, para_experiment_count);
% 计算
for iter_exp = 1:para_experiment_count
    % 准备
    temp_struct = result_original;
    analysis_original_ss = ss(temp_struct.A, temp_struct.B, temp_struct.C, temp_struct.D, para_sim_step);
    temp_struct = result_full_cell{iter_exp};
    analysis_full_idss = idss(temp_struct.A, temp_struct.B, temp_struct.C, temp_struct.D, temp_struct.K, ...
        zeros(temp_struct.x_size, 1), para_sim_step);%, 'NoiseVariance', temp_struct.cov(temp_struct.x_size+1:temp_struct.x_size+para_ysize, temp_struct.x_size+1:temp_struct.x_size+para_ysize));
    temp_struct = result_reduce_cell{iter_exp};
    analysis_reduce_idss = idss(temp_struct.A, temp_struct.B, temp_struct.C, temp_struct.D, temp_struct.K, ...
        zeros(temp_struct.x_size, 1), para_sim_step);%, 'NoiseVariance', temp_struct.cov(temp_struct.x_size+1:temp_struct.x_size+para_ysize, temp_struct.x_size+1:temp_struct.x_size+para_ysize));
    % xsize
    analysis_xsize(1, iter_exp) = para_xsize;
    analysis_xsize(2, iter_exp) = result_full_cell{iter_exp}.x_size;
    analysis_xsize(3, iter_exp) = result_reduce_cell{iter_exp}.x_size;
    % 误差系统的H范数 (仅关注第一个输入到第一个输出)
    [analysis_h2(1, iter_exp), analysis_hinf(1, iter_exp)] = anaNorm(analysis_original_ss);
    [analysis_h2(2, iter_exp), analysis_hinf(2, iter_exp)] = anaNorm(analysis_full_idss);
    [analysis_h2(3, iter_exp), analysis_hinf(3, iter_exp)] = anaNorm(analysis_reduce_idss);
end
% 求平均误差
analysis_criterion_h2 = mean(analysis_h2, 2); analysis_criterion_hinf = mean(analysis_hinf, 2);
disp(['R2 of H2 norm differences of determined paras: ', mat2str(analysis_criterion_h2(1), 2), ' / ', mat2str(analysis_criterion_h2(2), 2), ' / ', mat2str(analysis_criterion_h2(3), 2)]);
disp(['R2 of Hinf norm differences of determined paras: ', mat2str(analysis_criterion_hinf(1), 2), ' / ', mat2str(analysis_criterion_hinf(2), 2), ' / ', mat2str(analysis_criterion_hinf(3), 2)]);


%% 保存参数
if ~exist('save_valid'), save_valid = false; end
save_postfix = '_recursive';
if save_valid
    save(['results\mat_auto' save_postfix '.mat'], '-v7.3');
end

