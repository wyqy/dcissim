%#ok<*REPMAT>
%% 确定参数
% 仿真系统参数结构体
sim_plant_struct = load('parameter\paraPlant.mat', 'para');
sim_plant_struct = sim_plant_struct.para;
% 仿真系统
para_xsize = size(sim_plant_struct.A, 1);
para_ysize = size(sim_plant_struct.C, 1);
para_usize = size(sim_plant_struct.B, 2);

% 时间参数
para_sim_duration = 200;  % 总时间
para_sim_step = 1/1e3;  % 时间步长
% 输入信号 - 确定部分
para_sim_period = 20;  % 输入信号周期 (可选)
% 输入信号 - 随机部分
para_sim_x_snr = 30;  % 状态信号基础信噪比
para_sim_y_snr = 30;  % 输出信号基础信噪比
para_sim_u_snr = 30;  % 输入信号基础信噪比
para_sim_xyu_snr = [repmat(para_sim_x_snr, [para_xsize 1]); repmat(para_sim_y_snr, [para_ysize 1]); repmat(para_sim_u_snr, [para_usize 1])];
% 算法参数
para_iden_cutted_period = 1;  % 切除最开始未进入中心流形的部分
para_iden_x_size_bound = 10;   % X的变量数上界
para_iden_isim_excitation_type = 'reduced';  % 使用激励信号的数量
para_iden_sim_ss_bdx_type = 'analytical';  % BDX辨识方法: 'analytical', 'optimize'
para_iden_sim_ss_d_type = 'null';     % 是否假定D矩阵为0
para_iden_cov_cross_type = 'null';  % 是否假定存在互协方差

% 参数计算
para_sim_samples = fix(para_sim_duration/para_sim_step);  % 仿真采样数
para_sim_t = 0:para_sim_step:(para_sim_samples-1)*para_sim_step;  % 仿真时刻(可选)
para_sim_t_end = para_sim_t(end);  % 仿真终止时间

% 随机数种子
rng('shuffle', 'twister');
seed = randi(intmax('uint32'), 'uint32');
rng(seed, 'simdTwister');

%% 数据生成
% 激励信号生成
sim_excitation_struct = genExcitationSignal('signal_size', para_usize, 'type', 'outside', ...
    'duration', para_sim_duration, 'step', para_sim_step, 'period', para_sim_period);
% 导出信息
sim_excitation_signal = sim_excitation_struct.signal;
sim_excitation_samples = sim_excitation_struct.samples;
sim_period_samples = sim_excitation_struct.period_samples;  % 单周期仿真采样数

% 仿真状态生成
sim_result_struct = plantSimulation('period_samples', sim_period_samples, ...
    'excitation', sim_excitation_signal, 'snr', para_sim_xyu_snr);
un = sim_result_struct.un;
yn = sim_result_struct.yn;
period_samples = sim_result_struct.period_samples;

%% 辨识系统参数 - 非全量在线
% 初始化
iden_init_struct = idenDCISSIMLaucher(un, 'y_size', size(yn, 1), 'u_size', size(un, 1), 'period_samples', period_samples, 'cutted_periods', para_iden_cutted_period, ...
    'dcissim_type', 'online', 'isim_excitation_type', para_iden_isim_excitation_type, 'x_size_upbound', para_iden_x_size_bound, 'sim_ss_bdx_type', para_iden_sim_ss_bdx_type, 'sim_ss_d_type', para_iden_sim_ss_d_type, 'cov_cross_type', para_iden_cov_cross_type, ...
    'sim_x_size_type', 'fixed', 'cov_order_type', 'fixed', 'sim_x_size', 4, 'cov_order', 50);
% 在线辨识
iden_log_period = fix(0.1*period_samples);
iden_log_size = fix(sim_excitation_samples/iden_log_period);
iden_result_struct = cell(iden_log_size, 1);
iden_location = 1;
for iter_data = 1:sim_excitation_samples
    if mod(iter_data, iden_log_period) == 0
        iden_result_struct{iden_location} = idenDCISSIMRunner(iden_init_struct, yn(:, iter_data), un(:, iter_data));
        iden_location = iden_location + 1;
    else
        idenDCISSIMRunner(iden_init_struct, yn(:, iter_data), un(:, iter_data));
    end
end


%% 验证辨识结果
% 参数范数误差
analysis_ynorm2_identified = zeros(iden_log_size, 1);
for iter_ana = 1:iden_log_size
    [analysis_ynorm2_identified(iter_ana), ~] = anaNorm(iden_result_struct{iter_ana}.Y);
end
figure;
plot((1:iden_log_period:sim_excitation_samples)./period_samples, analysis_ynorm2_identified); hold on;
plot((1:iden_log_period:sim_excitation_samples)./period_samples, analysis_ynorm2_identified(end).*ones(iden_log_size, 1))
% 状态空间模型
analysis_identified_x_size = size(iden_result_struct{end}.A, 1);
analysis_y_size = size(iden_result_struct{end}.C, 1);
analysis_x0 = zeros(analysis_identified_x_size, 1);
analysis_original_ss = ss(sim_plant_struct.A, sim_plant_struct.B, sim_plant_struct.C, sim_plant_struct.D, para_sim_step);
analysis_indentified_idss = idss(iden_result_struct{end}.A, iden_result_struct{end}.B, iden_result_struct{end}.C, iden_result_struct{end}.D, iden_result_struct{end}.K, ...
    analysis_x0, para_sim_step, 'NoiseVariance', iden_result_struct{end}.cov_all(analysis_identified_x_size+1:analysis_identified_x_size+analysis_y_size, analysis_identified_x_size+1:analysis_identified_x_size+analysis_y_size));
% 特征值
disp(['Eig(A) of original system are: ', mat2str(eig(analysis_original_ss.A), 2)]);
disp(['Eig(A) of identified system are: ', mat2str(eig(analysis_indentified_idss.A), 2)]);
% Bode图
handle_bode = anaPlotBode(analysis_original_ss, analysis_indentified_idss);

