%#ok<*REPMAT>
%% 确定参数
% 仿真系统
para_xsize = 4; para_ysize = 1; para_usize = 1;

% 是否全量辨识
para_isfull = true; %#ok<*UNRCH>
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
para_iden_x_size_bound = 10;   % X的变量数上界
para_iden_sim_ss_bdx_type = 'analytical';  % BDX辨识方法: 'analytical', 'optimize'
para_iden_sim_ss_d_type = 'null';  % 是否假定D矩阵为0
para_iden_cov_cross_type = 'null';  % 是否假定存在互协方差

% 随机模型个数
para_experiment_count = 1;
result_original_cell = cell(para_experiment_count, 1);
result_identified_cell = cell(para_experiment_count, 1);

% 随机数种子
% rng('shuffle', 'twister'); seed = randi(intmax('uint32'), 'uint32'); % 非固定参数
seed = uint32(142857);  % 固定参数
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
    
    %% 辨识系统参数 - 离线
    % 初始化
    iden_offline_init_struct = idenDCISSIMLaucher(un, 'y_size', size(yn, 1), 'u_size', size(un, 1), 'period_samples', sim_period_samples, 'cutted_periods', para_iden_cutted_period, ...
        'dcissim_type', 'offline', 'isim_excitation_type', para_iden_isim_excitation_type, 'x_size_upbound', para_iden_x_size_bound, 'sim_ss_bdx_type', para_iden_sim_ss_bdx_type, 'sim_ss_d_type', para_iden_sim_ss_d_type, 'cov_cross_type', para_iden_cov_cross_type, ...
        'sim_x_size_type', 'fixed', 'cov_order_type', 'estimate', 'sim_x_size', para_xsize);
    % 离线辨识
    result_identified_cell{iter_exp} = idenDCISSIMRunner(iden_offline_init_struct, yn, un);
end

%% 验证确定性参数辨识结果
% 准备参数, 第一行为原系统, 第二行为辨识系统
analysis_xsize = zeros(2, para_experiment_count);
analysis_eig = zeros(2, para_xsize, para_experiment_count);
analysis_h2 = zeros(2, para_experiment_count);
analysis_hinf = zeros(2, para_experiment_count);
% 计算
for iter_exp = 1:para_experiment_count
    % 准备
    temp_struct = result_original_cell{iter_exp};
    analysis_original_ss = ss(temp_struct.A, temp_struct.B, temp_struct.C, temp_struct.D, para_sim_step);
    temp_struct = result_identified_cell{iter_exp};
    analysis_indentified_idss = idss(temp_struct.A, temp_struct.B, temp_struct.C, temp_struct.D, temp_struct.K, ...
        zeros(temp_struct.x_size, 1), para_sim_step, 'NoiseVariance', temp_struct.cov(temp_struct.x_size+1:temp_struct.x_size+para_ysize, temp_struct.x_size+1:temp_struct.x_size+para_ysize));
    % xsize
    analysis_xsize(1, iter_exp) = para_xsize;
    analysis_xsize(2, iter_exp) = result_identified_cell{iter_exp}.x_size;
    % 特征值
    analysis_eig(1, :, iter_exp) = sort(eig(analysis_original_ss.A), 'descend', 'ComparisonMethod', 'abs');
    temp_eig = sort(eig(analysis_indentified_idss.A), 'descend', 'ComparisonMethod', 'abs');
    if length(temp_eig) > para_xsize, temp_eig = temp_eig(1:para_xsize);
    else, temp_eig = [temp_eig(1:length(temp_eig)); zeros(para_xsize-length(temp_eig), 1)]; end
    analysis_eig(2, :, iter_exp) = temp_eig;
    % H范数
    [analysis_h2(1, iter_exp), analysis_hinf(1, iter_exp)] = anaNorm(analysis_original_ss);
    [analysis_h2(2, iter_exp), analysis_hinf(2, iter_exp)] = anaNorm(analysis_indentified_idss);
end
disp(['H2 norm differences of determined paras: ', mat2str(mean(power(diff(analysis_h2), 2)), 2)]);
disp(['Hinf norm differences of determined paras: ', mat2str(mean(power(diff(analysis_hinf), 2)), 2)]);

%% 验证随机参数辨识结果
% 准备参数
analysis_cov_count = 50;  % 验证阶数
analysis_cov_norm = zeros(2, analysis_cov_count, para_experiment_count);
% 计算
for iter_exp = 1:para_experiment_count
    % 方差
    if para_isfull
        temp_cov = result_original_cell{iter_exp}.cov(1:para_xsize+para_ysize, 1:para_xsize+para_ysize);
        [~, ~, ~, analysis_cov_norm(1, :, iter_exp)] = anaCov(temp_cov, result_original_cell{iter_exp}, analysis_cov_count);
        temp_cov = result_identified_cell{iter_exp}.cov(1:result_identified_cell{iter_exp}.x_size+para_ysize, 1:result_identified_cell{iter_exp}.x_size+para_ysize);
        [~, ~, ~, analysis_cov_norm(2, :, iter_exp)] = anaCov(temp_cov, result_identified_cell{iter_exp}, analysis_cov_count);
    end
end
if para_isfull, disp(['H2 norm differences of covariances: ', mat2str(mean(power(diff(analysis_cov_norm), 2), 'all'), 2)]); end

%% H2差异最大, 中位数和最小的Bode图
[~, analysis_hinf_sort] = sort(power(diff(analysis_h2), 2), 'descend');

% analysis_bode_location = analysis_hinf_sort(1); fig = anaPlotBode(result_original_cell{analysis_bode_location}, result_identified_cell{analysis_bode_location}, 'sample', para_sim_step); sgtitle(fig, 'Bode plot for maximum d\_{Hinf}');
% analysis_bode_location = analysis_hinf_sort(ceil(para_experiment_count/2)); fig = anaPlotBode(result_original_cell{analysis_bode_location}, result_identified_cell{analysis_bode_location}, 'sample', para_sim_step); sgtitle(fig, 'Bode plot for medium d\_{Hinf}');
analysis_bode_location = analysis_hinf_sort(end); fig = anaPlotBode(result_original_cell{analysis_bode_location}, result_identified_cell{analysis_bode_location}, 'sample', para_sim_step); sgtitle(fig, 'Bode plot for minimum d\_{Hinf}');

%% Covariance差异最小的矩估计图
if para_isfull
    [~, analysis_cov_sort] = sort(sum(shiftdim(power(diff(analysis_cov_norm), 2)), 1), 'descend');
    
    analysis_cov_location = analysis_cov_sort(end);
    temp_cov_ori = shiftdim(analysis_cov_norm(1, :, analysis_cov_location));
    temp_cov_ind = shiftdim(analysis_cov_norm(2, :, analysis_cov_location));
    ax = anaPlotCov(temp_cov_ori, temp_cov_ind); title(ax, 'Coavariance plot for minimum d\_{Hinf}');
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
