function ret_struct = idenDCISSIMLaucher(uk_test, varargin)
% IDENDCISSIMLAUCHER discrete-cISSIM系统辨识 - 初始化

    % 输入解析
    parser = inputParser;
    % 输入定义
    addParameter(parser, 'ysize', 1, @(i)(isnumeric(i)&&isscalar(i)));
    addParameter(parser, 'usize', 1, @(i)(isnumeric(i)&&isscalar(i)));
    addParameter(parser, 'xsize_upbound', 1, @(i)(isnumeric(i)&&isscalar(i)));
    addParameter(parser, 'samples_period', 1, @(i)(isnumeric(i)&&isscalar(i)));

    addParameter(parser, 'isstep', 'offline', @(i)(ischar(i)));
    addParameter(parser, 'use_freq', 'full', @(i)(ischar(i)));
    addParameter(parser, 'est_xsize', 'estimate', @(i)(ischar(i)));
    addParameter(parser, 'est_aorder', 'estimate', @(i)(ischar(i)));
    addParameter(parser, 'est_cov', 'simple', @(i)(ischar(i)));
    
    addParameter(parser, 'plant_d', 'null', @(i)(ischar(i)));
    addParameter(parser, 'cov_cross', 'null', @(i)(ischar(i)));
    addParameter(parser, 'prior', struct('xsize', 1, 'aorder', 1), @(i)(isstruct(i)));
    
    % 输入提取
    parse(parser, varargin{:});
    y_size = parser.Results.ysize;  % 输出信号维数
    u_size = parser.Results.usize;  % 输出信号维数
    x_size_upbound = parser.Results.xsize_upbound;  % 状态变量上界
    samples_period = parser.Results.samples_period;  % 单周期采样点数
    
    algo_type = parser.Results.isstep;  % 离线/在线运行
    freq_type = parser.Results.use_freq;  % ISIM激励方式
    xsize_est_type = parser.Results.est_xsize;  % 状态变量估计方式
    aorder_est_type = parser.Results.est_aorder;  % 相关函数计算量估计方式
    als_est_type = parser.Results.est_cov;  % 方差估计方法

    plant_d_type = parser.Results.plant_d;  % D矩阵是否为零(是否直通)
    cov_cross_type = parser.Results.cov_cross;  % 是否具有协方差
    prior = parser.Results.prior;  % 系统阶数和相关函数计算量的先验值
    

    % 运行时准备
    % 清除临时存储
    clear idenISIM idenDCISSIMRunner
    % 准备辨识频率
    freq_bound = persistentExcitationCondition(x_size_upbound, u_size);
    switch freq_type
        case 'full'  % 全激励频率
            freqs_list = 0:fix(samples_period/2);
        case 'reduced'  % 部分激励频率
            uk_test = uk_test(:, samples_period+1:2*samples_period);  % 从第二个周期开始
            freqs_list = reducedFreqList(uk_test, samples_period, freq_bound);
        otherwise, freqs_list = 0;
    end
    % 转换为角频率
    omega = (2*pi)/samples_period;
    freqs_list = omega.*freqs_list;
    % 初始化临界系统
    [mat_s, regressor] = invariantIniter(freqs_list);
    % 在线辨识 - 初始化参数
    if strcmp(algo_type, 'online') || strcmp(algo_type, 'online-test')
        v_size = size(mat_s, 1);
        idenISIM(zeros(y_size, 1), zeros(u_size, 1), zeros(v_size, 1), 'recursive');
    end

    % 返回值
    ret_struct = struct( 'mat_s', mat_s, 'regressor', regressor, 'prior', prior, ...
        'y_size', y_size, 'u_size', u_size, 'x_size_upbound', x_size_upbound, 'samples_period', samples_period,  ...
        'algo_type', algo_type, 'freq_type', freq_type, 'xsize_est_type', xsize_est_type, 'aorder_est_type', aorder_est_type, 'als_est_type', als_est_type, ...
        'plant_d_type', plant_d_type, 'cov_cross_type', cov_cross_type);

end

function freq_bound = persistentExcitationCondition(x_size_upbound, u_size)
% 持续激励条件(未知信号阶数, 用上界近似)
    freq_bound = u_size*2*x_size_upbound;
end

function freqs_list = reducedFreqList(uk_test, samples_period, freqs_bound)
% 选择频率
% 方法: 满足持续激励条件(再多几个频率点), 按照输入信号对应fft的幅值从高到低依次选择
% 最后总是加上零频率

    % 参数计算
    harmonic_size = fix(samples_period/2)+1;  % 频率上限
    freqs_bound = min(round(freqs_bound*3), harmonic_size);  % 升维度, 3为经验参数

    % FFT
    ufreq = fft(uk_test, samples_period, 2);
    ufreq_abs = abs(ufreq(:, 1:harmonic_size));

    % 提取频率点
    frequencies_selection = peakFinder(ufreq_abs, freqs_bound);
    
    % plot
    % figure; plot(ufreq_abs); hold on; stem(max(ufreq_abs, [], 'all')*frequencies_selection);
    % 返回值
    freqs_list = 0:harmonic_size-1;
    freqs_list = freqs_list(frequencies_selection == 1);
    % 加上零频率
    if freqs_list(1) ~= 0, freqs_list = [0 freqs_list]; end

end

function selected_freq = peakFinder(ufreq_abs, peak_number)

    % 从离群值点中取频率小且幅值大者
    % 参数计算
    freq_size = size(ufreq_abs, 2);
    selected_freq = zeros(1, freq_size);

    % 按通道取离群值
    outlier_thresold = max(100 * (1 - (peak_number/freq_size)) - 3, 3);  % -保证候选频率数量
    peak_outlier = isoutlier(ufreq_abs, 'percentiles', [0 outlier_thresold], 2);  % 按通道选取
    % 标记频率值点
    freq_order = 0:freq_size-1;
    freq_order = log10(freq_order + 1);
    % 按通道标记相对大小
    [~, maxval_order] = sort(ufreq_abs, 2, 'descend');  % maxval_order(i) = k: 表示第k个频率对应的幅值为第i大
    [~, maxval_order] = sort(maxval_order, 2, 'ascend');  % maxval_order(i) = k: 表示第i个频率对应的幅值为第k大
    maxval_order = min(maxval_order, [], 1);  % 取所有通道中的大小排名靠前者
    
    % 计算加权
    weights = (peak_outlier > 0).*(freq_order + min(maxval_order, [], 1)) + ...  % 离群点中频率小且幅值大, 则权重小
              (peak_outlier <= 0).*realmax;  % 不选择非离群点
    
    % 选择频率
    % 默认选择零频率
    selected_freq(1) = 1;
    weights(1) = realmax;
    % 低中高频均匀选择, 频段长度为对数递增, 每个频段选择权值最小者, 选不满则重复该过程
    count_select = 1;
    band_length = (log10(freq_size-1) - log10(1))/peak_number;
    while count_select < peak_number
        upperband = 0;
        while upperband < freq_size-1
            lowerband = upperband + 1;
            upperband = min(max(ceil(10^(log10(upperband) + band_length)), lowerband + 1), freq_size-1);
            [~, select_idx] = min(weights(lowerband:upperband));
            if weights(select_idx + lowerband - 1) < realmax
                selected_freq(select_idx + lowerband - 1) = 1;
                weights(select_idx + lowerband - 1) = realmax;
                count_select = count_select + 1;
            end
        end
    end

end

function [mat_s, regressor] = invariantIniter(raw_freq_list)
% S矩阵和回归元计算
    
    % 参数准备
    freq_size = length(raw_freq_list);
    raw_freq_list = unique(raw_freq_list, 'sorted');  % 默认升序排序

    % 分离常量和谐波
    if raw_freq_list(1) == 0
        v_size = 2*freq_size-1;
        harmonic_size = freq_size-1;
        harmonic_list = raw_freq_list(2:end);
    else
        v_size = 2*freq_size;
        harmonic_size = freq_size;
        harmonic_list = raw_freq_list;
    end
    % 初始化返回值
    mat_s = zeros(v_size, v_size);
    freq_list = zeros(v_size, 1);
    phi_list = zeros(v_size, 1);

    % 相角初始化
    harmonic_phi = harmonicPhaseIniter(harmonic_size);

    % 计算返回值
    % 直流部分(如有)
    if raw_freq_list(1) == 0
        mat_s(1, 1) = 1;
        freq_list(1) = 0;
        phi_list(1) = pi/2;
        loc_base = 1;
    else
        loc_base = 0;
    end
    % 谐波部分
    for iter = 1:harmonic_size
        iter_omega = harmonic_list(iter);
        iter_phi = harmonic_phi(iter);
        mat_s(loc_base+1:loc_base+2, loc_base+1:loc_base+2) = [cos(iter_omega) sin(iter_omega); -sin(iter_omega) cos(iter_omega)];
        freq_list(loc_base+1:loc_base+2) = iter_omega*ones(2, 1);  % [omega omega]
        phi_list(loc_base+1:loc_base+2) = iter_phi*ones(2, 1);  % [sin(phi) cos(phi)]
        phi_list(loc_base+2) = phi_list(loc_base+2) + pi/2;  % cos(phi)
        loc_base = loc_base + 2;
    end
    % 初值
    v0 = sin(phi_list);

    % 准备返回值
    regressor = struct('v0', v0, 'freq_list', freq_list, 'phi_list', phi_list);
    
end

function phi = harmonicPhaseIniter(harmonic_size)
% 谐波相角初始化, 不需要成对, 每个值即对应一组S^delta矩阵

    % 零初始化
    phi = zeros(harmonic_size, 1);

end
