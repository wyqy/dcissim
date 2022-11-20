function out = genMultisine(generator_info, signal_size, signal_step)
%GENMULTISINE multisine 信号生成
% 给定频率范围和输入范围, 返回multisine
% 输出size @ (size period_sample)

    % 参数提取
    continuous_frequencies = generator_info.multisine_freqs;
    signal_interval = generator_info.multisine_usat;
    has_zero_frequency = sum(mean(signal_interval, 2).^2) ~= 0;
    % 参数修正
    if signal_size > size(signal_interval, 1), signal_interval = repmat(signal_interval, [ceil(signal_size/size(signal_interval, 1)), 1]); end
    if signal_size < size(signal_interval, 1), signal_interval = signal_interval(1:signal_size, :); end

    % 计算周期
    [period, ~] = leastCommonPeriod(continuous_frequencies);
    period_sample = fix(2*pi*period/signal_step);
    omega = (2*pi)/period_sample;
    % 计算频率
    discrete_frequencies = frequenciesConverter(continuous_frequencies, signal_step, period_sample, has_zero_frequency);
    frequencies_number = length(discrete_frequencies);
    % 计算幅值
    amplitudes = amplitudesBuilder(signal_size, signal_interval, frequencies_number);
    % 计算相位
    phases = phaseBuilder(signal_size, frequencies_number, has_zero_frequency);

    % 生成信号 (一个周期)
    n = 0:period_sample-1;
    discrete_frequencies = reshape(discrete_frequencies, [], 1);
    phases = reshape(phases, [], 1);
    out = amplitudes*sin(omega.*discrete_frequencies*n + phases);

end

function discrete_frequencies = frequenciesConverter(continuous_frequencies, signal_step, period_sample, has_zero_frequency)
    
    % 将连续时域频率转换为对应采样时间下的离散频率值
    inverse_convertor = (signal_step*period_sample)/(2*pi);
    discrete_frequencies = inverse_convertor .* continuous_frequencies;
    % 取整
    discrete_frequencies = round(discrete_frequencies);
    % 去掉过大的频率点
    discrete_frequencies = discrete_frequencies(discrete_frequencies <= fix(period_sample/2));
    % 补充零频率
    if has_zero_frequency, discrete_frequencies = [0 discrete_frequencies]; end

end

function [retPeriod, retMultiply] = leastCommonPeriod(continuous_frequencies)
% 计算multisine的周期和对应各频率的倍数
% 使用符号表达(symbolic representation)
% V = [p/q s/t ...];
% t1 = 2*pi*(q/p) => q/p; t2 => t/s;
% T = a*t1 = b*t2;
% a/b = tp/qs;
    
    % 转为符号变量
    sym_freq = sym(continuous_frequencies, 'r');

    % 计算周期
    n = length(sym_freq);  % 输入长度
    [freq_num, freq_den] = numden(sym_freq);  % 输入分子, 分母
    if n > 1
        for iter_n = 1:n-1
            % 计算[a b]
            a = freq_num(iter_n)*freq_den(iter_n+1);
            b = freq_num(iter_n+1)*freq_den(iter_n);
            mul = gcd(a, b);
            % 计算中间量period
            a = a ./ mul;
            period = a*sym_freq(iter_n);
            % 替换N, D
            [nn, dd] = numden(period);
            freq_num(iter_n+1) = nn;
            freq_den(iter_n+1) = dd;
        end
    else
        period = freq_den/freq_num;
    end
    
    
    % 返回值
    retPeriod = double(period);
    retMultiply = repmat(period, [1 n]);
    retMultiply = retMultiply .* sym_freq;
    retMultiply = double(retMultiply);
end

function amplitudes = amplitudesBuilder(signal_size, signal_interval, frequencies_number)
% 方法:
% 将直流信号幅值定为(max+min)/2, 
% 余下正弦信号(正)幅值组合, 使得总和为(max-min)/2, 然后依据公式还原所需的Amp
% 由于是欠定方程, 且信号幅值不能过低(否则不易辨识), 采取如下算法:
% 1. 全部频率(除常值外)的幅值按顺序排成一个数组:[a1, a2, ..., ai, ..., aq], 设Σai = K
% 2. 为所有幅值ai分配最低值K/(2q): ai = K/(2q);
% 3. 为所有幅值ai分配剩下的一半: 将全部频率分为二半(若为奇数则少者在低端),
%    随机分配剩下的部分给低端, 剩下的给高端.
% 4. 递归运行.
% 注意: 由于幅值的正负可由相位决定, 因此幅值均假定为正
% 输出size @ (size frequencies_number)
    
    % 参数计算
    has_zero_frequency = sum(mean(signal_interval, 2).^2) ~= 0;
    harmonic_number = frequencies_number - double(has_zero_frequency);

    % 每行第一为下限, 第二为上限
    signal_interval = sort(signal_interval, 2, 'ascend');

    % 准备随机数种子
    rng('shuffle', 'twister'); seed = randi(intmax('uint32'));
    % seed = 142857;
    rand_stream = RandStream.create('mrg32k3a', 'NumStreams', signal_size, 'Seed', seed, 'CellOutput', true);
    
    % 确定常数幅值 (如果有)
    if has_zero_frequency, const_amplitudes = (signal_interval(:, 2)+signal_interval(:, 1)) ./ 2; end

    % 确定谐振幅值
    capacity = (signal_interval(:, 2) - signal_interval(:, 1)) ./ 2;
    % 确定谐振幅值本底值
    temp = capacity ./ (2*harmonic_number);
    harmonic_amplitudes = repmat(temp, [1 harmonic_number]);
    % 确定谐振幅值随机分配值
    % 递归运算
    harmonic_amplitudes = harmonic_amplitudes + recursiveDistributor(harmonic_number, capacity/2);
    % 递归函数
    function amplitudes = recursiveDistributor(harmonic_number, signal_capacity)
        % 递归到底
        if harmonic_number == 1, amplitudes = signal_capacity; return; end
        
        % 递归未到底
        % 计算分配的capacity
        lowcap = zeros(signal_size, 1);
        for iter = 1:signal_size, lowcap(iter) = signal_capacity(iter)*rand(rand_stream{iter}); end
        highcap = signal_capacity - lowcap;
        % 计算分配的lenth
        lowlen = fix(harmonic_number/2);
        highlen = harmonic_number - lowlen;
        % 递归运算
        amplitudes = [recursiveDistributor(lowlen, lowcap) recursiveDistributor(highlen, highcap)];
    end

    % 返回可能包含零频率的真实值
    if has_zero_frequency, amplitudes = [const_amplitudes harmonic_amplitudes];
    else, amplitudes = harmonic_amplitudes; end
    
end

function phases = phaseBuilder(signal_size, frequencies_number, has_zero_frequency)
% 按照给定的规则简单生成一个相位参数矩阵
% 输出size @ (size frequencies_number)

    % 计算参数
    harmonic_number = frequencies_number - double(has_zero_frequency);

    % 按照指定角度初始化
    harmonic_phases = 45.*ones(signal_size, harmonic_number);
    % 转换为弧度制
    harmonic_phases = deg2rad(harmonic_phases);
    
    % 返回可能包含零频率的真实值
    if has_zero_frequency, phases = [zeros(signal_size, 1) harmonic_phases];
    else, phases = harmonic_phases; end

end

