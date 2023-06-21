function out = excitationChirp(generator_info, samples_period, signal_size, ~)
% EXCITATIONCHIRP 扫频信号生成
% 为减少不同周期间的冲激, 逐元素减去y = (sig(N-1)-sig(0))/(N-1) * n的值
% 输出size @ (size period_sample)

    % 参数提取
    fmax_ratio = generator_info.fmax_ratio;
    % 参数修正
    if signal_size > size(fmax_ratio, 1), fmax_ratio = repmat(fmax_ratio, [ceil(signal_size/size(fmax_ratio, 1)), 1]); end
    if signal_size < size(fmax_ratio, 1), fmax_ratio = fmax_ratio(1:signal_size, :); end

    % 准备输出
    out = zeros(signal_size, samples_period);

    % 生成信号 (一个周期), 注意采样时间
    t = 0:samples_period-1;
    t1 = t(end);
    f0 = 0;
    for iter_signal = 1:signal_size
        % 线性扫频
        f1 = fmax_ratio(iter_signal)*0.5;
        out(iter_signal, :) = chirp(t, f0, t1, f1, 'linear');

        % 修正起终点
        out_fix_k = (out(iter_signal, samples_period)-out(iter_signal, 1))/(samples_period-1);
        out(iter_signal, :) = out(iter_signal, :) - out_fix_k.*t;
    end
end

