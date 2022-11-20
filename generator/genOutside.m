function out = genOutside(generator_info, period_sample, signal_size, signal_step)
%GENOUTSIDE outside 信号生成
% 给定外部信号, 返回[0, N-1]下的截断信号
% 为减少不同周期间的冲激, 逐元素减去y = (sig(N-1)-sig(0))/(N-1) * n的值
% 输出size @ (size period_sample)

    % 参数提取
    source_path = generator_info.outside_source;
    source_name = generator_info.outside_name;
    source_fsrate_name = generator_info.outside_fs_name;
    source_offest = generator_info.outside_offest;
    source_file = load(source_path, source_name, source_fsrate_name);
    source_signal = source_file.(source_name);
    source_fsrate = source_file.(source_fsrate_name);
    % 参数修正
    if signal_size > size(source_signal, 1), source_signal = repmat(source_signal, [ceil(signal_size/size(source_signal, 1)), 1]); end
    if signal_size < size(source_signal, 1), source_signal = source_signal(1:signal_size, :); end

    % 重采样
    target_fsrate = 1/signal_step;
    [resample_p, resample_q] = rat(target_fsrate/source_fsrate);
    source_signal = resample(source_signal, resample_p, resample_q, 'Dimension', 2);
    % 注意
    source_offest = round(size(source_signal, 2)*source_offest);
    
    % 生成信号 (一个周期), 忽略采样时间
    k = 0:period_sample-1;
    out = zeros(signal_size, period_sample);
    for iter_signal = 1:signal_size
        out(iter_signal, :) = source_signal(iter_signal, source_offest+1:source_offest+period_sample);
        % 修正起终点
        out_fix_k = (out(iter_signal, period_sample)-out(iter_signal, 1))/(period_sample-1);
        out(iter_signal, :) = out(iter_signal, :) - out_fix_k.*k;
    end

end

