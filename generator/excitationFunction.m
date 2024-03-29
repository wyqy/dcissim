function out = excitationFunction(generator_info, samples_period, signal_size, signal_step)
%EXCITATIONFUNCTION functional 信号生成
% 给定信号表达式, 返回[0, N-1]下的截断信号
% 为减少不同周期间的冲激, 逐元素减去y = (sig(N-1)-sig(0))/(N-1) * n的值
% 输出size @ (size period_sample)

    % 参数提取
    func_expression = generator_info.expression;
    % 参数修正
    if signal_size > size(func_expression, 1), func_expression = repmat(func_expression, [ceil(signal_size/size(func_expression, 1)), 1]); end
    if signal_size < size(func_expression, 1), func_expression = func_expression(1:signal_size, :); end

    % 准备输出
    out = zeros(signal_size, samples_period);

    % 转换表达式
    syms n  %#ok<NASGU>
    function_handles = cell(signal_size, 1);
    for iter_signal = 1:signal_size
        function_syms = eval(func_expression(iter_signal));
        function_handles{iter_signal} = matlabFunction(function_syms);
    end

    % 生成信号 (一个周期), 注意采样时间
    k_unit = 0:samples_period-1;
    k_step = signal_step.*k_unit;
    for iter_signal = 1:signal_size
        out(iter_signal, :) = function_handles{iter_signal}(k_step);
        % 修正起终点
        out_fix_k = (out(iter_signal, samples_period)-out(iter_signal, 1))/(samples_period-1);
        out(iter_signal, :) = out(iter_signal, :) - out_fix_k.*k_unit;
    end

end

