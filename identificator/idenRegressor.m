function regressor_struct = idenRegressor(period_sample, frequencies, signal_sample, regressor_type)
%IDENREGRESSOR 准备回归元
% 对v的初值无要求, 因此可以随机/指定初始化
    
    % 参数准备
    frequencies = unique(frequencies, 'sorted');  % 默认升序排序

    % 参数计算
    omega = (2*pi)/period_sample;
    frequencies_size = length(frequencies);
    if frequencies(1) == 0
        v_size = 2*frequencies_size-1;
        harmonic_size = frequencies_size-1;
        harmonic_frequencies = frequencies(2:end);
    else
        v_size = 2*frequencies_size;
        harmonic_size = frequencies_size;
        harmonic_frequencies = frequencies;
    end

    % 返回参数计算
    excite_freq = zeros(v_size, 1);
    excite_phi = zeros(v_size, 1);
    harmonic_phi = harmonicPhaseInit(harmonic_size);
    % 直流部分(如有)
    if frequencies(1) == 0
        excite_freq(1) = 0;
        excite_phi(1) = pi/2;
        location_base = 1;
    else
        location_base = 0;
    end
    % 谐波部分
    for iter = 1:harmonic_size
        excite_freq(location_base+1:location_base+2) = omega*harmonic_frequencies(iter)*ones(2, 1);
        excite_phi(location_base+1:location_base+2) = harmonic_phi(iter)*ones(2, 1);
        excite_phi(location_base+2) = excite_phi(location_base+2) + pi/2;
        location_base = location_base + 2;
    end
    % 初值
    v0 = sin(excite_phi);

    % 准备返回值
    switch regressor_type
        case 'ordinary'  % 返回计算结果
            n = 0:signal_sample-1;
            vn = sin((excite_freq*n) + excite_phi);
            regressor_struct = struct('vn', vn);
        case 'recursive'  % 返回构造参数
            regressor_struct = struct('v0', v0, 'excitation_frequencies', excite_freq, 'excitation_phi', excite_phi);
        otherwise, regressor_struct = 0;
    end

end

function phi = harmonicPhaseInit(harmonic_size)
% 谐波相角初始化, 不需要成对, 每个值即对应一组S^delta矩阵

    % 零初始化
    phi = zeros(harmonic_size, 1);
end

