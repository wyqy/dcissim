function ret_struct = genExcitationSignal(varargin)
% GENEXCITATIONSIGNAL 返回激励信号

    % 输入解析
    parser = inputParser;
    % 输入定义
    addParameter(parser, 'type', 'multisine', @(i)(ischar(i)));
    addParameter(parser, 'signal_size', 1, @(i)(isnumeric(i)&&isscalar(i)));
    addParameter(parser, 'duration', 1, @(i)(isnumeric(i)&&isscalar(i)));
    addParameter(parser, 'step', 1, @(i)(isnumeric(i)&&isscalar(i)));
    addParameter(parser, 'period', 1, @(i)(isnumeric(i)&&isscalar(i)));
    % 输入提取
    parse(parser, varargin{:});
    type = parser.Results.type;  % 输出信号种类
    signal_size = parser.Results.signal_size;  % 输出信号维数
    signal_duration = parser.Results.duration;  % 总时间
    signal_step = parser.Results.step;  % 时间步长
    signal_period = parser.Results.period;  % 单周期时间

    % 参数计算
    signal_sample = fix(signal_duration/signal_step);
    samples_period = fix(signal_period/signal_step);
    
    % 单周期信号生成
    switch type
        case 'multisine'
            signal = excitationMultisine(struct('freqs', [0.1,0.3,0.5,0.7,0.9,1,3,5,7,9], 'usat', [10, 20]), ...
                signal_size, signal_step);
        case 'chirp'
            signal = excitationChirp(struct('fmax_ratio', 0.8), samples_period, signal_size, signal_step);
        case 'function'
            signal = excitationFunction(struct('expression', "sin(n)+n^2+n+3"), ...
                samples_period, signal_size, signal_step);
        case 'outside'
            signal = excitationOutside(struct('source', 'parameter\source.mat', 'name', 'yt', 'fs_name', 'Fs', 'offest', 1/84.8), ...
                samples_period, signal_size, signal_step);
        case 'random'
            signal = excitationRandom(struct('mu', 10, 'covariance', 20), ...
                samples_period, signal_size, signal_step);
        otherwise, signal = 0;
    end
    % 修正周期采样点数
    samples_period = size(signal, 2);
    % 多周期信号复制
    signal = loopSignal(signal, signal_size, signal_sample, samples_period);
    
    % 返回值
    ret_struct = struct('signal', signal, 'samples', signal_sample, 'samples_period', samples_period);

end

function out = loopSignal(source, signal_size, signal_sample, samples_periods)

    % 参数计算
    loop_number = ceil(signal_sample/samples_periods);
    signal_length = loop_number*samples_periods;
    % 返回值初始化
    out = zeros(signal_size, signal_length);

    % 复制
    location_base = 0;
    for iter = 1:loop_number
        out(:, location_base+1:location_base+samples_periods) = source;
        location_base = location_base + samples_periods;
    end
    % 截断
    out = out(:, 1:signal_sample);

end

