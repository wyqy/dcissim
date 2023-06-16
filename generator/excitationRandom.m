function out = excitationRandom(generator_info, samples_period, signal_size, ~)
%EXCITATIONRANDOM random 信号生成
% 给定信号表达式, 返回[0, N-1]下的截断信号
% 为减少不同周期间的冲激, 逐元素减去y = (sig(N-1)-sig(0))/(N-1) * n的值
% 输出size @ (size period_sample)

    % 参数提取
    random_mu = generator_info.mu;
    random_covariance = generator_info.covariance;
    % 参数修正
    if signal_size > size(random_mu, 1)
        random_mu = repmat(random_mu, [ceil(signal_size/size(random_mu, 1)), 1]);
        random_covariance = repmat(random_covariance, [ceil(signal_size/size(random_covariance, 1)), 1]);
    end
    if signal_size < size(random_mu, 1)
        random_mu = random_mu(1:signal_size, :);
        random_covariance = random_covariance(1:signal_size, :);
    end

    % 准备输出
    out = zeros(signal_size, samples_period);
    % 准备种子
    seed = rng().Seed;
    noise_seed = RandStream.create('mrg32k3a', 'NumStreams', signal_size, 'Seed', seed, 'CellOutput', true);
    
    % 生成信号 (一个周期), 注意采样时间
    k_unit = 0:samples_period-1;
    random_covariance = diag(random_covariance);
    for iter_signal = 1:signal_size
        out(iter_signal, :) = randn(noise_seed{iter_signal}, 1, samples_period);
        % 修正起终点
        out_fix_k = (out(iter_signal, samples_period)-out(iter_signal, 1))/(samples_period-1);
        out(iter_signal, :) = out(iter_signal, :) - out_fix_k.*k_unit;
    end
    % 应用均值和方差
    out = sqrtm(random_covariance)*out + random_mu;

end

