function transferred_noise = genNoiser(plant_info, samples, xyu_snr, xyun)
%GENNOISER 根据返回满足指定信噪比和协方差矩阵的扰动信号

    % 参数提取
    seed = plant_info.Seed;
    covariance = plant_info.Cov;

    % 参数计算
    xyu_size = size(xyun, 1);

    % 计算初始正态分布
    rs = RandStream.create('mrg32k3a', 'NumStreams', xyu_size, 'Seed', seed, 'CellOutput', true);
    noise = zeros(xyu_size, samples);
    for iter = 1:xyu_size, noise(iter, :) = randn(rs{iter}, 1, samples); end

    % 计算放缩后的协方差矩阵
    signal_power = signalPowermeter(xyun);
    scale_covariance = noiseAutocovarianceCalculator(covariance, xyu_snr, signal_power);

    % 计算不含NaN的snr的正态分布
    scale_covariance_selection = ~isnan(diag(scale_covariance));
    scale_covariance_rational = scale_covariance(scale_covariance_selection, scale_covariance_selection);
    scale_transfer = zeros(xyu_size, xyu_size);
    scale_transfer(scale_covariance_selection, scale_covariance_selection) = sqrtm(scale_covariance_rational);
    transferred_noise = scale_transfer*noise;

end

function signal_power = signalPowermeter(signal)
% 得到指定信号的功率, 注意维度

    signal_power = bandpower(signal.');
    signal_power = signal_power.';

end

function scale_covariance = noiseAutocovarianceCalculator(covariance, snr, signal_power)
% 返回变换后的协方差矩阵

    % 计算噪声功率
    noise_power = signal_power./(10.^(snr./10));
    % 计算放缩系数
    autocovariance = diag(covariance);
    scale_vector = sqrt(noise_power./autocovariance);
    % 应用到协方差矩阵中
    scale_matrix = diag(scale_vector);
    scale_covariance = scale_matrix*covariance*scale_matrix;

end

