function ret_covariance = plantCovariance(plant_info, xyu_snr, xyuk)
%PLANTCOVARIANCE 根据返回满足指定信噪比的协方差矩阵
    % 参数提取
    covariance = plant_info.cov;

    % 参数计算
    xyu_size = size(xyuk, 1);

    % 计算放缩后的协方差矩阵
    signal_power = signalPowermeter(xyuk);
    scale_covariance = noiseAutocovarianceCalculator(covariance, xyu_snr, signal_power);

    % 将NaN替换为0
    nonan_selection = ~isnan(diag(scale_covariance));
    ret_covariance = zeros(xyu_size, xyu_size);
    ret_covariance(nonan_selection, nonan_selection) = scale_covariance(nonan_selection, nonan_selection);

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
    % 计算放缩系数的开方 -> 使得噪声的自协方差满足要求的SNR
    autocovariance = diag(covariance);
    scale_vector = sqrt(noise_power./autocovariance);
    % 应用到协方差矩阵中
    scale_matrix = diag(scale_vector);
    scale_covariance = scale_matrix*covariance*scale_matrix;

end

