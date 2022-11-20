function [fit_f, fit_vaf] = anaFit(ori_data, ind_data, datasize)
%ANAFIT 分别给定原系统和辨识系统对同一输入的响应数据, 输出指标
% 指标:
% 1. Normalized fit Metric: F = 100 (1 - ||y-\hat{y}||_2 / ||y-\bar{y}||_2)
%    [L. Ljung, System Identification Theory for the User, 2nd ed.]
% 2. Variance Account for Metric: VAF = 100 (1 - var(y-\hat{y} / var(y)));
% (Akaike's Final Prediction Error?)

    % 计算参数
    sampleSize = size(ori_data, 1);
    
    % 指标
    % F
    fit_f_mean_original = mean(ori_data, 1);
    fit_f_num = ori_data - ind_data;
    fit_f_den = ori_data - repmat(fit_f_mean_original, [sampleSize, 1]);
    fit_f = vecnorm(fit_f_num, 2, 1) ./ vecnorm(fit_f_den, 2, 1);
    fit_f = 100 .* (1 - fit_f);
    % VAF (equivalent to r^2)
    fit_vaf_num = ori_data - ind_data;
    fit_vaf = var(fit_vaf_num, 0, 1) ./ var(ori_data, 0, 1);
    fit_vaf = 100 .* (1 - fit_vaf);
    
    % 输出 (可选)
    fprintf('Model Fit Metric:\n');
    specs = specBuilder("F: ", datasize);
    fprintf(specs, fit_f);
    specs = specBuilder("VAF: ", datasize);
    fprintf(specs, fit_vaf);

end

function formatSpec = specBuilder(begin_str, vec_size)
    formatSpec = begin_str;
    for iter = 1:vec_size, formatSpec = formatSpec + "%.4f, "; end
    formatSpec = formatSpec + "\n";
end

