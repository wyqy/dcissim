function cov_covariance = anaInnoCov(inno_kalman, inno_ee)
%ANAINNOCOV 给定K, e, 和系统参数, 计算稳态下的输出矩阵

    % 参数计算
    cov_covariance = [inno_kalman*inno_ee*inno_kalman.' inno_kalman*inno_ee; inno_ee.'*inno_kalman.' inno_ee];

end

