function [moments, hnorms] = anaCov(covariance, ss_struct, out_size)
%ANACOV 根据协方差和(渐近稳定)系统矩阵计算稳态输出矩(输入的矩仅包含状态和输出)
% 同时返回计算结果和Hnorm-2对比
    
    % 参数提取
    mat_a = ss_struct.A;
    mat_c = ss_struct.C;
    % 参数计算
    z_size = size(mat_a, 1);
    cov_zz = covariance(1:z_size, 1:z_size);
    cov_zr = covariance(1:z_size, z_size+1:end);
    cov_rr = covariance(z_size+1:end, z_size+1:end);

    % 输出准备
    moments = cell(out_size, 1);
    hnorms = zeros(out_size, 1);

    % 计算
    stable_zz = dlyap(mat_a, cov_zz);
    stable_zr1 = mat_a*stable_zz*(mat_c.') + cov_zr;
    % 零阶矩
    moments{1} = mat_c*stable_zz*(mat_c.') + cov_rr;
    [hnorms(1), ~] = anaNorm(moments{1});
    % 非零阶矩
    for iter = 1:out_size-1
        moments{iter+1} = mat_c*mpower(mat_a, iter-1)*stable_zr1;
        [hnorms(iter+1), ~] = anaNorm(moments{iter+1});
    end

end

