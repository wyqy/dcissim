function [err_rrs, err_tt, ret_cor] = errCorrelation(sys_original, sys_compared, compare_size)
%ERRCORRELATION 根据协方差和(渐近稳定)系统参数计算稳态输出矩
% 同时返回计算结果和Hnorm-2对比
% 可以同时对新息形系统进行比较
    
    % 参数提取
    x_size = size(sys_original.A, 1);
    y_size = size(sys_original.C, 1);
    mat_a_original = sys_original.A;                                                        mat_a_compared = sys_compared.A;
    mat_c_original = sys_original.C;                                                        mat_c_compared = sys_compared.C;
    cov_ww_original = sys_original.cov(1:x_size, 1:x_size);                                 cov_ww_compared = sys_compared.cov(1:x_size, 1:x_size);
    cov_vv_original = sys_original.cov(x_size+1:x_size+y_size, x_size+1:x_size+y_size);     cov_vv_compared = sys_compared.cov(x_size+1:x_size+y_size, x_size+1:x_size+y_size);
    cov_tt_original = sys_original.cov(x_size+y_size+1:end, x_size+y_size+1:end);           cov_tt_compared = sys_compared.cov(x_size+y_size+1:end, x_size+y_size+1:end);

    
    % 计算相关
    % zz
    cor_zz_original = dlyap(mat_a_original, cov_ww_original);
    cor_zz_compared = dlyap(mat_a_compared, cov_ww_compared);
    % rr0
    cor_rr0_original = mat_c_original*cor_zz_original*mat_c_original.' + cov_vv_original;
    cor_rr0_compared = mat_c_compared*cor_zz_compared*mat_c_compared.' + cov_vv_compared;
    % rri
    cor_rri_original_cell = cell(compare_size, 1);
    cor_rri_compared_cell = cell(compare_size, 1);
    for iter = 1:compare_size
        cor_rri_original_cell{iter} = mat_c_original*mpower(mat_a_original, iter)*cor_zz_original*mat_c_original.';
        cor_rri_compared_cell{iter} = mat_c_compared*mpower(mat_a_compared, iter)*cor_zz_compared*mat_c_compared.';
    end
    
    % 计算误差
    % err_rrs
    err_rrs = norm(cor_rr0_compared - cor_rr0_original, 'fro') ./ norm(cor_rr0_original, 'fro');
    for iter = 1:compare_size
        err_rrs = err_rrs + norm(cor_rri_compared_cell{iter} - cor_rri_original_cell{iter}, 'fro') ./ norm(cor_rri_original_cell{iter}, 'fro');
    end
    err_rrs = err_rrs ./ compare_size;

    % err_tt
    err_tt = norm(cov_tt_compared - cov_tt_original, 'fro') ./ norm(cov_tt_original, 'fro');
    
    ret_cor = struct('original', cor_rri_original_cell, 'compared', cor_rri_compared_cell);

end

