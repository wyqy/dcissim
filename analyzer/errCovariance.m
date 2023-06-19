function [err_ww, err_vv, err_tt, ret_cov] = errCovariance(sys_original, sys_compared)
%ERRCOVARIANCE 直接计算方差误差 (同时进行相似变换)
% 同时返回计算结果和Hnorm-2对比

    % 参数提取
    x_size = size(sys_original.A, 1);
    y_size = size(sys_original.C, 1);

    % 相似变换1 - 使得两系统均为对角矩阵, 且特征值从小到大依次排列(默认可对角化)
    sys_original = transDiagnal(sys_original);
    sys_compared = transDiagnal(sys_compared);

    % 计算err_ww (相似变化2 - 求解优化问题)
    cov_ww_original = sys_original.cov(1:x_size, 1:x_size);
    cov_ww_compared = sys_compared.cov(1:x_size, 1:x_size);
    opt_set = optimoptions('fminunc', 'Display', 'none', 'Algorithm', 'trust-region', 'SpecifyObjectiveGradient', true);
    [kopt, err_ww] = fminunc(@optfunc, ones(x_size, 1), opt_set);
    function [fx, jacob] = optfunc(k)
        kmat = diag(k, 0);
        fx = norm(kmat*cov_ww_compared*kmat.' - cov_ww_original, 'fro')^2;
        jacob = diag(2*(kmat*cov_ww_compared*kmat - cov_ww_original)*(kmat*cov_ww_compared + cov_ww_compared*kmat.'));
    end
    err_ww = sqrt(err_ww) ./ norm(cov_ww_original, 'fro');

    % 计算err_vv
    cov_vv_original = sys_original.cov(x_size+1:x_size+y_size, x_size+1:x_size+y_size);
    cov_vv_compared = sys_compared.cov(x_size+1:x_size+y_size, x_size+1:x_size+y_size);
    err_vv = norm(cov_vv_compared - cov_vv_original, 'fro') ./ norm(cov_vv_original, 'fro');

    % 计算err_tt
    cov_tt_original = sys_original.cov(x_size+y_size+1:end, x_size+y_size+1:end);
    cov_tt_compared = sys_compared.cov(x_size+y_size+1:end, x_size+y_size+1:end);
    err_tt = norm(cov_tt_compared - cov_tt_original, 'fro') ./ norm(cov_tt_original, 'fro');
    
    ret_cov = struct('original', struct('cov_xx', cov_ww_original, 'cov_yy', cov_vv_original, 'cov_uu', cov_tt_original), ...
            'compared', struct('cov_ww', diag(kopt)*cov_ww_compared*diag(kopt).', 'cov_yy', cov_vv_compared, 'cov_uu', cov_tt_compared));

end

function [ret_sys] = transDiagnal(ori_sys)
%ANASIMILARTRANS 对角化系统
   
    % 参数计算
    x_size = size(ori_sys.A, 1);
    ret_covxx = ori_sys.cov(1:x_size, 1:x_size);
    
    % 有序对角化A矩阵
    [mat_t, eig_a] = eig(ori_sys.A, 'vector');
    [~, sort_ind] = sort(eig_a, 'ascend', 'ComparisonMethod', 'abs');
    mat_t = mat_t(:, sort_ind);
    mat_t = inv(mat_t);

    % 相似变换
    ret_sys = ori_sys;
    ret_sys.A = mat_t*ret_sys.A/mat_t; %#ok<*MINV>
    ret_sys.B = mat_t*ret_sys.B;
    ret_sys.C = ret_sys.C/mat_t;
    ret_sys.D = ret_sys.D;
    % 协方差变换
    ret_covxx = mat_t*ret_covxx*(mat_t.');
    ret_sys.cov(1:x_size, 1:x_size) = ret_covxx;
    if isfield(ori_sys, 'covariance_struct')
        ori_sys.covariance_struct.cov_all(1:x_size, 1:x_size) = ret_covxx;
        ori_sys.covariance_struct.cov_xx = ret_covxx;
    elseif isfield(ori_sys, 'noise_para')
        ori_sys.noise_para.cov_all(1:x_size, 1:x_size) = ret_covxx;
        ori_sys.noise_para.cov_xx = ret_covxx;
    end

end

