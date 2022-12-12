function [mat_a_sim, mat_c_sim, x_size_sim] = idenACN(y_isim, u_isim, mat_s, x_size_upbound, sim_x_size, sim_x_size_type)
%IDENACN 辨识A, C矩阵和XSize参数
% 应尽量满足v_size >= (u_size+y_size)*x_size_upbound!

    % 参数计算
    y_size = size(y_isim, 1);
    u_size = size(u_isim, 1);

    % 参数定义
    para_orthogal_bound = 2;  % 最少为2
    
    % 判断是否使用正交化, 从而使用不同的方法
    % 生成U1, Y1的数据矩阵
    % 强行使用经典方法
    % mat_u_observe = normalBuilder(u_isim, mat_s, x_size_upbound);
    % mat_y_observe = normalBuilder(y_isim, mat_s, x_size_upbound);
    if x_size_upbound <= para_orthogal_bound  % 若x_size_upbound不够大
        % 使用经典方法
        mat_u_observe = normalBuilder(u_isim, mat_s, x_size_upbound);
        mat_y_observe = normalBuilder(y_isim, mat_s, x_size_upbound);
    else  % 若InXSize足够大
        % 使用Forsyth多项式正交化方法
        [~, ~, mat_u_observe] = orthogalBuilder(u_isim, mat_s, x_size_upbound);
        [~, mat_z_orthogal, mat_y_observe] = orthogalBuilder(y_isim, mat_s, x_size_upbound);
    end

    % 聚合两矩阵
    mat_z_observe = [mat_u_observe.' mat_y_observe.'];
    % 计算QR分解的三角矩阵
    mat_z_qr_mat_r = qr(mat_z_observe);
    % 提取R矩阵的右下角部分, 注意需要修正至指定尺寸
    mat_z_qr_mat_r22 = mat_z_qr_mat_r(x_size_upbound*u_size+1:end, x_size_upbound*u_size+1:end);
    if size(mat_z_qr_mat_r22, 1) >= x_size_upbound*y_size
        mat_z_qr_mat_r22 = mat_z_qr_mat_r22(1:x_size_upbound*y_size, :);
    else
        mat_z_qr_mat_r22 = [mat_z_qr_mat_r22; zeros(x_size_upbound*y_size-size(mat_z_qr_mat_r22, 1), x_size_upbound*y_size)];
    end
    % 计算R的右下角矩阵的SVD分解
    [mat_r22_svd_mat_w, mat_r22_svd_mat_s, ~] = svd(mat_z_qr_mat_r22.', 'vector');

    % 估计系统状态阶数OutXSize
    % 询问(离线)或者直接估计/固定值(在线)
    switch sim_x_size_type
        case 'ask'
            disp('The eigenvalues of H are:')
            disp(mat_r22_svd_mat_s.')
            x_size_sim = input('Enter the estimated sizes of states: ');
            if ~(isscalar(x_size_sim) && isnumeric(x_size_sim)), x_size_sim = rank(mat_z_qr_mat_r22, 1e-2*max(mat_r22_svd_mat_s)); end
        case 'estimate'
            estimate_threshold = 1e-2*mean(mat_r22_svd_mat_s(mat_r22_svd_mat_s > 1e-1*max(mat_r22_svd_mat_s)));
            x_size_sim = rank(mat_z_qr_mat_r22, estimate_threshold);
        case 'fixed'
            x_size_sim = sim_x_size;
        otherwise, x_size_sim = 1;
    end
    % 得到矩阵$$\mathbf{O}$$
    mat_o_observe = mat_r22_svd_mat_w(:, 1:x_size_sim) * sqrt(diag(mat_r22_svd_mat_s(1:x_size_sim)));
    
    % 计算A, C矩阵
    % 强行使用经典方法
    % mat_a_sim = pinv(mat_o_observe(1:(x_size_upbound-1)*y_size, :)) * mat_o_observe(y_size+1:end, :);
    % mat_c_sim = mat_o_observe(1:y_size, :);
    if x_size_upbound <= para_orthogal_bound  % 若InXSize不够大
        % 使用经典方法
        mat_a_sim = pinv(mat_o_observe(1:(x_size_upbound-1)*y_size, :)) * mat_o_observe(y_size+1:end, :);
        mat_c_sim = mat_o_observe(1:y_size, :);
    else  % 若InXSize足够大
        % 使用正交方法
        [mat_d1_orthogal, mat_d2_orthogal] = alterShift(mat_z_orthogal);
        mat_a_sim = pinv(mat_d1_orthogal * mat_o_observe(y_size+1:(x_size_upbound-1)*y_size, :)) * (mat_o_observe(2*y_size+1:end, :) - mat_d2_orthogal * mat_o_observe(1:(x_size_upbound-2)*y_size, :));
        mat_c_sim = sqrt(mat_z_orthogal{1}) * mat_o_observe(1:y_size, :);
    end
    
end

function mat_x_observe = normalBuilder(mat_x, mat_d, order)
% 直接生成"可观测形"

    x_size = size(mat_x, 1); d_size = size(mat_d, 1);
    mat_x_observe = zeros(order*x_size, d_size);
    for iter_i = 1:order  % 0, ..., i-1
        location_base = (iter_i-1)*x_size;
        mat_x_observe(location_base+1:location_base+x_size, :) = mat_x * mpower(mat_d, iter_i-1);
    end

end

function [mat_r_orthogal, mat_z_orthogal, mat_x_observe] = orthogalBuilder(mat_x, mat_d, order)
% Forsyth多项式正交化下的"可观测形", 参见(Overschee & Moor, 1996)

    % 计算(实数意义下的)正交式
    % 准备, D需要为方阵!
    x_size = size(mat_x, 1); d_size = size(mat_d, 1);
    mat_r_orthogal = cell(order, 1);  % R集合, 每个矩阵shape: (XSize, DSize)
    mat_z_orthogal = cell(order, 1);  % Z集合, 每个矩阵shape: (DSize, DSize)
    % 初值准备
    mat_r_orthogal{1} = mat_x ;
    mat_z_orthogal{1} = diag(diag(mat_r_orthogal{1} * mat_r_orthogal{1}.'));  % R0, Z0
    mat_r_orthogal{2} = mat_r_orthogal{1} * mat_d;
    mat_z_orthogal{2} = diag(diag(mat_r_orthogal{2} * mat_r_orthogal{2}.'));  % R1, Z1
    % 迭代
    for iter_i = 3:order  % 2, ..., i-1
        mat_r_orthogal{iter_i} = mat_r_orthogal{iter_i-1} * mat_d + mat_z_orthogal{iter_i-1} / mat_z_orthogal{iter_i-2} * mat_r_orthogal{iter_i-2};
        mat_z_orthogal{iter_i} = diag(diag(mat_r_orthogal{iter_i} * mat_r_orthogal{iter_i}.'));
    end

    % 计算"可观测形"
    % 准备
    mat_x_observe = zeros(order*x_size, d_size);
    % 计算
    for iter_i = 1:order  % 0, ..., i-1
        location_base = (iter_i-1)*x_size;
        mat_x_observe(location_base+1:location_base+x_size, :) = sqrt(mat_z_orthogal{iter_i}) \ mat_r_orthogal{iter_i};
    end
end

function [mat_d1, mat_d2] = alterShift(mat_z)
% 计算alternative shift-structure

    % 参数计算
    i = size(mat_z, 1); d_size = size(mat_z{1}, 1);

    % 计算D1
    mat_d1 = zeros((i-2)*d_size, (i-2)*d_size);
    for iter_i = 1:i-2
        location_base = (iter_i-1)*d_size;
        mat_d1(location_base+1:location_base+d_size, location_base+1:location_base+d_size) = sqrt(mat_z{iter_i+1}) / sqrt(mat_z{iter_i+2});
    end

    % 计算D2
    mat_d2 = zeros((i-2)*d_size, (i-2)*d_size);
    for iter_i = 1:i-2
        location_base = (iter_i-1)*d_size;
        mat_d2(location_base+1:location_base+d_size, location_base+1:location_base+d_size) = mat_z{iter_i+1} / sqrt(mat_z{iter_i} * mat_z{iter_i+2});
    end
end
