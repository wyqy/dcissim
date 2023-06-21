function [mat_b_est, mat_d_est, xv_est] = idenBDX(yv_est, uv_est, mat_s, eig_mat, eig_vec, mat_a_est, mat_c_est, x_size_est, plant_d_type)
%IDENBDX 辨识B, D, X矩阵

    % 替换参数
    if strcmp(plant_d_type, 'null'), plant_d_type = 0;
    else, plant_d_type = 1; end

    % 参数计算
    y_size = size(yv_est, 1);
    u_size = size(uv_est, 1);
    v_size = size(mat_s, 1);
    
    % 参数矩阵变换
    y_isim_similar = yv_est * eig_mat;
    u_isim_similar = uv_est * eig_mat;
    

    % 计算参数矩阵
    loc_base = 0;
    if plant_d_type == 0  % 假定 D = 0
        mat_pinv = zeros(v_size*y_size, u_size*x_size_est);
        for iter_v = 1:v_size
            mat_pinv(loc_base+1:loc_base+y_size, :) = kron(u_isim_similar(:, iter_v).', mat_c_est / (eig_vec(iter_v)*eye(x_size_est)-mat_a_est));
            loc_base = loc_base + y_size;
        end
    else  % 不做假定
        mat_pinv = zeros(eig_size*y_size, u_size*(x_size_est+y_size));
        for iter_v = 1:v_size
            mat_pinv(loc_base+1:loc_base+y_size, 1:u_size*x_size_est) = kron(u_isim_similar(:, iter_v).', mat_c_est / (eig_vec(iter_v)*eye(x_size_est)-mat_a_est));
            mat_pinv(loc_base+1:loc_base+y_size, u_size*x_size_est+1:end) = kron(u_isim_similar(:, iter_v).', eye(y_size));
            loc_base = loc_base + y_size;
        end
    end

    % 计算最小二乘
    retVec = pinv(mat_pinv) * reshape(y_isim_similar, [v_size*y_size 1]);

    % 还原结果为矩阵
    if plant_d_type == 0  % 假定 D = 0
        mat_b_est = reshape(retVec, [x_size_est u_size]);
        mat_d_est = zeros(y_size, u_size);
    else  % 不做假定
        mat_b_est = reshape(retVec(1:x_size_est*u_size), [x_size_est u_size]);
        mat_d_est = reshape(retVec(x_size_est*u_size+1:end), [y_size u_size]);
    end

    % 计算X矩阵
    x_isim_similar = zeros(x_size_est, v_size);
    for iter_x = 1:v_size, x_isim_similar(:, iter_x) = (eig_vec(iter_x)*eye(x_size_est)-mat_a_est) \ mat_b_est * uv_est(:, iter_x); end
    xv_est = x_isim_similar / eig_mat;

    % 返回实数值
    mat_b_est = real(mat_b_est);
    mat_d_est = real(mat_d_est);
    xv_est = real(xv_est);

end

