function [mat_b_est, mat_d_est, xv_est] = idenBDX(yv_est, uv_est, mat_s, mat_a_est, mat_c_est, x_size_est, plant_d_type)
%IDENBDX 辨识B, D, X矩阵

    % 替换参数
    if strcmp(plant_d_type, 'null'), plant_d_type = 0;
    else, plant_d_type = 1; end

    % 通过对角化分解原方程, 并提取非负特征值
    [v_eig_mat, v_eig_vec] = eig_quick(mat_s);

    % 参数计算
    y_size = size(yv_est, 1);
    u_size = size(uv_est, 1);
    v_size = size(mat_s, 1);
    
    % 参数矩阵变换
    y_isim_similar = yv_est * v_eig_mat;
    u_isim_similar = uv_est * v_eig_mat;
    
    % 准备矩阵分解
    decompose_character_mat_a = cell(v_size, 1);
    for iter_ds = 1:v_size, decompose_character_mat_a{iter_ds} = decomposition(v_eig_vec(iter_ds)*eye(x_size_est)-mat_a_est); end

    % 计算参数矩阵
    loc_base = 0;
    if plant_d_type == 0  % 假定 D = 0
        mat_pinv = zeros(v_size*y_size, u_size*x_size_est);
        for iter_v = 1:v_size
            mat_pinv(loc_base+1:loc_base+y_size, :) = kron(u_isim_similar(:, iter_v).', mat_c_est / decompose_character_mat_a{iter_v});
            loc_base = loc_base + y_size;
        end
    else  % 不做假定
        mat_pinv = zeros(eig_size*y_size, u_size*(x_size_est+y_size));
        for iter_v = 1:v_size
            mat_pinv(loc_base+1:loc_base+y_size, 1:u_size*x_size_est) = kron(u_isim_similar(:, iter_v).', mat_c_est / decompose_character_mat_a{iter_v});
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
    for iter_x = 1:v_size, x_isim_similar(:, iter_x) = decompose_character_mat_a{iter_x} \ mat_b_est * uv_est(:, iter_x); end
    xv_est = x_isim_similar / v_eig_mat;

    % 返回实数值
    mat_b_est = real(mat_b_est);
    mat_d_est = real(mat_d_est);
    xv_est = real(xv_est);

end

function [eig_mat, eig_vec] = eig_quick(mat_s)
% 快速生成特征值和特征向量

    % 准备参数
    if mat_s(1, 1) == 1, freq_mat_s = mat_s(2:end, 2:end);
    else, freq_mat_s = mat_s; end
    freq_size = size(freq_mat_s, 1)/2;
    
    % 生成谐波特征值和特征向量
    eig_mat = zeros(freq_size*2, freq_size*2);
    eig_vec = zeros(freq_size*2, 1);
    for iter_freq = 1:freq_size
        temp_angle = atan2(freq_mat_s(iter_freq*2-1, iter_freq*2), freq_mat_s(iter_freq*2, iter_freq*2));
        eig_vec(iter_freq*2-1) = exp(-1i*temp_angle);
        eig_vec(iter_freq*2) = exp(1i*temp_angle);
        eig_mat(iter_freq*2-1:iter_freq*2, iter_freq*2-1:iter_freq*2) = [1i -1i; 1 1];
    end

    % 加直流分量
    if mat_s(1, 1) == 1
        eig_vec = [1; eig_vec];
        eig_mat = blkdiag(1, eig_mat);
    end

end
