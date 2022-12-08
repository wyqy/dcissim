function [mat_b_sim, mat_d_sim, mat_x_sim] = idenBDX(y_isim, u_isim, mat_s, mat_a_sim, mat_c_sim, x_size_sim, sim_lsqnonlin_type, sim_ss_d_type)
%IDENBDX 辨识B, D, X矩阵

    % 替换参数
    if strcmp(sim_ss_d_type, 'null'), sim_ss_d_type = 0;
    else, sim_ss_d_type = 1; end

    % 不同辨识方法 (选择一种即可)
    switch sim_lsqnonlin_type
        case 'analytical'  % 第一种方法: 使用解析方法求解
            [mat_b_sim, mat_d_sim, mat_x_sim] = decomposeSolver(y_isim, u_isim, mat_s, mat_a_sim, mat_c_sim, x_size_sim, sim_ss_d_type);
        case 'optimize'  % 第二种方法: 使用Optimization Toolbox求解(lsqnonlin)
            [mat_b_sim, mat_d_sim, mat_x_sim] = optimalSolver(y_isim, u_isim, mat_s, mat_a_sim, mat_c_sim, x_size_sim, sim_ss_d_type);
        otherwise, mat_b_sim = 0; mat_d_sim = 0; mat_x_sim = 0;
    end

end

function [mat_b_sim, mat_d_sim, mat_x_sim] = decomposeSolver(y_isim, u_isim, mat_s, mat_a_sim, mat_c_sim, x_size_sim, sim_ss_d_type)

    % 通过对角化分解原方程, 并提取非负特征值
    [v_eig_mat, v_eig_vec] = eig_quick(mat_s);

    % 参数计算
    y_size = size(y_isim, 1);
    u_size = size(u_isim, 1);
    v_size = size(mat_s, 1);
    
    % 参数矩阵变换
    y_isim_similar = y_isim * v_eig_mat;
    u_isim_similar = u_isim * v_eig_mat;
    
    % 准备矩阵分解
    decompose_character_mat_a = cell(v_size, 1);
    for iter_ds = 1:v_size, decompose_character_mat_a{iter_ds} = decomposition(v_eig_vec(iter_ds)*eye(x_size_sim)-mat_a_sim); end

    % 计算参数矩阵
    location_base = 0;
    if sim_ss_d_type == 0  % 假定 D = 0
        mat_pinv = zeros(v_size*y_size, u_size*x_size_sim);
        for iter_v = 1:v_size
            mat_pinv(location_base+1:location_base+y_size, :) = kron(u_isim_similar(:, iter_v).', mat_c_sim / decompose_character_mat_a{iter_v});
            location_base = location_base + y_size;
        end
    else  % 不做假定
        mat_pinv = zeros(eig_size*y_size, u_size*(x_size_sim+y_size));
        for iter_v = 1:v_size
            mat_pinv(location_base+1:location_base+y_size, 1:u_size*x_size_sim) = kron(u_isim_similar(:, iter_v).', mat_c_sim / decompose_character_mat_a{iter_v});
            mat_pinv(location_base+1:location_base+y_size, u_size*x_size_sim+1:end) = kron(u_isim_similar(:, iter_v).', eye(y_size));
            location_base = location_base + y_size;
        end
    end

    % 计算最小二乘
    retVec = pinv(mat_pinv) * reshape(y_isim_similar, [v_size*y_size 1]);

    % 还原结果为矩阵
    if sim_ss_d_type == 0  % 假定 D = 0
        mat_b_sim = reshape(retVec, [x_size_sim u_size]);
        mat_d_sim = zeros(y_size, u_size);
    else  % 不做假定
        mat_b_sim = reshape(retVec(1:x_size_sim*u_size), [x_size_sim u_size]);
        mat_d_sim = reshape(retVec(x_size_sim*u_size+1:end), [y_size u_size]);
    end

    % 计算X矩阵
    x_isim_similar = zeros(x_size_sim, v_size);
    for iter_x = 1:v_size, x_isim_similar(:, iter_x) = decompose_character_mat_a{iter_x} \ mat_b_sim * u_isim(:, iter_x); end
    mat_x_sim = x_isim_similar / v_eig_mat;

end

function [mat_b_sim, mat_d_sim, mat_x_sim] = optimalSolver(y_isim, u_isim, mat_s, mat_a_sim, mat_c_sim, x_size_sim, sim_ss_d_type)
    
    % 通过对角化分解原方程, 并提取非负特征值
    [v_eig_mat, v_eig_vec] = eig_quick(mat_s);
    v_eig_vec_angle = angle(v_eig_vec);
    v_eig_vec_nneg = v_eig_vec(v_eig_vec_angle >= 0);

    % 参数计算
    y_size = size(y_isim, 1);
    u_size = size(u_isim, 1);
    v_size = size(mat_s, 1);
    eig_size = length(v_eig_vec_nneg);
    
    % 参数矩阵变换
    y_isim_similar = y_isim * v_eig_mat;
    u_isim_similar = u_isim * v_eig_mat;
    y_isim_similar_nneg = y_isim_similar(:, v_eig_vec_angle >= 0);
    u_isim_similar_nneg = u_isim_similar(:, v_eig_vec_angle >= 0);
    
    % 准备矩阵分解
    decompose_character_mat_a = cell(v_size, 1);
    for iter_ds = 1:v_size, decompose_character_mat_a{iter_ds} = decomposition(v_eig_vec(iter_ds)*eye(x_size_sim)-mat_a_sim); end
    decompose_character_mat_a_nneg = decompose_character_mat_a(v_eig_vec_angle >= 0);
    
    % 用优化方法计算结果
    seed = rng().Seed;
    rs = RandStream('dsfmt19937', 'Seed', seed);
    if sim_ss_d_type == 0, xinit = rand(rs, x_size_sim, u_size);  % 假定 D = 0
    else, xinit = rand(rs, x_size_sim+y_size, u_size); end  % 不做假定
    op_option = optimoptions('lsqnonlin', 'Algorithm', 'trust-region-reflective', ...
        'SpecifyObjectiveGradient', true, 'Display', 'final');  % 定义优化选项
    [retMat, ~, ~, exitflag, ~] = lsqnonlin(@funcCaller, xinit, [], [], op_option);
    % 返回数值
    if exitflag >= 1 || exitflag <= 4  % 求解"成功"
        if sim_ss_d_type == 0  % 假定 D = 0
            mat_b_sim = retMat;
            mat_d_sim = zeros(y_size, u_size);
        else  % 不做假定
            mat_b_sim = retMat(1:x_size_sim, :);
            mat_d_sim = retMat(x_size_sim+1:end, :);
        end
    else  % 求解失败
        mat_b_sim = zeros(x_size_sim, u_size);
        mat_d_sim = zeros(y_size, u_size);
    end
    
    % 计算X矩阵
    x_isim_similar = zeros(x_size_sim, v_size);
    for iter_x = 1:v_size, x_isim_similar(:, iter_x) = decompose_character_mat_a{iter_x} \ mat_b_sim * u_isim(:, iter_x); end
    mat_x_sim = x_isim_similar / v_eig_mat;

    % 定义指标函数
    function [vec_fun, mat_jacob] = funcCaller(x)
        % 定义状态向量值函数
        vec_fun = zeros(eig_size*y_size, 1);
        if sim_ss_d_type == 0  % 假定 D = 0
            mat_jacob = zeros(eig_size*y_size, x_size_sim*u_size);
            mat_b_iter = x;
            mat_d_iter = zeros(y_size, u_size);
        else  % 不做假定
            mat_jacob = zeros(eig_size*y_size, (x_size_sim+y_size)*u_size);
            mat_b_iter = x(1:x_size_sim, :);
            mat_d_iter = x(x_size_sim+1:end, :);
            mat_i = eye(y_size, y_size);  % 辅助矩阵
        end
        
        for iter_eig = 1:eig_size
            % 辅助参数
            location_base = (iter_eig-1)*y_size;
            mat_aux_k = mat_c_sim / decompose_character_mat_a_nneg{iter_eig};
            mat_transfer_g = mat_aux_k * mat_b_iter + mat_d_iter;
            
            % 计算返回值F
            vec_fun(location_base+1:location_base+y_size) = y_isim_similar_nneg(:, iter_eig) - mat_transfer_g * u_isim_similar_nneg(:, iter_eig);
            
            % 计算Jacobian矩阵J
            if nargout > 1  % 检查输出参数数目
                for iter_y = 1:y_size
                    % 逐个输出计算
                    if sim_ss_d_type == 0  % 假定 D = 0
                        mat_jacob(location_base+iter_y, :) = reshape(-(mat_aux_k(iter_y, :).') * (u_isim_similar_nneg(:, iter_eig).'), [1 x_size_sim*u_size]);
                    else  % 不做假定
                        mat_jacob_slice = zeros(x_size_sim+y_size, u_size);
                        mat_jacob_slice(1:x_size_sim, :) = -(mat_aux_k(iter_y, :).') * (u_isim_similar_nneg(:, iter_eig).');
                        mat_jacob_slice(x_size_sim+1:end, :) = -mat_i(:, iter_y) * (u_isim_similar_nneg(:, iter_eig).');
                        mat_jacob(location_base+iter_y, :) = reshape(mat_jacob_slice, [1 (x_size_sim+y_size)*u_size]);
                    end
                end
            end
        end
    end

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
