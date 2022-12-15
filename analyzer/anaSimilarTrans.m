function sim_sys_b = anaSimilarTrans(sys_a, sys_b)
%ANASIMILARTRANS 对模型b进行相似变换, 从而和a的形式相似
% 要求A, B阶数相等, 且有类似的特征值结构
   
    % 参数计算
    x_size = size(sys_a.A, 1);
    
    % 对模型a有序对角化
    [eig_mat_a, eig_value_a] = eig(sys_a.A, 'vector');
    [~, sort_ind_a] = sort(eig_value_a, 'ascend', 'ComparisonMethod', 'abs');
    eig_mat_a = eig_mat_a(:, sort_ind_a);

    % 对模型b有序对角化
    [eig_mat_b, eig_value_b] = eig(sys_b.A, 'vector');
    [~, sort_ind_b] = sort(eig_value_b, 'ascend', 'ComparisonMethod', 'abs');
    eig_mat_b = eig_mat_b(:, sort_ind_b);

    % 计算变换矩阵
    mat_t = eig_mat_b\eig_mat_a;
    % 考虑B, C矩阵的参数变换
    norm_sys_a = norm(sys_a.C, 'fro');
    norm_sys_b = norm(sys_b.C, 'fro');
    mat_t = mat_t.*(norm_sys_a/norm_sys_b);

    % 变换
    sim_sys_b = sys_b;
    sim_sys_b.A = mat_t\sys_b.A*mat_t;
    sim_sys_b.B = mat_t\sys_b.B;
    sim_sys_b.C = sys_b.C*mat_t;
    sim_sys_b.cov(1:x_size, 1:x_size) = mat_t\(sys_b.cov(1:x_size, 1:x_size))/(mat_t.');

end


