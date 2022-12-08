function nominal_kalman = anaInnoK(covariance, inno_ee, ss_struct)
%ANAINNOK 求解innovation form中, 满足给定的输出特性和新息方差的C/A*K矩阵
% 仅用来验证! 要求A, inno_ee可逆!

    % 参数提取
    mat_a = ss_struct.A;
    mat_c = ss_struct.C;
    % 参数计算
    z_size = size(mat_a, 1); r_size = size(mat_c, 1);
    ori_zz = covariance(1:z_size, 1:z_size);
    ori_zr = covariance(1:z_size, z_size+1:end);
    ori_rr = covariance(z_size+1:end, z_size+1:end);

    % 求稳态特性
    stable_zz = dlyap(mat_a, ori_zz);
    stable_zr = mat_a*stable_zz*mat_c.' + ori_zr;  % E
    stable_rr = mat_c*stable_zz*mat_c.' + ori_rr;  % F
    % 求名义kalman增益
    nominal_kalman = (mat_c/mat_a*stable_zr - stable_rr)/inno_ee + eye(r_size);
   
end

